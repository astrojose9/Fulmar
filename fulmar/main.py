import os
from os import path
import numpy as np

from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.timeseries import (
    aggregate_downsample, TimeSeries)
import astropy.units as u

import lightkurve as lk
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from transitleastsquares import (
    transitleastsquares,
    cleaned_array,
    catalog_info,
    transit_mask)

import exoplanet as xo
import multiprocessing
from pathlib import Path

import warnings


##############################################################################

# FULMAR parts
import fulmar.fulmar_constants as fulmar_constants
from fulmar.func import (
    read_lc_from_file,
    normalize_lc,
    ts_binner,
    fbn,
    GP_fit
)
from fulmar.utils import (
    FulmarWarning,
    rjd_to_astropy_time
)

from fulmar.mission_dic_manager import read_json_dic

# # FONCTIONS A CLASSIFIER ######################################################################################

# utils?


# class FulmarWarning(Warning):
#     """ Class form warning to be displayed as
#     "FulmarWarning"
#     """

#     pass


# def warning_on_one_line(message, category, filename, lineno,
#                         file=None, line=None):
#     return ' %s:%s: %s: %s' % (filename, lineno, category.__name__, message)


# warnings.formatwarning = warning_on_one_line


# # FONCTIONS A CLASSIFIER ######################################################################################


class target:
    "A LightCurve class with its associated parameters"

    def __init__(self, targname, mission=None):

        targname = str(targname)

        if mission is not None:
            self.mission = fulmar_constants.MISSION_DIC[mission.lower()]

        else:
            if targname.isdigit():
                warnings.warn('Please add the catalog prefix and/or \
                    the mission parameter to your input.', FulmarWarning)
            else:  # Detect mission from the prefix
                self.target_identifier(targname, None)

        self.mission = fulmar_constants.MISSION_DIC[self.mission.lower()]

        # Identifiers and stellar parameters
        if self.mission == 'TESS':
            # TESS identifiers
            self.TIC, self.TOI, self.TIC_num = self.target_identifier(
                targname, self.mission)

            # Stellar parameters from TIC
            self.ab, self.M_star, self.M_star_min, self.M_star_max, self.R_star, self.R_star_min, self.R_star_max = catalog_info(
                TIC_ID=self.TIC_num)
        elif self.mission == 'Kepler':
            # Kepler identifiers
            self.KIC, self.kep, self.KIC_num = self.target_identifier(
                targname, self.mission)

            # Stellar parameters from KIC
            self.ab, self.M_star, self.M_star_min, self.M_star_max, self.R_star, self.R_star_min, self.R_star_max = catalog_info(
                KIC_ID=self.KIC_num)
        elif self.mission == 'K2':
            # K2 identifiers
            self.EPIC, self.K2, self.EPIC_num = self.target_identifier(
                targname, self.mission)

            # Stellar parameters from EPIC
            self.ab, self.M_star, self.M_star_min, self.M_star_max, self.R_star, self.R_star_min, self.R_star_max = catalog_info(
                EPIC_ID=self.EPIC_num)
        else:
            warnings.warn(self.mission, 'mission is not supported... yet(?)',
                          FulmarWarning)

        self.author = ('SPOC', 'TESS_SPOC', 'Kepler', 'K2', 'EVEREST')

        # Flux keywords
        self.flux_kw = 'flux'
        self.flux_err_kw = 'flux_err'

        # Lightcurve folder
        if self.mission == 'TESS':
            self.lc_folder = os.getcwd() + '/{}/'.format(self.TOI)
        elif self.mission == 'Kepler':
            self.lc_folder = os.getcwd() + '/{}/'.format(self.kep)
        elif self.mission == 'K2':
            self.lc_folder = os.getcwd() + '/{}/'.format(self.K2)
        else:
            warnings.warn(self.mission, 'mission is not supported... yet(?)',
                          FulmarWarning)

    def target_identifier(self, target, mission=None):
        """Translate the target identifiers between different catalogs
        such as TIC to TOI in the case of TESS or EPIC to K" for K2
        Updates the mission parameter in case it wasn't passed by the user.
        """
        if target[:3].upper() == 'TIC':
            inputCatalogID = 'TIC' + str(''.join(filter(str.isdigit, target)))
            tic2toi = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'TIC2TOI.json'))
            missionCatalogID = tic2toi[inputCatalogID]
            ICnum = int(inputCatalogID[3:])
            self.mission = 'TESS'

        elif target[:3].upper() == 'TOI':
            missionCatalogID = 'TOI-' + \
                str(''.join(filter(str.isdigit, target)))
            toi2tic = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'TOI2TIC.json'))
            inputCatalogID = toi2tic[missionCatalogID]
            ICnum = int(inputCatalogID[3:])
            self.mission = 'TESS'

        elif target[:3].upper() == 'KIC':
            inputCatalogID = 'KIC' + str(''.join(filter(str.isdigit, target)))
            kic2kepler = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'KIC2Kepler.json'))
            missionCatalogID = kic2kepler[inputCatalogID]
            ICnum = int(inputCatalogID[3:])
            self.mission = 'Kepler'

        elif target[:3].upper() == 'KEP':
            missionCatalogID = 'Kepler-' + \
                str(''.join(filter(str.isdigit, target)))
            kep2kic = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'Kepler2KIC.json'))
            inputCatalogID = kep2kic[missionCatalogID]
            ICnum = int(inputCatalogID[3:])
            self.mission = 'Kepler'

        elif target[:4].upper() == 'EPIC':
            inputCatalogID = 'EPIC' + str(''.join(filter(str.isdigit, target)))
            epic2k2 = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'EPIC2K2.json'))
            missionCatalogID = epic2k2[inputCatalogID]
            ICnum = int(inputCatalogID[4:])
            self.mission = 'K2'

        elif target[:2].upper() == 'K2':
            missionCatalogID = 'K2-' + \
                str(''.join(filter(str.isdigit, target[2:])))
            k22epic = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'K22EPIC.json'))
            inputCatalogID = k22epic[missionCatalogID]
            ICnum = int(inputCatalogID[4:])
            self.mission = 'K2'

        elif target.isdigit():
            if mission == 'TESS':
                inputCatalogID = 'TIC' + target
                tic2toi = read_json_dic(
                    path.join(fulmar_constants.fulmar_dir, 'TIC2TOI.json'))
                missionCatalogID = tic2toi[inputCatalogID]
                ICnum = int(inputCatalogID[3:])
                warnings.warn('No prefix was passed, target is assumed to be \
                    TIC {}'.format(ICnum), FulmarWarning)

            elif mission == 'Kepler':
                inputCatalogID = 'KIC' + target
                kic2kepler = read_json_dic(
                    path.join(fulmar_constants.fulmar_dir, 'KIC2Kepler.json'))
                missionCatalogID = kic2kepler[inputCatalogID]
                ICnum = int(inputCatalogID[3:])
                warnings.warn('No prefix was passed, target is assumed to be \
                    KIC {}'.format(ICnum), FulmarWarning)

            elif mission == 'K2':
                inputCatalogID = 'EPIC' + target
                epic2k2 = read_json_dic(
                    path.join(fulmar_constants.fulmar_dir, 'EPIC2K2.json'))
                missionCatalogID = epic2k2[inputCatalogID]
                ICnum = int(inputCatalogID[4:])
                warnings.warn('No prefix was passed, target is assumed to be \
                    EPIC {}'.format(ICnum), FulmarWarning)

        return inputCatalogID, missionCatalogID, ICnum

    # Setters

    def set_lc_folder(self, inptfolder):
        """Set the folder in which the lightcurves are downloaded"""
        self.lc_folder = inptfolder

    def set_flux_kw(self, flux_kw):
        if isinstance(flux_kw, str):
            self.flux_kw = flux_kw
        else:
            warnings.warn('flux_kw parameter should be string', TypeError)

    def set_flux_err_kw(self, flux_err_kw):
        if isinstance(flux_err_kw, str):
            self.flux_err_kw = flux_err_kw
        else:
            warnings.warn('flux_err_kw parameter should be string', TypeError)

    # Retrieving data

    def search_data(self, author=None, exptime=None, download=False):
        """Search for available lightcurves ofr the target"""
        self.exptime = exptime

        # Author for lightcurves ‘Kepler’, ‘K2’, and ‘SPOC’ are the officials
        if author is None:
            self.author = ('SPOC', 'TESS_SPOC', 'Kepler', 'K2', 'EVEREST')
        else:
            self.author = author

        if self.mission == 'TESS':
            self.srch = lk.search.search_lightcurve(
                self.TIC, exptime=self.exptime,
                author=self.author, mission=self.mission)

        elif self.mission == 'Kepler':
            self.srch = lk.search.search_lightcurve(
                self.KIC, exptime=self.exptime,
                author=self.author, mission=self.mission)

        elif self.mission == 'K2':
            self.srch = lk.search.search_lightcurve(
                self.EPIC, exptime=self.exptime,
                author=self.author, mission=self.mission)

        else:
            warnings.warn(self.mission, 'mission is not supported... yet(?)',
                          FulmarWarning)

        if download is True:
            self.srch.download_all(download_dir=self.lc_folder)
            self.lc_files = [path.as_posix() for path in Path(
                self.lc_folder).rglob('*lc.fits')]

        sctrs = np.unique(self.srch.mission)
        self.author = np.unique(self.srch.author)
        self.exptime = np.unique(self.srch.exptime)
        print(sctrs)
        return self.srch

    def build_lightcurve(self, filelist=None, author=None, exptime=None,
                         colnames=None):
        """Build the lightcurve by reading selected data and stitching
        Creates a LightCurve from a file.
        Parameters
        ----------
        filelist : list or str (optional)
            List of paths/path to the file containing the light curve data.
        author : str (optional)
            Name of the pipeline used to reduce the data.
        exptime : float (optional)
            Exposure time of the observation, in seconds.
        colnames : list (of str) (optional)
            Names of the columns. Should have the same number of
            items as the number of columns.
        Examples
        --------

        Returns
        -------
        ts_stitch : 'TimeSeries'
            TimeSeries object combining data from selected lightcurves.
        """
        lc_col = lk.LightCurveCollection([])

        if filelist is None:
            srch = self.search_data(author=author,
                                    exptime=exptime, download=True)
            lc_col = srch.download_all()

            if author is None:
                try:
                    if len(self.author) > 1:
                        warnings.warn("You are combining data from different "
                                      "pipelines ({})), which is probably not "
                                      "what you want. If such is your goal, "
                                      "please provide an 'author' "
                                      "parameter.".format(self.author),
                                      FulmarWarning)
                except TypeError:
                    pass

            if exptime is None:
                try:
                    if len(self.exptime) > 1:
                        warnings.warn("You are combining data with different "
                                      "exposure times ({}), which is probably "
                                      "not what you want. If such is your "
                                      "goal, please provide an 'exptime' "
                                      "parameter.".format(self.exptime),
                                      FulmarWarning)
                except TypeError:
                    pass

        else:
            for f in filelist:
                lc_col.append(read_lc_from_file(f, author=author,
                                                exptime=exptime,
                                                colnames=colnames))
        stitched_lc = lc_col.stitch(corrector_func=lambda x: normalize_lc(x))
        self.ts_stitch = TimeSeries(stitched_lc)
        return self.ts_stitch

    # Work the data

    def mask_outliers(self, timeseries=None, sigma=3,
                      sigma_lower=None, sigma_upper=None):
        """Creates a mask to remove outliers from the lightcurve
        with special care to avoid removing transits.
        Parameters
        ----------
        timeseries : astropy.TimeSeries or astropy.Table (optional)
            TimeSeries ot Table object containing the data to filter
        Examples
        --------

        Returns
        -------
        clean : np.array
            mask where outliers are marked as "false" """

        # # Select which flux to work on. Expected values are norm or corr
        # if flux_mode == 'norm':
        #     flux = self.ts_stitch[self.flux_kw + '_norm']
        # elif flux_mode == 'corr':
        #     flux = self.ts_stitch[self.ts_stitch[self.flux_kw + '_corr']]
        if timeseries is None:
            flux = self.ts_stitch[self.flux_kw]

        else:
            flux = timeseries[self.flux_kw]

        if len(flux) == 0:
            warnings.warn('flux is empty. Try updating the flux_kw',
                          FulmarWarning)

        robustmean, robustmedian, robustrms = sigma_clipped_stats(
            flux, sigma=sigma, maxiters=10)

        # remove bottom outliers with a sigma clipping
        # Indices of invalid data
        cleanlow = np.where((flux - robustmean) / robustrms < -sigma)[0]

        # finding 3 consecutive outliers (for transits)
        diff = -(cleanlow - np.roll(cleanlow, -1))
        diff2 = (cleanlow - np.roll(cleanlow, 1))
        flux4 = flux.copy()

        for i in range(0, len(cleanlow)):

            if np.logical_and(np.not_equal(diff[i], 1),
                              np.not_equal(diff2[i],
                                           1)):  # true outliers are set to 0
                flux4[cleanlow[i]] = 0

        # Creating the mask
        clean = np.logical_and(np.less((flux - robustmean) / robustrms, sigma),
                               np.not_equal(flux4, 0))

        # If sigma_upper or sigma_lower are passed, combines it with the
        # Previous mask, not caring for transits
        if sigma_upper is not None or sigma_lower is not None:
            clip = sigma_clip(flux, sigma_upper=sigma_upper,
                              sigma_lower=sigma_lower)
            mask = np.logical_or(~clean, clip.mask)  # Outliers flagged as True
            return ~mask

        return clean

    def clean_subt_activity_flatten(
            self,
            timeseries=None,
            sigma=3,
            wl=1501,
            polyorder=2,
            return_trend=False,
            remove_outliers=True,
            break_tolerance=5,
            niters=3,
            mask=None):
        """Removes the low frequency trend using scipy's Savitzky-Golay filter.
        This method wraps `scipy.signal.savgol_filter`.
        Parameters
        ----------
        sigma : int
            Number of sigma above which to remove outliers from the flatten
        timeseries : 'TimeSeries' (optional)
            TimeSeries ot Table object containing the data to filter
        wl : int
            Window_length
            The length of the filter window (i.e. the number of coefficients).
            ``window_length`` must be a positive odd integer.
        polyorder : int
            The order of the polynomial used to fit the samples. ``polyorder``
            must be less than window_length.
        return_trend : bool
            If `True`, the method will return a tuple of two elements
            (ts_clean, trend_ts) where trend_ts is the removed trend.
        remove_outliers : bool
            If 'True', the method uses mask_outliers to created a mask of valid
            datapoints to be applied to the products before returning them.
        break_tolerance : int
            If there are large gaps in time, flatten will split the flux into
            several sub-lightcurves and apply `savgol_filter` to each
            individually. A gap is defined as a period in time larger than
            `break_tolerance` times the median gap.  To disable this feature,
            set `break_tolerance` to None.
        niters : int
            Number of iterations to iteratively sigma clip and flatten. If more
            than one, will perform the flatten several times,
            removing outliers each time.

        mask : boolean array with length of self.time
            Boolean array to mask data with before flattening. Flux values
            where mask is True will not be used to flatten the data. An
            interpolated result will be provided for these points. Use this
            mask to remove data you want to preserve, e.g. transits.
        Returns
        -------
        ts_clean : 'TimeSeries'
            New TimeSeries object with long-term trends removed.
        If ``return_trend`` is set to ``True``, this method will also return:
        trend_ts : 'TimeSeries'
            New TimeSeries object containing the trend that was removed from
            the flux.
        """
        if timeseries is None:
            self.ts_clean = self.ts_stitch.copy()

        else:
            self.ts_clean = timeseries.copy()

        flux = self.ts_stitch[self.flux_kw]
        lc = lk.LightCurve(self.ts_stitch)

        # Inital robust stats
        robustmean1, robustmedian1, robustrms1 = sigma_clipped_stats(
            flux, sigma=sigma, maxiters=10, mask=mask)

        # prefiltering
        # lc = lk.LightCurve(time=self.ts_stitch.time.value,
        #                    flux=self.ts_stitch[self.flux_kw + '_norm'].value,
        #                    flux_err=self.ts_stitch[self.flux_err_kw].value)

        clc = lc.flatten(window_length=wl,
                         polyorder=polyorder,
                         break_tolerance=5,
                         sigma=sigma,
                         mask=mask)

        flux_filtered1 = clc.flux

        # sigma clipping for big transits
        clip = sigma_clip(flux_filtered1, sigma)

        finflat = lc.flatten(window_length=wl, polyorder=polyorder,
                             break_tolerance=5, sigma=sigma,
                             return_trend=True, mask=clip.mask)

        clc1 = finflat[0]
        trend_ts = TimeSeries(finflat[1])

        flux_filtered = clc1.flux

        # Final robust stats
        robustmean, robustmedian, robustrms = sigma_clipped_stats(
            flux_filtered, sigma=sigma, maxiters=10, mask=mask)

        # Warn the user if std_dev is bigger than the initial.
        if robustrms > robustrms1:
            warnings.warn('Standard deviaton of the flux did not decrease \
                after filtering (before : {}, after : {}). Try using a mask\
                for the transits' .format(robustrms1, robustrms),
                          FulmarWarning)

        self.ts_clean[self.flux_kw] = flux_filtered

        if remove_outliers is True:
            out_mask = self.mask_outliers(timeseries=self.ts_clean,
                                          sigma=sigma)
            self.ts_clean = self.ts_clean[out_mask]

            trend_ts = trend_ts[out_mask]

        if return_trend is True:
            trend_ts = trend_ts[out_mask]

            return self.ts_clean, trend_ts

        return self.ts_clean

    def clean_subt_activity_GP(
            self,
            timeseries=None,
            bin_duration=40 * u.min,
            period_min=0.2,
            period_max=100.0,
            tune=2500,
            draws=2500,
            chains=2,
            target_accept=0.95,
            ncores=None,
            return_trend=False,
            remove_outliers=True,
            sigma_out=3,
            mask=None):
        """Corrects the stellar rotation using GP
        Parameters
        ----------
        timeseries : 'TimeSeries' (optional)
            TimeSeries ot Table object containing the data to filter
        bin_duration : 'astropy.units.Quantity' or float
            Time interval for the binned time series.
            (Default is in units of days)
        period_min :
            ###########################################################################
        period_max :
            ###########################################################################
        ncores : int (optional)
            Number of cores to use for processing. (Default: all)
        return_trend : bool
            If `True`, the method will return a tuple of two elements
            (ts_clean, trend_ts) where trend_ts is the removed trend.
        remove_outliers : bool
            If 'True', the method uses mask_outliers to created a mask of valid
            datapoints to be applied to the products before returning them.
        sigma_out : int
            Number of sigma above which to remove outliers from the flatten
        mask : boolean array with length of self.time
            Boolean array to mask data with before flattening. Flux values
            where mask is True will not be used to flatten the data. An
            interpolated result will be provided for these points. Use this
            mask to remove data you want to preserve, e.g. transits.
        Returns
        -------
        ts_clean : 'TimeSeries'
            New TimeSeries object with long-term trends removed.
        If ``return_trend`` is set to ``True``, this method will also return:
        trend_ts : 'TimeSeries'
            New TimeSeries object containing the trend that was removed from
            the flux."""
        if ncores is None:
            ncores = multiprocessing.cpu_count()

        if timeseries is None:
            self.ts_clean = self.ts_stitch.copy()

        else:
            self.ts_clean = timeseries.copy()

        time = self.ts_clean.time.value
        flux = self.ts_clean[self.flux_kw]

        if isinstance(bin_duration, float):
            bin_duration = bin_duration * u.d

        # binning
        if mask is not None:
            self.ts_binned = ts_binner(self.ts_clean[mask], bin_duration)
        else:
            self.ts_binned = ts_binner(self.ts_clean, bin_duration)

        t_bin = self.ts_binned['time_bin_mid'].value
        y_bin = self.ts_binned[self.flux_kw].value
        y_err_bin = self.ts_binned[self.flux_err_kw].value
        # print(y_err_bin)

        # Case 1: flux_err only contains NaNs (common with EVEREST)
        if len(cleaned_array(t_bin, y_bin, y_err_bin)[0]) == 0:
            t_bin, y_bin = cleaned_array(t_bin, y_bin)
            y_err_bin = np.full_like(y_bin, np.std(y_bin))
        # Case 2: flux_err contains (at least some) valid data
        else:
            t_bin, y_bin, y_err_bin = cleaned_array(t_bin, y_bin, y_err_bin)

        t_bin, y_bin, y_err_bin = cleaned_array(t_bin, y_bin, y_err_bin)

        # Guess the rotation period
        time1, flux1 = cleaned_array(time, flux)
        ls = xo.estimators.lomb_scargle_estimator(
            time1,
            flux1,
            max_peaks=1,
            min_period=period_min,
            max_period=period_max,
            samples_per_peak=100)

        peak = ls["peaks"][0]

        print('guessed period is {} days'.format(peak['period']))

        # GP_fit
        trace, flat_samps = GP_fit(
            t_bin, y_bin, y_err_bin, per=peak['period'], ncores=ncores)
        gp_mod = np.median(flat_samps["pred"].values, axis=1)

        # interpolate the model to the data, as it is currently binned
        gp_int = np.interp(time, t_bin, gp_mod + 1)
        flux_filtered = flux / gp_int

        self.ts_clean[self.flux_kw] = flux_filtered

        trend_ts = self.ts_clean.copy()
        trend_ts[self.flux_kw] = gp_int

        if remove_outliers is True:
            mask = self.mask_outliers(timeseries=self.ts_clean,
                                      sigma=sigma_out)
            self.ts_clean = self.ts_clean[mask]

            trend_ts = trend_ts[mask]

        if return_trend is True:

            return self.ts_clean, trend_ts

        return self.ts_clean

    def time_flux_err(
            self,
            timeseries=None,
            cleaned=True,
            replace_nan_err=True):
        """ Returns 3 arrays with time, flux and flux_err
        Parameters
            ----------
            timeseries : 'TimeSeries'
                TimeSeries object
            cleaned : bool (optional)
                Whether the periodogram should be conducted on the cleaned
                or the stitched timeseries (default: True)
            replace_nan_err : bool
                To deal with the case where the flux_err column is full
                of NaNs. When True, creates an array filled with the standard
                deviation of the flux.
                When false, returns an array full of NaNs.
            Returns
            -------
            t : array
                Array of time where values of type NaN, None, inf, and negative
                have been removed, as well as masks.
            flux : array
                Array of flux where values of type NaN, None, inf, and negative
                have been removed, as well as masks.
            flux_err : array
                Array of flux_err where values of type NaN, None, inf,
                and negative have been removed, as well as masks.
                If the original was filled with NaNs: if replace_nan_err=True,
                returns an array filled with the standard deviation of flux.
                if replace_nan_err=False, returns an array filled with Nans.
        """
        t = timeseries.time.value
        flux = np.array(timeseries[self.flux_kw], dtype=np.float64)
        flux_err = np.array(timeseries[self.flux_err_kw], dtype=np.float64)

        # The case where flux_err is filled with Nan
        if len(cleaned_array(t, flux, flux_err)[0]) == 0:
            t, flux = cleaned_array(t, flux)

            if replace_nan_err is True:
                flux_err = np.full_like(flux, np.std(flux))
            else:
                flux_err = np.nan * flux
        # Case 2: y_err contains (at least some) valid data
        else:
            t, flux, flux_err = cleaned_array(t, flux, flux_err)

        return t, flux, flux_err

    def plot_transitcheck(self, best_period, epoch0, duration,
                          savefig=False, periodn=None):
        """
        Plots a transitcheck
        ts_fold : folded timeseries
        ts_fold_bin : binned folded timeseries
        """
        ts_fold, ts_fold_bin = fbn(
            self.ts_stitch, best_period, epoch0, duration)

        ts_fold['phase_norm'] = ts_fold.time / (best_period * u.d)
        ts_fold['phase_norm'][
            ts_fold['phase_norm'].value <
            -0.3] += 1 * ts_fold['phase_norm'].unit  # For the occultation
        ts_fold_bin['phase_norm'] = ts_fold_bin['time_bin_mid'] / \
            (best_period * u.d)
        ts_fold_bin['phase_norm'][ts_fold_bin['phase_norm'].value <
                                  -0.3] += 1 * ts_fold_bin['phase_norm'].unit
        # For the occultation

        # Plots the graphs
        fig = plt.figure(figsize=[9.6, 4.8], constrained_layout=True)
        gs = GridSpec(2, 3, figure=fig)

        ax1 = fig.add_subplot(gs[0:, :-1])
        ax1.plot(ts_fold['phase_norm'],
                 ts_fold[self.flux_kw],
                 '.k',
                 alpha=0.25,
                 # color='xkcd:charcoal',
                 marker='.',
                 linestyle='None',
                 ms=1.1)
        ax1.plot(ts_fold_bin['phase_norm'].value,
                 ts_fold_bin[self.flux_kw],
                 color='xkcd:green',
                 marker='.',
                 alpha=0.36,
                 linestyle='None',
                 ms=1.7)
        ax1.set_xlabel('Phase')
        ax1.set_xlim(-0.3, 0.7)
        # ax1.set_ylim(0.9965, 1.0024)
        ax1.set_ylabel('Flux')
        if self.mission == 'TESS':
            ax1.set_title(
                self.TOI + ' Phase folded at {0:.4f} d'.format(best_period))
        elif self.mission == 'Kepler':
            ax1.set_title(
                self.kep + ' Phase folded at {0:.4f} d'.format(best_period))
        elif self.mission == 'K2':
            ax1.set_title(
                self.K2 + ' Phase folded at {0:.4f} d'.format(best_period))
        ax2 = fig.add_subplot(gs[0, 2])
        ax2.plot(ts_fold_bin['phase_norm'].value,
                 ts_fold_bin[self.flux_kw],
                 color='xkcd:green',
                 marker='.',
                 linestyle='None',
                 ms=1.6)
        ax2.set_xlim(-2 * duration, 2 * duration)
        # ax2.set_xlim(-0.1,0.1)
        ax2.get_yaxis().get_major_formatter().set_useOffset(False)
        ax2.set_ylabel('Flux')
        ax2.set_title('Transit')

        ax3 = fig.add_subplot(gs[1, 2])
        ax3.plot(ts_fold_bin['phase_norm'].value,
                 ts_fold_bin[self.flux_kw],
                 color='xkcd:green',
                 marker='.',
                 linestyle='None',
                 ms=1.6)
        ax3.set_xlim(0.5 - 2 * duration, 0.5 + 2 * duration)
        # ax3.set_xlim(0.4,0.6)
        ax3.get_yaxis().get_major_formatter().set_useOffset(False)
        ax3.set_xlabel('Phase')
        ax3.set_ylabel('Flux')
        ax3.set_title('Occultation')
        if savefig is True:
            if periodn is None:
                periodn = 1
            plt.savefig(self.lc_folder + 'transitcheck' + str(periodn),
                        facecolor='white', dpi=240)
        plt.show()
        plt.close()

    # def tls_periodogram(
    #         self,
    #         timeseries=None,
    #         cleaned=True,
    #         period_min=None,
    #         period_max=None,
    #         n_transits_min=None,
    #         transit_depth_min=1e-5,
    #         mask=None):
    #     """Computes the tls periodogram of the selected lightcurve
    #     Parameters
    #     ----------
    #     timeseries : 'TimeSeries' (optional)
    #         TimeSeries ot Table object containing the data to filter
    #     cleaned : bool (optional)
    #         Whether the periodogram should be conducted on the cleaned or the
    #         stitched timeseries (default: True)
    #     period_min : float
    #         Minimum trial period (in units of days). If none is given,
    #         the limit is derived from the Roche limit
    #     period_max : float
    #         Maximum trial period (in units of days) (default: Half the duration
    #          of the time series)
    #     n_transits_min : int (optional)
    #         Minimum number of transits required. Overrules period_max.
    #         (default=2)
    #     mask : boolean array with length of time
    #         Boolean array to mask data, typically transits. Data where mask is
    #         "True" will not be taken into account for the periodogram.
    #     Returns
    #     -------
    #     tls_results : transitleastsquaresresults
    #     """

    #     if timeseries is None:

    #         if cleaned is True:
    #             t = self.ts_clean.time.value
    #             y = np.array(
    #                 self.ts_clean[self.flux_kw], dtype=np.float64)
    #             y_err = np.array(
    #                 self.ts_clean[self.flux_err_kw], dtype=np.float64)
    #         else:
    #             t = self.ts_stitch.time.value
    #             y = np.array(
    #                 self.ts_stitch[self.flux_kw], dtype=np.float64)
    #             y_err = np.array(
    #                 self.ts_stitch[self.flux_err_kw], dtype=np.float64)
    #     else:
    #         t = timeseries.time.value
    #         y = np.array(timeseries[self.flux_kw], dtype=np.float64)
    #         y_err = np.array(timeseries[self.flux_err_kw], dtype=np.float64)

    #     # Accounts for possible mask
    #     if mask is not None:  # intransit = True
    #         t = t[~mask]
    #         y = y[~mask]
    #         y_err = y_err[~mask]

    #     # Initialize the tls model

    #     # Case 1: y_err only contains NaNs (common with EVEREST)
    #     if len(cleaned_array(t, y, y_err)[0]) == 0:
    #         t, y = cleaned_array(t, y)
    #         model = transitleastsquares(t, y)
    #     # Case 2: y_err contains (at least some) valid data
    #     else:
    #         t, y, yerr = cleaned_array(t, y, y_err)
    #         model = transitleastsquares(t, y, y_err)

    #     if n_transits_min is None:
    #         if period_max is not None:
    #             n_transits_min = int((max(t) - min(t)) / period_max)
    #         else:
    #             n_transits_min = 2

    #     if period_max is None:
    #         period_max = (max(t) - min(t)) / n_transits_min

    #     if period_min is None:
    #         # Compute the periodogram
    #         self.tls_results = model.power(
    #             u=self.ab, R_star=self.R_star,
    #             R_star_min=self.R_star - self.R_star_min,
    #             R_star_max=self.R_star + self.R_star_max,
    #             M_star=self.M_star, M_star_min=self.M_star - self.M_star_min,
    #             M_star_max=self.M_star + self.M_star_max,
    #             period_max=period_max, n_transits_min=n_transits_min,
    #             transit_depth_min=transit_depth_min)
    #     else:
    #         # Compute the periodogram
    #         self.tls_results = model.power(
    #             u=self.ab, R_star=self.R_star,
    #             R_star_min=self.R_star - self.R_star_min,
    #             R_star_max=self.R_star + self.R_star_max,
    #             M_star=self.M_star, M_star_min=self.M_star - self.M_star_min,
    #             M_star_max=self.M_star + self.M_star_max,
    #             period_min=period_min, period_max=period_max,
    #             n_transits_min=n_transits_min,
    #             transit_depth_min=transit_depth_min)

    #     return self.tls_results
    def tls_periodogram(
            self,
            timeseries=None,
            cleaned=True,
            mask=None,
            **kwargs):
        """Computes the tls periodogram of the selected lightcurve
        Parameters
        ----------
        timeseries : 'TimeSeries' (optional)
            TimeSeries ot Table object containing the data to filter
        cleaned : bool (optional)
            Whether the periodogram should be conducted on the cleaned or the
            stitched timeseries (default: True)
        period_min : float
            Minimum trial period (in units of days). If none is given,
            the limit is derived from the Roche limit
        period_max : float
            Maximum trial period (in units of days) (default: Half the duration
             of the time series)
        n_transits_min : int (optional)
            Minimum number of transits required. Overrules period_max.
            (default=2)
        mask : boolean array with length of time
            Boolean array to mask data, typically transits. Data where mask is
            "True" will not be taken into account for the periodogram.
        Returns
        -------
        tls_results : transitleastsquaresresults
        """

        if timeseries is None:

            if cleaned is True:
                t = self.ts_clean.time.value
                y = np.array(
                    self.ts_clean[self.flux_kw], dtype=np.float64)
                y_err = np.array(
                    self.ts_clean[self.flux_err_kw], dtype=np.float64)
            else:
                t = self.ts_stitch.time.value
                y = np.array(
                    self.ts_stitch[self.flux_kw], dtype=np.float64)
                y_err = np.array(
                    self.ts_stitch[self.flux_err_kw], dtype=np.float64)
        else:
            t = timeseries.time.value
            y = np.array(timeseries[self.flux_kw], dtype=np.float64)
            y_err = np.array(timeseries[self.flux_err_kw], dtype=np.float64)

        # Accounts for possible mask
        if mask is not None:  # intransit = True
            t = t[~mask]
            y = y[~mask]
            y_err = y_err[~mask]

        # Initialize the tls model

        # Case 1: y_err only contains NaNs (common with EVEREST)
        if len(cleaned_array(t, y, y_err)[0]) == 0:
            t, y = cleaned_array(t, y)
            model = transitleastsquares(t, y)
        # Case 2: y_err contains (at least some) valid data
        else:
            t, y, y_err = cleaned_array(t, y, y_err)
            model = transitleastsquares(t, y, y_err)


        # Compute the periodogram
        self.tls_results = model.power(**kwargs)

        return self.tls_results

# # Create a transit/planet object?
# class transit:
#     "A transit object with its associated parameters"


def params_optimizer(time, flux, flux_err, period_guess, t0_guess, depth_guess, ab, r_star, target_id, tran_window=0.25, ncores=None, mask=None):
    if ncores is None:
        ncores = multiprocessing.cpu_count()
    print('running on {} cores'.format(ncores))
#     x = ts_stitch.time.value
#     y = ts_stitch[flux_kw + '_clean'].value
#     yerr = ts_stitch[flux_err_kw+'_clean'].value
    if r_star is None: #####################################################################
        r_star = catalog_info(TIC_ID=target_identifier(target_id)[0])[4]

    x = time.copy()
    y = flux.copy()
    yerr = flux_err.copy()

    if mask is not None:
        x = x[mask]
        y = y[mask]
        yerr = yerr[mask]

    transitMask = (np.abs(
        (x - t0_guess + 0.5 * period_guess) % period_guess - 0.5 * period_guess) < tran_window)
    x = np.ascontiguousarray(x[transitMask])
    y = np.ascontiguousarray(y[transitMask]) - 1
    yerr = np.ascontiguousarray(yerr[transitMask])


#     plt.figure(figsize=(8, 4))
#     x_fold = (
#         x - t0_guess + 0.5 * period_guess
#     ) % period_guess - 0.5 * period_guess
#     plt.scatter(x_fold, y, c=x, s=3)
#     plt.xlabel("time since transit [days]")
#     plt.ylabel("relative flux [ppt]")
#     plt.colorbar(label="time [days]")
#     _ = plt.xlim(-tran_window, tran_window)

    import pymc3 as pm
    import aesara_theano_fallback.tensor as tt

    import pymc3_ext as pmx
#     from celerite2.theano import terms, GaussianProcess

    with pm.Model() as model:

        # Stellar parameters
        mean = pm.Normal("mean", mu=0.0, sigma=10.0)
#         u = xo.distributions.QuadLimbDark("u", testval=np.array(ab))
#         star_params = [mean, u]
        u = ab
        star_params = [mean]

        # Planet parameters
        log_ror = pm.Normal(
            "log_ror", mu=0.5 * np.log(depth_guess), sigma=10.0
        )
        ror = pm.Deterministic("ror", tt.exp(log_ror))
        r_pl = pm.Deterministic("r_pl", ror * r_star)
        # Orbital parameters
        log_period = pm.Normal(
            "log_period", mu=np.log(period_guess), sigma=1.0)
        period = pm.Deterministic("period", tt.exp(log_period))
        t0 = pm.Normal("t0", mu=t0_guess, sigma=1.0)
        log_dur = pm.Normal("log_dur", mu=np.log(0.06), sigma=10.0)
        dur = pm.Deterministic("dur", tt.exp(log_dur))
        b = xo.distributions.ImpactParameter("b", ror=ror)

        # Set up the orbit
        orbit = xo.orbits.KeplerianOrbit(
            period=period, duration=dur, ror=ror, t0=t0, b=b)

        # We're going to track the implied density for reasons that will become clear later
        pm.Deterministic("rho_circ", orbit.rho_star)

        # Set up the mean transit model
        light_curves = xo.LimbDarkLightCurve(
            u).get_light_curve(orbit=orbit, r=ror, t=x)

    #     def lc_model(t):
    #         return mean + tt.sum(lc, axis=-1
    #         )
        light_curve = pm.math.sum(light_curves, axis=-1) + mean

        # Here we track the value of the model light curve for plotting
        # purposes
        pm.Deterministic("light_curves", light_curves)

        # Finally the GP observation model
    #     gp = GaussianProcess(
    #         kernel, t=x, diag=yerr ** 2 + sigma ** 2, mean=lc_model
    #     )
    #     gp.marginal("obs", observed=y)
    #     pm.Deterministic("gp_pred", gp.predict(y))

        pm.Normal("obs", mu=light_curve, sd=np.median(yerr), observed=y)

        # Double check that everything looks good - we shouldn't see any NaNs!
        print(model.check_test_point())

        # Optimize the model
        map_soln = model.test_point
        # map_soln = pmx.optimize(map_soln, vars=[sigma, log_sigma_gp, log_rho_gp]
        #                        )
        #map_soln = pmx.optimize(map_soln, [sigma])
        map_soln = pmx.optimize(map_soln, [ror, b, dur])
        #map_soln = pmx.optimize(map_soln, noise_params)
        map_soln = pmx.optimize(map_soln, star_params)
        map_soln = pmx.optimize(map_soln)
        #map_soln = pmx.optimize()


#         plt.figure(figsize=(9, 5))
#         x_fold = (x - map_soln["t0"] + 0.5 * map_soln["period"]) % map_soln[
#             "period"
#         ] - 0.5 * map_soln["period"]
#         inds = np.argsort(x_fold)
#         plt.scatter(x_fold, 1 + y - map_soln["mean"], c=x, s=3)
#         plt.plot(x_fold[inds], 1 + map_soln["light_curves"][inds] - map_soln["mean"], "k")
#         plt.xlabel("time since transit [days]")
#         plt.ylabel("relative flux [ppt]")
#         plt.colorbar(label="time [days]")
#         _ = plt.xlim(-tran_window, tran_window)
#         plt.show()

        np.random.seed()
        with model:
            trace = pmx.sample(
                tune=2500,
                draws=2000,
                start=map_soln,
                chains=2,
                cores=ncores,
                target_accept=0.96,
                return_inferencedata=True,
            )

        import arviz as az
        az.summary(trace,
                   var_names=[
                       "period",
                       "t0",
                       "ror",
                       'dur',
                       'b',
                       #                     "u",
                       "mean"
                   ],)

        flat_samps = trace.posterior.stack(sample=("chain", "draw"))
        p = np.median(flat_samps["period"])
        t0 = np.median(flat_samps["t0"])
        dur = np.median(flat_samps["dur"])
        depth = np.median(flat_samps['ror'])**2
#         ab = tuple((np.median(flat_samps['u'], axis=-1)))

        # Plot the folded data
        x_fold = (x - t0 + 0.5 * p) % p - 0.5 * p
        plt.plot(x_fold, 1 + y, ".k", alpha=0.4, label="data", zorder=-1000)

        # Overplot the phase binned light curve
        bins = np.linspace(-0.41, 0.41, 50)
        denom, _ = np.histogram(x_fold, bins)
        num, _ = np.histogram(x_fold, bins, weights=y)
        denom[num == 0] = 1.0
        plt.plot(
            0.5 * (bins[1:] + bins[:-1]), 1 + num / denom, "o", color="C1", label="binned", alpha=0.7
        )

        # Plot the folded model
        inds = np.argsort(x_fold)
        inds = inds[np.abs(x_fold)[inds] < 0.3]
        pred = np.percentile(
            flat_samps["light_curves"][inds, 0], [16, 50, 84], axis=-1
        )
        plt.plot(x_fold[inds], 1 + pred[1], color="xkcd:green", label="model")
        art = plt.fill_between(
            x_fold[inds], 1 + pred[0], 1 + pred[2], color="xkcd:green", alpha=0.2, zorder=1000
        )
        art.set_edgecolor("none")

        # Annotate the plot with the planet's period
        txt = "period = {0:.5f} +/- {1:.5f} d".format(
            np.mean(flat_samps["period"].values), np.std(
                flat_samps["period"].values)
        )
        plt.annotate(
            txt,
            (0, 0),
            xycoords="axes fraction",
            xytext=(5, 5),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=12,
        )

        plt.legend(fontsize=10, loc=4)
        plt.title(target_id)
        plt.xlim(-0.5 * p, 0.5 * p)
        plt.xlabel("time since transit [days]")
        plt.ylabel("de-trended flux")
        _ = plt.xlim(-tran_window, tran_window)
        plt.show()

        return p, t0, dur, depth, ab, flat_samps
