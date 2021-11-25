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
    time_flux_err,
    ts_binner,
    fbn,
    GP_fit
)
from fulmar.utils import (
    FulmarWarning,
    print_version,
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
    """
    A target object encompassing lightcurves and relevant parameters for their
    analysis.


    Parameters
    ----------
    targname : str or int
        Name of the target as a string, e.g. "TOI-175" or, if mission
    mission : str, optional
        'Kepler', 'K2', or 'TESS'

    Attributes
    ----------
    ab : tuple of floats
        Quadratic limb darkening parameters a, b.
    M_star : float
        Stellar mass (in units of solar masses)
    M_star_min : float
        1-sigma upper confidence interval on stellar mass
        (in units of solar mass)
    M_star_max : float
        1-sigma lower confidence interval on stellar mass
        (in units of solar mass)
    R_star : float
        Stellar radius (in units of solar radii).
    R_star_min : float
        1-sigma lower confidence interval on stellar radius
        (in units of solar radii)
    R_star_max : float
        1-sigma upper confidence interval on stellar radius
        (in units of solar radii)
    flux_kw : str
        Keyword for the column containing the flux values
        (Default: 'flux')
    flux_err_kw : str
        Keyword for the column containing the flux uncertainty values
        (Default: 'flux_err')

    Notes
    -----

    Examples
    --------
    >>> import fulmar
    >>> lc_targ = fulmar.target('TOI-175')
    >>> lc_targ.R_star
    array(0.31416)
    >>> lc_targ.mission
    'TESS'
    """

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
        """Sets the folder in which the lightcurves are downloaded"""
        self.lc_folder = inptfolder

    def set_flux_kw(self, flux_kw):
        """
        Sets the keyword for the column containing the flux values.
        """
        if isinstance(flux_kw, str):
            self.flux_kw = flux_kw
        else:
            warnings.warn('flux_kw parameter should be string', TypeError)

    def set_flux_err_kw(self, flux_err_kw):
        """
        Sets the keyword for the column containing the flux uncertainty values.
        """
        if isinstance(flux_err_kw, str):
            self.flux_err_kw = flux_err_kw
        else:
            warnings.warn('flux_err_kw parameter should be string', TypeError)

    # Retrieving data

    def search_data(self, author=None, exptime=None, download=False):
        """Search for available lightcurves of the target

        Parameters
        ----------
        author : str, tuple of str, or "any"
            Author of the data product (`provenance_name` in the MAST API).
            Official Kepler, K2, and TESS pipeline products have author names
            'Kepler', 'K2', and 'SPOC'.
        exptime : 'long', 'short', 'fast', or float
            'long' selects 10-min and 30-min cadence products;
            'short' selects 1-min and 2-min products;
            'fast' selects 20-sec products.
            Alternatively, you can pass the exact exposure time in seconds as
            an int or a float, e.g., ``exptime=600`` selects 10-minute cadence.
            By default, all cadence modes are returned.
        download : bool
            Whether the data should be downloaded or not.

        Returns
        -------
        srch : `SearchResult` object
            Object detailing the data products found.
        """
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
        """Build the timeseries by reading and stitching selected data.

        Parameters
        ----------
        filelist : list or str, optional
            List of paths/path to the file containing the light curve data.
        author : str, optional
            Name of the pipeline used to reduce the data.
        exptime : float, optional
            Exposure time of the observation, in seconds.
        colnames : list (of str), optional
            Names of the columns. Should have the same number of
            items as the number of columns.

        Returns
        -------
        ts_stitch : `~astropy.timeseries.TimeSeries`
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
        timeseries : `~astropy.timeseries.TimeSeries` or `~astropy.table.Table`, optional
            TimeSeries ot Table object containing the data to filter
        sigma : float, optional
            The number of standard deviations to use for the clipping limit
        sigma_lower : float, optional
            The number of standard deviations to use as the lower bound
            for the clipping limit.
        sigma_upper : float, optional
            The number of standard deviations to use as the upper bound
            for the clipping limit.

        Returns
        -------
        clean : np.array
            mask where outliers are marked as "False"
        """

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
            time_window=None,
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
        sigma : float, optional
            Number of standard deviations to use for the clipping limit.
        timeseries : `~astropy.timeseries.TimeSeries`, optional
            TimeSeries ot Table object containing the data to filter
        wl : int
            Window_length
            The length of the filter window (i.e. the number of coefficients).
            ``window_length`` must be a positive odd integer.
        time_window : '~astropy.units.Quantity' or float, optional
            Time length of the filter window. Window_lenght will be set to the
            closest odd integer taking exposition time into account.
            Overrules wl.
        polyorder : int, optional
            The order of the polynomial used to fit the samples. ``polyorder``
            must be less than window_length.
        return_trend : bool, optional
            If `True`, the method will return a tuple of two elements
            (ts_clean, trend_ts) where trend_ts is the removed trend.
        remove_outliers : bool
            If 'True', the method uses mask_outliers to created a mask of valid
            datapoints to be applied to the products before returning them.
        break_tolerance : int, optional
            If there are large gaps in time, flatten will split the flux into
            several sub-lightcurves and apply `savgol_filter` to each
            individually. A gap is defined as a period in time larger than
            `break_tolerance` times the median gap.  To disable this feature,
            set `break_tolerance` to None.
        niters : int, optional
            Number of iterations to iteratively sigma clip and flatten. If more
            than one, will perform the flatten several times,
            removing outliers each time.

        mask : boolean array with length of time, optional
            Boolean array to mask data with before flattening. Flux values
            where mask is True will not be used to flatten the data. An
            interpolated result will be provided for these points. Use this
            mask to remove data you want to preserve, e.g. transits.

        Returns
        -------
        ts_clean : `~astropy.timeseries.TimeSeries`
            New `TimeSeries` object with long-term trends removed.
        If ``return_trend`` is set to ``True``, this method will also return:
        trend_ts : `~astropy.timeseries.TimeSeries`
            New `TimeSeries` object containing the trend that was removed from
            the flux.
        """
        if timeseries is None:
            self.ts_clean = self.ts_stitch.copy()

        else:
            self.ts_clean = timeseries.copy()

        flux = self.ts_stitch[self.flux_kw]
        lc = lk.LightCurve(self.ts_stitch)

        if mask is None:
            mask = np.full_like(self.ts_clean.time.value, False, dtype=bool)

        # Inital robust stats
        robustmean1, robustmedian1, robustrms1 = sigma_clipped_stats(
            flux, sigma=sigma, maxiters=10, mask=mask)

        if time_window is not None:
            if isinstance(time_window, float):
                time_window = time_window * u.d
            # observation time interval
            dt = (self.ts_clean.time[1] - self.ts_clean.time[0]).to(
                time_window.unit)
            nobs = round((time_window / dt).value)
            # make sure the number is odd
            wl = nobs + 1 - (nobs % 2)

        # prefiltering
        # lc = lk.LightCurve(time=self.ts_stitch.time.value,
        #                    flux=self.ts_stitch[self.flux_kw + '_norm'].value,
        #                    flux_err=self.ts_stitch[self.flux_err_kw].value)

        clc = lc.flatten(window_length=wl,
                         polyorder=polyorder,
                         break_tolerance=break_tolerance,
                         sigma=sigma,
                         mask=mask)

        flux_filtered1 = clc.flux

        # sigma clipping for big transits
        clip = sigma_clip(flux_filtered1, sigma)

        finflat = lc.flatten(window_length=wl, polyorder=polyorder,
                             break_tolerance=break_tolerance, sigma=sigma,
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
            mask=None,
            store_trace=False):
        """Corrects the stellar rotation using GP

        Parameters
        ----------
        timeseries : `~astropy.timeseries.TimeSeries`, optional
            TimeSeries ot Table object containing the data to filter
        bin_duration : 'astropy.units.Quantity' or float
            Time interval for the binned time series.
            (Default is in units of days)
        period_min : float, optional
            Minimum value for the rotation period of the star. (In days)
        period_max : float, optional
            Maximum value for the rotation period of the star. (In days)
        ncores : int, optional
            Number of cores to use for processing. (Default: all)
        return_trend : bool, optional
            If `True`, the method will return a tuple of two elements
            (ts_clean, trend_ts) where trend_ts is the removed trend.
        remove_outliers : bool
            If 'True', the method uses mask_outliers to created a mask of valid
            datapoints to be applied to the products before returning them.
        sigma_out : int, optional
            Number of sigma above which to remove outliers from the flatten
        mask : boolean array with length of self.time
            Boolean array to mask data with before flattening. Flux values
            where mask is True will not be used to flatten the data. An
            interpolated result will be provided for these points. Use this
            mask to remove data you want to preserve, e.g. transits.
        store_trace : bool, optional
            If True, the posterior sampling of the GP model will be stored as
            an attribute called "activity_gp_trace" (target.activity_gp_trace).
            It can be useful to run convergence checks.

        Returns
        -------
        ts_clean : `~astropy.timeseries.TimeSeries`
            New TimeSeries object with long-term trends removed.
        If ``return_trend`` is set to ``True``, this method will also return:
        trend_ts : `~astropy.timeseries.TimeSeries`
            New TimeSeries object containing the trend that was removed from
            the flux.
        """
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

        t_bin = self.ts_binned.time.value
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

        if store_trace is True:
            self.activity_gp_trace = trace
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

    def plot_transitcheck(self, period, epoch0, duration=3 * u.h, nbin=40,
                          timeseries=None, savefig=False, fig_id=None):
        """
        Plots a transitcheck image. A visual check at a given period and epoch
        useful to probe signals detected in RV.

        Parameters
        ----------
        period : '~astropy.time.Time' or float
            The period to use for folding.
        epoch0 : '~astropy.units.Quantity' or float
            The time to use as the reference epoch.
        duration : '~astropy.units.Quantity' or float
            Duration of the transit.
        nbin : int
            Number of bins in the transit window.
        timeseries : `~astropy.timeseries.TimeSeries`, optional
            TimeSeries object
        savefig : bool
            If True, saves the resulting figure on the disk.
        fig_id : str or int
            Suffix for the filename when the figure is exported.

        """
        if timeseries is not None:
            ts_fold, ts_fold_bin = fbn(
                timeseries, period, epoch0, duration)
        else:
            ts_fold, ts_fold_bin = fbn(
                self.ts_stitch, period, epoch0, duration)

        # Wrap the phase for the occultation
        ts_fold['phase_norm'][
            ts_fold['phase_norm'].value <
            -0.3] += 1 * ts_fold['phase_norm'].unit

        ts_fold_bin['phase_norm'][ts_fold_bin['phase_norm'].value <
                                  -0.3] += 1 * ts_fold_bin['phase_norm'].unit

        # For the occultation

        # Plots the graphs
        fig = plt.figure(figsize=[9.6, 4.8], constrained_layout=True)
        gs = GridSpec(2, 3, figure=fig)

        ax1 = fig.add_subplot(gs[0:, :-1])
        ax1.plot(ts_fold['phase_norm'].value,
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
                self.TOI + ' Phase folded at {0:.4f}'.format(period))
        elif self.mission == 'Kepler':
            ax1.set_title(
                self.kep + ' Phase folded at {0:.4f}'.format(period))
        elif self.mission == 'K2':
            ax1.set_title(
                self.K2 + ' Phase folded at {0:.4f}'.format(period))
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
            if fig_id is None:
                fig_id = 1
            plt.savefig(self.lc_folder + 'transitcheck' + str(fig_id),
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
    #     timeseries : `~astropy.timeseries.TimeSeries`, optional
    #         TimeSeries ot Table object containing the data to filter
    #     cleaned : bool, optional
    #         Whether the periodogram should be conducted on the cleaned or the
    #         stitched timeseries (default: True)
    #     period_min : float
    #         Minimum trial period (in units of days). If none is given,
    #         the limit is derived from the Roche limit
    #     period_max : float
    #         Maximum trial period (in units of days) (default: Half the duration
    #          of the time series)
    #     n_transits_min : int, optional
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
        timeseries : `~astropy.timeseries.TimeSeries`, optional
            TimeSeries ot Table object containing the data to filter
        cleaned : bool, optional
            Whether the periodogram should be conducted on the cleaned or the
            stitched timeseries (default: True)
        period_min : float, optional
            Minimum trial period (in units of days). If none is given,
            the limit is derived from the Roche limit
        period_max : float, optional
            Maximum trial period (in units of days) (default: Half the duration
             of the time series)
        n_transits_min : int, optional
            Minimum number of transits required. Overrules period_max.
            (default=2)
        mask : boolean array with length of time
            Boolean array to mask data, typically transits. Data where mask is
            "True" will not be taken into account for the periodogram.

        Returns
        -------
        tls_results : `transitleastsquaresresults`
        """

        if timeseries is None:

            if cleaned is True:
                t, y, y_err = time_flux_err(
                    self.ts_clean,
                    flux_kw=self.flux_kw,
                    flux_err_kw=self.flux_err_kw)

                # t = self.ts_clean.time.value
                # y = np.array(
                #     self.ts_clean[self.flux_kw], dtype=np.float64)
                # y_err = np.array(
                #     self.ts_clean[self.flux_err_kw], dtype=np.float64)
            else:
                t, y, y_err = time_flux_err(
                    self.ts_stitch,
                    flux_kw=self.flux_kw,
                    flux_err_kw=self.flux_err_kw)
                # t = self.ts_stitch.time.value
                # y = np.array(
                #     self.ts_stitch[self.flux_kw], dtype=np.float64)
                # y_err = np.array(
                #     self.ts_stitch[self.flux_err_kw], dtype=np.float64)
        else:
            t, y, y_err = time_flux_err(
                timeseries,
                flux_kw=self.flux_kw,
                flux_err_kw=self.flux_err_kw)
            # t = timeseries.time.value
            # y = np.array(timeseries[self.flux_kw], dtype=np.float64)
            # y_err = np.array(timeseries[self.flux_err_kw], dtype=np.float64)

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


