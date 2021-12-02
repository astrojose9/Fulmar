#!/usr/bin/python

import arviz as az

from astropy.io.registry import IORegistryError
import astropy.units as u
from astropy.table import Table
from astropy.time import Time
from astropy.timeseries import (
    aggregate_downsample, TimeSeries)
from astropy.stats import sigma_clipped_stats
import exoplanet as xo
import lightkurve as lk
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import multiprocessing
import numpy as np
import os
from os import path
from transitleastsquares import cleaned_array

import warnings

import fulmar.fulmar_constants as fulmar_constants
from fulmar.utils import (
    FulmarError,
    FulmarWarning,
    print_version,
    rjd_to_astropy_time
)

from fulmar.mission_dic_manager import read_json_dic

##############################################################################


# class FulmarWarning(Warning):
#     """ Class form warning to be displayed as
#     "FulmarWarning"
#     """

#     pass

##############################################################################


def mission_identifier(target):
    """Identifies the mission ('Kepler', 'K2' or 'TESS') from the identifier
    of the target.

    Parameters
    ----------
    target : str or int
        Name of the target as a string, e.g. "TOI-175" or, if mission is
        passed as the numerical identifier of the input catalog.

    Returns
    -------
    mission : str
        'Kepler', 'K2' or 'TESS'


    Raises
    ------
    TargError
        If the target was not resolved or linked to a supported mission
    """
    if not isinstance(target, str):
        target = str(target)

    # Deal with dot as a separator. e.g "TOI-175.01"
    target = target.split('.')[0]

    if target[:3].upper() == 'TIC':
        mission = 'TESS'

    elif target[:3].upper() == 'TOI':
        mission = 'TESS'

    elif target[:3].upper() == 'KIC':
        mission = 'Kepler'

    elif target[:3].upper() == 'KEP':
        mission = 'Kepler'

    elif target[:4].upper() == 'EPIC':
        mission = 'K2'

    elif target[:2].upper() == 'K2':
        mission = 'K2'

    else:
        raise ValueError(
            "targname {} could not be linked to a supported mission "
            "('Kepler', 'K2' or 'TESS')".format(str(target)))

    return mission


def target_identifier(target, mission=None):
    """Translates the target identifiers between different catalogs
    such as TIC to TOI in the case of TESS or EPIC to KIC" for K2

    Parameters
    ----------
    target : str or int
        Name of the target as a string, e.g. "TOI-175" or, if mission is
        passed as the numerical identifier of the input catalog.
    mission : str, optional
        'Kepler', 'K2' or 'TESS'

    Returns
    -------
    inputCatalogID : str
        Identifier in the format of the input catalog, e.g. "TIC307210830"
    missionCatalogID : str or None
        Identifier in the format of the mission catalog, e.g. "TOI-175".
        None if no ID was found.
    ICnum : int
        Number of the input catalog, e.g. "307210830"

    Raises
    ------
    TargError
        If the target was not resolved or linked to a supported mission
    """
    if isinstance(target, int):
        target = str(target)

        if mission == 'TESS':
            inputCatalogID = 'TIC' + target
            tic2toi = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'TIC2TOI.json'))
            try:
                missionCatalogID = tic2toi[inputCatalogID]
            except KeyError:
                warnings.warn('Could not find a TOI identifier for {}.'.format(
                    inputCatalogID) + ' TOI set to None', FulmarWarning)
                missionCatalogID = None

            ICnum = int(inputCatalogID[3:])
            warnings.warn('No prefix was passed, target is assumed to be '
                          'TIC {}'.format(ICnum), FulmarWarning)

        elif mission == 'Kepler':
            inputCatalogID = 'KIC' + target
            kic2kepler = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'KIC2Kepler.json'))

            try:
                missionCatalogID = kic2kepler[inputCatalogID]
            except KeyError:
                warnings.warn('Could not find a Kepler identifier '
                              'for {}.'.format(
                                  inputCatalogID) + ' Kepler set to None',
                              FulmarWarning)
                missionCatalogID = None
            ICnum = int(inputCatalogID[3:])
            if (ICnum < 1) or (ICnum > 13161029):
                raise ValueError("KIC ID must be in range 1 to 13161029")
            warnings.warn('No prefix was passed, target is assumed to be '
                          'KIC {}'.format(ICnum), FulmarWarning)

        elif mission == 'K2':
            inputCatalogID = 'EPIC' + target
            epic2k2 = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'EPIC2K2.json'))

            try:
                missionCatalogID = epic2k2[inputCatalogID]
            except KeyError:
                warnings.warn('Could not find a K2 identifier '
                              'for {}.'.format(
                                  inputCatalogID) + ' K2 set to None',
                              FulmarWarning)
                missionCatalogID = None

            ICnum = int(inputCatalogID[4:])
            if (ICnum < 201000001) or (ICnum > 251813738):
                raise ValueError(
                    "EPIC ID must be in range 201000001 to 251813738")
            warnings.warn('No prefix was passed, target is assumed to be '
                          'EPIC {}'.format(ICnum), FulmarWarning)

        elif mission is None:
            raise ValueError('mission parameter should be passed '
                             'when target is int.')
        else:
            raise ValueError("mission {} could not be linked to a supported "
                             "mission ('Kepler', 'K2' or 'TESS')".format(
                                 mission))

    elif isinstance(target, str):

        # Deal with dot as a separator. e.g "TOI-175.01"
        target = target.split('.')[0]

        if target[:3].upper() == 'TIC':
            inputCatalogID = 'TIC' + str(''.join(filter(str.isdigit, target)))
            tic2toi = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'TIC2TOI.json'))
            try:
                missionCatalogID = tic2toi[inputCatalogID]
            except KeyError:
                warnings.warn('Could not find a TOI identifier for {}.'.format(
                    inputCatalogID) + ' TOI set to None', FulmarWarning)
                missionCatalogID = None

            ICnum = int(inputCatalogID[3:])

        elif target[:3].upper() == 'TOI':
            missionCatalogID = 'TOI-' + \
                str(''.join(filter(str.isdigit, target)))
            toi2tic = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'TOI2TIC.json'))
            inputCatalogID = toi2tic[missionCatalogID]
            ICnum = int(inputCatalogID[3:])

        elif target[:3].upper() == 'KIC':
            inputCatalogID = 'KIC' + str(''.join(filter(str.isdigit, target)))
            kic2kepler = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'KIC2Kepler.json'))
            try:
                missionCatalogID = kic2kepler[inputCatalogID]
            except KeyError:
                warnings.warn('Could not find a Kepler identifier '
                              'for {}.'.format(
                                  inputCatalogID) + ' Kepler set to None',
                              FulmarWarning)
                missionCatalogID = None

            ICnum = int(inputCatalogID[3:])
            if (ICnum < 1) or (ICnum > 13161029):
                raise ValueError("KIC ID must be in range 1 to 13161029")

        elif target[:3].upper() == 'KEP':
            missionCatalogID = 'Kepler-' + \
                str(''.join(filter(str.isdigit, target)))
            kep2kic = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'Kepler2KIC.json'))
            inputCatalogID = kep2kic[missionCatalogID]
            ICnum = int(inputCatalogID[3:])
            if (ICnum < 1) or (ICnum > 13161029):
                raise ValueError("KIC ID must be in range 1 to 13161029")

        elif target[:4].upper() == 'EPIC':
            inputCatalogID = 'EPIC' + str(''.join(filter(str.isdigit, target)))
            epic2k2 = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'EPIC2K2.json'))

            try:
                missionCatalogID = epic2k2[inputCatalogID]
            except KeyError:
                warnings.warn('Could not find a K2 identifier '
                              'for {}.'.format(
                                  inputCatalogID) + ' K2 set to None',
                              FulmarWarning)
                missionCatalogID = None
            ICnum = int(inputCatalogID[4:])
            if (ICnum < 201000001) or (ICnum > 251813738):
                raise ValueError(
                    "EPIC ID must be in range 201000001 to 251813738")

        elif target[:2].upper() == 'K2':
            missionCatalogID = 'K2-' + \
                str(''.join(filter(str.isdigit, target[2:])))
            k22epic = read_json_dic(
                path.join(fulmar_constants.fulmar_dir, 'K22EPIC.json'))
            inputCatalogID = k22epic[missionCatalogID]
            ICnum = int(inputCatalogID[4:])
            if (ICnum < 201000001) or (ICnum > 251813738):
                raise ValueError(
                    "EPIC ID must be in range 201000001 to 251813738")
        else:
            raise ValueError(
                'targname {} could not be linked to '
                'a supported mission'.format(target))

    else:
        raise ValueError('target should be str or int')

    return inputCatalogID, missionCatalogID, ICnum


def read_lc_from_file(
        file,
        author=None,
        exptime=None,
        timeformat=None,
        colnames=None):
    """Creates a LightCurve from a file.

    Parameters
    ----------
    file : str
        Path to the file containing the light curve data.
    author : str, optional.
        Name of the pipeline used to reduce the data.
    exptime : `~astropy.units.Quantity` or float, optional.
        Exposure time of the observation, in seconds.
    timeformat : str, optional.
        Format of the Time values. Should be 'rjd', 'bkjd', 'btjd', or a valid
        `~astropy.time.Time` format. Refer to the docs here:
        (https://docs.astropy.org/en/stable/time/index.html#time-format)
    colnames : list (of str), optional.
        Names of the columns. Should have the same number of
        items as the number of columns.

    Returns
    -------
    lc : '~lightkurve.lightcurve.LightCurve'
        LightCurve object with data from the file.
    """
    if str(file).split('.')[-1] == 'fits':
        lc = lk.read(file)
    else:
        try:  # First try using astropy.table's autodetect
            t_1 = Table.read(file)
        except IORegistryError:  # Helps the astropy reader
            t_1 = Table.read(file, format='ascii', comment='#')

        if colnames is None:
            if t_1.colnames[0] == 'col1':
                t_1.rename_column('col1', 'time')
                t_1.rename_column('col2', 'flux')
                if len(t_1.colnames) > 2:
                    t_1.rename_column('col3', 'flux_err')
            colnames = t_1.colnames
        else:
            try:
                t_1 = Table(t_1, names=colnames)
            except ValueError:
                warnings.warn('number of items in colnames should match ' +
                              'the number of columns in the data',
                              FulmarWarning)

        if 'time' not in colnames:
            raise ValueError("A 'time' column is required")
        if 'flux' not in colnames:
            raise ValueError("A 'flux' column is required")
        if len(colnames) > 2 and 'flux_err' not in colnames:
            warnings.warn(
                "No 'flux_err' column was passed, 'flux_err' set to Nan",
                FulmarWarning)

        # elif t_1.colnames[0] == 'col1':
        #     t_1.rename_column('col1', 'time')
        #     t_1.rename_column('col2', 'flux')
        #     if len(t_1.colnames) > 2:
        #         t_1.rename_column('col3', 'flux_err')

        if timeformat is not None:
            t_1.rename_column('time', 'no_unit_time')
            if timeformat == 'rjd':
                ts = TimeSeries(
                    t_1, time=rjd_to_astropy_time(t_1['no_unit_time']))
            else:
                ts = TimeSeries(t_1, time=Time(
                    t_1['no_unit_time'], format=timeformat))
            ts.remove_column('no_unit_time')

            lc = lk.LightCurve(ts)
        else:
            lc = lk.LightCurve(t_1)

    if author is not None:
        try:
            lc.meta['AUTHOR']
            if lc.meta['AUTHOR'] != author:
                warnings.warn('author parameter ({}) does not match the \
                    metadata of the LightCurve file ({})'.format(
                    author, lc.meta['AUTHOR']), FulmarWarning)
        except KeyError:
            lc.meta['AUTHOR'] = author

    if exptime is not None:
        if isinstance(exptime, float):
            exptime = exptime * u.s
    else:
        dt = (lc.time[1] - lc.time[0]).to(u.s)
        exptime = round(dt.value) * u.s

    lc.meta['EXPTIME'] = exptime
    lc['exptime'] = np.ones(len(lc)) * exptime

    return lc


def normalize_lc(lc_in, unit='unscaled'):
    """Returns a normalized version of the light curve using robust stats.

    Parameters
    ----------
    lc_in : '~lightkurve.lightcurve.LightCurve'
        LightCurve object.
    unit : 'unscaled', 'percent', 'ppt', 'ppm'
        The desired relative units of the normalized light curve;
        'ppt' means 'parts per thousand', 'ppm' means 'parts per million'.

    Returns
    -------
    lc : '~lightkurve.lightcurve.LightCurve'
        A new light curve object in which ``flux`` and ``flux_err`` have
        been divided by the median flux.

    Warns
    -----
    LightkurveWarning
        If the median flux is negative or within half a standard deviation
        from zero.
    """
    lk.utils.validate_method(unit, ["unscaled", "percent", "ppt", "ppm"])
    # median_flux = np.nanmedian(lc.flux)
    # std_flux = np.nanstd(lc.flux)
    mean_flux, median_flux, std_flux = sigma_clipped_stats(lc_in.flux)

    # If the median flux is within half a standard deviation from zero, the
    # light curve is likely zero-centered and normalization makes no sense.
    if (median_flux == 0) or (
        np.isfinite(std_flux) and (np.abs(median_flux) < 0.5 * std_flux)
    ):
        warnings.warn(
            "The light curve appears to be zero-centered "
            "(median={:.2e} +/- {:.2e}); `normalize()` will divide "
            "the light curve by a value close to zero, which is "
            "probably not what you want."
            "".format(median_flux, std_flux),
            lk.LightkurveWarning,
        )
    # If the median flux is negative, normalization will invert the light
    # curve and makes no sense.
    if median_flux < 0:
        warnings.warn(
            "The light curve has a negative median flux ({:.2e});"
            " `normalize()` will therefore divide by a negative "
            "number and invert the light curve, which is probably"
            "not what you want".format(median_flux),
            lk.LightkurveWarning,
        )

    # Create a new light curve instance and normalize its values
    lc = lc_in.copy()
    lc.flux = lc.flux / median_flux
    lc.flux_err = lc.flux_err / median_flux
    if not lc.flux.unit:
        lc.flux *= u.dimensionless_unscaled
    if not lc.flux_err.unit:
        lc.flux_err *= u.dimensionless_unscaled

    # Set the desired relative (dimensionless) units
    if unit == "percent":
        lc.flux = lc.flux.to(u.percent)
        lc.flux_err = lc.flux_err.to(u.percent)
    elif unit in ("ppt", "ppm"):
        lc.flux = lc.flux.to(unit)
        lc.flux_err = lc.flux_err.to(unit)

    lc.meta["NORMALIZED"] = True
    return lc


def time_flux_err(
        timeseries,
        flux_kw='flux',
        flux_err_kw='flux_err',
        replace_nan_err=True):
    """Extracts 3 arrays with time, flux and flux_err from a TimeSeries.

    Parameters
    ----------
    timeseries : `~astropy.timeseries.TimeSeries`
        TimeSeries object
    flux_kw : str, optional
        Keyword for the column containing the flux values
        (Default: 'flux')
    flux_err_kw : str, optional
        Keyword for the column containing the flux uncertainty values
        (Default: 'flux_err')
    replace_nan_err : bool, optional
        To deal with the case where the flux_err column is full
        of NaNs. When True, creates an array filled with the standard
        deviation of the flux. When false, returns an array full of NaNs.
        (default : True)

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
    flux = np.array(timeseries[flux_kw], dtype=np.float64)
    flux_err = np.array(timeseries[flux_err_kw], dtype=np.float64)

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


def ts_binner(timeseries, bin_duration):
    """Wrap around for astropy's aggregate_downsample returning a
    `~astropy.timeseries.TimeSeries`object with Time at the center of the bins.

    Parameters
    ----------
    timeseries : `~astropy.timeseries.TimeSeries`
        TimeSeries object

    bin_duration : `~astropy.units.Quantity` or float
        Time interval for the binned time series.
        (Default is in units of days)

    Returns
    -------
    ts_binned : `~astropy.timeseries.TimeSeries`
        TimeSeries which has been binned.
        With time corresponding to the center time of all time bins
    """
    if isinstance(bin_duration, float):
        bin_duration = bin_duration * u.d

    ts_bin = aggregate_downsample(timeseries, time_bin_size=bin_duration)
    # ts_bin['time_bin_mid'] = ts_bin['time_bin_start'] + \
    #     ts_bin['time_bin_size'].to(timeseries.time.unit)
    ts_binned = TimeSeries(ts_bin,
                           time=ts_bin.time_bin_center)
    return ts_binned


def fbn(timeseries, period, epoch0, duration=3 * u.h, nbin=40):
    """fbn for "fold, bin, norm"

    Parameters
    ----------
    timeseries : `~astropy.timeseries.TimeSeries`
        TimeSeries object
    period : `~astropy.time.Time` or float
        The period to use for folding.
    epoch0 : `~astropy.units.Quantity` or float
        The time to use as the reference epoch.
    duration : `~astropy.units.Quantity` or float, optional
        Duration of the transit. (Default is 3 hours)
    nbin : int, optional
        Number of bins in the transit window. (Default is 40)

    Returns
    -------
    ts_fold : `~astropy.timeseries.TimeSeries`
        The folded time series with an extra column 'phase_norm'.
    ts_fold_bin : `~astropy.timeseries.TimeSeries`
        The folded binned time series with an extra column 'phase_norm'.

    """
    if isinstance(period, float):
        period = period * u.d

    if isinstance(epoch0, float):
        epoch0 = Time(epoch0, format=timeseries.time.format)

    if isinstance(duration, float):
        duration = duration * u.d

    ts_fold = timeseries.fold(period=period,
                              epoch_time=epoch0)
    # ts_fold_bin = aggregate_downsample(
    #     ts_fold, time_bin_size=duration * period * u.day / nbin)
    # ts_fold_bin['time_bin_mid'] = ts_fold_bin['time_bin_start'] + \
    #     ts_fold_bin['time_bin_size']
    # ts_fold_bin = aggregate_downsample(ts_fold, time_bin_size=duration.to(
    #     period.unit) * period / nbin)
    ts_fold_bin = ts_binner(ts_fold, duration.to(
        period.unit) * period / nbin)
    # Normalize the phase

    ts_fold['phase_norm'] = ts_fold.time / (period)
    # ts_fold['phase_norm'][
    #     ts_fold['phase_norm'].value <
    #     -0.3] += 1 * ts_fold['phase_norm'].unit  # For the occultation
    ts_fold_bin['phase_norm'] = ts_fold_bin['time_bin_center'] / \
        (period)
    # ts_fold_bin['phase_norm'][
    #     ts_fold_bin['phase_norm'].value <
    #     -0.3] += 1 * ts_fold_bin['phase_norm'].unit  # For the occultation

    return ts_fold, ts_fold_bin


def GP_fit(time, flux, flux_err=None, mode='rotation',
           period_min=0.2, period_max=100,
           tune=2500, draws=2500, chains=2, target_accept=0.95,
           per=None, ncores=None):
    """Uses Gaussian Processes to model stellar activity.

    Parameters
    ----------
    time : array
        array of times at which data were taken
    flux : array
        array of flux at corresponding time
    flux_err : array, optional
        array of measurment errors of the flux data.
        Defaults to np.std(flux)
    mode : 'rotation', others to be implemented, optional
        Type of stellar variablity to correct.
    period_min : float, optional
        Minimum value for the rotation period of the star. (In days)
    period_max : float, optional
        Maximum value for the rotation period of the star. (In days)
    tune : int, optional
        number of tune iterations
    draws : int, optional
        number of draws iterations
    chains : int, optional
        number of chains to sample
    target_accept : float, optional
        number should be between 0 and 1
    per : float, optional
        Estimation of the variability period.
    ncores : int, optional
        Number of cores to use for processing. (Default: all)

    Returns
    -------
    trace:

    flat_samps :

    """

    if ncores is None:
        ncores = multiprocessing.cpu_count()

    # time flux and flux_err need to be arrays with no Nans
    if not isinstance(time, np.ndarray):
        try:
            time = time.value
        except AttributeError:
            time = np.array(time)

    if not isinstance(flux, np.ndarray):
        try:
            flux = flux.value
        except AttributeError:
            flux = np.array(flux)

    if flux_err is None:
        flux_err = np.full_like(flux, np.std(flux))

    elif not isinstance(flux_err, np.ndarray):
        try:
            flux_err = flux_err.value
        except AttributeError:
            flux_err = np.array(flux_err)

    # Case 1: flux_err only contains NaNs (common with EVEREST)
    if len(cleaned_array(time, flux, flux_err)[0]) == 0:
        time, flux = cleaned_array(time, flux)
        flux_err = np.full_like(flux, np.std(flux))
    # Case 2: flux_err contains (at least some) valid data
    else:
        time, flux, flux_err = cleaned_array(time, flux, flux_err)

    time, flux, flux_err = cleaned_array(time, flux, flux_err)

    # flux should be centered around 0 for the GP
    if np.median(flux) > 0.5:
        flux = flux - 1

    if per is None:
        ls = xo.estimators.lomb_scargle_estimator(
            time, flux, max_peaks=1,
            min_period=period_min, max_period=period_max, samples_per_peak=100)
        peak = ls["peaks"][0]
        per = peak["period"]

    # Initialize the GP model
    import pymc3 as pm
    import pymc3_ext as pmx
    import aesara_theano_fallback.tensor as tt
    from celerite2.theano import terms, GaussianProcess

    with pm.Model() as model:

        # The mean flux of the time series
        mean = pm.Normal("mean", mu=0.0, sd=10.0)

        # A jitter term describing excess white noise
        log_jitter = pm.Normal(
            "log_jitter", mu=np.log(np.mean(flux_err)), sd=2.0)

        # A term to describe the non-periodic variability
        sigma = pm.InverseGamma(
            "sigma", **pmx.estimate_inverse_gamma_parameters(1.0, 5.0)
        )
        rho = pm.InverseGamma(
            "rho", **pmx.estimate_inverse_gamma_parameters(0.5, 3.0)
        )

        # The parameters of the RotationTerm kernel
        sigma_rot = pm.InverseGamma(
            "sigma_rot", **pmx.estimate_inverse_gamma_parameters(1.0, 5.0)
        )
        log_period = pm.Normal("log_period", mu=np.log(per), sd=3.0)
        period = pm.Deterministic("period", tt.exp(log_period))
        log_Q0 = pm.HalfNormal("log_Q0", sd=2.0)
        log_dQ = pm.Normal("log_dQ", mu=0.0, sd=2.0)
        f = pm.Uniform("f", lower=0.1, upper=1.0)

        # Set up the Gaussian Process model
        kernel = terms.SHOTerm(sigma=sigma, rho=rho, Q=1 / 3.0)
        kernel += terms.RotationTerm(
            sigma=sigma_rot,
            period=period,
            Q0=tt.exp(log_Q0),
            dQ=tt.exp(log_dQ),
            f=f,
        )
        gp = GaussianProcess(
            kernel,
            t=time,
            diag=flux_err ** 2 + tt.exp(2 * log_jitter),
            mean=mean,
            quiet=True,
        )

        # Compute the Gaussian Process likelihood and add it into the
        # the PyMC3 model as a "potential"
        gp.marginal("gp", observed=flux)

        # Compute the mean model prediction for plotting purposes
        pm.Deterministic("pred", gp.predict(flux))

        # Optimize to find the maximum a posteriori parameters
        map_soln = pmx.optimize()

    # plt.plot(t, y, "k", label="data")
    # plt.plot(t_bin, map_soln["pred"], color="C1", label="model")
    # plt.xlim(t.min(), t.max())
    # plt.xlim(59149,59160)
    # plt.legend(fontsize=10)
    # plt.xlabel("time [days]")
    # plt.ylabel("relative flux [ppt]")
    # _ = plt.title(target_TOI+" map model")

    # Sampling the model
    np.random.seed()
    with model:
        trace = pmx.sample(
            tune=tune,
            draws=draws,
            start=map_soln,
            chains=chains,
            cores=ncores,
            target_accept=target_accept,
            return_inferencedata=True,
        )

    az.summary(
        trace,
        var_names=[
            "f",
            "log_dQ",
            "log_Q0",
            "log_period",
            "sigma_rot",
            "rho",
            "sigma",
            "log_jitter",
            "mean",
            "pred",
        ],
    )

    flat_samps = trace.posterior.stack(sample=("chain", "draw"))
    return trace, flat_samps


def params_optimizer(timeseries, period_guess, t0_guess, depth_guess, ab, r_star, target_id, tran_window=0.25, ncores=None, mask=None):
    if ncores is None:
        ncores = multiprocessing.cpu_count()
    print('running on {} cores'.format(ncores))
#     x = ts_stitch.time.value
#     y = ts_stitch[flux_kw + '_clean'].value
#     yerr = ts_stitch[flux_err_kw+'_clean'].value

    if mask is not None:
        x, y, yerr = time_flux_err(timeseries[mask])

    x, y, yerr = time_flux_err(timeseries)

    # x = time.copy()
    # y = flux.copy()
    # yerr = flux_err.copy()

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

        # We're going to track the implied density
        pm.Deterministic("rho_circ", orbit.rho_star)

        pm.Deterministic("incl", orbit.incl)

        # Set up the mean transit model
        light_curves = xo.LimbDarkLightCurve(
            u).get_light_curve(orbit=orbit, r=ror, t=x)

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

        map_soln = pmx.optimize(map_soln, [ror, b, dur])

        map_soln = pmx.optimize(map_soln, star_params)
        map_soln = pmx.optimize(map_soln)
        map_soln = pmx.optimize()


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


def estimate_planet_mass(
        R_p,
        rho_p):
    """
    Estimates the mass of an exoplanet from its radius and density.

    Parameters
    ----------
    R_p : `~astropy.units.Quantity` or float
        Radius of the exolanet. (defaults to units of Earth radii)
    rho_p : `~astropy.units.Quantity`, float or str
        Density of the exoplanet in kg * m^-3. Can be "Earth" or "Neptune".

    Returns
    -------
    M_planet : `~astropy.units.Quantity`
        Estimated mass of the exoplanet.
    """
    dens_dic = {'earth': 5514 * (u.kg / u.m**3),
                'neptune': 1638 * (u.kg / u.m**3)}

    if isinstance(R_p, (int, float)):
        R_p = R_p * u.earthRad

    elif isinstance(R_p, u.Quantity):
        R_p = R_p.to(u.earthRad)

    else:
        raise TypeError('R_p should be `astropy.units.Quantity` or float')

    if isinstance(rho_p, (int, float)):
        rho_p = rho_p * (u.kg / u.m**3)

    elif isinstance(rho_p, str):
        if rho_p.lower() in dens_dic.keys():
            rho_p = dens_dic[rho_p.lower()]
        else:
            raise ValueError(
                'Accepted str values for rho_p are "Earth" and "Neptune.')
    else:
        raise TypeError(
            'rho_p should be `astropy.units.Quantity`, float or str.')

    M_planet = (R_p.value ** 3 * rho_p / dens_dic['earth']) * u.earthMass
    return M_planet


def estimate_semi_amplitude(
        period,
        M_star,
        M_planet=None,
        R_planet=None,
        rho_planet=None,
        inc=90 * u.deg,
        ecc=0):
    """
    Estimates the radial velocity semi-amplitude corresponding to a planet of
    given parameters.

    Parameters
    ----------
    period : `~astropy.units.Quantity` or float
        The period to use for folding. (defaults to units of days)
    M_star : `~astropy.units.Quantity` or float
        Stellar mass (defaults to units of solar masses)
    M_planet : `~astropy.units.Quantity` or float
        Mass of the exolanet. (defaults to units of Earth masses)
    R_p : `~astropy.units.Quantity` or float
        Radius of the exolanet. (defaults to units of Earth radii)
    rho_p : `~astropy.units.Quantity`, float or str
        Density of the exoplanet in kg * m^-3. Can be "Earth" or "Neptune".
    inc : `~astropy.units.Quantity` or float
        Orbital inclination (in degrees). Default: 90.
    ecc : float
        Orbital eccentricity. Default: 0.

    Returns
    -------
    K : `~astropy.units.Quantity`
        Estimated semi-amplitude of the RV.
    """
    if isinstance(period, (int, float)):
        period = period * u.d
    elif not isinstance(period, u.Quantity):
        raise TypeError(
            'shoud be ```astropy.units.Quantity` or float')

    if isinstance(M_star, (int, float)):
        M_star = M_star * u.solMass
    elif not isinstance(M_star, u.Quantity):
        raise TypeError(
            'shoud be ```astropy.units.Quantity` or float')

    if isinstance(inc, (int, float)):
        inc = inc * u.deg
    elif not isinstance(inc, u.Quantity):
        raise TypeError(
            'shoud be ```astropy.units.Quantity` or float')

    if M_planet is None:
        if R_planet is None or rho_planet is None:
            raise ValueError('R_planet and rho_planet are both required ' +
                             'when M_planet is not given')
        else:
            M_planet = estimate_planet_mass(R_planet, rho_planet)
    else:
        if R_planet is not None or rho_planet is not None:
            warnings.warn(
                'M_planet overrides R_planet and rho_planet', FulmarWarning)

        if isinstance(M_planet, (int, float)):
            M_planet = M_planet * u.earthMass

        elif not isinstance(M_planet, u.Quantity):
            raise TypeError(
                'M_planet should be `astropy.units.Quantity` or float')

    inc = inc.to(u.deg)

    M_p_jovian = M_planet.to(u.jupiterMass).value
    M_tot_solar = (M_star + M_planet).to(u.solMass).value

    # Formula (14) from Lovis & Fischer 2010
    K = 28.4329 * (u.m / u.s) * \
        M_p_jovian * \
        np.sin(inc) * np.power(M_tot_solar, -2 / 3) * \
        np.power(period.to(u.year).value, -1 / 3) / np.sqrt(1 - ecc)

    return K
