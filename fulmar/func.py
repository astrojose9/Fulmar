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
import multiprocessing
import numpy as np
from transitleastsquares import cleaned_array

import warnings

from fulmar.utils import (
    FulmarWarning,
    rjd_to_astropy_time
)
##############################################################################


# class FulmarWarning(Warning):
#     """ Class form warning to be displayed as
#     "FulmarWarning"
#     """

#     pass

##############################################################################


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
    author : str (optional)
        Name of the pipeline used to reduce the data.
    exptime : float (optional)
        Exposure time of the observation, in seconds.
    timeformat : str
        Format of the Time values. Should be 'rjd', 'bkjd', 'btjd', or a valid
        astropy.time format. Refer to the docs here:
        (https://docs.astropy.org/en/stable/time/index.html#time-format)
    colnames : list (of str) (optional)
        Names of the columns. Should have the same number of
        items as the number of columns.
    Examples
    --------

    Returns
    -------
    lc : 'LightCurve'
        LightCurve object with data from the file.
    """
    if str(file).split('.')[-1] == 'fits':
        lc = lk.read(file)
    else:
        try:  # First try using astropy.table's autodetect
            t_1 = Table.read(file)
        except IORegistryError:  # Helps the astropy reader
            t_1 = Table.read(file, format='ascii', comment='#')
        if colnames is not None:
            try:
                t_1 = Table(t_1, names=colnames)
            except ValueError:
                warnings.warn('number of items in colnames should match \
                    the number of columns in the data', FulmarWarning)
        elif t_1.colnames[0] == 'col1':
            t_1.rename_column('col1', 'time')
            t_1.rename_column('col2', 'flux')
            if len(t_1.colnames) > 2:
                t_1.rename_column('col3', 'flux_err')

        if timeformat is not None:
            t_1.rename_column('time', 'no_unit_time')
            if timeformat == 'rjd':
                ts = TimeSeries(
                    t_1, time=rjd_to_astropy_time(t_1['no_unit_time']))
            else:
                ts = TimeSeries(t_1, time=Time(
                    t_1['no_unit_time'], format='timeformat'))
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
        lc.meta['EXPTIME'] = exptime

    return lc


def normalize_lc(lc_in, unit='unscaled'):
    """Returns a normalized version of the light curve.
    Using robust stats.
    Parameters
    ----------
    lc_in : LightCurve
        LightCurve object.
    unit : 'unscaled', 'percent', 'ppt', 'ppm'
        The desired relative units of the normalized light curve;
        'ppt' means 'parts per thousand', 'ppm' means 'parts per million'.
    Examples
    --------
        >>> import lightkurve as lk

        >>> lc = lk.LightCurve(time=[1, 2, 3],
                               flux=[25945.7, 25901.5, 25931.2],
                               flux_err=[6.8, 4.6, 6.2])
        >>> normalized_lc = normalize_lc(lc)
        >>> normalized_lc.flux
        <Quantity [1.00055917, 0.99885466, 1.        ]>
        >>> normalized_lc.flux_err
        <Quantity [0.00026223, 0.00017739, 0.00023909]>
    Returns
    -------
    normalized_lightcurve : `LightCurve`
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


def ts_binner(ts, bin_duration):
    """
    Wrap around for astropy's aggregate_downsample with centered time
    Parameters
        ----------
        ts : 'TimeSeries'
            TimeSeries object
        bin_duration : 'astropy.units.Quantity' or float
            Time interval for the binned time series.
            (Default is in units of days)
        Returns
        -------
        ts_binned : 'TimeSeries'
            TimeSeries which has been binned.
    """
    if isinstance(bin_duration, float):
        bin_duration = bin_duration * u.d

    ts_binned = aggregate_downsample(ts, time_bin_size=bin_duration)
    ts_binned['time_bin_mid'] = ts_binned['time_bin_start'] + \
        ts_binned['time_bin_size'].to(u.d)
    ts_binned
    return ts_binned


def fbn(ts, best_period, epoch0, duration=3 * u.h, nbin=40):
    """
    fbn for "fold, bin, norm"

    epoch0 : 'astropy.time.Time' or float
        The time to use as the reference epoch
    """
    if isinstance(best_period, float):
        best_period = best_period * u.d
    if isinstance(epoch0, float):
        epoch0 = Time(epoch0, format=ts.time.format)

    if isinstance(duration, float):
        duration = duration * u.d

    ts_fold = ts.fold(period=best_period,
                      epoch_time=epoch0)
    # ts_fold_bin = aggregate_downsample(
    #     ts_fold, time_bin_size=duration * best_period * u.day / nbin)
    # ts_fold_bin['time_bin_mid'] = ts_fold_bin['time_bin_start'] + \
    #     ts_fold_bin['time_bin_size']
    ts_fold_bin = ts_binner(ts_fold, duration.to(
        best_period.unit) * best_period / nbin)
    # Normalize the phase

    ts_fold['phase_norm'] = ts_fold.time / (best_period)
    # ts_fold['phase_norm'][
    #     ts_fold['phase_norm'].value <
    #     -0.3] += 1 * ts_fold['phase_norm'].unit  # For the occultation
    ts_fold_bin['phase_norm'] = ts_fold_bin['time_bin_mid'] / \
        (best_period)
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
        flux_err : array (optional)
            array of measurment errors of the flux data.
            Defaults to np.std(flux)
        mode : 'rotation', others to be implemented
            Type of stellar variablity to correct.
            Defaults to 'rotation'
        period_min :
            ###########################################################################
        period_max :
            ###########################################################################            
        tune : int
            number of tune iterations
        draws : int
            number of draws iterations
        chains : int
            number of chains to sample
        target_accept : float
            number should be between 0 and 1
        per : float (optional)
            Estimation of the variability period.
        ncores : int (optional)
            Number of cores to use for processing. (Default: all)
        Returns
        -------
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
