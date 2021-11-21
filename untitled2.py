#!/usr/bin/python

import arviz as az

import astropy.units as u
from astropy.time import Time
from astropy.timeseries import (
    aggregate_downsample, TimeSeries)
from astropy.stats import sigma_clipped_stats
import exoplanet as xo
from fulmar.func import read_lc_from_file
import lightkurve as lk
import multiprocessing
import numpy as np

import warnings

############################################################


class FulmarWarning(Warning):
    """ Class form warning to be displayed as 
    "FulmarWarning"
    """

    pass


def normalize_flux(self, flux_kw, flux_err_kw):
    'Normalize flux using robust stats'
    mean_flux, median_flux, stddev = sigma_clipped_stats(self.lc[flux_kw])

    # If the median flux is negative, normalization will invert the light
    # curve and makes no sense.
    if median_flux < 0:
        warnings.warn(
            "The light curve has a negative median flux ({:.2e});"
            " `normalize_flux()` will therefore divide by a negative "
            "number and invert the light curve, which is probably"
            "not what you want".format(median_flux),
            FulmarWarning,
        )

    # Create a new light curve instance and normalize its values
    lc = self.lc.copy()
    lc.flux = lc[flux_kw] / median_flux
    lc.flux_err = lc[flux_err_kw] / median_flux

    if not lc.flux.unit:
        lc.flux *= u.dimensionless_unscaled
    if not lc.flux_err.unit:
        lc.flux_err *= u.dimensionless_unscaled

    lc.meta["NORMALIZED"] = True
    return lc

    def build_lightcurve(self):

        self.lc_files = [path.as_posix() for path in Path(
            self.lc_folder).rglob('*lc.fits')]  # PATH FROM MISSION
        self.lc_files = [f.replace(os.getcwd(), '')[1:] for f in self.lc_files]
        if (len(self.lc_files) == 0):  # Download LC
            self.download_data()

        if self.mission == 'TESS':
            lc_reader = read_TESS_LC
        elif (self.mission == 'Kepler' or self.mission == 'K2'):
            lc_reader = read_Kepler_LC
        else:
            warnings.warn(self.mission, 'mission is not supported... yet(?)',
                          FulmarWarning)

        for i, file in enumerate(self.lc_files):
            print(str(i), '/', str(len(self.lc_files)))
            self.ts = lc_reader(file)

            self.ts[self.flux_kw +
                    '_norm'] = normalize_flux(self.ts, self.flux_kw)

            # LC_plot_BJD(ts.time.value, ts[flux_kw + '_norm'].value)

            # plt.savefig(self.folder + str(i))

            try:
                ts_stitch = astropy.table.vstack(
                    [ts_stitch, self.ts], metadata_conflicts='silent')
            except NameError:
                self.ts_stitch = self.ts

            del self.ts
            self.ts_stitch = ts_stitch
            self.ts_stitch.sort('time')

            return self.ts_stitch

    def build_lightcurve(self, file_list=None, author=None, exptime=None):

        if file_list is not None:  #
            lc_col = lk.LightCurveCollection([])  # creates an empty collection
            for f in file_list:
                lc_col.append(read_lc_from_file(f, author, exptime))

        else:
            lc_col = self.download_data(author=author, exptime=exptime)
            if len(np.unique(self.srch.author)) > 1:
                if author is None:
                    warnings.warn("Data comes from different pipelines. You \
                        probably don't want to combine them. If it's the case,\
                        provide an explicit 'author' parameter", FulmarWarning)
            if len(np.unique(self.srch.exptime)) > 1:
                if exptime is None:
                    warnings.warn("Lightcurves have different exposure times. \
                        You probably don't want to combine them. If it's the \
                        case, provide an explicit 'exptime' parameter",
                                  FulmarWarning)

            lc_col.stitch


def GP_fit(time, flux, flux_err=None, mode='rotation',
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

    if flux_err is None:
        flux_err = np.std(flux)

    # flux should be centered around 0 for the GP
    if np.median(flux) > 0.5:
        flux = flux - 1

    if per is None:
        ls = xo.estimators.lomb_scargle_estimator(
            time, flux, max_peaks=1,
            min_period=0.2, max_period=100.0, samples_per_peak=100)
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
        log_Q0 = pm.HalfNormal("log_Q0", mu=0.0, sd=2.0)
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
    return flat_samps
