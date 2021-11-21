    def plot_transitcheck(self, best_period, epoch0, duration, periodn):
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
                 ts_fold[self.flux_kw + '_clean'],
                 '.k',
                 alpha=0.25,
                 # color='xkcd:charcoal',
                 marker='.',
                 linestyle='None',
                 ms=1.1)
        ax1.plot(ts_fold_bin['phase_norm'].value,
                 ts_fold_bin[self.flux_kw + '_clean'],
                 color='xkcd:green',
                 marker='.',
                 alpha=0.36,
                 linestyle='None',
                 ms=1.7)
        ax1.set_xlabel('Phase')
        ax1.set_xlim(-0.3, 0.7)
        ax1.set_ylim(0.9965, 1.0024)
        ax1.set_ylabel('Flux')
        ax1.set_title(
            self.TOI + ' Phase folded at {0:.4f} d'.format(best_period))

        ax2 = fig.add_subplot(gs[0, 2])
        ax2.plot(ts_fold_bin['phase_norm'].value,
                 ts_fold_bin[self.flux_kw + '_clean'],
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
                 ts_fold_bin[self.flux_kw + '_clean'],
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

        plt.savefig(self.folder + 'transitcheck' + periodn,
                    facecolor='white', dpi=240)
        plt.show()
        plt.close()


def target_identifier(target, mission=self.mission):
    """Translate the target identifiers between different catalogs
    such as TIC to TOI in the case of TESS or EPIC to K" for K2
    Updates the mission parameter in case it wasn't passed by the user.
    """
    if target[:3].upper() == 'TIC':
        inputCatalogID = 'TIC' + str(''.join(filter(str.isdigit, target)))
        tic2toi = read_json_dic('TIC2TOI.json')
        missionCatalogID = tic2toi[inputCatalogID]
        ICnum = int(inputCatalogID[3:])
        self.mission = 'TESS'

    elif target[:3].upper() == 'TOI':
        missionCatalogID = 'TOI-' + str(''.join(filter(str.isdigit, target)))
        toi2tic = read_json_dic('TOI2TIC.json')
        inputCatalogID = toi2tic[missionCatalogID]
        ICnum = int(inputCatalogID[3:])
        self.mission = 'TESS'

    elif target[:3].upper() == 'KIC':
        inputCatalogID = 'KIC' + str(''.join(filter(str.isdigit, target)))
        kic2kepler = read_json_dic('KIC2Kepler.json')
        missionCatalogID = kic2kepler[inputCatalogID]
        ICnum = int(inputCatalogID[3:])
        self.mission = 'Kepler'

    elif target[:3].upper() == 'KEP':
        missionCatalogID = 'Kepler-' + \
            str(''.join(filter(str.isdigit, target)))
        kep2kic = read_json_dic('Kep2KIC.json')
        inputCatalogID = kep2kic[missionCatalogID]
        ICnum = int(inputCatalogID[3:])
        self.mission = 'Kepler'

    elif target[:4].upper() == 'EPIC':
        inputCatalogID = 'EPIC' + str(''.join(filter(str.isdigit, target)))
        epic2k2 = read_json_dic('EPIC2K2.json')
        missionCatalogID = epic2k2[inputCatalogID]
        ICnum = int(inputCatalogID[4:])
        self.mission = 'K2'

    elif target[:2].upper() == 'K2':
        missionCatalogID = 'K2-' + str(''.join(filter(str.isdigit, target)))
        k22epic = read_json_dic('K22EPIC.json')
        inputCatalogID = k22epic[missionCatalogID]
        ICnum = int(inputCatalogID[4:])
        self.mission = 'K2'

    elif target.isdigit():
        if mission == 'TESS':
            inputCatalogID = 'TIC' + target
            tic2toi = read_json_dic('TIC2TOI.json')
            missionCatalogID = tic2toi[inputCatalogID]
            ICnum = int(inputCatalogID[3:])
            print('As no prefix was passed, target was assumed to be \
                TIC {}'.format(ICnum))

        elif mission == 'Kepler':
            inputCatalogID = 'KIC' + target
            kic2kepler = read_json_dic('KIC2Kepler.json')
            missionCatalogID = kic2kepler[inputCatalogID]
            ICnum = int(inputCatalogID[3:])
            print('As no prefix was passed, target was assumed to be \
                KIC {}'.format(ICnum))

        elif mission == 'K2':
            inputCatalogID = 'EPIC' + target
            epic2k2 = read_json_dic('EPIC2K2.json')
            missionCatalogID = epic2k2[inputCatalogID]
            ICnum = int(inputCatalogID[4:])
            print('As no prefix was passed, target was assumed to be \
                EPIC {}'.format(ICnum))

    return inputCatalogID, missionCatalogID, ICnum

    def mask_outliers(self, timeseries=None, sigma=3):
        """Creates a mask to remove outliers from the lightcurve
        with special care to avoid removing transits.
        Parameters
        ----------
        timeseries : TimeSeries or Table (optional)
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
            flux = self.ts_stitch['flux']

        else:
            flux = timeseries[self.flux_kw]

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
        mask=None,
    ):
        """Removes the low frequency trend using scipy's Savitzky-Golay filter.
        This method wraps `scipy.signal.savgol_filter`.
        Parameters
        ----------
        sigma : int
            Number of sigma above which to remove outliers from the flatten
        timeseries : TimeSeries (optional)
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
        flatten_lc : `LightCurve`
            New light curve object with long-term trends removed.
        If ``return_trend`` is set to ``True``, this method will also return:
        trend_lc : `LightCurve`
            New light curve object containing the trend that was removed.
        """
        if timeseries is None:
            flux = self.ts_stitch['flux']
            lc = lk.LightCurve(timeseries)
            self.ts_clean = timeseries.copy()

        else:
            flux = timeseries[self.flux_kw]
            lc = lk.LightCurve(self.ts_stitch)
            self.ts_clean = self.ts_stitch.copy()

        robustmean1, robustmedian1, robustrms1 = sigma_clipped_stats(
            flux, sigma=sigma, maxiters=10)

        # prefiltering
        # lc = lk.LightCurve(time=self.ts_stitch.time.value,
        #                    flux=self.ts_stitch[self.flux_kw + '_norm'].value,
        #                    flux_err=self.ts_stitch[self.flux_err_kw].value)

        clc = lc.flatten(window_length=wl,
                         polyorder=polyorder,
                         break_tolerance=5,
                         sigma=sigma,
                         mask=masking)

        flux_filtered1 = clc.flux

        # sigma clipping for big transits
        clip = sigma_clip(flux_filtered1, sigma)

        if return_trend is True:
            finflat = lc.flatten(window_length=wl, polyorder=polyorder,
                                 break_tolerance=5, sigma=sigma,
                                 mask=clip.mask)
            clc1 = finflat[0]
            trend_ts = TimeSeries(finflat[1])

        flux_filtered = clc1.flux

        robustmean, robustmedian, robustrms = sigma_clipped_stats(
            flux_filtered, sigma=sigma, maxiters=10)
        print(robustrms)
        if robustrms > robustrms1:
            warnings.warn('error in the SGfit, flattening should not be used,\
             try using a mask with the transits', FulmarWarning)

        self.ts_clean[self.flux_kw] = flux_filtered

        if remove_outliers is True:
            out_mask = self.mask_outliers(sigma=sigma)
            self.ts_clean = self.ts_clean[out_mask]

            if return_trend is True:
                trend_ts = trend[out_mask]

                return self.ts_clean, trend_ts
            else:
                return self.ts_clean

        else:

            if return_trend is True:
                trend_ts = trend[out_mask]

                return self.ts_clean, trend_ts

            return self.ts_clean


def params_optimizer(
        timeseries=None,
        period_guess,
        t0_guess,
        depth_guess,
        ab,
        r_star,
        target_id,
        tran_window=0.25,
        tune=2500,
        draws=2500,
        chains=2,
        target_accept=0.95,
        ncores=None,
        mask=None):
    if ncores is None:
        ncores = multiprocessing.cpu_count()
    print('running on {} cores'.format(ncores))
#     x = ts_stitch.time.value
#     y = ts_stitch[flux_kw + '_clean'].value
#     yerr = ts_stitch[flux_err_kw+'_clean'].value
    if r_star is None:
        r_star = catalog_info(TIC_ID=target_identifier(target_id)[0])[4]

    x = time.copy()
    y = flux.copy()
    yerr = flux_err.copy()

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
                tune=tune,
                draws=draws,
                start=map_soln,
                chains=chains,
                cores=ncores,
                target_accept=target_accept,
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
            x_fold[inds], 1 + pred[0], 1 + pred[2], color="xkcd:green", alpha=0.1, zorder=1000
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
