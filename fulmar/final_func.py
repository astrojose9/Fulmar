# Fulmar functions which are ready to ship


# main.py
import os
from os import path
import numpy as np

from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.timeseries import aggregate_downsample
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

#############################################################

    def target_identifier(self, target, mission):
        """Translate the target identifiers between different catalogs
        such as TIC to TOI in the case of TESS or EPIC to K" for K2
        Updates the mission parameter in case it wasn't passed by the user.
        """
        if target[:3].upper() == 'TIC':
            inputCatalogID = 'TIC' + str(''.join(filter(str.isdigit, target)))
            tic2toi = read_json_dic(
                os.path.join(fulmar_constants.fulmar_dir, 'TIC2TOI.json'))
            missionCatalogID = tic2toi[inputCatalogID]
            ICnum = int(inputCatalogID[3:])
            self.mission = 'TESS'

        elif target[:3].upper() == 'TOI':
            missionCatalogID = 'TOI-' + \
                str(''.join(filter(str.isdigit, target)))
            toi2tic = read_json_dic(
                os.path.join(fulmar_constants.fulmar_dir, 'TOI2TIC.json'))
            inputCatalogID = toi2tic[missionCatalogID]
            ICnum = int(inputCatalogID[3:])
            self.mission = 'TESS'

        elif target[:3].upper() == 'KIC':
            inputCatalogID = 'KIC' + str(''.join(filter(str.isdigit, target)))
            kic2kepler = read_json_dic(
                os.path.join(fulmar_constants.fulmar_dir, 'KIC2Kepler.json'))
            missionCatalogID = kic2kepler[inputCatalogID]
            ICnum = int(inputCatalogID[3:])
            self.mission = 'Kepler'

        elif target[:3].upper() == 'KEP':
            missionCatalogID = 'Kepler-' + \
                str(''.join(filter(str.isdigit, target)))
            kep2kic = read_json_dic(
                os.path.join(fulmar_constants.fulmar_dir, 'Kepler2KIC.json'))
            inputCatalogID = kep2kic[missionCatalogID]
            ICnum = int(inputCatalogID[3:])
            self.mission = 'Kepler'

        elif target[:4].upper() == 'EPIC':
            inputCatalogID = 'EPIC' + str(''.join(filter(str.isdigit, target)))
            epic2k2 = read_json_dic(
                os.path.join(fulmar_constants.fulmar_dir, 'EPIC2K2.json'))
            missionCatalogID = epic2k2[inputCatalogID]
            ICnum = int(inputCatalogID[4:])
            self.mission = 'K2'

        elif target[:2].upper() == 'K2':
            missionCatalogID = 'K2-' + \
                str(''.join(filter(str.isdigit, target[2:])))
            k22epic = read_json_dic(
                os.path.join(fulmar_constants.fulmar_dir, 'K22EPIC.json'))
            inputCatalogID = k22epic[missionCatalogID]
            ICnum = int(inputCatalogID[4:])
            self.mission = 'K2'

        elif target.isdigit():
            if mission == 'TESS':
                inputCatalogID = 'TIC' + target
                tic2toi = read_json_dic(
                    os.path.join(fulmar_constants.fulmar_dir, 'TIC2TOI.json'))
                missionCatalogID = tic2toi[inputCatalogID]
                ICnum = int(inputCatalogID[3:])
                warnings.warn('As no prefix was passed, target was assumed to be \
                    TIC {}'.format(ICnum), FulmarWarning)

            elif mission == 'Kepler':
                inputCatalogID = 'KIC' + target
                kic2kepler = read_json_dic(
                    os.path.join(fulmar_constants.fulmar_dir, 'KIC2Kepler.json'))
                missionCatalogID = kic2kepler[inputCatalogID]
                ICnum = int(inputCatalogID[3:])
                warnings.warn('As no prefix was passed, target was assumed to be \
                    KIC {}'.format(ICnum), FulmarWarning)

            elif mission == 'K2':
                inputCatalogID = 'EPIC' + target
                epic2k2 = read_json_dic(
                    os.path.join(fulmar_constants.fulmar_dir, 'EPIC2K2.json'))
                missionCatalogID = epic2k2[inputCatalogID]
                ICnum = int(inputCatalogID[4:])
                warnings.warn('As no prefix was passed, target was assumed to be \
                    EPIC {}'.format(ICnum), FulmarWarning)

        return inputCatalogID, missionCatalogID, ICnum

# func.py

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

import warnings


def read_lc_from_file(file, author=None, exptime=None, colnames=None):
    """Creates a LightCurve from a file.
    Parameters
    ----------
    file : str
        path to the file containing the light curve data.
    author : str (optional)
        Name of the pipeline used to reduce the data.
    exptime : float
        Exposure time of the observation, in seconds.
    colnames : list (of str)
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
            t = Table.read(file)
        except IORegistryError:  # Helps the astropy reader
            t = Table.read(file, format='ascii', comment='#')
        if colnames is not None:
            try:
                t = Table(t, names=colnames)
            except ValueError:
                warnings.warn('number of items in colnames should match \
                    the number of columns in the data', FulmarWarning)
        elif t.colnames[0] == 'col1':
            t.rename_column('col1', 'time')
            t.rename_column('col2', 'flux')
            if len(t.colnames) > 2:
                t.rename_column('col3', 'flux_err')
        lc = lk.LightCurve(t)

    if author is not None:
        try lc.meta['AUTHOR']:
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
