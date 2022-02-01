"""
Adds the RJD time format for use by Astropy's `Time` object.
Caution: AstroPy time objects make a distinction between a time's format
(e.g. ISO, JD, MJD) and its scale (e.g. UTC, TDB).
"""
from astropy.time import Time
from astropy.time.formats import TimeFromEpoch
import numpy as np


class TimeRJD(TimeFromEpoch):
    """
    Reduced Julian Date (RJD) time format.
    This represents the number of days since noon on November 16, 1858.
    For example, 51544.5 in RJD is midnight on January 1, 2000.
    Typically used in DACE.
    RJD = JD − 2400000

    """
    name = 'rjd'
    unit = 1.0
    epoch_val = 2400000
    epoch_val2 = None
    epoch_scale = 'tdb'
    epoch_format = 'jd'


def rjd_to_astropy_time(rjd) -> Time:
    """Converts Reduced Julian Day (RJD) time values to an
    `~astropy.time.Time` object.
    Reduced Julian Day (RJD) is a Julian day minus 2400000.0
    (UTC=January 1, 2000 12:00:00)..
    The time is in the Barycentric Dynamical Time frame (TDB), which is a
    time system that is not affected by leap seconds.

    Parameters
    ----------
    rjd : float or array of floats
        Reduced Julian Day.

    Returns
    -------
    time : `astropy.time.Time` object
        Resulting time object.
    """
    rjd = np.atleast_1d(rjd)
    # Some data products have missing time values;
    # we need to set these to zero or `Time` cannot be instantiated.
    rjd[~np.isfinite(rjd)] = 0
    return Time(rjd, format="rjd", scale="tdb")
