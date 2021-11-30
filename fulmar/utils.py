from astropy.time import Time
from fulmar.fulmar_constants import FULMAR_VERSION_STR
import numpy as np

##############################################################################


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


class FulmarError(Exception):
    """To raise exceptions related to FULMAR
    """
    pass


class FulmarWarning(Warning):
    """Class from warning to be displayed as
    "FulmarWarning"
    """

    pass


def warning_on_one_line(message, category, filename, lineno,
                        file=None, line=None):
    """Function to display warnings on one line, as to avoid displaying
    'warnings.warn('warning message')' under the 'warning message',
    for clarity."""
    return ' %s:%s: %s: %s' % (filename, lineno, category.__name__, message)


# warnings.formatwarning = warning_on_one_line
def print_version():
    """prints the version of fulmar used."""
    print(FULMAR_VERSION_STR)
    return
