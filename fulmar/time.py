"""
Adds the RJD time format for use by Astropy's `Time` object.
Caution: AstroPy time objects make a distinction between a time's format
(e.g. ISO, JD, MJD) and its scale (e.g. UTC, TDB).
Note: the classes below derive from an AstroPy meta class which
will automatically register the formats for use in AstroPy Time objects.
"""
from astropy.time.formats import TimeFromEpoch


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
