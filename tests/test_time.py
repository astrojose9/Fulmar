from astropy.time import Time
from fulmar.time import TimeRJD
import pytest


def test_rjd():
    """Tests for the Reduced Julian Date (RJD) time format."""
    # Sanity checks
    t0 = Time(0, format="rjd")
    assert t0.format == "rjd"
    assert t0.scale == "tdb"
    assert t0.iso == "1858-11-16 12:00:00.000"
