import os
import sys
#sys.path.insert(0, os.path.abspath('/home/jrodrigues/Documents/PhD/fulmar'))

from astropy.time import Time
from fulmar.time import (
    TimeRJD,
    rjd_to_astropy_time
)
import pytest


def test_rjd():
    """Tests for the Reduced Julian Date (RJD) time format."""
    # Sanity checks
    t0 = Time(0, format="rjd")
    assert t0.format == "rjd"
    assert t0.scale == "tdb"
    assert t0.iso == "1858-11-16 12:00:00.000"


def test_rjd_to_astropy_time():
    """Tests for the Reduced Julian Date (RJD) time format."""
    # Sanity checks
    t0 = rjd_to_astropy_time(0)
    assert t0.format == "rjd"
    assert t0.scale == "tdb"
    assert t0.iso == "1858-11-16 12:00:00.000"
