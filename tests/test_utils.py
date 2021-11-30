from astropy.time import Time
from fulmar.time import TimeRJD
from fulmar.utils import (
    rjd_to_astropy_time,
    FulmarError,
    FulmarWarning)

import pytest
import warnings


def test_rjd_to_astropy_time():
    """Tests for the Reduced Julian Date (RJD) time format."""
    # Sanity checks
    t0 = rjd_to_astropy_time(0)
    assert t0.format == "rjd"
    assert t0.scale == "tdb"
    assert t0.iso == "1858-11-16 12:00:00.000"


def test_FulmarError():
    """Tests if FulmarError can be raised."""
    with pytest.raises(FulmarError):
        raise FulmarError


def test_FulmarWarning():
    """Tests if FulmarWarning can be warned"""
    with pytest.warns(FulmarWarning):
        warnings.warn('FulmarWarning was raised', FulmarWarning)
