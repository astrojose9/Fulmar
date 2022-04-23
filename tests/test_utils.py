import os
import sys
#sys.path.insert(0, os.path.abspath('/home/jrodrigues/Documents/PhD/fulmar'))

from fulmar.utils import (
    FulmarError,
    FulmarWarning)

import pytest
import warnings


def test_FulmarError():
    """Tests if FulmarError can be raised."""
    with pytest.raises(FulmarError):
        raise FulmarError


def test_FulmarWarning():
    """Tests if FulmarWarning can be warned"""
    with pytest.warns(FulmarWarning):
        warnings.warn('FulmarWarning was raised', FulmarWarning)
