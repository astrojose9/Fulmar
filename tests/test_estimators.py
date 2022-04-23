import os
import sys
# sys.path.insert(0, os.path.abspath('/home/jrodrigues/Documents/PhD/fulmar'))

from fulmar.estimators import (
    estimate_planet_mass,
    estimate_semi_amplitude
)
from fulmar.utils import (
    FulmarWarning
)

import astropy.units as u
import numpy as np
import numpy.testing as npt
from astropy.units import UnitConversionError

import pytest


def test_estimate_planet_mass():
    """test if estimate_planet_mass behaves as expected"""
    npt.assert_equal(estimate_planet_mass(
        1, 'Earth').value, 1)  # * u.earthMass)
    npt.assert_almost_equal(estimate_planet_mass(
        1, 'Neptune').value, 0.29706202)  # * u.earthMass)
    npt.assert_almost_equal(estimate_planet_mass(
        1, 5514).value, 1)  # * u.earthMass)

    with pytest.raises(TypeError, match='`astropy.units.Quantity` or float'):
        estimate_planet_mass('string', 'Earth')

    with pytest.raises(ValueError, match="Accepted str values for rho_p"):
        estimate_planet_mass(1, 'Uranus')

    with pytest.raises(UnitConversionError):
        estimate_planet_mass(1 * u.s, 'neptune')


def test_estimate_semi_amplitude():
    """test if estime_semi_amplitude behaves as exected"""
    npt.assert_almost_equal(estimate_semi_amplitude(
        365, 1, 1).value, 0.08948015)  # * u.m / u.s)
    npt.assert_almost_equal(estimate_semi_amplitude(
        365, 1, R_planet=1, rho_planet='earth').value, 0.08948015)  # * u.m / u.s)
    npt.assert_equal(estimate_semi_amplitude(
        365, 1, 1, inc=0).value, 0)
    npt.assert_almost_equal(estimate_semi_amplitude(
        365, 1, 1, ecc=0.5).value, 0.12654404)  # * u.m / u.s)

    with pytest.raises(TypeError, match='`astropy.units.Quantity` or float'):
        estimate_planet_mass('1 earthRad', 'earth')
        estimate_semi_amplitude('1 year', 1, 1)
        estimate_semi_amplitude(365, '1 solMass', 1)
        estimate_semi_amplitude(365, 1, M_planet='1 earthMass')
        estimate_semi_amplitude(365, 1, M_planet=1 * u.earthMass, inc='90 deg')

    with pytest.raises(ValueError, match='required when M_planet is not'):
        estimate_semi_amplitude(365, 1)
        estimate_semi_amplitude(365, 1, R_planet=1)
        estimate_semi_amplitude(365, 1, rho_planet=1)

    with pytest.warns(FulmarWarning, match='overrides'):
        estimate_semi_amplitude(365, 1, 1, R_planet=1, rho_planet='earth')
