from astopy.table import Table
from astropy.time import Time
from astropy.timeseries import TimeSeries
import astropy.units as u
import UnitConversionError
from fulmar.func import (
    mission_identifier,
    target_identifier,
    read_lc_from_file,
    normalize_lc,
    time_flux_err,
    ts_binner,
    fbn,
    estimate_planet_mass,
    estimate_semi_amplitude
)
from fulmar.time import TimeRJD
from fulmar.utils import FulmarWarning
import lightkurve as lk
import numpy as np
import numpy.testing as npt

import pytest


def test_mission_identifier():
    """Test the mission_identifier function""
    """
    # Start by checking if ValueError is raised for a bad target
    with pytest.raises(ValueError):
        mission_identifier('Hello')

    # Check for the Kepler mission
    assert mission_identifier('Kepler-10') == 'Kepler'
    assert mission_identifier('KIC11904151') == 'Kepler'

    # Check for K2
    assert mission_identifier('K2-109') == 'K2'
    assert mission_identifier('EPIC201437844') == 'K2'

    # Check for TESS
    assert mission_identifier('TOI-175') == 'TESS'
    assert mission_identifier('TIC307210830') == 'TESS'


def test_target_identifier():
    """Test the target_identifier function
    """
    # Start by checking if ValueError is raised for a bad target
    with pytest.raises(ValueError):
        target_identifier('World')

    # target is int
    # Check if the correct message is passed for Kepler
    with pytest.raises(ValueError, match='range 1 to 13161029'):
        target_identifier(0, mission='Kepler')
        target_identifier(123456789, mission='Kepler')

    # Check if the correct message is passed for K2
    with pytest.raises(ValueError, match='range 201000001 to 251813738'):
        target_identifier(0, mission='K2')
        target_identifier(260000000, mission='Kepler')

    # Check if the correct message is passed for no mission
    with pytest.raises(ValueError, match='supported mission'):
        target_identifier(42, mission='H2G2')

    # Check for the No prefix warning
    with pytest.warns(FulmarWarning, match='target is assumed'):
        target_identifier(11904151, mission='Kepler')
        target_identifier(201537844, mission='K2')
        target_identifier(307210830, mission='TESS')


def test_read_lc_from_file():
    """Test the read_lc_from_file function
    """
    # First, the simple case of a 3 column file with no headers
    lc = read_lc_from_file('simple_table.txt')

    assert lc.colnames == ['time', 'flux', 'flux_err']
    assert lc.time.format == 'jd'
    npt.assert_array_equal(lc.time.value, np.array([1, 2, 3]))
    npt.assert_array_equal(lc.flux.value, np.array([1, 1, 1]))
    npt.assert_array_equal(lc.flux_err.value, np.array([0.1, 0.1, 0.1]))

    # Test if the exptime column behaves as expected
    npt.assert_array_equal(lc['exptime'], np.array([1, 1, 1]))

    # Now let's try naming the columns
    # Using expected names:
    lc = read_lc_from_file('simple_table.txt',
                           colnames=['time', 'flux', 'flux_err'])
    npt.assert_array_equal(lc.time.value, np.array([1, 2, 3]))
    npt.assert_array_equal(lc.flux.value, np.array([1, 1, 1]))
    npt.assert_array_equal(lc.flux_err.value, np.array([0.1, 0.1, 0.1]))

    # Not using 'time'
    with pytest.raises(ValueError, match='time'):
        read_lc_from_file('simple_table.txt',
                          colnames=['t', 'flux', 'flux_err'])
        read_lc_from_file('simple_table_wrong_t.txt')

    # Not using 'flux'
    with pytest.raises(ValueError, match='flux'):
        read_lc_from_file('simple_table.txt',
                          colnames=['time', 'y', 'flux_err'])
        read_lc_from_file('simple_table_wrong_flux.txt')

    # Using 'time' anf 'flux' but not 'flux_err'
    # Correct number of columns
    with pytest.warns(FulmarWarning):
        lc = read_lc_from_file('simple_table.txt',
                               colnames=['time', 'flux', 'err'])
        # Check if a 'flux_err' column was created with Nans
        assert len(lc.colnames) == 5
        npt.assert_array_equal(
            np.array([True, True, True]), np.isnan(lc.flux_err))
    # Wrong number of columns
    with pytest.warns(FulmarWarning, match='should match'):
        read_lc_from_file('simple_table.txt',
                          colnames=['time', 'flux', 'flux_err', '4th_column'])

    # Then, the same file but with author, exptime and timeformat
    lc = read_lc_from_file(
        'simple_table.txt', author='TELESCOPE', exptime=120, timeformat='mjd')

    assert lc.author == 'TELESCOPE'
    assert lc.meta['EXPTIME'] == 120
    assert lc.time.format == 'mjd'
    npt.assert_array_equal(lc['exptime'], np.array([120, 120, 120]))

    # Correct headers in the file
    lc = read_lc_from_file('simple_table_correct_cols.txt')
    assert lc.colnames == ['time', 'flux', 'flux_err']
    assert lc.time.format == 'jd'
    npt.assert_array_equal(lc.time.value, np.array([1, 2, 3]))
    npt.assert_array_equal(lc.flux.value, np.array([1, 1, 1]))
    npt.assert_array_equal(lc.flux_err.value, np.array([0.1, 0.1, 0.1]))


def test_normalize_lc():
    """Test if normalize_lc actually normalizes the lightcurve"""
    lc = lk.lightkurve(time=np.arange(5), flux=3 * np.ones(5),
                       flux_err=0.1 * np.ones(5))

    npt.assert_allclose(np.median(normalize_lc(lc).flux), 1)
    npt.assert_allclose(np.median(normalize_lc(lc).flux_err), 0.1 / 4)


def test_time_flux_err():
    """Test if time_flux_err returns the expected values"""
    dirty_array = np.ones(10, dtype=object)
    time_array = np.linspace(1, 10, 10)
    dy_array = np.ones(10, dtype=object)
    dirty_array[1] = None
    dirty_array[2] = np.inf
    dirty_array[3] = -np.inf
    dirty_array[4] = np.nan
    dirty_array[5] = -99
    time_array[8] = np.nan
    dy_array[9] = np.inf
    ts = TimeSeries([Time(time_array, format='jd'), dirty_array, dy_array],
                    names=['time', 'flux', 'flux_err'])

    t, y, dy = time_flux_err(ts)
    npt.assert_equal(t, [1, 7, 8])
    npt.assert_equal(y, [1, 1, 1])
    npt.assert_equal(dy, [1, 1, 1])

    # Now with flux_err as only nan
    dy_array = time_array * np.nan
    ts = TimeSeries([Time(time_array, format='jd'), dirty_array, dy_array],
                    names=['time', 'flux', 'flux_err'])
    # replace_nan_err is True
    t, y, dy = time_flux_err(ts)
    npt.assert_equal(dy, [0, 0, 0])
    # replace_nan_err is False
    t, y, dy = time_flux_err(ts, replace_nan_err=False)
    npt.assert_array_equal(np.array([True, True, True]), np.isnan(dy))


def test_ts_binner():
    """test if ts_binner behaves as expected"""

    testtime = np.linspace(1, 10, 100)
    testflux = np.linspace(0.99, 1.01, 100)
    ts = TimeSeries([Time(testtime * u.d, format='jd'),
                     testflux], names=['time', 'flux'])
    ts_bin = ts_binner(ts, bin_duration=1 * u.d)

    expected_t = np.array([1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
    expected_flux = np.array([0.99111111, 0.99343434, 0.99565657, 0.99787879,
                              1.00010101, 1.00232323, 1.00454545, 1.00676768,
                              1.00888889])
    npt.assert_almost_equal(ts_bin.time.value, expected_t)
    npt.assert_almost_equal(ts_bin['flux'], expected_flux)
    npt.assert_equal(len(ts_bin), 9)


def test_fbn():
    """test if fbn behaves as expected"""

    testtime = np.linspace(1, 10, 1000)

    testflux = np.sin(np.pi * testtime)

    ts = TimeSeries([Time(testtime * u.d, format='jd'),
                     testflux], names=['time', 'flux'])
    ts_fold, ts_fold_bin = fbn(ts, 2., 5., 2 * u.d)

    npt.assert_equal(len(ts_fold), 1000)
    npt.assert_equal(len(ts_fold_bin), 20)
    expected_flux_bin = np.array([-0.95, -0.85, -0.75, -0.65, -0.55, -0.45,
                                  -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25,
                                  0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95])

    npt.assert_almost_equal(ts_fold['phase_norm'], ts_fold.time.value / 2)
    npt.assert_almost_equal(ts_fold_bin['flux'], expected_flux_bin)
    npt.assert_almost_equal(ts_fold['phase_norm'].max(), 0.4954955)
    npt.assert_almost_equal(ts_fold['phase_norm'].min(), -.5)


def test_estimate_planet_mass():
    """test if estimate_planet_mass behaves as expected"""
    npt.assert_equal(estimate_planet_mass(
        1, 'Earth'), 1 * u.earthMass)
    npt.assert_equal(estimate_planet_mass(
        1, 'Neptune'), 0.29706202 * u.earthMass)
    npt.assert_equal(estimate_planet_mass(
        1, 5514), 1 * u.earthMass)

    with pytest.raises(TypeError, match='`astropy.units.Quantity` or float'):
        estimate_planet_mass('string', 'Earth')

    with pytest.raises(ValueError, match="Accepted str values for rho_p"):
        estimate_planet_mass(1, 'Uranus')

    with pytest.raises(UnitConversionError):
        estimate_planet_mass(1 * u.s, 'neptune')


def test_estimate_semi_amplitude():
    """test if estime_semi_amplitude behaves as exected"""
    npt.assert_equal(estimate_semi_amplitude(
        365, 1, 1), 0.08948015 * u.m / u.s)
    npt.assert_equal(estimate_semi_amplitude(
        365, 1, R_planet=1, rho_planet='earth'), 0.08948015 * u.m / u.s)
    npt.assert_equal(estimate_semi_amplitude(365, 1, 1, inc=0).value, 0)
    npt.assert_equal(estimate_semi_amplitude(
        365, 1, 1, ecc=0.5), 0.12654404 * u.m / u.s)

    with pytest.raises(TypeError, match='`astropy.units.Quantity` or float'):
        estimate_planet_mass('1 year', 1, 1)
        estimate_planet_mass(365, '1 solMass', 1)
        estimate_planet_mass(365, 1, '1 earthMass')
        estimate_planet_mass(365, 1, 1, inc='90 deg')

    with pytest.raises(ValueError, match='required when M_planet is not'):
        estimate_planet_mass(365, 1)
        estimate_planet_mass(365, 1, R_planet=1)
        estimate_planet_mass(365, 1, rho_planet=1)

    with pytest.warns(FulmarWarning, match='overrides'):
        estimate_planet_mass(365, 1, 1, R_planet=1, rho_planet='earth')
