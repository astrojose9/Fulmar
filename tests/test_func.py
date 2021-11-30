from astopy.table import Table
from astropy.time import Time
from astropy.timeseries import TimeSeries
from fulmar.func import (
    mission_identifier,
    target_identifier,
    read_lc_from_file,
    normalize_lc,
    time_flux_err
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
        assert len(lc.colnames) == 4
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
    assert lc.exptime == 120
    assert lc.time.format == 'mjd'

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
    t = np.range(5)
    flux = np.ones(5)
    flux_err = 0.05 * np.ones(5)
    tbl = Table([t, flux, flux_err])
    ts = TimeSeries(tbl, time=Time(t, format='jd'))
