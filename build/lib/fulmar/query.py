import numpy as np
from astropy.table import Table


def teff_logg_TIC(TIC_ID):
    """Takes TIC_ID, returns Teff, logg, from online catalog using Vizier

    Parameters
    ----------
    TIC_ID : int
        For TESS targets. Number of the input catalog, e.g. "307210830"

    Returns
    -------
    Teff : float
        Effective temperature
    Teff_err : float
        Error on `Teff`
    logg : float
        Spectroscopic surface gravity
    logg_err : float
        Error on `logg`

    Raises
    ------
    ImportError
        If astroquery package failed to import
    """
    if type(TIC_ID) is not int:
        raise TypeError('TIC_ID ID must be of type "int"')
    try:
        from astroquery.mast import Catalogs
    except ModuleNotFoundError:
        raise ImportError("Package astroquery required but failed to import")
    Teff, Teff_err, logg, logg_err = Catalogs.query_criteria(
        catalog="Tic", ID=TIC_ID)[
        'Teff', 'e_Teff', 'logg', 'e_logg'].as_array()[0]

    return Teff, Teff_err, logg, logg_err


def teff_logg_KIC(KIC_ID):
    """Takes KIC_ID, returns Teff, logg, from online catalog using Vizier

    Parameters
    ----------
    KIC_ID : int
        For Kepler targets. Number of the input catalog, e.g. "11904151"

    Returns
    -------
    Teff : float
        Effective temperature
    Teff_err : float
        Error on `Teff`
    logg : float
        Spectroscopic surface gravity
    logg_err : float
        Error on `logg`

    Raises
    ------
    ImportError
        If astroquery package failed to import
    """
    if type(KIC_ID) is not int:
        raise TypeError('KIC_ID ID must be of type "int"')
    try:
        from astroquery.vizier import Vizier
    except ModuleNotFoundError:
        raise ImportError("Package astroquery required but failed to import")

    columns = ["Teff", 'e_Teff', 'log(g)', 'e_log(g)']
    catalog = "J/ApJS/229/30/catalog"
    Teff, Teff_err, logg, logg_err = Vizier(columns=columns).query_constraints(
        KIC=KIC_ID, catalog=catalog)[0].as_array()[0]

    return Teff, Teff_err, logg, logg_err


def teff_logg_EPIC(EPIC_ID):
    """Takes EPIC_ID, returns Teff, logg, from online catalog using Vizier

    Parameters
    ----------
    EPIC_ID : int
        For K2 targets. Number of the input catalog, e.g. "201437844"

    Returns
    -------
    Teff : float
        Effective temperature
    Teff_err : float
        Error on `Teff`
    logg : float
        Spectroscopic surface gravity
    logg_err : float
        Error on `logg`

    Raises
    ------
    ImportError
        If astroquery package failed to import
    """
    if type(EPIC_ID) is not int:
        raise TypeError('EPIC_ID ID must be of type "int"')
    try:
        from astroquery.vizier import Vizier
    except ModuleNotFoundError:
        raise ImportError("Package astroquery required but failed to import")

    columns = ["Teff", 'e_Teff', 'logg', 'e_logg']
    catalog = "IV/34/epic"
    Teff, Teff_err, logg, logg_err = Vizier(columns=columns).query_constraints(
        ID=EPIC_ID, catalog=catalog)[0].as_array()[0]

    return Teff, Teff_err, logg, logg_err
