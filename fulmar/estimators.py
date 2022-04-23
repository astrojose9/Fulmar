import astropy.units as u
import numpy as np

import warnings
from fulmar.utils import (
    FulmarWarning
)


def estimate_planet_mass(
        R_p,
        rho_p):
    """
    Estimates the mass of an exoplanet from its radius and density.

    Parameters
    ----------
    R_p : `~astropy.units.Quantity` or float
        Radius of the exolanet. (defaults to units of Earth radii)
    rho_p : `~astropy.units.Quantity`, float or str
        Density of the exoplanet in kg * m^-3. Can be "Earth" or "Neptune".

    Returns
    -------
    M_planet : `~astropy.units.Quantity`
        Estimated mass of the exoplanet.
    """
    dens_dic = {'earth': 5514 * (u.kg / u.m**3),
                'neptune': 1638 * (u.kg / u.m**3)}

    if isinstance(R_p, (int, float)):
        R_p = R_p * u.earthRad

    elif isinstance(R_p, u.Quantity):
        R_p = R_p.to(u.earthRad)

    else:
        raise TypeError('R_p should be `astropy.units.Quantity` or float')

    if isinstance(rho_p, (int, float)):
        rho_p = rho_p * (u.kg / u.m**3)

    elif isinstance(rho_p, str):
        if rho_p.lower() in dens_dic.keys():
            rho_p = dens_dic[rho_p.lower()]
        else:
            raise ValueError(
                'Accepted str values for rho_p are "Earth" and "Neptune".')
    else:
        raise TypeError(
            'rho_p should be `astropy.units.Quantity`, float or str.')

    M_planet = (R_p.value ** 3 * rho_p / dens_dic['earth']) * u.earthMass
    return M_planet


def estimate_semi_amplitude(
        period,
        M_star,
        M_planet=None,
        R_planet=None,
        rho_planet=None,
        inc=90 * u.deg,
        ecc=0):
    """
    Estimates the radial velocity semi-amplitude corresponding to a planet of
    given parameters.

    Parameters
    ----------
    period : `~astropy.units.Quantity` or float
        The period to use for folding. (defaults to units of days)
    M_star : `~astropy.units.Quantity` or float
        Stellar mass (defaults to units of solar masses)
    M_planet : `~astropy.units.Quantity` or float
        Mass of the exolanet. (defaults to units of Earth masses)
    R_planet : `~astropy.units.Quantity` or float
        Radius of the exolanet. (defaults to units of Earth radii)
    rho_planet : `~astropy.units.Quantity`, float or str
        Density of the exoplanet in kg * m^-3. Can be "Earth" or "Neptune".
    inc : `~astropy.units.Quantity` or float
        Orbital inclination (in degrees). Default: 90.
    ecc : float
        Orbital eccentricity. Default: 0.

    Returns
    -------
    K : `~astropy.units.Quantity`
        Estimated semi-amplitude of the RV.
    """
    if isinstance(period, (int, float)):
        period = period * u.d
    elif not isinstance(period, u.Quantity):
        raise TypeError(
            'period shoud be `astropy.units.Quantity` or float')

    if isinstance(M_star, (int, float)):
        M_star = M_star * u.solMass
    elif not isinstance(M_star, u.Quantity):
        raise TypeError(
            'M_star shoud be `astropy.units.Quantity` or float')

    if isinstance(inc, (int, float)):
        inc = inc * u.deg
    elif not isinstance(inc, u.Quantity):
        raise TypeError(
            'inc shoud be `astropy.units.Quantity` or float')

    if M_planet is None:

        if R_planet is None or rho_planet is None:
            raise ValueError('R_planet and rho_planet are both required '
                             'when M_planet is not given')
        else:
            M_planet = estimate_planet_mass(R_planet, rho_planet)
    else:
        if R_planet is not None or rho_planet is not None:
            warnings.warn(
                'M_planet overrides R_planet and rho_planet', FulmarWarning)

        if isinstance(M_planet, (int, float)):
            M_planet = M_planet * u.earthMass

        elif not isinstance(M_planet, u.Quantity):
            raise TypeError(
                'M_planet should be `astropy.units.Quantity` or float')

    inc = inc.to(u.deg)

    M_p_jovian = M_planet.to(u.jupiterMass).value
    M_tot_solar = (M_star + M_planet).to(u.solMass).value

    # Equation (14) from Lovis & Fischer 2010
    K = 28.4329 * (u.m / u.s) * \
        M_p_jovian * \
        np.sin(inc) * np.power(M_tot_solar, -2 / 3) * \
        np.power(period.to(u.year).value, -1 / 3) / np.sqrt(1 - ecc)

    return K
