FULMAR
======

Target
------
.. py:class:: fulmar.target(targname, mission=None)

   :param targname: Name of the target as a string, e.g. "TOI-175" or, if mission is passed
            as the numerical identifier of the input catalog.
   :type targname: str or int
   :param mission: 'Kepler', 'K2', or 'TESS'
   :type mission: str, optional

   When the object is created, some parameters are automatically created/updated:

   - :ab: *(tuple of floats)* Quadratic limb darkening parameters a, b.
   - :M_star: *(float)* Stellar mass (in units of solar masses)
   - :M_star_min: *(float)* 1-sigma lower confidence interval on stellar mass (in units of solar mass)
   - :M_star_max: *(float)* 1-sigma upper confidence interval on stellar mass (in units of solar mass)
   - :R_star: *(float)* Stellar radius (in units of solar radii)
   - :R_star_min: *(float)* 1-sigma upper confidence interval on stellar radius (in units of solar radii)
   - :R_star_max: *(float)* 1-sigma lower confidence interval on stellar radius (in units of solar radii)

   Once the target is defined, one can look for data or directly build the lightcurve.

   .. automethod:: search_data
   .. automethod:: build_lightcurve

   Then one can work on the data by cleaning it.

   .. automethod:: mask_outliers
   .. automethod:: clean_subt_activity_flatten
   .. automethod:: clean_subt_activity_GP

   Having a look at a known promising signal.

   .. automethod:: plot_transitcheck

   Or looking for signals.

   .. automethod:: tls_periodogram


Functions
---------
.. autofunction:: fulmar.mission_identifier
.. autofunction:: fulmar.target_identifier
.. autofunction:: fulmar.read_lc_from_file
.. autofunction:: fulmar.normalize_lc
.. autofunction:: fulmar.time_flux_err
.. autofunction:: fulmar.ts_binner
.. autofunction:: fulmar.fbn
.. autofunction:: fulmar.GP_fit
.. autofunction:: fulmar.params_optimizer
.. autofunction:: fulmar.estimate_planet_mass
.. autofunction:: fulmar.estimate_semi_amplitude


Time
----
.. automodule:: fulmar.time
   :members:

Utils
-----
.. automodule:: fulmar.utils
   :members: