API documentation
=================

Target object
-------------
.. py:class:: fulmar.target(targname, mission=None)

   :param targname: Name of the target as a string, e.g. "TOI-175" or, if mission is passed
            as the numerical identifier of the input catalog.
   :type targname: str or int
   :param mission: 'Kepler', 'K2', or 'TESS'
   :type mission: str, optional

   When the object is created, some parameters are automatically retrieved:

   - **ab** *(tuple of floats)* – Quadratic limb darkening parameters a, b.
   - **M_star** *(float)* Stellar mass (in units of solar masses)
   - **M_star_min** *(float)* – 1-sigma lower confidence interval on stellar mass (in units of solar mass)
   - **M_star_max** *(float)* – 1-sigma upper confidence interval on stellar mass (in units of solar mass)
   - **R_star** *(float)* Stellar radius (in units of solar radii)
   - **R_star_min** *(float)* – 1-sigma upper confidence interval on stellar radius (in units of solar radii)
   - **R_star_max** *(float)* – 1-sigma lower confidence interval on stellar radius (in units of solar radii)
   - **flux_kw** *(str)* – Keyword for the column containing the flux values
        (Default: 'flux')
   - **flux_err_kw** *(str)* – Keyword for the column containing the flux uncertainty values
        (Default: 'flux_err')
   - **lc_folder** *(str)* – Path to the folder where data is produced.
   - For Kepler targets:
      - **KIC** *(str)* – Identifier in the format of the Kepler input
           catalog, e.g. "KIC11904151"
      - **kep** *(str)* – Identifier in the format of the Kepler mission
           catalog, e.g. "Kepler-10"
      - **KIC_num** *(int)* – Number of the input catalog, e.g. "11904151"
   - For K2 targets:
      - **EPIC** *(str)* – Identifier in the format of the Ecliptic Plane input
           catalog, e.g. "EPIC201437844"
      - **K2** *(str)* – Identifier in the format of the K2 mission
           catalog, e.g. "K2-109"
      - **EPIC_num** *(int)* – Number of the input catalog, e.g. "201437844"
   - For TESS targets:
      - **TIC** *(str)* – Identifier in the format of the TESS input catalog,
           e.g. "TIC307210830"
      - **TOI** *(str)* – Identifier in the format of the TESS mission catalog,
           e.g. "TOI-175"
      - **TIC_num** *(int)* – Number of the input catalog,
           e.g. "307210830"  


   .. automethod:: set_lc_folder
   .. automethod:: set_flux_kw
   .. automethod:: set_flux_err_kw
   .. automethod:: search_data
   .. automethod:: build_lightcurve

   .. automethod:: mask_outliers
   .. automethod:: clean_subt_activity_flatten
   .. automethod:: clean_subt_activity_GP

   .. automethod:: plot_transitcheck

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
.. autofunction:: fulmar.perioplot
.. autofunction:: fulmar.modelplot


Time
----
.. automodule:: fulmar.time
   :members:

Utils
-----
.. automodule:: fulmar.utils
   :members: