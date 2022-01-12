---
title: "FULMAR: Follow-Up Lightcurves Multitool Assisting Radial velocities"
tags:
  - Python
  - astronomy
authors:
 - name: J. Rodrigues
   orcid: 0000-0001-5164-3602
   affiliation: "1, 2" # (Multiple affiliations must be quoted)
 - name: S. C. C. Barros
   orcid: 0000-0003-2434-3625
   affiliation: "1, 2"
 - name: N. C. Santos
   orcid: 0000-0003-4422-2919
   affiliation: "1, 2"


affiliations:
 - name: Instituto de Astrofísica e Ciências do Espaço, Universidade do Porto, CAUP, Rua das Estrelas, 4150-762, Porto, Portugal
   index: 1
 - name: Departamento de Física e Astronomia, Faculdade de Ciências, Universidade do Porto, Rua do Campo Alegre, 4169-007, Porto, Portugal
   index: 2

date: 27 December 2021
bibliography: paper.bib
---


# Summary
<!-- With the advent of all sky photometric surveys looking for exoplanets, a big number of interesting systems are detected in a short amount of time. Given the high instrumental time required for an efficient radial velocity follow-up, most of these targets cannot realistically be observed using ground facilities. An efficient follow-up strategy is thus very important. -->
There are two main ways of detecting exoplanets: _transit photometry_ and _radial velocities_ (RV). If a planet crosses in front of the star it orbits around, the observed brightness of the host star drops a bit. The depth of this drop (called a _transit_) is directly related to the relative sizes of the star and the planet and can therefore be used to determine the planet's _radius_.

An orbiting planet gravitationnaly attracts its host star, and both orbit around their common center of mass. This motion creates a periodic variation of the velocity of the star projected on the obesrver's line of sight. By measuring the doppler shifts in stellar spectra, one can infer the relative masses of the star and the planet and therefore determine the planet's _mass_.

Combining both measurements allow to estimate the _bulk density_ of the exoplanet, and thus probe its internal composition. Such results are necessary for astronomers to bring new understanding to formation and evolution models of exoplanets.

**FULMAR** is an open source Python package that was created to assist astronomers involved in radial velocity follow-up programs of transiting exoplanets by making the analysis of the light curves easier.
It provides tools to download lightcurves, correct stellar activity, to look for transits, to refine transit parameters, to estimate the amplitude of the corresponding RV signal, and to visually probe signals detected in RV.
It was build in a modular way, making new features easier to implement. 
**FULMAR** is currently used in the working group 3 ("RV follow-up of K2 and TESS transiting planets") ESPRESSO [@pepe21] for the selection and preliminary analysis of candidates.


# Statement of Need

With the advent of all sky photometric surveys looking for exoplanets, a big number of interesting systems are detected in a short amount of time. Given the high instrumental time required for an efficient radial velocity follow-up, most of these targets cannot realistically be observed using ground facilities. <!-- An efficient follow-up strategy is thus very important. -->
Our tool aims at selecting suitable RV follow-up targets more effectively and making their analysis easier and more convenient. By combining information coming from photometric observations with spectroscopic ones, one might reduce the necessary observation time per target, allowing for the characterization of more systems with a given instrument.
We also aimed for practicity, by using a single tool to retrieve the lightcurves and stellar parameters, correct the activity, look for exoplanets and estimate their fundamental parameters. 


# Documentation 

The documentation for **FULMAR** is available at [fulmar-astro.readthedocs.io](https://fulmar-astro.readthedocs.io/en/latest/) where it is hosted on [ReadTheDocs](https://readthedocs.org). It contains installation instructions, API and basic tutorials.


# Similar tools

The field of exoplanet research is very thriving. Many tools are arose with the ever-growing need of detecting and characterizing exoplanets. 

Some of the most popular tools in the field are (not exhaustive) `EXOFAST` [@eastman13; @eastman19], `PYTRANSIT` [@parviainen15], `Lightkurve` [@lightkurve18],`juliet` [@espinoza19], `exostriker` [@trifonov19], `PYANETI`[@barragan19], `allesfitter` [@guenther20], and `LATTE` [@eisner20].

These tools mostly focus on either fitting planetary parameters to the data, getting the data ready for further analysis, or vetting candidates using a GUI. **FULMAR** does not disregard the (excellent) existing tools, it was actually built upon them, namely `Astropy`[@astropy13; @astropy18], `exoplanet` [@foremanmackey21], `Lightkurve` [@lightkurve18], and `transitleastsquares` [@hippke19] to provide an "all-in one" solution.


# Acknowledgements

This work was supported by FCT through national funds (PTDC/FIS-AST/28953/2017) and by FEDER - Fundo Europeu de Desenvolvimento Regional through COMPETE2020 - Programa Operacional Competitividade e Internacionalização (POCI-01-0145-FEDER-028953) and through national funds (PIDDAC) by the grant UID/FIS/04434/2019. 


# References
