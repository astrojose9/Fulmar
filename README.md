![Logo](https://raw.githubusercontent.com/astrojose9/fulmar/main/docs/source/FULMAR_logo_title.png)
### A modular tool for analyzing light curves in support of RV follow-up programs.
[![PyPI - License](https://img.shields.io/pypi/l/fulmar-astro?color=brightgreen)](https://github.com/astrojose9/fulmar/blob/main/LICENSE) [![PyPI](https://img.shields.io/pypi/v/fulmar-astro?color=brightgreen)](https://pypi.org/project/fulmar-astro/) [![readthedocs](https://img.shields.io/badge/read%20the-docs-brightgreen)](https://fulmar-astro.readthedocs.io/en/latest/)

**FULMAR** is an open source Python package that was created to assist radial velocity follow-up programs by making the analysis of the light curves easier.
It provides tools to download lightcurves, correct stellar activity, to look for transits, to refine transit parameters, to estimate the amplitude of the corresponding RV signal, and to visually probe signals detected in RV.
Our tool aims at selecting suitable RV follow-up targets more effectively and making their analysis easier. It was build in a modular way, making new features easier to implement.



## Installation

FULMAR can be installed using: `pip install fulmar-astro`

If you have multiple versions of Python and pip on your machine, try: `pip3 install fulmar-astro`

The latest version can be pulled from github::
```
git clone https://github.com/astrojose9/fulmar.git
cd fulmar
python setup.py install
```

If the command `python` does not point to Python 3 on your machine, you can try to replace the last line with `python3 setup.py install`. If you don't have `git` on your machine, you can find installation instructions [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).



**Dependencies**:
Python 3,
[Arviz](https://arviz-devs.github.io/arviz/),
[Astropy](https://www.astropy.org/),
[celerite2](https://celerite2.readthedocs.io/en/latest/)
[corner](https://github.com/dfm/corner.py),
[exoplanet](https://docs.exoplanet.codes/en/latest/),
[Lightkurve](https://docs.lightkurve.org/),
[Matplotlib](https://matplotlib.org/),
[NumPy](https://www.numpy.org/),
[pymc3-ext](https://github.com/exoplanet-dev/pymc3-ext),
[TransitLeastSquares](https://github.com/hippke/tls)


If you have trouble installing, please [open an issue](https://github.com/astrojose9/fulmar/issues).


## Documentation
The [documentation](https://fulmar-astro.readthedocs.io/en/latest/) is hosted on [ReadTheDocs](https://readthedocs.org).



## Contributing Code, Bugfixes, or Feedback
We welcome and encourage contributions. If you have any trouble, [open an issue](https://github.com/astrojose9/fulmar/issues).



## License
FULMAR is distributed under MIT License.



Copyright 2021, José Rodrigues.
