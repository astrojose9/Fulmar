![Logo](https://raw.githubusercontent.com/astrojose9/fulmar/main/docs/source/FULMAR_logo_title.png)
### A modular tool for analyzing light curves in support of RV follow-up programs.
[![Image](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/astrojose9/fulmar/blob/main/LICENSE)

**FULMAR** is an open source Python package that was created to assist RV follow-up programs by making the analysis of the light curves easier. It provides tools to correct stellar activity, to look for transits, to refine transit parameters and to visually probe signals detected in RV.
Our tool aims at selecting suitable RV follow-up targets more effectively and making their analysis easier. It was build in a modular way, making new features easier to implement.



## Installation

FULMAR can be installed using: `pip install fulmar-astro`

If you have multiple versions of Python and pip on your machine, try: `pip3 install fulmar-astro`

pip might output the following error message:
```
ERROR: Could not find a version that satisfies the requirement argparse (from transitleastsquares) (from versions: none)
ERROR: No matching distribution found for argparse
```
It comes from the new behaviour of pip. More info [here](https://github.com/pypa/pip/issues/9035#issuecomment-714595232). ~~A pull request with a patch fixing the issue was made on [TransitLeastSquares](https://github.com/hippke/tls)' repo.~~ The issue was fixed in transitleastsquares 1.0.31.


The latest version can be pulled from github::
```
git clone https://github.com/astrojose9/fulmar.git
cd fulmar
python setup.py install
```

If the command `python` does not point to Python 3 on your machine, you can try to replace the last line with `python3 setup.py install`. If you don't have `git` on your machine, you can find installation instructions [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git). TLS also runs on Python 2, but without multi-threading.



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
Read the [documentation](https://fulmar-astro.readthedocs.io/en/latest/).



## Contributing Code, Bugfixes, or Feedback
We welcome and encourage contributions. If you have any trouble, [open an issue](https://github.com/astrojose9/fulmar/issues).



## License
FULMAR is distributed under MIT License.



Copyright 2021, José Rodrigues.
