#!/usr/bin/python

from os import path
import fulmar.version as versinfo

# Mission dictionary

MISSION_DIC = {
    'tess': 'TESS',
    'kepler': 'Kepler',
    'k2': 'K2'
}

# Path for the module
fulmar_dir = path.join(path.dirname(__file__))

FULMAR_VERSION_STR = (
    "Fulmar " + versinfo.FULMAR_VERSION + " (" + versinfo.FULMAR_DATE + ")"
)
