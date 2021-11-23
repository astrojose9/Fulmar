#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __init__.py
from __future__ import division

from .main import *
from .time import *
from .utils import *

import warnings

# Change default behaviour of printing warnings
warnings.formatwarning = warning_on_one_line
