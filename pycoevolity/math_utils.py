#! /usr/bin/env python

import math

def diff_almost_zero(x, y, abs_tol = 1e-8):
    return math.isclose(x - y, 0.0, abs_tol = abs_tol, rel_tol = 0)
