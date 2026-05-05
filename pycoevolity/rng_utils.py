#!/usr/bin/env python

import os
import math
import random
import numpy as np

def get_safe_seed(rng = None):
    """
    Get a random seed (int) between 0 and 2^31, which is safe to write and read
    across systems and is safe for seeding numpy (which must be between 0 and
    2^32-1).

    Parameters
    ----------
    rng : `random.Random` object
        An instance of a `random.Random` object

    Returns
    -------
    int
        A random integer from the range 1 to 2^31-1 (inclusive)
    """
    if not rng:
        return random.randint(1, (2**31)-1)
    return rng.randint(1, (2**31)-1)

def get_safe_seeds(rng, n):
    return (get_safe_seed(rng) for _ in range(n))

def get_numpy_rng(self, seed = None):
    if seed is None:
        return np.random.default_rng()
    return np.random.default_rng(seed)
