#! /usr/bin/env python

import sys
import os
import random
import string

def random_str(length=8,
        char_pool=string.ascii_letters + string.digits):
    return ''.join(random.choice(char_pool) for i in range(length))
