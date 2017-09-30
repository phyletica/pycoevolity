#! /usr/bin/env python

class PycoevolityError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class ArgumentError(PycoevolityError):
    def __init__(self, *args, **kwargs):
        PycoevolityError.__init__(self, *args, **kwargs)

class TempFSError(PycoevolityError):
    def __init__(self, *args, **kwargs):
        PycoevolityError.__init__(self, *args, **kwargs)
