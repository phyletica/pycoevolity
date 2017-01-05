#! /usr/bin/env python

class SumcoevolityError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class ArgumentError(SumcoevolityError):
    def __init__(self, *args, **kwargs):
        SumcoevolityError.__init__(self, *args, **kwargs)

class TempFSError(SumcoevolityError):
    def __init__(self, *args, **kwargs):
        SumcoevolityError.__init__(self, *args, **kwargs)
