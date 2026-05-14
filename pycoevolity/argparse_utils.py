#! /usr/bin/env python

import os
import argparse


class SmartHelpFormatter(argparse.HelpFormatter):
    '''
    A class to allow customizable line breaks for an argument help message on a
    per argument basis.
    '''

    def _split_lines(self, text, width):
        ret = split_arg_help_lines(text)
        if ret is None:
            return argparse.HelpFormatter._split_lines(self, text, width)
        return ret


class SmartDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    '''
    A class to allow customizable line breaks for an argument help message on a
    per argument basis, and to include argument defaults in the message.
    '''

    def _split_lines(self, text, width):
        ret = split_arg_help_lines(text)
        if ret is None:
            return argparse.ArgumentDefaultsHelpFormatter._split_lines(self, text, width)
        return ret

def split_arg_help_lines(text):
    if text.startswith('RAW|'):
        return text[4:].splitlines()
    return None

def arg_is_path(path):
    try:
        if not os.path.exists(path):
            raise
    except:
        msg = 'path {0!r} does not exist'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def arg_is_file(path):
    try:
        if not os.path.isfile(path):
            raise
    except:
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def arg_is_dir(path):
    try:
        if not os.path.isdir(path):
            raise
    except:
        msg = '{0!r} is not a directory'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def arg_is_nonnegative_int(i):
    try:
        if int(i) < 0:
            raise
    except:
        msg = '{0!r} is not a non-negative integer'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)

def arg_is_positive_int(i):
    try:
        if int(i) < 1:
            raise
    except:
        msg = '{0!r} is not a positive integer'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)

def arg_is_positive_float(i):
    try:
        if float(i) <= 0.0:
            raise
    except:
        msg = '{0!r} is not a positive real number'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return float(i)

def arg_is_nonnegative_float(i):
    try:
        if float(i) < 0.0:
            raise
    except:
        msg = '{0!r} is not a non-negative real number'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return float(i)

def arg_is_dir_or_new_dir(path):
    """
    Returns the passed string if it is a valid path to a directory, or its
    parent is a valid directory. Otherwise raises an `ArgumentTypeError`.

    Examples
    --------
    >>> d = os.path.abspath(os.path.dirname(__file__))
    >>> returned = arg_is_dir_or_new_dir(d)
    >>> returned == d
    True
    >>> new_dir = os.path.join(d, "probably-not-a-dir-in-this-dir")
    >>> returned = arg_is_dir_or_new_dir(new_dir)
    >>> returned == new_dir
    True
    """
    if os.path.isdir(path):
        return path
    elif os.path.exists(path):
        msg = 'path {0!r} exists but is not a directory'.format(path)
    elif os.path.sep not in path:
        # just dir name which can be created in working dir with mkdir
        return path
    elif os.path.isdir(os.path.dirname(path)):
        # path doesn't exist, but is in an existing parent directory
        return path
    else:
        msg = '{0!r} is not a directory nor is its parent'.format(path)
    raise argparse.ArgumentTypeError(msg)

def process_output_dir_arg(output_dir):
    if not output_dir:
        output_dir = os.curdir
    else:
        if not os.path.exists(output_dir):
            try:
                os.mkdir(output_dir)
            except Exception as e:
                sys.stderr.write(
                    f"ERROR: Could not create output directory '{output_dir}'\n"
                )
                raise e
    return output_dir
