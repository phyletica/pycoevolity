#! /usr/bin/env python

import sys
import os

from .metadata import *


PACKAGE_DIR = os.path.abspath(__file__)
BASE_DIR = os.path.dirname(PACKAGE_DIR)

# NOTE: Imports to populate the namespace can break the scripts' control of the
# logging level, because imported modules will initiate their loggers before
# the CLI scripts can update LoggingControl.

import pycoevolity.stats
import pycoevolity.parsing
import pycoevolity.posterior
import pycoevolity.argparse_utils
import pycoevolity.partition

def _get_git_data(repo_path):
    try:
        import subprocess
        import datetime

        p = subprocess.Popen(
                ["git", "rev-parse", "HEAD"],
                shell = False,
                cwd = repo_path,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                universal_newlines = True)
        stdout, stderr = p.communicate()
        exit_code = p.wait()
        commit = stdout.strip()[0:7]

        p = subprocess.Popen(
                ["git", "name-rev", "--name-only", "HEAD"],
                shell = False,
                cwd = repo_path,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                universal_newlines = True)
        stdout, stderr = p.communicate()
        exit_code = p.wait()
        branch = stdout.strip()

        p = subprocess.Popen(
                ["git", "show", "--quiet", "--pretty=format:'%at'", "HEAD"],
                shell = False,
                cwd = repo_path,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                universal_newlines = True)
        stdout, stderr = p.communicate()
        exit_code = p.wait()
        t = stdout.strip().replace("'", "").replace('"', '')
        commit_time = datetime.datetime.fromtimestamp(float(t))

        return branch, commit, commit_time
    except:
        return None, None, None

__homedir__ = None
try:
    try:
        __homedir__ = __path__[0]
    except AttributeError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except IndexError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
except:
    pass

__gitinfo__ = ""

__branch__, __commit__, __committime__ = _get_git_data(__homedir__)

def get_description():
    d = "{0} version {1}".format(__project__, __version__)
    if __branch__:
        d += " {0}".format(__branch__)
    if __commit__:
        d += " {0}".format(__commit__)
    if __committime__:
        d += " {0}".format(__committime__)
    return d

def write_splash(stream, console_width = 72):
    w = console_width
    stream.write("{0}\n".format("=" * w))
    stream.write("{0:^{1}}\n".format(__project__.capitalize(), w))
    stream.write("{0:^{1}}\n\n".format(
            "Summarizing evolutionary coevality",
            w))
    stream.write("{0:^{1}}\n\n".format(
            "A Python package for analyzing the output of Ecoevolity",
            w))
    stream.write("{0:^{1}}\n\n".format(
            "Version {v} ({b} {c}: {t})".format(
                    v = __version__,
                    b = __branch__,
                    c = __commit__,
                    t = __committime__),
            w))
    stream.write("{0:^{1}}\n".format(
            "License: {}".format(__license_short__),
            w))
    # for line in __license__.strip().split("\n"):
    #     stream.write("{0:^{1}}\n".format(line, w))
    stream.write("{0}\n".format("=" * w))
