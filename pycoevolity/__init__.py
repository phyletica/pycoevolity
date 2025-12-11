#! /usr/bin/env python

import sys
import os

PACKAGE_DIR = os.path.abspath(__file__)
BASE_DIR = os.path.dirname(PACKAGE_DIR)

__project__ = "pycoevolity"

# NOTE: Imports to populate the namespace can break the scripts' control of the
# logging level, because imported modules will initiate their loggers before
# the CLI scripts can update LoggingControl.

import pycoevolity.stats
import pycoevolity.parsing
import pycoevolity.posterior
import pycoevolity.argparse_utils
import pycoevolity.partition

def _get_git_data(repo_path):
    version = "unknown"
    branch = "unknown"
    commit = "unknown"
    commit_time = "unknown"
    if repo_path is not None:
        try:
            import subprocess
            import datetime

            try:
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
            except:
                pass

            try:
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
            except:
                pass

            try:
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
            except:
                pass

            try:
                p = subprocess.Popen(
                        ["git", "describe", "--abbrev=0"],
                        shell = False,
                        cwd = repo_path,
                        stdin = subprocess.PIPE,
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE,
                        universal_newlines = True)
                stdout, stderr = p.communicate()
                exit_code = p.wait()
                version = stdout.strip()
            except:
                pass
        except:
            pass

    return version, branch, commit, commit_time

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

# Get version and git data by calling git
__version__, __branch__, __commit__, __committime__ = _get_git_data(__homedir__)

# If importlib is available, use it to get the version and other metadata
__license__ = "unknown"
__description__ = "Package for summarizing output of ecoevolity"
try:
    import importlib.metadata
    __version__ = importlib.metadata.version(__project__)
    pkg_metadata = importlib.metadata.metadata(__project__)
    __license__ = pkg_metadata.get("License", __license__)
    __description__ = pkg_metadata.get("Description", __description__)
except:
    pass

def get_description():
    d = "{0} version {1}".format(__project__, __version__)

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
            "Version {v}".format(
                    v = __version__),
            w))
    stream.write("{0:^{1}}\n".format(
            "License: {}".format(__license__),
            w))
    stream.write("{0}\n".format("=" * w))

def write_dummy_biallelic_data_file(
        nspecies = 2,
        ngenomes = 10,
        ncharacters = 1000,
        prefix = "sp",
        out = sys.stdout):
    nspecies_padding = len(str(nspecies))
    out.write("---\n");
    out.write("markers_are_dominant: false\n")
    out.write("population_labels: [")
    for sp_idx in range(nspecies):
        if sp_idx > 0:
            out.write(", ")
        sp_label = prefix + "{n:0{padding}d}".format(
                n = sp_idx + 1,
                padding = nspecies_padding)
        out.write("{0}".format(sp_label))
    out.write("]\n")
    out.write("allele_count_patterns:\n")
    out.write("    - [")
    for sp_idx in range(nspecies):
        if sp_idx > 0:
            out.write(", ")
        out.write("[1, {0}]".format(ngenomes))
    out.write("]\n")
    out.write("pattern_weights: [{0}]\n".format(ncharacters))
