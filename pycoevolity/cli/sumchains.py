#! /usr/bin/env python

"""
CLI program for diagnosing MCMC convergence from ecoevolity state log files.
"""

import os
import sys
import argparse

import pycoevolity

def main(argv = sys.argv):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('log_paths',
            metavar = 'IPYRAD-LOCI-FILE-PATH',
            nargs = '+',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Paths to ecoevolity state log files.'))
    parser.add_argument('-s', '--burnin-step-size',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            default = 100,
            help = ('The step size (in the number of samples) to evaluate '
                    'convergence.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    cs = pycoevolity.posterior.ChainConvergenceSummary(
            paths = args.log_paths,
            burnin_step_size = args.burnin_step_size)

    cs.write_summary(out = sys.stdout, err = sys.stderr)


if __name__ == "__main__":
    main()
