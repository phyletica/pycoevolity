#! /usr/bin/env python

"""
CLI program for summarizing divergence times from ecoevolity state log files.
"""

import os
import sys
import argparse
import subprocess

import sumcoevolity

def main(argv = sys.argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('log_paths',
            metavar = 'IPYRAD-LOCI-FILE-PATH',
            nargs = '+',
            type = sumcoevolity.argparse_utils.arg_is_file,
            help = ('Paths to ecoevolity state log files.'))
    parser.add_argument('-b', '--burnin',
            action = 'store',
            type = sumcoevolity.argparse_utils.arg_is_nonnegative_int,
            default = 0,
            help = ('The number of samples to remove from the beginning of '
                    'each log file as burn in.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    posterior = sumcoevolity.posterior.PosteriorSample(paths = args.log_paths,
            burnin = args.burnin)
    labels, heights = posterior.get_heights_2d(
            label_map = None,
            include_model_indices = True)

    rscript = """#! /usr/bin/Rscript

library(ggplot2)
library(ggjoy)

divergence_time = c({heights})
comparison = c({labels})

data <- data.frame(divergence_time = divergence_time, comparison = comparison) 

ggplot(data, aes(x = divergence_time, y = comparison, height = ..density..)) +
    geom_joy(stat = "density", scale = 8) +
    theme_minimal(base_size = 14) +
    theme(axis.text.y = element_text(vjust = 0)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_discrete(expand = c(0.01, 0)) +
    labs(x = "Population divergence time") +
    labs(y = "Comparison")

ggsave("divtimes.pdf")
""".format(
            heights = ", ".join(str(h) for h in heights),
            labels = ", ".join(repr(l) for l in labels))

    r_path = "plot-div-times.R"
    with os.fdopen(os.open(r_path, os.O_WRONLY | os.O_CREAT, 0o744), "w") as out:
        out.write("{0}\n".format(rscript))

    sout = subprocess.PIPE
    serr = subprocess.PIPE
    process = subprocess.Popen(["Rscript", "--vanilla", r_path],
            stdout = sout,
            stderr = serr,
            shell = False)
    stdout, stderr = process.communicate()
    exit_code = process.wait()
    if exit_code != 0:
        sys.stderr.write(
                "ERROR: The R plotting script exited with an error code.\n"
                "However, the script is available at\n"
                "{r_script_path!r}.\n"
                "You may need to install the R packages ggplot2 and ggjoy and "
                "re-run the R script.\n"
                "Here is the stderr from R:\n{stderr}".format(
                    r_script_path = r_path,
                    stderr = stderr))


if __name__ == "__main__":
    main()
