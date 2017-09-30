#! /usr/bin/env python

"""
CLI program for summarizing event times from ecoevolity state log files.
"""

import os
import sys
import stat
import argparse
import subprocess
import logging

import pycoevolity

def main(argv = sys.argv):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('log_paths',
            metavar = 'IPYRAD-LOCI-FILE-PATH',
            nargs = '+',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Paths to ecoevolity state log files.'))
    parser.add_argument('-b', '--burnin',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
            default = 0,
            help = ('The number of samples to remove from the beginning of '
                    'each log file as burn in.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "",
            help = ('A prefix to prepend to all output files.'))
    parser.add_argument('-f', '--force',
            action = 'store_true',
            help = ('Overwrite any existing output files. By default, an error '
                    'is thrown if an output path exists.'))
    parser.add_argument('-l', '--label',
            action = 'append',
            type = str,
            nargs = 2,
            metavar = ("COMPARISON-LABEL", "REPLACEMENT-LABEL"),
            help = ('This option takes two arguments: (1) The original label '
                    'used to identify a population (i.e., the prefix or suffix '
                    'in the original alignment analyzed by Ecoevolity), and '
                    '(2) A label to use in the plots in place of the '
                    'original population label. '
                    'Note, you need to quote labels with spaces. '
                    'This option can be used multiple times to specify label '
                    'replacements for multiple taxa.'))
    parser.add_argument('-z', '--include-zero',
            action = 'store_true',
            help = ('By default, ggplot2 auto-magically determines the limits '
                    'of the time axis, which often excludes zero (present). '
                    'This option ensures that the time axis starts from zero.'))
    parser.add_argument('--include-map-model',
            action = 'store_true',
            help = ('Include event indices associated with the maximum a '
                    'posteriori (MAP) model in the comparison labels.'))
    parser.add_argument('-x', '--x-label',
            action = 'store',
            type = str,
            default = "Time",
            help = ('Label for the X-axis. Default: \'Time\'.'))
    parser.add_argument('-y', '--y-label',
            action = 'store',
            type = str,
            default = "Comparison",
            help = ('Label for the Y-axis. Default: \'Comparison\'.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    prefix = args.prefix
    if not prefix:
        prefix = os.path.join(os.curdir, "")

    r_path = prefix + "pycoevolity-plot-times.R"
    pdf_path = prefix + "pycoevolity-times.pdf"
    png_path = prefix + "pycoevolity-times.png"
    svg_path = prefix + "pycoevolity-times.svg"
    output_dir = os.path.dirname(r_path)
    if not output_dir:
        output_dir = os.curdir
    if not args.force:
        for p in [r_path, pdf_path, png_path, svg_path]:
            if os.path.exists(p):
                raise Exception(
                        "\nERROR: File {0!r} already exists.\n"
                        "Use \'-p/--prefix\' option to specify a different prefix,\n"
                        "or the \'-f/--force\' option to overwrite existing "
                        "files.".format(p))
    label_map = {}
    if args.label:
        for label, replacement in args.label:
            label_map[label] = replacement

    sys.stderr.write("Parsing log files...\n")
    posterior = pycoevolity.posterior.PosteriorSample(paths = args.log_paths,
            burnin = args.burnin)
    sys.stderr.write("Parsed {0} samples from {1} files.\n".format(
            posterior.number_of_samples,
            len(posterior.paths)))
    labels, heights = posterior.get_heights_2d(
            label_map = label_map,
            include_model_indices = args.include_map_model)

    plot_width = 7.0
    plot_height = plot_width / 1.618034
    plot_units = "in"
    plot_scale = 8
    plot_base_size = 14
    scale_x_continuous_args = ["expand = c(0.05, 0)"]
    if args.include_zero:
        scale_x_continuous_args.append("limits = c(0, NA)")

    rscript = """#! /usr/bin/env Rscript

library(ggplot2)
library(ggjoy)

time = c({heights})
comparison = c(\"{labels}\")

data <- data.frame(time = time, comparison = comparison)

ggplot(data, aes(x = time, y = comparison, height = ..density..)) +
    geom_joy(stat = \"density\", scale = {plot_scale}, rel_min_height = 0.001) +
    theme_minimal(base_size = {plot_base_size}) +
    theme(axis.text.y = element_text(vjust = 0)) +
    scale_x_continuous({scale_x_continuous_args}) +
    scale_y_discrete(expand = c(0.01, 0)) +
    labs(x = \"{x_label}\") +
    labs(y = \"{y_label}\")

ggsave(\"{pdf_path}\", width = {plot_width}, height = {plot_height}, units = \"{plot_units}\")
ggsave(\"{png_path}\", width = {plot_width}, height = {plot_height}, units = \"{plot_units}\")
r <- tryCatch(
    {{
        ggsave(\"{svg_path}\", width = {plot_width}, height = {plot_height}, units = \"{plot_units}\")
    }},
    error = function(cond) {{
        message(\"An error occurred while trying to save plot as SVG.\")
        message(\"The plot has been saved in PDF and PNG format.\")
        message(\"If you want the SVG file, you may need to install additional R packages.\")
        message(\"Here's the original error message for details:\")
        message(cond)
    }},
    warning = function(cond) {{
        message(\"A warning occurred while trying to save the plot in SVG format.\")
        message(\"The plot has been saved in PDF and PNG format.\")
        message(\"If you want the SVG file, you may need to install additional R packages.\")
        message(\"Here's the original warning message for details:\")
        message(cond)
    }},
    finally =  {{}})
""".format(
            heights = ", ".join(str(h) for h in heights),
            labels = "\", \"".join(labels),
            plot_scale = plot_scale,
            plot_base_size= plot_base_size,
            scale_x_continuous_args = ", ".join(scale_x_continuous_args),
            plot_width = plot_width,
            plot_height = plot_height,
            plot_units = plot_units,
            x_label = args.x_label,
            y_label = args.y_label,
            pdf_path = os.path.basename(pdf_path),
            png_path = os.path.basename(png_path),
            svg_path = os.path.basename(svg_path))

    with open(r_path, "w") as out:
        out.write("{0}".format(rscript))
    file_stat = os.stat(r_path)
    os.chmod(r_path, file_stat.st_mode | stat.S_IEXEC)

    sys.stderr.write("Running R script to generate plots...\n")
    sout = subprocess.PIPE
    serr = subprocess.PIPE
    process = subprocess.Popen([r_path],
            cwd = output_dir,
            stdout = sout,
            stderr = serr,
            shell = False,
            universal_newlines = True)
    stdout, stderr = process.communicate()
    exit_code = process.wait()
    if exit_code != 0:
        sys.stderr.write(
                "The R plotting script exited with an error code.\n"
                "However, the script is available at\n"
                "{r_script_path!r}.\n"
                "You may need to install the R packages ggplot2 and ggjoy and "
                "re-run the R script.\n"
                "Here is the stderr from R:\n{stderr}\n".format(
                    r_script_path = r_path,
                    stderr = stderr))
    else:
        if stderr:
            sys.stderr.write("Here is the stderr returned by R:\n")
            sys.stderr.write("{0}\n".format("-" * 72))
            sys.stderr.write("{0}\n".format(stderr))
            sys.stderr.write("{0}\n".format("-" * 72))
    if os.path.exists(r_path):
        sys.stderr.write("Here are the outputs:\n")
        sys.stderr.write("    R script: {0!r}\n".format(r_path))
        if os.path.exists(pdf_path):
            sys.stderr.write("    PDF plot: {0!r}\n".format(pdf_path))
        if os.path.exists(png_path):
            sys.stderr.write("    PNG plot: {0!r}\n".format(png_path))
        if os.path.exists(svg_path):
            sys.stderr.write("    SVG plot: {0!r}\n".format(svg_path))


if __name__ == "__main__":
    main()
