#! /usr/bin/env python

"""
CLI program for summarizing population sizes from ecoevolity state log files.
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
            metavar = 'ECOEVOLITY-STATE-LOG-PATH',
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
    parser.add_argument('-i', '--ignore',
            action = 'append',
            type = str,
            metavar = ("COMPARISON-LABEL"),
            help = ('Ignore (i.e., do not plot) the specified comparison. '
                    'This option takes one argument: The original label '
                    'used to identify a population (i.e., the prefix or suffix '
                    'in the original alignment analyzed by Ecoevolity) '
                    'Note, you need to quote labels with spaces. '
                    'This option can be used multiple times to ignore '
                    'multiple comparisons.'))
    parser.add_argument('--violin',
            action = 'store_true',
            help = ('Produce a violin plot, rather than a ridge plot.'))
    parser.add_argument('--relative',
            action = 'store_true',
            help = ('Only plot the relative sizes of ancestral populations.'))
    parser.add_argument('--x-limits',
            action = 'store',
            nargs = 2,
            type = float,
            metavar = ("LOWER-LIMIT", "UPPER-LIMIT"),
            help = ('Lower and upper limits for the X-axis. By default, '
                    'ggplot2 determines the limits.'))
    parser.add_argument('-x', '--x-label',
            action = 'store',
            type = str,
            help = ('Label for the X-axis. Default: \'Size\'.'))
    parser.add_argument('-y', '--y-label',
            action = 'store',
            type = str,
            default = "Population",
            help = ('Label for the Y-axis. Default: \'Population\'.'))
    parser.add_argument('-w', '--width',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_float,
            default = 7.0,
            help = ('The width of the plot. Default: 7.0.'))
    parser.add_argument('--base-font-size',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_float,
            default = 14.0,
            help = ('The base font size. Default: 14.0.'))
    parser.add_argument('--no-plot',
            action = 'store_true',
            help = ('Skip plotting; only report summary table.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    prefix = args.prefix
    if len(prefix.split(os.path.sep)) < 2:
        prefix = os.path.join(os.curdir, prefix)

    r_path = prefix + "pycoevolity-plot-sizes.R"
    pdf_path = prefix + "pycoevolity-sizes.pdf"
    png_path = prefix + "pycoevolity-sizes.png"
    svg_path = prefix + "pycoevolity-sizes.svg"
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

    if args.x_label is None:
        args.x_label = "Size"
        if args.relative:
            args.x_label = "Relative size of ancestral population"

    sys.stderr.write("Parsing log files...\n")
    posterior = pycoevolity.posterior.PosteriorSample(paths = args.log_paths,
            burnin = args.burnin)
    sys.stderr.write("Parsed {0} samples from {1} files.\n".format(
            posterior.number_of_samples,
            len(posterior.paths)))
    sys.stderr.write("Summary of sizes:\n")
    posterior.write_size_summary(out = sys.stdout)

    if args.no_plot:
        sys.exit(0)

    populations_to_ignore = []
    if args.ignore:
        populations_to_ignore = args.ignore

    plot_width = args.width 
    plot_height = plot_width / 1.618034

    if args.violin:
        try:
            import matplotlib as mpl
            import matplotlib.pyplot as plt
            from matplotlib import gridspec
        except ImportError as e:
            sys.stderr.write('ERROR: matplotlib could not be imported; it is needed for violin plots')
            raise e

        # Use TrueType (42) fonts rather than Type 3 fonts
        mpl.rcParams["pdf.fonttype"] = 42
        mpl.rcParams["ps.fonttype"] = 42

        ## These settings are nice, but require latex to be installed
        # tex_font_settings = {
        #         "text.usetex": True,
        #         "font.family": "sans-serif",
        #         "text.latex.preamble" : [
        #                 "\\usepackage[T1]{fontenc}",
        #                 "\\usepackage[cm]{sfmath}",
        #                 ]
        # }
        # mpl.rcParams.update(tex_font_settings)

        if args.relative:
            labels, sizes = posterior.get_labels_and_relative_sizes(
                    label_map = label_map,
                    populations_to_ignore = populations_to_ignore)
        else:
            labels, sizes = posterior.get_labels_and_sizes(
                    label_map = label_map,
                    populations_to_ignore = populations_to_ignore)

        fig = plt.figure(figsize = (plot_width, plot_height))
        gs = gridspec.GridSpec(1, 1,
                wspace = 0.0,
                hspace = 0.0)
        ax = plt.subplot(gs[0, 0])
        positions = range(1, len(sizes) + 1)
        v = ax.violinplot(sizes,
                positions = positions,
                vert = False,
                widths = 0.9,
                showmeans = False,
                showextrema = False,
                showmedians = False,
                points = 100,
                bw_method = None,
                )

        colors = ["gray"] * len(labels)
        for i in range(len(v["bodies"])):
            v["bodies"][i].set_alpha(1)
            v["bodies"][i].set_facecolor(colors[i])
            v["bodies"][i].set_edgecolor(colors[i])

        means = []
        ci_lower = []
        ci_upper = []
        for sample in sizes:
            summary = pycoevolity.stats.get_summary(sample)
            means.append(summary["mean"])
            ci_lower.append(summary["qi_95"][0])
            ci_upper.append(summary["qi_95"][1])
        ax.hlines(positions, ci_lower, ci_upper,
                colors = "black",
                linestyle = "solid",
                zorder = 100)
        ax.scatter(ci_lower, positions,
                marker = "|",
                color = "black",
                s = 120,
                zorder = 200,
                )
        ax.scatter(ci_upper, positions,
                marker = "|",
                color = "black",
                s = 120,
                zorder = 200,
                )
        ax.scatter(means, positions,
                marker = ".",
                color = "white",
                s = 50,
                zorder = 300,
                )

        if args.x_limits:
            xlims = sorted(args.x_limits)
            ax.set_xlim(xlims[0], xlims[1])

        ax.yaxis.set_ticks(range(1, len(labels) + 1))
        ytick_labels = [item for item in ax.get_yticklabels()]
        assert(len(ytick_labels) == len(labels))
        for i in range(len(ytick_labels)):
            ytick_labels[i].set_text(labels[i])
        ax.set_yticklabels(ytick_labels)

        ax.set_xlabel(
                args.x_label)
        ax.set_ylabel(
                args.y_label)

        fig.tight_layout()
        plt.savefig(pdf_path)
        sys.stderr.write("Here are the outputs:\n")
        sys.stderr.write("    PDF plot: {0!r}\n".format(pdf_path))
        sys.exit(0)

    if args.relative:
        labels, sizes = posterior.get_relative_population_sizes_2d(
                label_map = label_map)
    else:
        labels, sizes = posterior.get_population_sizes_2d(
                label_map = label_map)

    plot_units = "in"
    plot_scale = 8
    plot_base_size = args.base_font_size
    scale_x_continuous_args = ["expand = c(0.05, 0)"]
    if args.x_limits:
        scale_x_continuous_args.append("limits = c({0}, {1})".format(
                *sorted(args.x_limits)))

    rscript = """#! /usr/bin/env Rscript

library(ggplot2)
library(ggridges)

size = c({sizes})
comparison = c(\"{labels}\")

data <- data.frame(size = size, comparison = comparison)
data$comparison = factor(data$comparison, levels = rev(unique(as.character(data$comparison))))

ggplot(data, aes(x = size, y = comparison, height = ..density..)) +
    geom_density_ridges(stat = \"density\", scale = {plot_scale}, rel_min_height = 0.001) +
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
            sizes = ", ".join(str(s) for s in sizes),
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
                "You may need to install the R packages ggplot2 and ggridges and "
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
