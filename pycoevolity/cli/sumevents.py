#! /usr/bin/env python

"""
CLI program for summarizing the number of events from sumcoevolity output.
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

    parser.add_argument('nevents_path',
            metavar = 'SUMCOEVOLITY-NEVENTS-FILE-PATH',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Path to the \'sumcoevolity-results-nevents.txt\' file '
                    'output by sumcoevolity.'))
    parser.add_argument('--no-legend',
            action = 'store_true',
            help = ('Do not include legend in the plot.'))
    parser.add_argument('--legend-in-plot',
            action = 'store_true',
            help = ('Place legend within plot (rather than above).'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "",
            help = ('A prefix to prepend to all output files.'))
    parser.add_argument('-f', '--force',
            action = 'store_true',
            help = ('Overwrite any existing output files. By default, an error '
                    'is thrown if an output path exists.'))
    parser.add_argument('-x', '--x-label',
            action = 'store',
            type = str,
            default = "Number of events",
            help = ('Label for the X-axis. Default: \'Number of events\'.'))
    parser.add_argument('-y', '--y-label',
            action = 'store',
            type = str,
            default = "Probability",
            help = ('Label for the Y-axis. Default: \'Probability\'.'))
    parser.add_argument('--bf-font-size',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_nonnegative_float,
            default = 4.0,
            help = ('Font size for the Bayes factor labels. Default: 4.0.'))
    parser.add_argument('-w', '--width',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_float,
            default = 7.0,
            help = ('The width of the plot. Default: 7.0.'))
    parser.add_argument('--full-prob-axis',
            action = 'store_true',
            help = ('Force probability (Y) axis to range from 0-1.'))
    parser.add_argument('--use-r',
            action = 'store_true',
            help = ('Create and execute an R script for plotting with ggplot2. '
                    'Default is to use matplotlib for plotting.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    prefix = args.prefix
    if len(prefix.split(os.path.sep)) < 2:
        prefix = os.path.join(os.curdir, prefix)

    pdf_path = prefix + "pycoevolity-nevents.pdf"
    output_dir = os.path.dirname(pdf_path)
    if not output_dir:
        output_dir = os.curdir
    if not args.force:
        if os.path.exists(pdf_path):
            raise Exception(
                    "\nERROR: File {0!r} already exists.\n"
                    "Use \'-p/--prefix\' option to specify a different prefix,\n"
                    "or the \'-f/--force\' option to overwrite existing "
                    "files.".format(pdf_path))

    sys.stderr.write("Parsing nevents file...\n")
    nevents = pycoevolity.posterior.SumcoevolityNeventsTable(args.nevents_path)

    max_prob = max(nevents.posterior_probs)
    prior_probs = [0.0] * nevents.number_of_elements
    bfs = ["0"] * nevents.number_of_elements
    if not nevents.no_prior:
        max_prob = max(nevents.posterior_probs + nevents.prior_probs)
        prior_probs = nevents.prior_probs 
        bfs = []
        for bayes_factor in nevents.bayes_factors:
            if bayes_factor is None:
                bfs.append("")
                continue
            bfs.append("{:.3g}".format(bayes_factor))
        for i, a in enumerate(nevents.bayes_factors_annotations):
            if a:
                bfs[i] = a + bfs[i]

    plot_width = args.width
    plot_height = plot_width / 1.618034

    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib import gridspec
    except ImportError as e:
        sys.stderr.write('WARNING: matplotlib could not be imported; '
                'it is needed for Python plotting. Will try using '
                'R instead.')
        args.use_r = True

    add_legend = (not nevents.no_prior) and (not args.no_legend)
    bump_bfs = add_legend and (args.legend_in_plot)

    if not args.use_r:
        #######################################################################
        # Using matplotlib for plotting
        #######################################################################
        # Use TrueType (42) fonts rather than Type 3 fonts
        mpl.rcParams["pdf.fonttype"] = 42
        mpl.rcParams["ps.fonttype"] = 42

        fig = plt.figure(figsize = (plot_width, plot_height))
        gs = gridspec.GridSpec(1, 1,
                wspace = 0.0,
                hspace = 0.0)
        ax = plt.subplot(gs[0, 0])

        nevents_indices = [float(x) for x in range(nevents.number_of_elements)]
        bar_width = 0.45
        if nevents.no_prior:
            bar_width *= 2.0
        posterior_color = "0.3"
        prior_color = "0.85"

        bars_posterior = ax.bar(
                nevents_indices,
                nevents.posterior_probs,
                bar_width,
                color = posterior_color,
                label = "Posterior")
        if not nevents.no_prior:
            bars_prior = ax.bar(
                    [x + bar_width for x in nevents_indices],
                    prior_probs,
                    bar_width,
                    color = prior_color,
                    label = "Prior")

        ax.set_xlabel("Number of events")
        ax.set_ylabel("Probability")
        x_tick_labels = [str(i + 1) for i in range(nevents.number_of_elements)]
        if nevents.no_prior:
            ax.set_ylabel("Posterior probability")
            plt.xticks(
                    nevents_indices,
                    x_tick_labels
                    )
            if args.full_prob_axis:
                ax.set_ylim(0.0, 1.0)
        else:
            y_min, y_max = ax.get_ylim()
            y_max *= 1.08
            if args.full_prob_axis:
                y_min = 0.0
                if y_max < 1.0:
                    y_max = 1.0
            ax.set_ylim(y_min, y_max)
            bar_midpoints = [x + (bar_width / 2.0) for x in nevents_indices]
            plt.xticks(
                    bar_midpoints,
                    x_tick_labels
                    )
            for i, x in enumerate(bar_midpoints):
                upper_loc = 0.99
                lower_loc = 0.94
                if bump_bfs:
                    upper_loc = 0.89
                    lower_loc = 0.84
                bf = bfs[i]
                y = y_max * upper_loc
                if ((i + 1) % 2) == 0:
                    y = y_max * lower_loc
                ax.text(x, y, bf,
                        horizontalalignment = "center",
                        verticalalignment = "top",
                        size = 8.0,
                        zorder = 300)

        if add_legend:
            loc = 'lower center'
            if args.legend_in_plot:
                loc = 'upper center'
            ax.legend(loc=loc, bbox_to_anchor=(0.5, 1.0), ncol=2)
        else:
            l = ax.legend()
            if l:
                l.remove()
        fig.tight_layout()
        plt.savefig(pdf_path)

        sys.exit()


    ###########################################################################
    # Creating and executing R script for plotting
    ###########################################################################
    r_path = prefix + "pycoevolity-plot-nevents.R"
    png_path = prefix + "pycoevolity-nevents.png"
    svg_path = prefix + "pycoevolity-nevents.svg"
    if not args.force:
        for p in [r_path, png_path, svg_path]:
            if os.path.exists(p):
                raise Exception(
                        "\nERROR: File {0!r} already exists.\n"
                        "Use \'-p/--prefix\' option to specify a different prefix,\n"
                        "or the \'-f/--force\' option to overwrite existing "
                        "files.".format(p))

    plot_prior = "FALSE"
    if not nevents.no_prior:
        plot_prior = "TRUE"
    plot_units = "in"
    plot_scale = 8
    plot_base_size = 14

    legend_arg = "legend.title = element_blank()"
    if args.no_legend:
        legend_arg = "legend.position = \"none\""

    rscript = """#! /usr/bin/env Rscript

library(ggplot2)
library(ggridges)

plot_prior <- {plot_prior}
nevents <- c(\"{nevents}\")
posterior_probs <- c({posterior_probs})
prior_probs <- c({prior_probs})
bf_labels <- c(\"{bayes_factor_strings}\")
max_prob <- {max_prob}
bf_position_bottom <- max_prob + (max_prob * 0.04)
bf_position_top <- bf_position_bottom + 0.04
bf_positions <- c()
for (i in 1:length(nevents)) {{
    if (i %% 2 == 0) {{
        bf_positions = c(bf_positions, bf_position_bottom)
    }} else {{
        bf_positions = c(bf_positions, bf_position_top)
    }}
}}

posterior_df <- data.frame(nevents = nevents, probability = posterior_probs, label = rep("posterior", {number_of_comparisons}))
prior_df <- data.frame(nevents = nevents, probability = prior_probs, label = rep("prior", {number_of_comparisons}))
if (plot_prior) {{
    data <- rbind(posterior_df, prior_df)
    bar_colors <- c("posterior" = "gray30", "prior" = "gray85")
}} else {{
    data <- posterior_df
    bar_colors <- c("posterior" = "gray30")
}}

if (plot_prior) {{
    ggplot(data, aes(x = nevents, y = probability, fill = label)) +
        geom_col(position = "dodge") +
        scale_x_discrete(limits = nevents) +
        theme_minimal(base_size = {plot_base_size}) +
        theme({legend_arg}) +
        scale_colour_manual(values = bar_colors) +
        scale_fill_manual(values = bar_colors) +
        labs(x = \"{x_label}\") +
        labs(y = \"{y_label}\") +
        annotate("text", x = nevents, y = bf_positions, label = bf_labels, size = {bf_font_size})
}} else {{
    ggplot(data, aes(x = nevents, y = probability, fill = label)) +
        geom_col(position = "dodge") +
        scale_x_discrete(limits = nevents) +
        theme_minimal(base_size = {plot_base_size}) +
        theme({legend_arg}) +
        scale_colour_manual(values = bar_colors) +
        scale_fill_manual(values = bar_colors) +
        labs(x = \"{x_label}\") +
        labs(y = \"{y_label}\")
}}

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
            plot_prior = plot_prior,
            nevents = "\", \"".join(str(i + 1) for i in range(nevents.number_of_elements)),
            posterior_probs = ", ".join(str(p) for p in nevents.posterior_probs),
            prior_probs = ", ".join(str(p) for p in prior_probs),
            bayes_factor_strings = "\", \"".join(bfs),
            max_prob = max_prob,
            number_of_comparisons = nevents.number_of_elements,
            legend_arg = legend_arg,
            x_label = args.x_label,
            y_label = args.y_label,
            bf_font_size = args.bf_font_size,
            plot_base_size= plot_base_size,
            plot_width = plot_width,
            plot_height = plot_height,
            plot_units = plot_units,
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
