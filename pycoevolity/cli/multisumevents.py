#! /usr/bin/env python

"""
CLI program for summarizing the number of events from multiple sumcoevolity outputs.
"""

import os
import sys
import stat
import argparse
import subprocess
import logging
import math

import pycoevolity

def main(argv = sys.argv):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('--paths',
            metavar = 'SUMCOEVOLITY-NEVENTS-FILE-PATH',
            type = pycoevolity.argparse_utils.arg_is_file,
            nargs = '+',
            required = True,
            help = ('Paths to the \'sumcoevolity-results-nevents.txt\' files '
                    'output by sumcoevolity. Priors among files must match.'))
    parser.add_argument('--labels',
            metavar = 'LABEL-FOR-NEVENTS-FILE',
            type = str,
            nargs = '+',
            help = ('Label for each sumcoevolity nevents file.'))
    parser.add_argument('--colors',
            metavar = 'BAR-COLOR-FOR-NEVENTS-FILE',
            type = str,
            nargs = '+',
            help = ('Bar color for each sumcoevolity nevents file.'))
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
            default = None,
            help = ('Label for the Y-axis. Default: \'Probability\' '
                    'if there are posterior and prior probabilities to plot; '
                    '\'Posterior probability\' if there are only posterior '
                    'probabilities to plot.'))
    parser.add_argument('-w', '--width',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_float,
            default = 7.0,
            help = ('The width of the plot. Default: 7.0.'))
    parser.add_argument('--full-prob-axis',
            action = 'store_true',
            help = ('Force probability (Y) axis to range from 0-1.'))
    parser.add_argument('--log-bf',
            action = 'store_true',
            help = ('Log transform Bayes factors.'))
    parser.add_argument('--bf-marker-size',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_float,
            default = 8.0,
            help = ('Size BF markers. Default: 8.0.'))
    parser.add_argument('--bf-marker-border-color',
            action = 'store',
            type = str,
            default = "white",
            help = ('Border color for BF markers. Default: \'white\'.'))
    parser.add_argument('--no-legend',
            action = 'store_true',
            help = ('Do not include legend in the plot.'))
    parser.add_argument('--legend-in-plot',
            action = 'store_true',
            help = ('Place legend within plot (rather than above).'))
    parser.add_argument('--legend-font-size',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_float,
            default = 9.0,
            help = ('Font size used in legend. Default: 9.0.'))
    parser.add_argument('--prior-label',
            action = 'store',
            type = str,
            default = "Prior",
            help = ('Label for prior probability in legend. '
                    'Default: \'Prior\'.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    if len(args.paths) < 2:
        raise Exception(
            "ERROR: Multiple sumcoevolity nevents files are required"
        )

    if not args.labels:
        args.labels = [f'file-{i+1}' for i in range(len(args.paths))]

    if not args.colors:
        args.colors = [f'C{i}' for i in range(len(args.paths))]

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

    sys.stderr.write("Parsing nevents files...\n")
    data = pycoevolity.posterior.MultiSumcoevolityNeventsTable(args.paths)

    if args.y_label is None:
        args.y_label = "Probability"
        if data.no_prior:
            args.y_label = "Posterior probability"

    max_prob = data.max_prob
    if not data.no_prior:
        for nevents_data in data.nevents_tables:
            pretty_bfs = []
            for bf in nevents_data.bayes_factors:
                if bf is None:
                    pretty_bfs.append("")
                    continue
                pretty_bfs.append("{:.3g}".format(bf))
            for i, a in enumerate(nevents_data.bayes_factors_annotations):
                if a:
                    pretty_bfs[i] = a + pretty_bfs[i]
            nevents_data.pretty_bayes_factors = pretty_bfs

    plot_width = args.width
    plot_height = plot_width / 1.618034

    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib import gridspec
    except ImportError as e:
        raise Exception('ERROR: matplotlib could not be imported; '
                'it is needed for Python plotting.')

    add_legend = (not data.no_prior) and (not args.no_legend)

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

    nevents_indices = [float(x) for x in range(data.number_of_elements)]
    category_buffer = 0.08
    category_width = 1.0 - category_buffer 
    bar_width = category_width / (len(args.paths) + 1)
    if data.no_prior:
        bar_width = category_width / len(args.paths)

    prior_color = "0.65"

    bars_posterior = []
    for i in range(len(args.paths)):
        bars_post = ax.bar(
                [x + (bar_width * i) for x in nevents_indices],
                data.nevents_tables[i].posterior_probs,
                bar_width,
                color = args.colors[i],
                label = f'{args.labels[i]} posterior')
    if not data.no_prior:
        bars_prior = ax.bar(
                [x + (bar_width * len(args.paths)) for x in nevents_indices],
                data.prior_probs,
                bar_width,
                color = prior_color,
                label = args.prior_label)

    ax.set_xlabel(args.x_label)
    ax.set_ylabel(args.y_label)

    inner_legend_bump = 1.0
    if add_legend and args.legend_in_plot:
        inner_legend_bump = 1.08

    y_min, y_max = ax.get_ylim()
    y_max *= inner_legend_bump
    if args.full_prob_axis:
        y_min = 0.0
        if y_max < 1.0:
            y_max = 1.0
    ax.set_ylim(y_min, y_max)

    x_tick_labels = [str(i + 1) for i in range(data.number_of_elements)]

    # x values in nevents_indices (and x ticks) are at the center of the
    # first bar in each nevent category. Need to shift x ticks to the
    # center of each nevent category
    n_bars_per_nevent = len(args.paths)
    if not data.no_prior:
        n_bars_per_nevent += 1
    even_n_bars = (n_bars_per_nevent % 2) == 0
    middle_bar_idx = math.floor(n_bars_per_nevent / 2.0)
    x_tick_shift = middle_bar_idx * bar_width
    if even_n_bars:
        x_tick_shift -= (bar_width / 2.0)
    category_midpoints = [x + x_tick_shift for x in nevents_indices]
    plt.xticks(
            category_midpoints,
            x_tick_labels)

    partition_left_shift = (bar_width / 2.0) + (category_buffer / 2.0)
    category_partitions = [x - partition_left_shift for x in nevents_indices[1:]]
    for cat_part in category_partitions:
        ax.axvline(
            x = cat_part,
            color = '0.8',
            linestyle = '-',
            linewidth = 0.5)
    partition_span = category_partitions[1] - category_partitions[0]
    x_min = category_partitions[0] - partition_span
    x_max = category_partitions[-1] + partition_span
    ax.set_xlim(x_min, x_max)

    if not data.no_prior:
        bf_ax = ax.twinx()
        if args.log_bf:
            bf_ax.set_ylabel("Log$_{10}$ Bayes factor")
        else:
            bf_ax.set_ylabel("Bayes factor")
        for post_idx in range(len(args.paths)):
            # x values in nevents_indices are at center of first bar in each
            # nevent category. Centers of other bars are shifted bar_width to
            # the right
            bar_midpoints = [x + (bar_width * post_idx) for x in nevents_indices]
            bfs = data.nevents_tables[post_idx].bayes_factors
            if args.log_bf:
                bfs = [math.log10(f) for f in bfs]
            line, = bf_ax.plot(
                bar_midpoints,
                bfs,
                label = f'{args.labels[post_idx]} BF',
            )
            plt.setp(line,
                marker = 'd', # thin diamond
                markerfacecolor = args.colors[post_idx],
                markeredgecolor = args.bf_marker_border_color,
                markeredgewidth = 1.0,
                linestyle = '',
                rasterized = False,
                markersize = args.bf_marker_size, # mpl default = 6
                zorder = 500,
            )

        if add_legend and args.legend_in_plot:
            bf_y_min, bf_y_max = bf_ax.get_ylim()
            bf_y_max *= inner_legend_bump
            bf_ax.set_ylim(bf_y_min, bf_y_max)

    legend_handles = []
    legend_labels = []
    for current_ax in (ax, bf_ax):
        handles, labels = current_ax.get_legend_handles_labels()
        if handles:
            legend_handles.extend(handles)
        if labels:
            legend_labels.extend(labels)
        if current_ax.get_legend():
            current_ax.get_legend().remove()

    if add_legend:
        ncol = len(legend_labels)
        loc = 'lower center'
        if args.legend_in_plot:
            loc = 'upper center'
        leg = ax.legend(
            legend_handles,
            legend_labels,
            bbox_to_anchor = (0.5, 1.0),
            fontsize = args.legend_font_size, # mpl default = 10
            ncol = ncol,
            loc = loc)
        # Remove legend border, but keep the background
        leg.get_frame().set_linewidth(0.0)
    fig.tight_layout()
    plt.savefig(pdf_path)


if __name__ == "__main__":
    main()
