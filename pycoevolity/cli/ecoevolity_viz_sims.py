#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import seaborn as sns

import pycoevolity


def parse_config_label_order_arg(arg, sep):
    labels = [x.strip() for x in arg.split(sep)]
    return tuple(labels)

def parse_cli_args():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'results_summary_path',
        metavar = 'RESULTS-SUMMARY-TSV-FILE',
        type = pycoevolity.argparse_utils.arg_is_file,
        help = (
            'Path to results-summary.tsv.gz file created by pyco-eco-sum-sims.'
        ),
    )
    parser.add_argument(
        '-o', '--config-label-order',
        action = 'store',
        type = str,
        help = (
            'An optional list of config labels that determines the order of '
            'rows and columns of grids of plots. By default, the config labels '
            'should be separated by semi-colons (\';\'). '
            'If your config labels contain semi-colons, you can specify a '
            'different separator using the '
            '\'-s\'/\'--label-sep\' argument. '
            'If you used the \'-l\'/\'--config-label-file\' option with '
            'pyco-eco-sum-sims, the labels you provide with this argument '
            'should match the labels in the config label file. '
            'Otherwise the labels should match the name '
            'of the config files (NOT the paths). '
            'If you are not sure, you can check the values in the '
            '\'simulation_config\' and \'inference_config\' columns of the '
            'results-summary.tsv.gz file; your labels should match the labels '
            'in these columns.'
            'If your labels, or separations between them, include spaces, you '
            'will need to quote the argument; e.g., '
            '-o \'Independent model; Shared model; Uniform model; DPP model\'.'
        ),
    )
    parser.add_argument(
        '-s', '--label-sep',
        action = 'store',
        type = str,
        default = ';',
        help = (
            'The character used to separate config labels in the '
            '\'-o\'/\'--config-label-order\' argument.'
        ),
    )

    args = parser.parse_args()
    return args

def main_cli():
    args = parse_cli_args()

    df = pd.read_csv(
        args.results_summary_path,
        sep = "\t",
    )

    ordered_labels = None
    if args.config_label_order:
        ordered_labels = parse_config_label_order_arg(
            args.config_label_order,
            args.label_sep,
        )
    plot_dir = os.path.dirname(args.results_summary_path)

    grid = pycoevolity.plotting.plot_nevents_heatmap_grid(
        df,
        sim_config_col = "simulation_config",
        inference_config_col = "inference_config",
        ordered_labels = ordered_labels,
        height = 6.5,
        annotate_counts = False,
        include_cbar = True,
        outline_identity = True,
        annotate_stats = True,
        cred_level = 0.95,
    )
    plot_path = os.path.join(
        plot_dir,
        "nevents-heatmap-grid.pdf",
    )
    grid.savefig(plot_path)

    grid = pycoevolity.plotting.plot_abs_error_grid(
        df,
        error_col = "map_num_events_distance",
        error_lower_col = "hpdi_95_lower_num_events_distance",
        error_upper_col = "hpdi_95_upper_num_events_distance",
        true_val_col = "true_num_events",
        sim_config_col = "simulation_config",
        inference_config_col = "inference_config",
        ordered_labels = ordered_labels,
        height = 4.5,
    )
    grid.set_axis_labels(
        "Simulation replicate",
        "Number of events error",
    )
    plot_path = os.path.join(
        plot_dir,
        "nevents-abs-error-grid.pdf",
    )
    grid.savefig(plot_path)

    grid = pycoevolity.plotting.process_error_scatter_grid(
        df,
        parameters = ["concentration"],
        parameter_root = None,
        use_mean = True,
        use_hpdi = True,
        xlabel = "True concentration",
        ylabel = "Mean concentration",
        ess_min = 200,
        psrf_max = 1.2,
        bad_sampling_color = "C1",
        ordered_labels = ordered_labels,
        height = 6.5,
        annotate_stats = True,
        stat_label = r"\alpha",
        annot_x_position = 0.02,
        annot_y_position = 0.98,
        cred_level = 0.95,
    )
    plot_path = os.path.join(
        plot_dir,
        "concentration-scatter-grid.pdf",
    )
    grid.savefig(plot_path)

    parameter_root = "root_height"
    parameters = pycoevolity.plotting.get_all_parameters(
        parameter_prefix = parameter_root,
        column_headers = df.columns,
    )
    grid = pycoevolity.plotting.process_error_scatter_grid(
        df,
        parameters = parameters,
        parameter_root = parameter_root,
        use_mean = True,
        use_hpdi = True,
        xlabel = "True divergence time",
        ylabel = "Mean divergence time",
        ess_min = 200,
        psrf_max = 1.2,
        bad_sampling_color = "C1",
        ordered_labels = ordered_labels,
        height = 6.5,
        annotate_stats = True,
        stat_label = r"\tau",
        annot_x_position = 0.02,
        annot_y_position = 0.98,
        cred_level = 0.95,
    )
    plot_path = os.path.join(
        plot_dir,
        "div-time-scatter-grid.pdf",
    )
    grid.savefig(plot_path)

if __name__ == "__main__":
    sns.set_theme(context = "talk", style = "ticks", palette = "colorblind")
    main_cli()
