#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
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
        metavar = "\'CONFIG-LABEL-1;CONFIG-LABEL-2;...\'",
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
        '-c', '--comparison',
        action = 'append',
        type = str,
        nargs = 2,
        metavar = ("INFERENCE-CONFIG-LABEL-1", "INFERENCE-CONFIG-LABEL-2"),
        help = (
            'This option takes two arguments: '
            'Two inference config labels. '
            'If provided, statistical tests comparing the results under these '
            'two configs will be performed, and the results will be annotated '
            'on some of the plots. '
            'If you used the \'-l\'/\'--config-label-file\' option with '
            'pyco-eco-sum-sims, the labels you provide with this argument '
            'should match the labels in the config label file. '
            'Otherwise the labels should match the name '
            'of the config files (NOT the paths). '
            'If you are not sure, you can check the values in the '
            '\'simulation_config\' and \'inference_config\' columns of the '
            'results-summary.tsv.gz file; your labels should match the labels '
            'in these columns.'
            'Note, you need to quote config labels that have spaces. '
            'This option can be used multiple times to specify multiple '
            'statistical comparisons.'
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
    parser.add_argument(
        '-p', '--plot-height',
        type = pycoevolity.argparse_utils.arg_is_positive_float,
        default = 4.5,
        help = (
            'The height of each plot. Adjusting this is useful for changing '
            'the relative size of text on the plots (e.g., increase the plot '
            'height to decrease the size of the text.'
        ),
    )
    parser.add_argument(
        '--use-median-model-distance',
        action = 'store_true',
        help = (
            'Use the posterior median model distance. Default: Use posterior '
            'mean model distance.'
        ),
    )
    parser.add_argument(
        '-v', '--violin-plot-height',
        type = pycoevolity.argparse_utils.arg_is_positive_float,
        default = 5.0,
        help = (
            'The height of each violin plot. Adjusting this is useful for '
            'changing the relative size of text on violin plots (e.g., '
            'increase the plot height to decrease the size of the text.'
        ),
    )
    parser.add_argument(
        '--violin-label-size',
        type = pycoevolity.argparse_utils.arg_is_positive_float,
        help = (
            'Font size for inference config labels on violin plots.'
        ),
    )
    parser.add_argument(
        '--context',
        type = str,
        default = 'talk',
        help = (
            'Value of \'context\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    parser.add_argument(
        '--style',
        type = str,
        default = 'ticks',
        help = (
            'Value of \'style\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    parser.add_argument(
        '--palette',
        type = str,
        default = 'colorblind',
        help = (
            'Value of \'palette\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    parser.add_argument(
        '--font',
        type = str,
        default = 'sans-serif',
        help = (
            'Value of \'font\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    parser.add_argument(
        '--font-scale',
        type = pycoevolity.argparse_utils.arg_is_positive_float,
        default = 1.0,
        help = (
            'Value of \'font-scale\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    args = parser.parse_args()
    return args

def main_cli(args = None):
    if args is None:
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

    ##################################################################
    # Plot number of events heat maps
    ##################################################################
    grid = pycoevolity.plotting.plot_nevents_heatmap_grid(
        df,
        sim_config_col = "simulation_config",
        inference_config_col = "inference_config",
        ordered_labels = ordered_labels,
        height = args.plot_height,
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

    ##################################################################
    # Plot absolute number of events error
    ##################################################################

    grid = pycoevolity.plotting.plot_abs_error_grid(
        df,
        true_val_col = "true_num_events",
        est_col = "map_num_events",
        est_lower_col = "hpdi_95_lower_num_events",
        est_upper_col = "hpdi_95_upper_num_events",
        sim_config_col = "simulation_config",
        inference_config_col = "inference_config",
        ordered_labels = ordered_labels,
        annotate_true_values = True,
        height = args.plot_height,
        scatter_kwargs = {},
        annotate_kwargs = {},
    )
    grid.set_axis_labels(
        "True number of events",
        "Number of events error",
    )
    plot_path = os.path.join(
        plot_dir,
        "nevents-abs-error-grid.pdf",
    )
    grid.savefig(plot_path)

    ##################################################################
    # Parameter scatter plots
    ##################################################################

    # scatter_params = (
    #     (
    #         "concentration",
    #         r"\alpha",
    #         "concentration",
    #     ),
    #     (
    #         "root_height",

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
        height = args.plot_height,
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
        height = args.plot_height,
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

    model_dist = "mean_model_distance"
    model_dist_label = "Mean model distance"
    if args.use_median_model_distance:
        model_dist = "median_model_distance"
        model_dist_label = "Median model distance"
    grid = pycoevolity.plotting.plot_violin_grid(
        data = df,
        value_col = model_dist,
        sim_config_col = "simulation_config",
        inference_config_col = "inference_config",
        spaghettify_col = "simulation_id",
        spaghettify = True,
        value_label = model_dist_label,
        inference_label = "Inference model",
        sim_label_template = "True model = {col_name}",
        ordered_labels = ordered_labels,
        comparisons = args.comparison,
        height = args.violin_plot_height,
        inference_label_size = args.violin_label_size,
        violin_kwargs = {},
        spaghetti_kwargs = {},
    )
    plot_path = os.path.join(
        plot_dir,
        "model-error-grid.pdf",
    )
    grid.savefig(plot_path)


if __name__ == "__main__":
    # Calling sns.set_theme within main_cli has wonky results, so we need to
    # parse args here so we can call sns.set_theme outside of the main_cli
    # function.
    args = parse_cli_args()

    sns.set_theme(
        context = args.context,
        style = args.style,
        palette = args.palette,
        font = args.font,
        font_scale = args.font_scale,
    )
    main_cli(args)
