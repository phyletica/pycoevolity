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
        formatter_class = pycoevolity.argparse_utils.SmartDefaultsHelpFormatter,
    )

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
            'in these columns. '
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
        '--use-median',
        action = 'store_true',
        help = (
            'Use the posterior median when plotting. Default: Use posterior '
            'mean.'
        ),
    )
    parser.add_argument(
        '--use-eti',
        action = 'store_true',
        help = (
            'Use the equal-tailed credible intervals. Default: Use highest '
            'posterior density intervals.'
        ),
    )
    parser.add_argument(
        '--nevents-cred-level',
        type = pycoevolity.argparse_utils.arg_is_proportion,
        default = 0.95,
        help = (
            'Credibility level to use when plotting then number of events.'
        ),
    )
    parser.add_argument(
        '--psrf-max',
        type = pycoevolity.argparse_utils.arg_is_positive_float,
        default = 1.2,
        help = (
            'The maximum value for the potential scale reduction factor. '
            'Any parameters with a value greater than this for a simulation '
            'replicate will be highlighted to indicate poor MCMC mixing.'
            'Note, this argument is only used if multiple MCMC chains were '
            'run on each simulated data set.'
        ),
    )
    parser.add_argument(
        '--ess-min',
        type = pycoevolity.argparse_utils.arg_is_positive_float,
        default = 200,
        help = (
            'The minimum value for the effective sample size. '
            'Any parameters with a value less than this for a simulation '
            'replicate will be highlighted to indicate poor MCMC mixing.'
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
        default = 'notebook',
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
    parser.add_argument(
        '--output-dir',
        type = pycoevolity.argparse_utils.arg_is_dir_or_new_dir,
        help = (
            'The directory to which to write output files. By default, the '
            'directory of the results summary input file will be used.'
        ),
    )
    parser.add_argument(
        '--prefix',
        type = str,
        help = (
            'Prefix to add to the file name of every output file.'
        ),
    )
    parser.add_argument(
        '--plot-ext',
        type = str,
        default = 'svg',
        help = (
            'The file extension (and format) to use for output plotting files. '
            'Examples: '
            '\'--plot-ext svg\' (default), '
            '\'--plot-ext pdf\', '
            '\'--plot-ext png\', '
            '\'--plot-ext jpg\', etc. '
            'Any file formats supported by matplotlib should work.'
        ),
    )
    parser.add_argument(
        '--force',
        action = 'store_true',
        help = (
            'Overwrite output files if they already exist.'
        ),
    )
    parser.add_argument(
        '-f', '--parameter-file',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_file,
        help = (
            'RAW|'
            'Path to a YAML-formatted file that contains information about '
            'parameters to plot.\n'
            'The general format should be:\n'
            '    ---\n'
            '    parameter_name:\n'
            '        symbol: \'x\'\n'
            '        label: \"text to use for axis labels\"\n'
            'A specific example:\n'
            '    ---\n'
            '    time_prior_parameter_0:\n'
            '        symbol: \'\\mu\'\n'
            '        label: \"time prior mean\"\n'
        ),
    )
    parser.add_argument(
        '--exclude-column-titles',
        action = 'store_true',
        help = (
            'Exclude column titles from grid plots.'
        ),
    )
    parser.add_argument(
        '--exclude-row-titles',
        action = 'store_true',
        help = (
            'Exclude row titles from grid plots.'
        ),
    )
    args = parser.parse_args()
    return args

def write_existing_path_warning(path, out_stream = sys.stderr):
    out_stream.write(
        f"\n"
        f"WARNING: Skipping plotting of '{path}' because it "
        f"already exists. If you wish to overwrite this file, please "
        f"use the '--force' option.\n\n"
    )

def parse_parameter_yaml(path):
    d = pycoevolity.fileio.load_yaml(path)
    ret = {k : {} for k in d}
    for key in d:
        for sub_key in ('label', 'symbol'):
            ret[key][sub_key] = d[key].get(sub_key, None)
    return ret

def main_cli():
    args = parse_cli_args()

    sns.set_theme(
        context = args.context,
        style = args.style,
        palette = args.palette,
        font = args.font,
        font_scale = args.font_scale,
    )

    parameters_to_plot = {}
    if args.parameter_file:
        parameters_to_plot = parse_parameter_yaml(args.parameter_file)

    df = pd.read_csv(
        args.results_summary_path,
        sep = "\t",
    )

    # This copy is needed to avoid Panda's "PerformanceWarning: DataFrame is
    # highly fragmented"
    df = df.copy()
    df['abs_num_events_error'] = (df['true_num_events'] - df['map_num_events']).abs()

    cred_interval_percent = pycoevolity.plotting.get_cred_interval_percent(
        df.columns,
    )
    ci_str = str(cred_interval_percent)
    ci_prefix = f"hpdi_{ci_str}"
    if args.use_eti:
        ci_prefix = f"eti_{ci_str}"
    cred_level = cred_interval_percent / 100.0

    ordered_labels = None
    if args.config_label_order:
        ordered_labels = parse_config_label_order_arg(
            args.config_label_order,
            args.label_sep,
        )
    plot_dir = os.path.dirname(args.results_summary_path)
    if args.output_dir:
        plot_dir = pycoevolity.argparse_utils.process_output_dir_arg(args.output_dir)

    prefix = ''
    if args.prefix:
        prefix = args.prefix
    plot_prefix = os.path.join(
        plot_dir,
        prefix,
    )

    col_title_template = "{col_name}"
    row_title_template = "{row_name}"

    if args.exclude_column_titles:
        col_title_template = ""
    if args.exclude_row_titles:
        row_title_template = ""

    ##################################################################
    # Plot number of events heat maps
    ##################################################################
    plot_path = f"{plot_prefix}nevents-heatmap-grid.{args.plot_ext}"

    if (not args.force) and os.path.exists(plot_path):
        write_existing_path_warning(plot_path, sys.stderr)
    else:
        grid = pycoevolity.plotting.plot_nevents_heatmap_grid(
            df,
            row_col = "simulation_config",
            column_col = "inference_config",
            ordered_labels = ordered_labels,
            height = args.plot_height,
            annotate_counts = False,
            include_cbar = True,
            outline_identity = True,
            annotate_stats = True,
            cred_level = args.nevents_cred_level,
            col_title_template = col_title_template,
            row_title_template = row_title_template,
        )
        grid.savefig(plot_path)

    ##################################################################
    # Plot number of events error scatter
    ##################################################################
    plot_path = f"{plot_prefix}nevents-error-scatter-grid.{args.plot_ext}"

    if (not args.force) and os.path.exists(plot_path):
        write_existing_path_warning(plot_path, sys.stderr)
    else:
        grid = pycoevolity.plotting.plot_error_scatter_grid(
            df,
            true_val_col = "true_num_events",
            est_col = "map_num_events",
            row_col = "simulation_config",
            column_col = "inference_config",
            id_col = "simulation_id",
            est_lower_col = f"hpdi_{ci_str}_lower_num_events",
            est_upper_col = f"hpdi_{ci_str}_upper_num_events",
            ess_col = None,
            psrf_max = None,
            bad_sampling_color = "C1",
            ordered_labels = ordered_labels,
            annotate_true_values = True,
            annotate_stats = False,
            cred_level = cred_level,
            height = args.plot_height,
            scatter_kwargs = {},
            annot_true_vals_kwargs = {},
            col_title_template = col_title_template,
            row_title_template = row_title_template,
        )
        grid.set_axis_labels(
            "True number of events",
            "Number of events error",
        )
        grid.savefig(plot_path)

    ##################################################################
    # Violin plots of model performance
    ##################################################################
    violin_plot_stats = {
        "true_model_p" :
        {
            "label" : "True model posterior probability",
            "plot_path" : f"{plot_prefix}true-model-prob-violin-grid.{args.plot_ext}",
        },
        "true_model_cred_level" :
        {
            "label" : "True model posterior rank sum",
            "plot_path" : f"{plot_prefix}true-model-rank-sum-violin-grid.{args.plot_ext}",
        },
        "abs_num_events_error" :
        {
            "label" : "MAP $k$ absolute error",
            "plot_path" : f"{plot_prefix}map-nevents-abs-error-violin-grid.{args.plot_ext}",
        },
        "true_num_events_p" :
        {
            "label" : "True $k$ posterior probability",
            "plot_path" : f"{plot_prefix}true-nevents-prob-violin-grid.{args.plot_ext}",
        },
        "true_num_events_cred_level" :
        {
            "label" : "True $k$ posterior rank sum",
            "plot_path" : f"{plot_prefix}true-nevents-rank-sum-violin-grid.{args.plot_ext}",
        },
    }
    if args.use_median:
        violin_plot_stats["median_map_model_distance"] = {
            "label" : "MAP model distance",
            "plot_path" : f"{plot_prefix}map-model-distance-violin-grid.{args.plot_ext}",
        }
    else:
        violin_plot_stats["mean_map_model_distance"] = {
            "label" : "MAP model distance",
            "plot_path" : f"{plot_prefix}map-model-distance-violin-grid.{args.plot_ext}",
        }

    for stat_key, plotting_args in violin_plot_stats.items():
        plot_path = plotting_args['plot_path']
        stat_label = plotting_args['label']

        if (not args.force) and os.path.exists(plot_path):
            write_existing_path_warning(plot_path, sys.stderr)
        else:
            grid = pycoevolity.plotting.plot_violin_grid(
                data = df,
                value_col = stat_key,
                plot_col = "simulation_config",
                categorical_col = "inference_config",
                spaghettify_col = "simulation_id",
                spaghettify = True,
                value_label = stat_label,
                categorical_label = "Inference model",
                plot_label_template = "True model = {col_name}",
                ordered_labels = ordered_labels,
                comparisons = args.comparison,
                height = args.violin_plot_height,
                categorical_label_size = args.violin_label_size,
                violin_kwargs = {},
                spaghetti_kwargs = {},
            )
            grid.savefig(plot_path)

    ##################################################################
    # Parameter scatter plots
    ##################################################################

    scatter_params = {
        "concentration" : {
            "symbol" : r"\alpha",
            "label" : "concentration",
        },
        "root_height" : {
            "symbol" : r"\tau",
            "label" : "divergence time",
        },
        "pop_size_root" : {
            "symbol" : r"N_e",
            "label" : "ancestral $N_e$",
        },
        "time_prior_parameter_0" : {
            "symbol" : None,
            "label" : "time prior parameter 0",
        },
        "time_prior_parameter_1" : {
            "symbol" : None,
            "label" : "time prior parameter 1",
        },
    }
    scatter_params.update(parameters_to_plot)

    for param_key, param_info in scatter_params.items():
        parameter_root = param_key
        parameters = pycoevolity.plotting.get_all_comparison_parameters(
            parameter_prefix = param_key,
            column_headers = df.columns,
        )
        if not parameters:
            # param_key was not a prefix for comparison parameters, so check to
            # see if it's a global parameter
            if f"mean_{param_key}" in df.columns:
                parameters = [param_key]
                parameter_root = None
            else:
                # This parameter is not in the data frame; we don't throw an
                # error, because some parameters don't end up in the summary
                # table if they were constrained/fixed
                continue
        xlabel = None
        ylabel = None
        stat_label = None
        if param_info['label']:
            xlabel = f"True {param_info['label']}"
            ylabel = f"Mean {param_info['label']}"
        if param_info['symbol']:
            stat_label = param_info['symbol']
        plot_path = f"{plot_prefix}{param_key}-scatter-grid.{args.plot_ext}"
        if (not args.force) and os.path.exists(plot_path):
            write_existing_path_warning(plot_path, sys.stderr)
        else:
            grid = pycoevolity.plotting.process_scatter_grid(
                df,
                parameters = parameters,
                row_col = "simulation_config",
                column_col = "inference_config",
                parameter_root = parameter_root,
                use_mean = (not args.use_median),
                use_hpdi = (not args.use_eti),
                xlabel = xlabel,
                ylabel = ylabel,
                ess_min = args.ess_min,
                psrf_max = args.psrf_max,
                bad_sampling_color = "C1",
                ordered_labels = ordered_labels,
                height = args.plot_height,
                annotate_stats = True,
                stat_label = stat_label,
                annot_position = (0.02, 0.98),
                cred_percent = cred_interval_percent,
                scatter_kwargs = {},
                annotate_kwargs = {},
                col_title_template = col_title_template,
                row_title_template = row_title_template,
            )
            if grid:
                grid.savefig(plot_path)

        plot_path = f"{plot_prefix}{param_key}-error-scatter-grid.{args.plot_ext}"
        if (not args.force) and os.path.exists(plot_path):
            write_existing_path_warning(plot_path, sys.stderr)
        else:
            grid = pycoevolity.plotting.process_error_scatter_grid(
                df,
                parameters = parameters,
                row_col = "simulation_config",
                column_col = "inference_config",
                id_col = "simulation_id",
                parameter_root = parameter_root,
                use_mean = (not args.use_median),
                use_hpdi = (not args.use_eti),
                ess_min = args.ess_min,
                psrf_max = args.psrf_max,
                bad_sampling_color = "C1",
                ordered_labels = ordered_labels,
                annotate_true_values = False,
                annotate_stats = True,
                stat_label = stat_label,
                annot_stats_position = (0.02, 0.98),
                cred_percent = cred_interval_percent,
                height = 4.5,
                scatter_kwargs = {},
                annot_true_vals_kwargs = {},
                annot_stats_kwargs = {},
                col_title_template = col_title_template,
                row_title_template = row_title_template,
            )
            if grid:
                if param_info['label']:
                    label_list = param_info['label'].split()
                    label_list[0] = label_list[0].capitalize()
                    label_list.append("error")
                    ylabel = " ".join(label_list)
                    grid.set_ylabels(ylabel)
                grid.savefig(plot_path)


if __name__ == "__main__":
    main_cli()
