#!/usr/bin/env python

import os
import sys
import math
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd

import pycoevolity


def get_errors(values, lowers, uppers):
    values = tuple(values)
    lowers = tuple(lowers)
    uppers = tuple(uppers)
    n = len(values)
    assert(n == len(lowers))
    assert(n == len(uppers))
    return [[values[i] - lowers[i] for i in range(n)],
            [uppers[i] - values[i] for i in range(n)]]

def get_true_v_map_nevents_data_frame(
    data_frame,
    sim_config_col = "simulation_config",
    inference_config_col = "inference_config",
):
    num_comparisons = len(data_frame["map_model"][0].split(","))
    nevent_labels = tuple(range(1, num_comparisons + 1))
    sim_configs = data_frame[sim_config_col].unique()
    inf_configs = data_frame[inference_config_col].unique()
    rows = []
    for sim_conf in sim_configs:
        for inf_conf in inf_configs:
            counts = [ [0 for i in range(num_comparisons)] for j in range(num_comparisons) ]
            sub_df = data_frame.loc[
                (data_frame[sim_config_col] == sim_conf)
                & (data_frame[inference_config_col] == inf_conf)
            ]
            for true_n, map_n in zip(
                sub_df["true_num_events"],
                sub_df["map_num_events"],
            ):
                counts[true_n - 1][map_n - 1] += 1
            for true_idx in range(num_comparisons):
                for map_idx in range(num_comparisons):
                    count = counts[true_idx][map_idx]
                    rows.append({
                        'simulation_config' : sim_conf,
                        'inference_config' : inf_conf,
                        'true_num_events' : true_idx + 1,
                        'map_num_events' : map_idx + 1,
                        'count' : count,
                    })
    return pd.DataFrame(rows)

def get_all_population_labels(column_headers):
    pop_labels = set()
    for col in column_headers:
        if col.startswith("mean_root_height_"):
            pop_lab = col[len("mean_root_height_"):]
            pop_labels.add(pop_lab)
        elif col.startswith("mean_pop_size_root_"):
            pop_lab = col[len("mean_pop_size_root_"):]
            pop_labels.add(pop_lab)
        elif col.startswith("mean_pop_size_"):
            pop_lab = col[len("mean_pop_size_"):]
            pop_labels.add(pop_lab)
    return tuple(sorted(pop_labels))

def get_comparison_labels(column_headers):
    comp_labels = []
    for col in column_headers:
        if col.startswith("mean_root_height_"):
            comp_lab = col[len("mean_root_height_"):]
            comp_labels.append(comp_lab)
    assert len(comp_labels) == len(set(comp_labels))
    return tuple(comp_labels)

def get_all_parameters(parameter_prefix, column_headers):
    pop_labels = get_all_population_labels(column_headers)
    prefix = parameter_prefix.rstrip("_")
    possible_params = [f"{prefix}_{l}" for l in pop_labels]
    params = [p for p in possible_params if f"mean_{p}" in column_headers]
    return tuple(params)

def get_stacked_parameter_data_frame(data_frame, parameters, parameter_root):
    param_root = parameter_root.rstrip("_")
    param_keys = [
        f"psrf_{param_root}",
        f"true_{param_root}",
        f"true_{param_root}_rank",
        f"mean_{param_root}",
        f"median_{param_root}",
        f"stddev_{param_root}",
        f"hpdi_95_lower_{param_root}",
        f"hpdi_95_upper_{param_root}",
        f"eti_95_lower_{param_root}",
        f"eti_95_upper_{param_root}",
        f"ess_{param_root}",
        f"ess_sum_{param_root}",
    ]
    config_keys = [
        "simulation_config",
        "inference_config",
    ]
    columns = {k : [] for k in param_keys + config_keys}
    for param in parameters:
        for k in param_keys:
            param_label = k.replace(param_root, param)
            columns[k].extend(data_frame[param_label])
        for conf_key in config_keys:
            columns[conf_key].extend(data_frame[conf_key])
    n = None
    for k, vals in columns.items():
        if n is None:
            n = len(vals)
        else:
            assert len(vals) == n
    return pd.DataFrame(columns)

def process_error_scatter_grid(
    data_frame,
    parameters,
    parameter_root = None,
    use_mean = True,
    use_hpdi = True,
    xlabel = None,
    ylabel = None,
    ess_min = 200,
    psrf_max = 1.2,
    bad_sampling_color = "C1",
    ordered_labels = None,
    height = 6.5,
    annotate_stats = True,
    stat_label = None,
    annot_x_position = 0.02,
    annot_y_position = 0.98,
    cred_level = 0.95,
    **kwargs,
):
    if not parameters:
        raise Exception(
            "parameters are empty"
        )
    elif len(parameters) > 1:
        if not parameter_root:
            raise Exception(
                "parameter_root is required when processing multiple parameters"
            )
        df = get_stacked_parameter_data_frame(data_frame, parameters, parameter_root)
        parameter = parameter_root
    else:
        df = data_frame
        parameter = parameters[0]
    true_col = f"true_{parameter}"
    true_val_rank_col = f"true_{parameter}_rank"
    est_col = f"median_{parameter}"
    if use_mean:
        est_col = f"mean_{parameter}"
    est_error_lower_col = f"eti_95_lower_{parameter}"
    est_error_upper_col = f"eti_95_upper_{parameter}"
    if use_hpdi:
        est_error_lower_col = f"hpdi_95_lower_{parameter}"
        est_error_upper_col = f"hpdi_95_upper_{parameter}"
    ess_col = f"ess_{parameter}"
    psrf_col = f"psrf_{parameter}"
    grid = plot_error_scatter_grid(
        data_frame = df,
        true_col = true_col,
        est_col = est_col,
        est_error_lower_col = est_error_lower_col,
        est_error_upper_col = est_error_upper_col,
        sim_config_col = "simulation_config",
        inference_config_col = "inference_config",
        true_val_rank_col = true_val_rank_col,
        xlabel = xlabel,
        ylabel = ylabel,
        ess_col = ess_col,
        psrf_col = ess_col,
        ess_min = ess_min,
        psrf_max = psrf_max,
        bad_sampling_color = bad_sampling_color,
        ordered_labels = ordered_labels,
        height = height,
        annotate_stats = annotate_stats,
        stat_label = stat_label,
        annot_x_position = annot_x_position,
        annot_y_position = annot_y_position,
        cred_level = cred_level,
        **kwargs,
    )
    return grid

def plot_error_scatter_grid(
    data_frame,
    true_col,
    est_col,
    est_error_lower_col,
    est_error_upper_col,
    sim_config_col = "simulation_config",
    inference_config_col = "inference_config",
    true_val_rank_col = None,
    xlabel = None,
    ylabel = None,
    ess_col = None,
    psrf_col = None,
    ess_min = 200,
    psrf_max = 1.2,
    bad_sampling_color = "C1",
    ordered_labels = None,
    height = 6.5,
    annotate_stats = True,
    stat_label = None,
    annot_x_position = 0.02,
    annot_y_position = 0.98,
    cred_level = 0.95,
    **kwargs,
):
    col_order = None
    row_order = None
    if ordered_labels:
        sim_labels = data_frame[sim_config_col].unique()
        inf_labels = data_frame[inference_config_col].unique()
        row_order = [l for l in ordered_labels if l in sim_labels]
        col_order = [l for l in ordered_labels if l in inf_labels]
    grid = sns.FacetGrid(
        data_frame,
        row = sim_config_col,
        col = inference_config_col,
        margin_titles = True,
        height = height,
        row_order = row_order,
        col_order = col_order,
        sharey = True,
        sharex = True,
    )
    grid.map_dataframe(
        plot_error_scatter,
        x = true_col,
        y = est_col,
        y_error_lower = est_error_lower_col,
        y_error_upper = est_error_upper_col,
        ess_col = ess_col,
        psrf_col = psrf_col,
        ess_min = ess_min,
        psrf_max = psrf_max,
        **kwargs,
    )
    mn = min(min(data_frame[true_col]), min(data_frame[est_col]))
    mx = max(max(data_frame[true_col]), max(data_frame[est_col]))
    for ax in grid.axes_dict.values():
        ax_id_line(
            ax = ax,
            mn = mn, 
            mx = mx,
            color = "0.8",
            linestyle = "-",
            linewidth = 1.0,
        )
    if annotate_stats:
        grid.map_dataframe(
            annotate_error_scatter,
            x = true_col,
            y = est_col,
            y_error_lower = est_error_lower_col,
            y_error_upper = est_error_upper_col,
            x_position = annot_x_position,
            y_position = annot_y_position,
            cred_level = cred_level,
            stat_label = stat_label,
            **kwargs,
        )
    if xlabel:
        grid.set_xlabels(xlabel)
    if ylabel:
        grid.set_ylabels(ylabel)
    grid.set_titles(
        col_template = "{col_name}",
        row_template = "{row_name}",
    )
    return grid

def annotate_error_scatter(
    data,
    x,
    y,
    y_error_lower,
    y_error_upper,
    x_position = 0.02,
    y_position = 0.98,
    cred_level = 0.95,
    stat_label = None,
    **kwargs,
):
    ax = plt.gca()
    num_within_ci = (
        (data[x] >= data[y_error_lower])
        & (data[x] <= data[y_error_upper])
    ).sum()
    prop_within_ci = num_within_ci / len(data[x])
    prop_est_under = (
        sum(data[x] > data[y])
        / len(data[x])
    )
    sum_sq_err = ((data[x] - data[y]) ** 2).sum()
    mean_sq_err = sum_sq_err / len(data)
    root_mean_sq_err = math.sqrt(mean_sq_err)
    if not stat_label:
        stat_label = "x"
    annot_str = (
        r"$p(\hat{{{stat_label}}} < {stat_label}) = {prop_under:.2g}$"
        "\n"
        r"$p({stat_label} \in {cred_level:.2f}\,\text{{CI}}) = {coverage:.2g}$"
        "\n"
        r"$\text{{RMSE}} = {rmse:.2g}$".format(
            stat_label = stat_label,
            prop_under = prop_est_under,
            cred_level = cred_level,
            coverage = prop_within_ci,
            rmse = root_mean_sq_err,
        )
    )
    default_args = {
        'horizontalalignment' : "left",
        'verticalalignment' : "top",
        'transform' : ax.transAxes,
        'zorder' : 200,
        'fontsize' : 'small',
        'bbox' : {
            'facecolor': 'white',
            # 'edgecolor': 'white',
            'pad': 2,
            'alpha': 0.5,
        },
    }
    default_args.update(kwargs)
    # seaborn passes in color keyword arg set to the same color used for
    # plotting points; overriding that here
    default_args["color"] = "black"
    ax.text(
        x_position, y_position,
        annot_str,
        **default_args,
    )

def plot_error_scatter(
    data,
    x,
    y,
    y_error_lower,
    y_error_upper,
    ess_col = None,
    psrf_col = None,
    ess_min = 200,
    psrf_max = 1.2,
    bad_sampling_color = "C1",
    **kwargs,
):
    df = data.copy()
    if ess_col and psrf_col:
        df["Poor MCMC sampling"] = (
            (df[ess_col] < ess_min)
            & (df[psrf_col] > psrf_max)
        )
    elif ess_col:
        df["Poor MCMC sampling"] = df[ess_col] < ess_min
    elif psrf_col:
        df["Poor MCMC sampling"] = df[psrf_col] > psrf_max
    ax = plt.gca()
    shared_args = {
        'elinewidth' : 1.0,
        'capsize' : 1.5,
        'barsabove' : False,
        'marker' : 'o',
        'linestyle' : '',
        'markeredgewidth' : 0.0,
        'markersize' : 6.5,
        'rasterized' : False,
        'alpha' : 0.5,
    }
    shared_args.update(kwargs)
    line = ax.errorbar(
        x = df[x],
        y = df[y],
        yerr = get_errors(df[y], df[y_error_lower], df[y_error_upper]),
        ecolor = 'C0',
        markerfacecolor = 'C0',
        markeredgecolor = 'C0',
        zorder = 100,
        **shared_args,
    )
    if "Poor MCMC sampling" in df.columns:
        d = df[df["Poor MCMC sampling"]].copy()
        if len(d) > 0:
            bad_line = ax.errorbar(
                x = d[x],
                y = d[y],
                yerr = get_errors(d[y], d[y_error_lower], d[y_error_upper]),
                ecolor = bad_sampling_color,
                markerfacecolor = bad_sampling_color,
                markeredgecolor = bad_sampling_color,
                zorder = 200,
                **shared_args,
            )

def plot_nevents_heatmap_grid(
    data_frame,
    sim_config_col = "simulation_config",
    inference_config_col = "inference_config",
    ordered_labels = None,
    height = 6.5,
    annotate_counts = True,
    include_cbar = True,
    outline_identity = True,
    annotate_stats = True,
    cred_level = 0.95,
):
    data = get_true_v_map_nevents_data_frame(data_frame)
    vmin = min(data["count"])
    vmax = max(data["count"])
    col_order = None
    row_order = None
    if ordered_labels:
        sim_labels = data_frame[sim_config_col].unique()
        inf_labels = data_frame[inference_config_col].unique()
        row_order = [l for l in ordered_labels if l in sim_labels]
        col_order = [l for l in ordered_labels if l in inf_labels]
    grid = sns.FacetGrid(
        data,
        row = sim_config_col,
        col = inference_config_col,
        margin_titles = True,
        height = height,
        row_order = row_order,
        col_order = col_order,
        sharey = True,
        sharex = True,
    )
    cbar = False
    cbar_ax = None
    if include_cbar:
        cbar = True
        cbar_ax = grid.fig.add_axes([0.92, 0.2, 0.02, 0.6])
    grid.map_dataframe(
        plot_heatmap,
        index = 'map_num_events',
        columns = 'true_num_events',
        values = 'count',
        annot = annotate_counts,
        # cmap = sns.cubehelix_palette(as_cmap = True),
        cmap = sns.color_palette("Blues", as_cmap=True),
        cbar = cbar,
        cbar_ax = cbar_ax,
        vmin = vmin,
        vmax = vmax,
        outline_identity = outline_identity,
    )
    if annotate_stats:
        grid.map_dataframe(
            annotate_heatmap,
            full_data = data_frame,
            sim_config_col = sim_config_col,
            inference_config_col = inference_config_col,
            true_val_col = "true_num_events",
            est_val_col = "map_num_events",
            true_val_cred_col = 'true_num_events_cred_level',
            true_val_prob_col = 'true_num_events_p',
            cred_level = cred_level,
        )
    grid.set_axis_labels(
        "True number of events",
        "MAP number of events")
    grid.set_titles(
        col_template = "{col_name}",
        row_template = "{row_name}",
    )
    if include_cbar:
        grid.fig.tight_layout(rect = [0, 0, 0.92, 1])
    # grid.figure.subplots_adjust(wspace = 0.05, hspace = 0.05)
    return grid

def annotate_heatmap(
    data,
    full_data,
    sim_config_col,
    inference_config_col,
    true_val_col,
    est_val_col,
    true_val_cred_col,
    true_val_prob_col,
    cred_level = 0.95,
    x_position = 0.02,
    y_position = 0.98,
    **kwargs,
):
    ax = plt.gca()
    sim_configs = data[sim_config_col].unique()
    assert len(sim_configs) == 1
    sim_config = sim_configs[0]
    inf_configs = data[inference_config_col].unique()
    assert len(inf_configs) == 1
    inf_config = inf_configs[0]
    df = full_data.loc[
        (full_data[sim_config_col] == sim_config)
        & (full_data[inference_config_col] == inf_config)
    ]
    num_within_cs = (
        (df[true_val_cred_col] <= cred_level)
    ).sum()
    prop_within_cs = num_within_cs / len(df[true_val_cred_col])
    prop_map_under = (
        sum(df[true_val_col] > df[est_val_col])
        / len(df[true_val_col])
    )
    median_prob_true = np.median(df[true_val_prob_col])
    annot_str = (
        r"$p(\hat{{k}} < k) = {prop_under:.2g}$"
        "\n"
        r"$p(k \in {cred_level:.2f}\,\text{{CS}}) = {coverage:.2g}$"
        "\n"
        r"median $p(k|D) = {med_prob:.2g}$".format(
            cred_level = cred_level,
            coverage = prop_within_cs,
            prop_under = prop_map_under,
            med_prob = median_prob_true,
        )
    )
    ax.text(
        x_position, y_position,
        annot_str,
        horizontalalignment = "left",
        verticalalignment = "top",
        transform = ax.transAxes,
        zorder = 200,
        fontsize = 'small',
        # bbox = {
        #     'facecolor': 'white',
        #     'edgecolor': 'white',
        #     'pad': 2},
    )

def plot_heatmap(
    data,
    index,
    columns,
    values,
    **kwargs):
    outline_identity = kwargs.pop('outline_identity', True)
    d = data.pivot(index = index, columns = columns, values = values)
    ax = sns.heatmap(d, **kwargs)
    ax.invert_yaxis()
    if outline_identity:
        for i in range(len(d)):
            ax.add_patch(matplotlib.patches.Rectangle(
                (i, i), 1, 1,
                fill = False,
                edgecolor = '0.7',
                lw = 2,
            ))
    if len(d) > 10:
        for i, label in enumerate(ax.get_xticklabels()):
            if i % 2 != 0:
                label.set_visible(False)
        for i, label in enumerate(ax.get_yticklabels()):
            if i % 2 != 0:
                label.set_visible(False)

def ax_id_line(
    ax,
    mn,
    mx,
    color = "0.8",
    linestyle = "-",
    linewidth = 1.0,
):
    identity_line, = ax.plot(
            [mn, mx],
            [mn, mx])
    plt.setp(identity_line,
            color = color,
            linestyle = linestyle,
            linewidth = linewidth,
            marker = '',
            zorder = 0)

def ax_compare_samples_to_cdf(
    ax,
    samples,
    prob_dist,
    include_ks_test = True,
):
    samples = sorted(samples)
    sample_cd = np.arange(len(samples)) / float(len(samples))
    emp_line = sns.lineplot(
        x = samples,
        y = sample_cd,
        ax = ax,
        label = "Empirical",
    )
    if hasattr(prob_dist, "ppf"):
        x = np.linspace(prob_dist.ppf(0.001), prob_dist.ppf(0.999), 100)
        model_cd = prob_dist.cdf(x)
    else:
        prob_dist = sorted(prob_dist)
        x = prob_dist
        model_cd = np.arange(len(prob_dist)) / float(len(prob_dist))

    mod_line = sns.lineplot(
        x = x,
        y = model_cd,
        ax = ax,
        label = "Model",
    )
    ax.set(
        xlabel = "X",
        ylabel = "Cumulative density",
    )
    if include_ks_test:
        prob_arg = prob_dist
        if hasattr(prob_dist, "cdf"):
            prob_arg = prob_dist.cdf
        res = st.kstest(
            samples,
            prob_arg,
            alternative = "two-sided",
        )
        ks_str = f"KS D = {res.statistic:.2g}\np = {res.pvalue:.2g}"
        ax.text(
            0.99, 0.02,
            ks_str,
            horizontalalignment = "right",
            verticalalignment = "bottom",
            transform = ax.transAxes,
            zorder = 500,
            fontsize = 'small',
            # bbox = {
            #     'facecolor': 'white',
            #     'edgecolor': 'white',
            #     'pad': 2},
        )
    # ax.legend(bbox_to_anchor = (1.01, 0.5), loc = "center left")
    ax.legend(loc = "center right")
    return emp_line, mod_line

def ax_compare_nevents_samples_to_cdf(
    ax,
    samples,
    prior_samples,
    number_of_comparisons,
    include_ks_test = True,
):
    emp_line = sns.ecdfplot(
        x = samples,
        ax = ax,
        label = "Empirical",
    )
    mod_line = sns.ecdfplot(
        x = prior_samples,
        ax = ax,
        label = "Model",
    )
    ax.set(
        xlabel = "Number of events",
        ylabel = "Cumulative probability",
        xlim = (1, number_of_comparisons),
    )
    if include_ks_test:
        res = st.kstest(
            samples,
            prior_samples,
            alternative = "two-sided",
        )
        ks_str = f"KS D = {res.statistic:.2g}\np = {res.pvalue:.2g}"
        ax.text(
            0.99, 0.02,
            ks_str,
            horizontalalignment = "right",
            verticalalignment = "bottom",
            transform = ax.transAxes,
            zorder = 500,
            fontsize = 'small',
            # bbox = {
            #     'facecolor': 'white',
            #     'edgecolor': 'white',
            #     'pad': 2},
        )
    # ax.legend(bbox_to_anchor = (1.01, 0.5), loc = "center left")
    ax.legend(loc = "center right")
    return emp_line, mod_line

def ax_qq(ax, samples, prob_dist):
    probs = np.linspace(0.01, 0.99, num = 99)
    if hasattr(prob_dist, "ppf"):
        q = prob_dist.ppf(probs)
    else:
        q = np.quantile(prob_dist, probs)
    sample_q = np.quantile(samples, probs)
    mn = min(min(q), min(sample_q))
    mx = max(max(q), max(sample_q))
    identity_line, = ax.plot(
        [mn, mx],
        [mn, mx],
    )
    plt.setp(
        identity_line,
        color = '0.7',
        linestyle = '-',
        linewidth = 1.0,
        marker = '',
        zorder = 0,
    )
    ax.set_xlim(mn, mx)
    ax.set_ylim(mn, mx)
    line, = ax.plot(q, sample_q)
    plt.setp(
        line,
        marker = 'o',
        linestyle = '',
        markerfacecolor = 'none',
        markeredgecolor = '0.35',
        markeredgewidth = 1.0,
        # markersize = 2.5,
        zorder = 100,
    )
    ax.set(
        xlabel = "Model quantiles",
        ylabel = "Sample quantiles",
    )
    return line

def plot_abs_error_scatter(
    data,
    error_col,
    error_lower_col,
    error_upper_col,
    true_val_col,
    true_val_annot_pos = None,
    true_val_annot_vert_aln = "bottom",
    **kwargs,
):
    df = data.copy()
    df["x"] = range(1, len(df) + 1)
    ax = plt.gca()
    shared_args = {
        'elinewidth' : 1.0,
        'capsize' : 1.5,
        'barsabove' : False,
        'marker' : 'o',
        'linestyle' : '',
        'markeredgewidth' : 0.0,
        'markersize' : 6.5,
        'rasterized' : False,
        'alpha' : 0.5,
    }
    shared_args.update(kwargs)
    line = ax.errorbar(
        x = df["x"],
        y = df[error_col],
        yerr = get_errors(df[error_col], df[error_lower_col], df[error_upper_col]),
        ecolor = 'C0',
        markerfacecolor = 'C0',
        markeredgecolor = 'C0',
        zorder = 100,
        **shared_args,
    )
    ax.axhline(
        y = 0.0,
        color = "0.8",
        linestyle = "-",
        linewidth = 1.0,
        zorder = 0,
    )
    # ax.xaxis.set_visible(False)
    ax.set_xticks([])
    uniq_true_values = df[true_val_col].unique()
    if len(uniq_true_values) < 30:
        annot_args = {
            "horizontalalignment" : "center",
            "verticalalignment" : true_val_annot_vert_aln,
            "zorder" : 200,
            "fontsize" : 'small',
            # "bbox" : {
            #     'facecolor': 'white',
            #     'edgecolor': 'white',
            #     'pad': 2,
            # },
        }
        prev_x_sep = 1.0
        prev_true_val = uniq_true_values[0]
        y_pos = true_val_annot_pos
        for i, true_val in enumerate(uniq_true_values[1:]):
            first_row = df.loc[df[true_val_col] == true_val].iloc[0]
            x_sep = first_row["x"] + 0.5
            ax.axvline(
                x = x_sep,
                color = "0.8",
                linestyle = "-",
                linewidth = 1.0,
                zorder = 0,
            )
            annot_str = f"{prev_true_val}"
            x_pos = prev_x_sep + ((x_sep - prev_x_sep) / 2.0)
            ax.text(
                x_pos, y_pos,
                annot_str,
                **annot_args,
            )
            prev_x_sep = x_sep
            prev_true_val = true_val
        x_sep = len(df) + 0.0
        annot_str = f"{prev_true_val}"
        x_pos = prev_x_sep + ((x_sep - prev_x_sep) / 2.0)
        ax.text(
            x_pos, y_pos,
            annot_str,
            **annot_args,
        )

def plot_abs_error_grid(
    data_frame,
    error_col,
    error_lower_col,
    error_upper_col,
    true_val_col,
    sim_config_col = "simulation_config",
    inference_config_col = "inference_config",
    sim_id_col = "simulation_id",
    ordered_labels = None,
    height = 6.5,
    **kwargs,
):
    col_order = None
    row_order = None
    if ordered_labels:
        sim_labels = data_frame[sim_config_col].unique()
        inf_labels = data_frame[inference_config_col].unique()
        row_order = [l for l in ordered_labels if l in sim_labels]
        col_order = [l for l in ordered_labels if l in inf_labels]
    df = data_frame[[
        sim_config_col,
        inference_config_col,
        true_val_col,
        error_col,
        error_lower_col,
        error_upper_col,
        sim_id_col,
    ]].copy()

    df.sort_values(
        by = ["true_num_events", "simulation_id"],
        ascending = [True, True],
        inplace = True,
    )

    max_y = max(df[error_upper_col])
    grid = sns.FacetGrid(
        df,
        row = sim_config_col,
        col = inference_config_col,
        margin_titles = True,
        height = height,
        row_order = row_order,
        col_order = col_order,
        sharey = True,
        sharex = True,
    )
    grid.map_dataframe(
        plot_abs_error_scatter,
        error_col = error_col,
        error_lower_col = error_lower_col,
        error_upper_col = error_upper_col,
        true_val_col = true_val_col,
        true_val_annot_pos = max_y,
        true_val_annot_vert_aln = "bottom",
    )
    grid.set_axis_labels(
        "Simulation replicate",
        "Absolute error")
    grid.set_titles(
        col_template = "{col_name}",
        row_template = "{row_name}",
    )
    # grid.figure.subplots_adjust(wspace = 0.05, hspace = 0.05)
    return grid

def plot_violin_grid(
    data,
    value_col,
    sim_config_col = "simulation_config",
    inference_config_col = "inference_config",
    spaghettify_col = "simulation_id",
    value_label = None,
    inference_label = "Inference model",
    sim_label_template = "True model = {col_name}",
    ordered_labels = None,
    height = 6.5,
    violin_kwargs = {},
    spaghetti_kwargs = {},
):
    sim_order = None
    inf_order = None
    if ordered_labels:
        sim_labels = data[sim_config_col].unique()
        inf_labels = data[inference_config_col].unique()
        sim_order = [l for l in ordered_labels if l in sim_labels]
        inf_order = [l for l in ordered_labels if l in inf_labels]
    grid = sns.FacetGrid(
        data,
        row = None,
        col = sim_config_col,
        margin_titles = True,
        height = height,
        row_order = None,
        col_order = sim_order,
        sharey = True,
        sharex = True,
    )
    grid.map_dataframe(
        plot_violin,
        value_col = value_col,
        categorical_col = inference_config_col,
        spaghettify_col = spaghettify_col,
        ordered_labels = ordered_labels,
        spaghetti_kwargs = spaghetti_kwargs,
        **violin_kwargs,
    )
    if value_label is not None:
        grid.set_ylabels(value_label)
    if inference_label is not None:
        grid.set_xlabels(inference_label)
    if sim_label_template is not None:
        grid.set_titles(
            col_template = sim_label_template,
        )
    # grid.figure.subplots_adjust(wspace = 0.05, hspace = 0.05)
    return grid

def plot_violin(
    data,
    value_col,
    categorical_col,
    spaghettify_col = "simulation_id",
    ordered_labels = None,
    spaghetti_kwargs = {},
    **violin_kwargs,
):
    df = data.copy()
    vio_kwargs = {
        'inner' : 'point',
        # 'inner' : None,
    }
    vio_kwargs.update(violin_kwargs)
    cat_labels = tuple(df[categorical_col].unique())
    cat_order = sorted(cat_labels)
    if ordered_labels:
        cat_order = [l for l in ordered_labels if l in cat_labels]
    vio_kwargs['order'] = cat_order
    ax = sns.violinplot(
        data = df,
        x = categorical_col,
        y = value_col,
        hue = None,
        **vio_kwargs,
    )
    if spaghettify_col:
        spag_kwargs = {
            'color' : '0.8',
            'marker' : '',
            'linestyle' : '-',
            'linewidth' : 1.0,
            'zorder' : 100,
            'alpha' : 0.5,
        }
        spag_kwargs.update(spaghetti_kwargs)

        x_positions = tuple(range(len(cat_labels)))

        shared_units = tuple(df[spaghettify_col].unique())
        for u in shared_units:
            sub_df = df[df[spaghettify_col] == u].copy()
            y_positions = []
            for cat in cat_order:
               ss_df = sub_df[sub_df[categorical_col] == cat]
               y_positions.append(ss_df[value_col])
            line, = ax.plot(
                x_positions,
                y_positions,
            )
            plt.setp(
                line,
                **spag_kwargs,
            )
