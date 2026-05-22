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


def get_errors(values, lowers = None, uppers = None):
    values = tuple(values)
    if lowers is None:
        lowers = tuple(values)
    else:
        lowers = tuple(lowers)
    if uppers is None:
        uppers = tuple(values)
    else:
        uppers = tuple(uppers)
    n = len(values)
    assert(n == len(lowers))
    assert(n == len(uppers))
    return [[values[i] - lowers[i] for i in range(n)],
            [uppers[i] - values[i] for i in range(n)]]

def get_nonneg_errors(values, lowers = None, uppers = None):
    values = tuple(values)
    if lowers is None:
        lowers = tuple(values)
    else:
        lowers = tuple(lowers)
    if uppers is None:
        uppers = tuple(values)
    else:
        uppers = tuple(uppers)
    n = len(values)
    assert(n == len(lowers))
    assert(n == len(uppers))
    return [[max(0.0, values[i] - lowers[i]) for i in range(n)],
            [max(0.0, uppers[i] - values[i]) for i in range(n)]]

def get_true_v_map_nevents_data_frame(
    data_frame,
    row_col,
    column_col,
):
    num_comparisons = len(data_frame["map_model"][0].split(","))
    nevent_labels = tuple(range(1, num_comparisons + 1))
    row_labels = data_frame[row_col].unique()
    col_labels = data_frame[column_col].unique()
    rows = []
    for row_label in row_labels:
        for col_label in col_labels:
            counts = [ [0 for i in range(num_comparisons)] for j in range(num_comparisons) ]
            sub_df = data_frame.loc[
                (data_frame[row_col] == row_label)
                & (data_frame[column_col] == col_label)
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
                        row_col : row_label,
                        column_col : col_label,
                        'true_num_events' : true_idx + 1,
                        'map_num_events' : map_idx + 1,
                        'count' : count,
                    })
    return pd.DataFrame(rows)

def get_cred_interval_percent(column_headers):
    perc = None
    for col in column_headers:
        if (col.startswith("hpdi_")) or (col.startswith("eti_")):
            return int(col.split("_")[1])

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

def get_all_comparison_parameters(parameter_prefix, column_headers):
    pop_labels = get_all_population_labels(column_headers)
    prefix = parameter_prefix.rstrip("_")
    possible_params = [f"{prefix}_{l}" for l in pop_labels]
    params = [p for p in possible_params if f"mean_{p}" in column_headers]
    return tuple(params)

def get_stacked_parameter_data_frame(
    data_frame,
    parameters,
    parameter_root,
    extra_cols_to_keep = [],
):
    ci = get_cred_interval_percent(data_frame.columns)
    param_root = parameter_root.rstrip("_")
    param_keys = [
        f"psrf_{param_root}",
        f"true_{param_root}",
        f"true_{param_root}_rank",
        f"mean_{param_root}",
        f"median_{param_root}",
        f"stddev_{param_root}",
        f"hpdi_{ci}_lower_{param_root}",
        f"hpdi_{ci}_upper_{param_root}",
        f"eti_{ci}_lower_{param_root}",
        f"eti_{ci}_upper_{param_root}",
        f"ess_{param_root}",
        f"ess_sum_{param_root}",
    ]
    columns = {k : [] for k in param_keys + extra_cols_to_keep}
    for param in parameters:
        for k in param_keys:
            param_label = k.replace(param_root, param)
            columns[k].extend(data_frame[param_label])
        for col_key in extra_cols_to_keep:
            columns[col_key].extend(data_frame[col_key])
    n = None
    for k, vals in columns.items():
        if n is None:
            n = len(vals)
        else:
            assert len(vals) == n
    return pd.DataFrame(columns)

def process_scatter_grid(
    data_frame,
    parameters,
    row_col,
    column_col,
    parameter_root = None,
    use_mean = True,
    use_hpdi = True,
    xlabel = None,
    ylabel = None,
    ess_min = 200,
    psrf_max = 1.2,
    bad_sampling_color = "C1",
    ordered_labels = None,
    height = 4.5,
    annotate_stats = True,
    stat_label = None,
    annot_position = (0.02, 0.98),
    cred_percent = 95,
    scatter_kwargs = {},
    annotate_kwargs = {},
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
        df = get_stacked_parameter_data_frame(
            data_frame,
            parameters,
            parameter_root,
            extra_cols_to_keep = [row_col, column_col],
        )
        parameter = parameter_root
    else:
        df = data_frame
        parameter = parameters[0]

    cred_percent = int(cred_percent)
    cred_level = cred_percent / 100.0
    ci_prefix = f"eti_{cred_percent}"
    if use_hpdi:
        ci_prefix = f"hpdi_{cred_percent}"

    # Remove rows for which the parameter was not estimated (standard deviation
    # is zero)
    std_dev_col = f"stddev_{parameter}"
    if std_dev_col in df.columns:
        df = df[df[std_dev_col] > 0.0]
    true_col = f"true_{parameter}"
    true_val_rank_col = f"true_{parameter}_rank"
    est_col = f"median_{parameter}"
    if use_mean:
        est_col = f"mean_{parameter}"
    est_lower_col = None
    est_upper_col = None
    if f"{ci_prefix}_lower_{parameter}" in df.columns:
        est_lower_col = f"{ci_prefix}_lower_{parameter}"
    if f"{ci_prefix}_upper_{parameter}" in df.columns:
        est_upper_col = f"{ci_prefix}_upper_{parameter}"
    ess_col = None
    if f"ess_{parameter}" in df.columns:
        ess_col = f"ess_{parameter}"
    psrf_col = None
    if f"psrf_{parameter}" in df.columns:
        psrf_col = f"psrf_{parameter}"
    grid = None
    if len(df) > 0:
        grid = plot_scatter_grid(
            data_frame = df,
            true_col = true_col,
            est_col = est_col,
            row_col = row_col,
            column_col = column_col,
            est_lower_col = est_lower_col,
            est_upper_col = est_upper_col,
            true_val_rank_col = true_val_rank_col,
            xlabel = xlabel,
            ylabel = ylabel,
            ess_col = ess_col,
            psrf_col = psrf_col,
            ess_min = ess_min,
            psrf_max = psrf_max,
            bad_sampling_color = bad_sampling_color,
            ordered_labels = ordered_labels,
            height = height,
            annotate_stats = annotate_stats,
            stat_label = stat_label,
            annot_position = annot_position,
            cred_level = cred_level,
            scatter_kwargs = scatter_kwargs,
            annotate_kwargs = annotate_kwargs,
        )
    return grid

def plot_scatter_grid(
    data_frame,
    true_col,
    est_col,
    row_col,
    column_col,
    est_lower_col = None,
    est_upper_col = None,
    true_val_rank_col = None,
    xlabel = None,
    ylabel = None,
    ess_col = None,
    psrf_col = None,
    ess_min = 200,
    psrf_max = 1.2,
    bad_sampling_color = "C1",
    ordered_labels = None,
    height = 4.5,
    annotate_stats = True,
    stat_label = None,
    annot_position = (0.02, 0.98),
    cred_level = 0.95,
    scatter_kwargs = {},
    annotate_kwargs = {},
):
    col_order = None
    row_order = None
    if ordered_labels:
        row_labels = data_frame[row_col].unique()
        col_labels = data_frame[column_col].unique()
        row_order = [l for l in ordered_labels if l in row_labels]
        col_order = [l for l in ordered_labels if l in col_labels]
    grid = sns.FacetGrid(
        data_frame,
        row = row_col,
        col = column_col,
        margin_titles = True,
        height = height,
        row_order = row_order,
        col_order = col_order,
        sharey = True,
        sharex = True,
    )
    grid.map_dataframe(
        plot_scatter,
        x = true_col,
        y = est_col,
        y_error_lower = est_lower_col,
        y_error_upper = est_upper_col,
        ess_col = ess_col,
        psrf_col = psrf_col,
        ess_min = ess_min,
        psrf_max = psrf_max,
        **scatter_kwargs,
    )
    y_limits = grid.axes.flat[0].get_ylim()
    x_limits = grid.axes.flat[0].get_xlim()
    mn = min(min(y_limits), min(x_limits))
    mx = max(max(y_limits), max(x_limits))
    # mn = min(min(data_frame[true_col]), min(data_frame[est_col]))
    # mx = max(max(data_frame[true_col]), max(data_frame[est_col]))
    for ax in grid.axes_dict.values():
        ax_id_line(
            ax = ax,
            mn = mn, 
            mx = mx,
            color = "0.7",
            linestyle = "-",
            linewidth = 1.0,
        )
    if annotate_stats:
        grid.map_dataframe(
            annotate_scatter,
            x = true_col,
            y = est_col,
            y_error_lower = est_lower_col,
            y_error_upper = est_upper_col,
            position = annot_position,
            cred_level = cred_level,
            stat_label = stat_label,
            **annotate_kwargs,
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

def annotate_scatter(
    data,
    x,
    y,
    y_error_lower = None,
    y_error_upper = None,
    position = (0.02, 0.98),
    cred_level = 0.95,
    stat_label = None,
    **kwargs,
):
    ax = plt.gca()
    prop_est_under = (
        sum(data[x] > data[y])
        / len(data[x])
    )
    sum_sq_err = ((data[x].values - data[y].values) ** 2).sum()
    mean_sq_err = sum_sq_err / len(data)
    root_mean_sq_err = math.sqrt(mean_sq_err)
    if not stat_label:
        stat_label = "x"
    wtest = st.wilcoxon(data[x].values - data[y].values)
    annot_str = (
        r"$p(\hat{{{stat_label}}} < {stat_label}) = {prop_under:.2g}$; $p = {pval:.2g}$"
        "\n"
        r"$\text{{RMSE}} = {rmse:.2g}$".format(
            stat_label = stat_label,
            prop_under = prop_est_under,
            pval = wtest.pvalue,
            rmse = root_mean_sq_err,
        )
    )
    if y_error_lower and y_error_upper:
        num_within_ci = (
            (data[x] >= data[y_error_lower])
            & (data[x] <= data[y_error_upper])
        ).sum()
        prop_within_ci = num_within_ci / len(data[x])
        annot_str = (
            r"$p(\hat{{{stat_label}}} < {stat_label}) = {prop_under:.2g}$; $p = {pval:.2g}$"
            "\n"
            r"$p({stat_label} \in {cred_level:.2f}\,\text{{CI}}) = {coverage:.2g}$"
            "\n"
            r"$\text{{RMSE}} = {rmse:.2g}$".format(
                stat_label = stat_label,
                prop_under = prop_est_under,
                pval = wtest.pvalue,
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
        position[0], position[1],
        annot_str,
        **default_args,
    )

def plot_scatter(
    data,
    x,
    y,
    y_error_lower = None,
    y_error_upper = None,
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
            | (df[psrf_col] > psrf_max)
        )
    elif ess_col:
        df["Poor MCMC sampling"] = df[ess_col] < ess_min
    elif psrf_col:
        df["Poor MCMC sampling"] = df[psrf_col] > psrf_max
    d = df
    if "Poor MCMC sampling" in df.columns:
        d = df[df["Poor MCMC sampling"] == False].copy()

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
    yerr = get_errors(d[y])
    if y_error_lower and y_error_upper:
        yerr = get_nonneg_errors(d[y], d[y_error_lower], d[y_error_upper])
    # TODO: ax.errorbar does not allow negative error values (i.e., the value
    # falling outside the error bar). However, with Bayesian cred intervals,
    # this can happen. E.g., the mean/mode can be outside the equal-tailed
    # credible interval. For now, we will set any negative values to zero for
    # plotting purposes. This does not effect coverage stats that are annotated
    # on plots (i.e., they use the un-fudged credible interval).
    line = ax.errorbar(
        x = d[x],
        y = d[y],
        yerr = yerr,
        ecolor = 'C0',
        markerfacecolor = 'C0',
        markeredgecolor = 'C0',
        zorder = 100,
        **shared_args,
    )
    if "Poor MCMC sampling" in df.columns:
        d = df[df["Poor MCMC sampling"]].copy()
        bad_yerr = get_errors(d[y])
        if y_error_lower and y_error_upper:
            bad_yerr = get_nonneg_errors(d[y], d[y_error_lower], d[y_error_upper])
        if len(d) > 0:
            bad_line = ax.errorbar(
                x = d[x],
                y = d[y],
                yerr = bad_yerr,
                ecolor = bad_sampling_color,
                markerfacecolor = bad_sampling_color,
                markeredgecolor = bad_sampling_color,
                zorder = 200,
                **shared_args,
            )

def plot_nevents_heatmap_grid(
    data_frame,
    row_col,
    column_col,
    ordered_labels = None,
    height = 4.5,
    annotate_counts = True,
    include_cbar = True,
    outline_identity = True,
    annotate_stats = True,
    cred_level = 0.95,
):
    data = get_true_v_map_nevents_data_frame(data_frame, row_col, column_col)
    vmin = min(data["count"])
    vmax = max(data["count"])
    col_order = None
    row_order = None
    if ordered_labels:
        row_labels = data_frame[row_col].unique()
        col_labels = data_frame[column_col].unique()
        row_order = [l for l in ordered_labels if l in row_labels]
        col_order = [l for l in ordered_labels if l in col_labels]
    grid = sns.FacetGrid(
        data,
        row = row_col,
        col = column_col,
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
            row_col = row_col,
            column_col = column_col,
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
    row_col,
    column_col,
    true_val_col,
    est_val_col,
    true_val_cred_col,
    true_val_prob_col,
    cred_level = 0.95,
    position = (0.02, 0.98),
    **kwargs,
):
    ax = plt.gca()
    row_labels = data[row_col].unique()
    assert len(row_labels) == 1
    row_label = row_labels[0]
    col_labels = data[column_col].unique()
    assert len(col_labels) == 1
    col_label = col_labels[0]
    df = full_data.loc[
        (full_data[row_col] == row_label)
        & (full_data[column_col] == col_label)
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
        position[0], position[1],
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
    color = "0.7",
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

def plot_error_scatter(
    data,
    true_val_col,
    est_col,
    est_lower_col = None,
    est_upper_col = None,
    ess_col = None,
    psrf_col = None,
    ess_min = 200,
    psrf_max = 1.2,
    bad_sampling_color = "C1",
    **kwargs,
):
    df = data.copy()
    df["x"] = range(1, len(df) + 1)
    error_col = f"{est_col}_distance"
    df[error_col] = df[est_col].values - df[true_val_col].values
    if ess_col and psrf_col:
        df["Poor MCMC sampling"] = (
            (df[ess_col] < ess_min)
            | (df[psrf_col] > psrf_max)
        )
    elif ess_col:
        df["Poor MCMC sampling"] = df[ess_col] < ess_min
    elif psrf_col:
        df["Poor MCMC sampling"] = df[psrf_col] > psrf_max
    d = df
    if "Poor MCMC sampling" in df.columns:
        d = df[df["Poor MCMC sampling"] == False].copy()

    # Errors will be all zeros
    yerr = get_errors(d[est_col])

    if est_lower_col and est_upper_col:
        yerr = get_nonneg_errors(d[est_col], d[est_lower_col], d[est_upper_col])

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
    # TODO: ax.errorbar does not allow negative error values (i.e., the value
    # falling outside the error bar). However, with Bayesian cred intervals,
    # this can happen. E.g., the mean/mode can be outside the equal-tailed
    # credible interval. For now, we will set any negative values to zero for
    # plotting purposes. This does not effect coverage stats that are annotated
    # on plots (i.e., they use the un-fudged credible interval).
    line = ax.errorbar(
        x = d["x"],
        y = d[error_col],
        yerr = yerr,
        ecolor = 'C0',
        markerfacecolor = 'C0',
        markeredgecolor = 'C0',
        zorder = 100,
        **shared_args,
    )
    ax.axhline(
        y = 0.0,
        color = "0.7",
        linestyle = "-",
        linewidth = 1.0,
        zorder = 0,
    )
    # ax.xaxis.set_visible(False)
    ax.set_xticks([])
    if "Poor MCMC sampling" in df.columns:
        d = df[df["Poor MCMC sampling"]].copy()
        bad_yerr = get_errors(d[est_col])
        if est_lower_col and est_upper_col:
            bad_yerr = get_nonneg_errors(d[est_col], d[est_lower_col], d[est_upper_col])
        if len(d) > 0:
            bad_line = ax.errorbar(
                x = d["x"],
                y = d[error_col],
                yerr = bad_yerr,
                ecolor = bad_sampling_color,
                markerfacecolor = bad_sampling_color,
                markeredgecolor = bad_sampling_color,
                zorder = 200,
                **shared_args,
            )

def annotate_true_values_on_error_scatter(
    data,
    true_val_col,
    annot_y_position,
    **kwargs,
):
    df = data.copy()
    df["x"] = range(1, len(df) + 1)
    ax = plt.gca()
    uniq_true_values = df[true_val_col].unique()
    if len(uniq_true_values) < 30:
        annot_args = {
            "horizontalalignment" : "center",
            "verticalalignment" : "bottom",
            "zorder" : 200,
            "fontsize" : "small",
            # "transform" : ax.transAxes,
            # "bbox" : {
            #     'facecolor': 'white',
            #     'edgecolor': 'white',
            #     'pad': 2,
            # },
        }
        annot_args.update(kwargs)
        annot_args["color"] = "black"
        first_last_buffer = len(df) * 0.05
        prev_x_sep = 1.0 - first_last_buffer
        prev_true_val = uniq_true_values[0]
        for i, true_val in enumerate(uniq_true_values[1:]):
            first_row = df.loc[df[true_val_col] == true_val].iloc[0]
            x_sep = first_row["x"] + 0.5
            ax.axvline(
                x = x_sep,
                color = "0.7",
                linestyle = "-",
                linewidth = 1.0,
                zorder = 0,
            )
            annot_str = f"{prev_true_val}"
            x_pos = prev_x_sep + ((x_sep - prev_x_sep) / 2.0)
            ax.text(
                x_pos, annot_y_position,
                annot_str,
                **annot_args,
            )
            prev_x_sep = x_sep
            prev_true_val = true_val
        x_sep = len(df) + first_last_buffer
        annot_str = f"{prev_true_val}"
        x_pos = prev_x_sep + ((x_sep - prev_x_sep) / 2.0)
        ax.text(
            x_pos, annot_y_position,
            annot_str,
            **annot_args,
        )

def plot_error_scatter_grid(
    data_frame,
    true_val_col,
    est_col,
    row_col,
    column_col,
    id_col,
    est_lower_col = None,
    est_upper_col = None,
    ess_col = None,
    psrf_col = None,
    ess_min = 200,
    psrf_max = 1.2,
    bad_sampling_color = "C1",
    ordered_labels = None,
    annotate_true_values = False,
    annotate_stats = True,
    stat_label = None,
    annot_stats_position = (0.02, 0.98),
    cred_level = 0.95,
    height = 4.5,
    scatter_kwargs = {},
    annot_true_vals_kwargs = {},
    annot_stats_kwargs = {},
):
    col_order = None
    row_order = None
    if ordered_labels:
        row_labels = data_frame[row_col].unique()
        col_labels = data_frame[column_col].unique()
        row_order = [l for l in ordered_labels if l in row_labels]
        col_order = [l for l in ordered_labels if l in col_labels]
    cols_to_keep = [
        row_col,
        column_col,
        true_val_col,
        est_col,
        id_col,
    ]
    if ess_col:
        cols_to_keep.append(ess_col)
    if psrf_col:
        cols_to_keep.append(psrf_col)
    if est_lower_col:
        cols_to_keep.append(est_lower_col)
    if est_upper_col:
        cols_to_keep.append(est_upper_col)
    df = data_frame[cols_to_keep].copy()

    df.sort_values(
        by = [true_val_col, id_col],
        ascending = [True, True],
        inplace = True,
    )

    grid = sns.FacetGrid(
        df,
        row = row_col,
        col = column_col,
        margin_titles = True,
        height = height,
        row_order = row_order,
        col_order = col_order,
        sharey = True,
        sharex = True,
    )
    grid.map_dataframe(
        plot_error_scatter,
        true_val_col = true_val_col,
        est_col = est_col,
        est_lower_col = est_lower_col,
        est_upper_col = est_upper_col,
        ess_col = ess_col,
        psrf_col = psrf_col,
        ess_min = ess_min,
        psrf_max = psrf_max,
        bad_sampling_color = bad_sampling_color,
        **scatter_kwargs
    )
    if annotate_stats:
        grid.map_dataframe(
            annotate_scatter,
            x = true_val_col,
            y = est_col,
            y_error_lower = est_lower_col,
            y_error_upper = est_upper_col,
            position = annot_stats_position,
            cred_level = cred_level,
            stat_label = stat_label,
            **annot_stats_kwargs,
        )
    grid.set_ylabels("Error")
    grid.set_xlabels("")
    if annotate_true_values:
        annot_y_position = grid.axes.flat[0].get_ylim()[0]
        grid.map_dataframe(
            annotate_true_values_on_error_scatter,
            true_val_col = true_val_col,
            annot_y_position = annot_y_position,
            **annot_true_vals_kwargs,
        )
        grid.set_xlabels("True value")
    grid.set_titles(
        col_template = "{col_name}",
        row_template = "{row_name}",
    )
    # grid.figure.subplots_adjust(wspace = 0.05, hspace = 0.05)
    return grid

def process_error_scatter_grid(
    data_frame,
    parameters,
    row_col,
    column_col,
    id_col,
    parameter_root = None,
    use_mean = True,
    use_hpdi = True,
    ess_min = 200,
    psrf_max = 1.2,
    bad_sampling_color = "C1",
    ordered_labels = None,
    annotate_true_values = False,
    annotate_stats = True,
    stat_label = None,
    annot_stats_position = (0.02, 0.98),
    cred_percent = 95,
    height = 4.5,
    scatter_kwargs = {},
    annot_true_vals_kwargs = {},
    annot_stats_kwargs = {},
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
        df = get_stacked_parameter_data_frame(
            data_frame,
            parameters,
            parameter_root,
            extra_cols_to_keep = [row_col, column_col, id_col],
        )
        parameter = parameter_root
    else:
        df = data_frame
        parameter = parameters[0]

    cred_percent = int(cred_percent)
    cred_level = cred_percent / 100.0
    ci_prefix = f"eti_{cred_percent}"
    if use_hpdi:
        ci_prefix = f"hpdi_{cred_percent}"

    # Remove rows for which the parameter was not estimated (standard deviation
    # is zero)
    std_dev_col = f"stddev_{parameter}"
    if std_dev_col in df.columns:
        df = df[df[std_dev_col] > 0.0]
    true_col = f"true_{parameter}"
    true_val_rank_col = f"true_{parameter}_rank"
    est_col = f"median_{parameter}"
    if use_mean:
        est_col = f"mean_{parameter}"
    est_lower_col = None
    est_upper_col = None
    if f"{ci_prefix}_lower_{parameter}" in df.columns:
        est_lower_col = f"{ci_prefix}_lower_{parameter}"
    if f"{ci_prefix}_upper_{parameter}" in df.columns:
        est_upper_col = f"{ci_prefix}_upper_{parameter}"
    ess_col = None
    if f"ess_{parameter}" in df.columns:
        ess_col = f"ess_{parameter}"
    psrf_col = None
    if f"psrf_{parameter}" in df.columns:
        psrf_col = f"psrf_{parameter}"
    grid = None
    if len(df) > 0:
        grid = plot_error_scatter_grid(
            data_frame = df,
            true_val_col = true_col,
            est_col = est_col,
            row_col = row_col,
            column_col = column_col,
            id_col = id_col,
            est_lower_col = est_lower_col,
            est_upper_col = est_upper_col,
            ess_col = ess_col,
            psrf_col = psrf_col,
            ess_min = ess_min,
            psrf_max = psrf_max,
            bad_sampling_color = bad_sampling_color,
            ordered_labels = ordered_labels,
            annotate_true_values = annotate_true_values,
            annotate_stats = annotate_stats,
            stat_label = stat_label,
            annot_stats_position = annot_stats_position,
            cred_level = cred_level,
            height = height,
            scatter_kwargs = scatter_kwargs,
            annot_true_vals_kwargs = annot_true_vals_kwargs,
            annot_stats_kwargs = annot_stats_kwargs,
        )
    return grid

def plot_violin_grid(
    data,
    value_col,
    plot_col,
    categorical_col,
    spaghettify_col,
    spaghettify = True,
    value_label = None,
    categorical_label = "Model",
    plot_label_template = "Model = {col_name}",
    ordered_labels = None,
    comparisons = None,
    height = 4.5,
    categorical_label_size = None,
    violin_kwargs = {},
    spaghetti_kwargs = {},
):
    plot_order = None
    cat_order = None
    if ordered_labels:
        plot_labels = data[plot_col].unique()
        cat_labels = data[categorical_col].unique()
        plot_order = [l for l in ordered_labels if l in plot_labels]
        cat_order = [l for l in ordered_labels if l in cat_labels]
    grid = sns.FacetGrid(
        data,
        row = None,
        col = plot_col,
        margin_titles = True,
        height = height,
        row_order = None,
        col_order = plot_order,
        sharey = True,
        sharex = True,
    )
    min_max_values = (min(data[value_col]), max(data[value_col]))
    grid.map_dataframe(
        plot_violin,
        value_col = value_col,
        categorical_col = categorical_col,
        spaghettify_col = spaghettify_col,
        spaghettify = spaghettify,
        ordered_labels = ordered_labels,
        comparisons = comparisons,
        min_max_values = min_max_values,
        categorical_label_size = categorical_label_size,
        spaghetti_kwargs = spaghetti_kwargs,
        **violin_kwargs,
    )
    if value_label is not None:
        grid.set_ylabels(value_label)
    if categorical_label is not None:
        grid.set_xlabels(categorical_label)
    if plot_label_template is not None:
        grid.set_titles(
            col_template = plot_label_template,
        )
    # grid.figure.subplots_adjust(wspace = 0.05, hspace = 0.05)
    return grid

def get_bracket_level(end_points, existing_end_points):
    level = 0
    p1, p2 = sorted(end_points)
    for ends in existing_end_points:
        e1, e2 = sorted(ends)
        if max(p1, e1) <= min(p2, e2):
            level += 1
    return level

def plot_violin(
    data,
    value_col,
    categorical_col,
    spaghettify_col,
    spaghettify = True,
    ordered_labels = None,
    comparisons = None,
    min_max_values = None,
    categorical_label_size = None,
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
    if categorical_label_size is not None:
        ax.tick_params(
            axis = 'x',
            labelsize = categorical_label_size,
        )
    
    x_positions = tuple(range(len(cat_labels)))
    if comparisons:
        if not spaghettify_col:
            sys.stderr.write(
                f"WARNING: Wilcoxon tests between comparisons requested "
                f"without the spaghettify column.\n"
            )
        else:
            bracket_kwargs = {
                'color' : '0.0',
                'marker' : '',
                'linestyle' : '-',
                'linewidth' : 1.0,
                'zorder' : 200,
                'alpha' : 1.0,
            }
            bracket_label_args = {
                'horizontalalignment' : "center",
                'verticalalignment' : "bottom",
                'zorder' : 200,
                'fontsize' : 'small',
            }
            if min_max_values:
                limits = min_max_values
            else:
                limits = (min(df[value_col]), max(df[value_col]))
            bracket_depth = (abs(limits[1] - limits[0]) * 0.01)
            dodge =  (abs(limits[1] - limits[0]) * 0.04)
            bracket_bottom = limits[1] + dodge
            bracket_top = bracket_bottom + bracket_depth
            existing_end_points = []
            for cat1, cat2 in comparisons:
                if cat1 not in cat_labels:
                    sys.stderr.write(
                        f"WARNING: category '{cat1}' not present; skipping "
                        f"Wilcoxon test of '{cat1}' vs '{cat2}'.\n"
                    )
                    continue
                if cat2 not in cat_labels:
                    sys.stderr.write(
                        f"WARNING: category '{cat2}' not present; skipping "
                        f"Wilcoxon test of '{cat1}' vs '{cat2}'.\n"
                    )
                    continue
                df1 = df[df[categorical_col] == cat1].sort_values(by = spaghettify_col)
                df2 = df[df[categorical_col] == cat2].sort_values(by = spaghettify_col)
                assert tuple(df1[spaghettify_col]) == tuple(df2[spaghettify_col])
                val_diff = df1[value_col].values - df2[value_col].values
                wtest = st.wilcoxon(val_diff)
                x1 = cat_order.index(cat1)
                x2 = cat_order.index(cat2)
                level = get_bracket_level((x1, x2), existing_end_points)
                existing_end_points.append((x1, x2))
                level_bump = level * (2.5 * dodge)
                bracket_y_pos = [
                    bracket_bottom + level_bump,
                    bracket_top + level_bump,
                    bracket_top + level_bump,
                    bracket_bottom + level_bump,
                ]
                bracket_x_pos = [x1, x1, x2, x2]
                bracket, = ax.plot(
                    bracket_x_pos,
                    bracket_y_pos,
                )
                plt.setp(
                    bracket,
                    **bracket_kwargs,
                )
                bracket_label_y = bracket_top + level_bump + bracket_depth
                bracket_label_x = sum((x1, x2)) / 2.0
                bracket_label = f"$p = {wtest.pvalue:.2g}$"
                ax.text(
                    bracket_label_x, bracket_label_y,
                    bracket_label,
                    **bracket_label_args,
                )
    if spaghettify and spaghettify_col:
        spag_kwargs = {
            'color' : '0.7',
            'marker' : '',
            'linestyle' : '-',
            'linewidth' : 1.0,
            'zorder' : 100,
            'alpha' : 0.5,
        }
        spag_kwargs.update(spaghetti_kwargs)

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
