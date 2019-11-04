#! /usr/bin/env python

"""
CLI program for summarizing ecoevolity output from analyzing datasets simulated
with simcoevolity.
"""

import os
import sys
import re
import math
import argparse
import glob

import pycoevolity

MATPLOTLIB_AVAILABLE = False
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    sys.stderr.write('WARNING: matplotlib could not be imported; '
            'plotting functionality not supported\n')

if MATPLOTLIB_AVAILABLE:
    # Use TrueType (42) fonts rather than Type 3 fonts
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    tex_font_settings = {
            "text.usetex": True,
            "font.family": "sans-serif",
            "text.latex.preamble" : [
                    "\\usepackage[T1]{fontenc}",
                    "\\usepackage[cm]{sfmath}",
                    ]
    }
    mpl.rcParams.update(tex_font_settings)



class ScatterData(object):
    def __init__(self,
            x = [],
            y = [],
            y_lower = [],
            y_upper = [],
            highlight_values = [],
            highlight_threshold = None,
            highlight_greater_than = True):
        self._x = x
        self._y = y
        self._y_lower = y_lower
        self._y_upper = y_upper
        self._highlight_values = highlight_values
        self._highlight_threshold = highlight_threshold
        self._highlight_greater_than = highlight_greater_than
        self._vet_data()
        self._highlight_indices = []
        self._populate_highlight_indices()
        self.highlight_color = (184 / 255.0, 90 / 255.0, 13 / 255.0) # pauburn

    @classmethod
    def init(cls, results, parameters,
            highlight_parameter_prefix = None,
            highlight_threshold = None,
            highlight_greater_than = True):
        d = cls()
        d._x = []
        d._y = []
        d._y_lower = []
        d._y_upper = []
        d._highlight_threshold = highlight_threshold
        d._highlight_values = []
        d._highlight_indices = []
        for parameter_str in parameters:
            d._x.extend(float(x) for x in results["true_{0}".format(parameter_str)])
            d._y.extend(float(x) for x in results["mean_{0}".format(parameter_str)])
            d._y_lower.extend(float(x) for x in results["eti_95_lower_{0}".format(parameter_str)])
            d._y_upper.extend(float(x) for x in results["eti_95_upper_{0}".format(parameter_str)])
            if highlight_parameter_prefix:
                d._highlight_values.extend(float(x) for x in results["{0}_{1}".format(
                        highlight_parameter_prefix,
                        parameter_str)])
        d._vet_data()
        d._populate_highlight_indices()
        return d

    def _vet_data(self):
        assert len(self._x) == len(self._y)
        if self._y_lower:
            assert len(self._x) == len(self._y_lower)
        if self._y_upper:
            assert len(self._x) == len(self._y_upper)
        if self._highlight_values:
            assert len(self._x) == len(self._highlight_values)

    def _populate_highlight_indices(self):
        if (self._highlight_values) and (self._highlight_threshold is not None):
            for i in range(len(self._x)):
                if self.highlight(i):
                    self._highlight_indices.append(i)

    def has_y_ci(self):
        return bool(self.y_lower) and bool(self.y_upper)

    def has_highlights(self):
        return bool(self._highlight_indices)

    def _get_x(self):
        return self._x
    x = property(_get_x)

    def _get_y(self):
        return self._y
    y = property(_get_y)

    def _get_y_lower(self):
        return self._y_lower
    y_lower = property(_get_y_lower)

    def _get_y_upper(self):
        return self._y_upper
    y_upper = property(_get_y_upper)

    def _get_highlight_indices(self):
        return self._highlight_indices
    highlight_indices = property(_get_highlight_indices)

    def _get_highlight_x(self):
        return [self._x[i] for i in self._highlight_indices]
    highlight_x = property(_get_highlight_x)

    def _get_highlight_y(self):
        return [self._y[i] for i in self._highlight_indices]
    highlight_y = property(_get_highlight_y)

    def _get_highlight_y_lower(self):
        return [self._y_lower[i] for i in self._highlight_indices]
    highlight_y_lower = property(_get_highlight_y_lower)

    def _get_highlight_y_upper(self):
        return [self._y_upper[i] for i in self._highlight_indices]
    highlight_y_upper = property(_get_highlight_y_upper)

    def highlight(self, index):
        if (not self._highlight_values) or (self._highlight_threshold is None):
            return False

        if self._highlight_greater_than:
            if self._highlight_values[index] > self._highlight_threshold:
                return True
            else:
                return False
        else:
            if self._highlight_values[index] < self._highlight_threshold:
                return True
            else:
                return False
        return False

def get_parameter_names(header):
    p = []
    for h in header:
        if h == "generation":
            continue
        if h.startswith("ln_likelihood"):
            continue
        if h.startswith("ln_prior"):
            continue
        if h.startswith("root_height_index"):
            continue
        if h.startswith("freq_1"):
            continue
        if h.startswith("mutation_rate"):
            continue
        p.append(h)
    return p

def get_results_header(header, number_of_comparisons):
    h = [
            "sim",
            "sample_size",
            "true_model",
            "map_model",
            "true_model_cred_level",
            "map_model_p",
            "true_model_p",
            "true_num_events",
            "map_num_events",
            "true_num_events_cred_level",
        ]
    for i in range(number_of_comparisons):
        h.append("num_events_{0}_p".format(i+1))

    for p in get_parameter_names(header):
        h.append("true_{0}".format(p))
        h.append("true_{0}_rank".format(p))
        h.append("mean_{0}".format(p))
        h.append("median_{0}".format(p))
        h.append("stddev_{0}".format(p))
        h.append("hpdi_95_lower_{0}".format(p))
        h.append("hpdi_95_upper_{0}".format(p))
        h.append("eti_95_lower_{0}".format(p))
        h.append("eti_95_upper_{0}".format(p))
        h.append("ess_{0}".format(p))
        h.append("ess_sum_{0}".format(p))
        h.append("psrf_{0}".format(p))
    return h

def append_results_from_sim_rep(
        sim_number,
        true_values,
        post_sample,
        number_of_runs,
        results):
    results["sample_size"].append(post_sample.number_of_samples)
    results["sim"].append(sim_number)

    true_model = tuple(int(true_values[h][0]) for h in post_sample.height_index_keys)
    true_model_p = post_sample.get_model_probability(true_model)
    true_model_cred = post_sample.get_model_credibility_level(true_model)
    map_models = post_sample.get_map_models()
    map_model = map_models[0]
    if len(map_models) > 1:
        if true_model in map_models:
            map_model = true_model
    map_model_p = post_sample.get_model_probability(map_model)
    results["true_model"].append(",".join((str(i) for i in true_model)))
    results["map_model"].append(",".join((str(i) for i in map_model)))
    results["true_model_cred_level"].append(true_model_cred)
    results["map_model_p"].append(map_model_p)
    results["true_model_p"].append(true_model_p)
    
    true_nevents = int(true_values["number_of_events"][0])
    true_nevents_p = post_sample.get_number_of_events_probability(true_nevents)
    true_nevents_cred = post_sample.get_number_of_events_credibility_level(true_nevents)
    map_numbers_of_events = post_sample.get_map_numbers_of_events()
    map_nevents = map_numbers_of_events[0]
    if len(map_numbers_of_events) > 1:
        if true_nevents in map_numbers_of_events:
            map_nevents = true_nevents
    results["true_num_events"].append(true_nevents)
    results["map_num_events"].append(map_nevents)
    results["true_num_events_cred_level"].append(true_nevents_cred)
    for i in range(post_sample.number_of_comparisons):
        results["num_events_{0}_p".format(i + 1)].append(post_sample.get_number_of_events_probability(i + 1))

    for p in get_parameter_names(post_sample.header):
        true_val = float(true_values[p][0])
        true_val_rank = post_sample.get_rank(p, true_val)
        ss = pycoevolity.stats.get_summary(post_sample.parameter_samples[p])
        ess = pycoevolity.stats.effective_sample_size(
                post_sample.parameter_samples[p])
        nsamples_per_chain = int(post_sample.number_of_samples / number_of_runs)
        ess_sum = 0.0
        samples_by_chain = []
        for i in range(number_of_runs):
            chain_samples = post_sample.parameter_samples[p][i * nsamples_per_chain : (i + 1) * nsamples_per_chain]
            assert(len(chain_samples) == nsamples_per_chain)
            ess_sum += pycoevolity.stats.effective_sample_size(chain_samples)
            if number_of_runs > 1:
                samples_by_chain.append(chain_samples)
        psrf = float("nan")
        if number_of_runs > 1:
            psrf = pycoevolity.stats.potential_scale_reduction_factor(samples_by_chain)
        results["true_{0}".format(p)].append(true_val)
        results["true_{0}_rank".format(p)].append(true_val_rank)
        results["mean_{0}".format(p)].append(ss["mean"])
        results["median_{0}".format(p)].append(ss["median"])
        results["stddev_{0}".format(p)].append(math.sqrt(ss["variance"]))
        results["hpdi_95_lower_{0}".format(p)].append(ss["hpdi_95"][0])
        results["hpdi_95_upper_{0}".format(p)].append(ss["hpdi_95"][1])
        results["eti_95_lower_{0}".format(p)].append(ss["qi_95"][0])
        results["eti_95_upper_{0}".format(p)].append(ss["qi_95"][1])
        results["ess_{0}".format(p)].append(ess)
        results["ess_sum_{0}".format(p)].append(ess_sum)
        results["psrf_{0}".format(p)].append(psrf)

def parse_simulation_results(true_value_paths,
        burnin = 0):
    true_value_path_pattern_str = r'^(?P<prefix>.*)(?P<sim_label>simcoevolity-sim)-(?P<sim_number>\d+)-true-values.txt$'
    true_value_path_pattern = re.compile(true_value_path_pattern_str)

    result_keys = None
    results = None
    number_of_comparisons = None
    number_of_runs = None
    number_of_samples = None
    number_of_samples_per_chain = None
    expected_number_of_samples = None
    height_keys = []
    ancestral_pop_size_keys = []
    descendant_pop_size_keys = []
    for i, true_val_path in enumerate(true_value_paths):
        true_val_path_match = true_value_path_pattern.match(
                os.path.basename(true_val_path))
        sim_number_str = true_val_path_match.group("sim_number")
        sys.stderr.write("Parsing simulation {0}...\n".format(sim_number_str))
        sim_base_label = "{0}-{1}".format(
                true_val_path_match.group("sim_label"),
                sim_number_str)
        sim_full_label = "{0}{1}".format(
                true_val_path_match.group("prefix"),
                sim_base_label)
        state_sample_paths =  sorted(glob.glob(os.path.join(
                os.path.dirname(true_val_path),
                "*{0}-config-state-run-*.log*".format(sim_full_label))))
        post_sample = pycoevolity.posterior.PosteriorSample(
                paths = state_sample_paths,
                burnin = burnin)
        true_header = tuple(pycoevolity.parsing.parse_header_from_path(true_val_path))
        assert true_header == post_sample.header, (
                "true value header does not match state log header for sim "
                "{0}".format(sim_number_str))
        true_values = pycoevolity.parsing.get_dict_from_spreadsheets(
                [true_val_path],
                sep = "\t",
                header = None)
        for v in true_values.values():
            assert(len(v) == 1)
        if i == 0:
            number_of_runs = len(state_sample_paths)
            number_of_samples = post_sample.number_of_samples
            number_of_samples_per_chain = pycoevolity.parsing.data_line_count(
                    state_sample_paths[0])
            number_of_comparisons = post_sample.number_of_comparisons
            result_keys = get_results_header(true_header, number_of_comparisons)
            results = dict(zip(
                    result_keys,
                    ([] for i in range(len(result_keys)))
                    ))
            height_keys = post_sample.get_height_keys()
            ancestral_pop_size_keys = post_sample.get_ancestral_pop_size_keys()
            descendant_pop_size_keys = post_sample.get_descendant_pop_size_keys()
        if number_of_runs != len(state_sample_paths):
            sys.stderr.write(
                "WARNING: found output from {nruns} runs for sim {sim_num};\n"
                "found {expected_nruns} runs for previous sims.\n".format(
                    nruns = len(state_sample_paths),
                    sim_num = sim_number_str,
                    expected_nruns = number_of_runs))
        if number_of_samples != post_sample.number_of_samples:
            sys.stderr.write(
                "WARNING: parsed {nsamples} samples for sim {sim_num};\n"
                "parsed {expected_nsamples} samples for previous sims.\n".format(
                    nsamples = post_sample.number_of_samples,
                    sim_num = sim_number_str,
                    expected_nsamples = number_of_samples))
        if number_of_comparisons != post_sample.number_of_comparisons:
            sys.stderr.write(
                "ERROR: found {ncomps} comparisons for sim {sim_num};\n"
                "parsed {expected_ncomps} comparisons for previous sims.\n".format(
                    ncomps = post_sample.number_of_comparisons,
                    sim_num = sim_number_str,
                    expected_ncomps = number_of_comparisons))
            sys.exit(1)
        expected_number_of_samples = (len(state_sample_paths) *
                (number_of_samples_per_chain - burnin))
        if post_sample.number_of_samples != expected_number_of_samples:
            sys.stderr.write(
                "WARNING: parsed {nsamples} samples for sim {sim_num};\n"
                "expected to find {expected_nsamples} samples.\n".format(
                    nsamples = post_sample.number_of_samples,
                    sim_num = sim_number_str,
                    expected_nsamples = expected_number_of_samples))
            sys.exit(1)

        append_results_from_sim_rep(
                sim_number = int(sim_number_str),
                true_values = true_values,
                post_sample = post_sample,
                number_of_runs = number_of_runs,
                results = results)
    return (results,
            result_keys,
            number_of_comparisons,
            height_keys,
            ancestral_pop_size_keys,
            descendant_pop_size_keys)

def get_errors(values, lowers, uppers):
    n = len(values)
    assert(n == len(lowers))
    assert(n == len(uppers))
    return [[values[i] - lowers[i] for i in range(n)],
            [uppers[i] - values[i] for i in range(n)]]

def truncate_color_map(cmap, min_val = 0.0, max_val = 10, n = 100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(
                    n = cmap.name,
                    a = min_val,
                    b = max_val),
            cmap(list(get_sequence_iter(min_val, max_val, n))))
    return new_cmap

def get_sequence_iter(start = 0.0, stop = 1.0, n = 10):
    assert(stop > start)
    step = (stop - start) / float(n - 1)
    return ((start + (i * step)) for i in range(n))

def generate_scatter_plot(
        data,
        plot_path,
        parameter_symbol = "t",
        title = None,
        title_size = 16.0,
        x_label = None,
        x_label_size = 16.0,
        y_label = None,
        y_label_size = 16.0,
        plot_width = 3.5,
        plot_height = 3.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        force_shared_xy_ranges = True,
        xy_limits = None,
        include_coverage = True,
        include_rmse = True,
        include_identity_line = True,
        include_error_bars = True):
    if xy_limits:
        x_axis_min, x_axis_max, y_axis_min, y_axis_max = xy_limits
    else:
        x_min = min(data.x)
        x_max = max(data.x)
        y_min = min(data.y)
        y_max = max(data.y)
        if force_shared_xy_ranges:
            mn = min(x_min, y_min)
            mx = max(x_max, y_max)
            x_min = mn
            y_min = mn
            x_max = mx
            y_max = mx
        x_buffer = math.fabs(x_max - x_min) * 0.05
        x_axis_min = x_min - x_buffer
        x_axis_max = x_max + x_buffer
        y_buffer = math.fabs(y_max - y_min) * 0.05
        y_axis_min = y_min - y_buffer
        y_axis_max = y_max + y_buffer

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)

    proportion_within_ci = 0.0
    if include_coverage and data.has_y_ci():
        proportion_within_ci = pycoevolity.stats.get_proportion_of_values_within_intervals(
                data.x,
                data.y_lower,
                data.y_upper)
    rmse = 0.0
    if include_rmse:
        rmse = pycoevolity.stats.root_mean_square_error(data.x, data.y)
    ax = plt.subplot(gs[0, 0])
    if include_error_bars and data.has_y_ci():
        line = ax.errorbar(
                x = data.x,
                y = data.y,
                yerr = get_errors(data.y, data.y_lower, data.y_upper),
                ecolor = '0.65',
                elinewidth = 0.5,
                capsize = 0.8,
                barsabove = False,
                marker = 'o',
                linestyle = '',
                markerfacecolor = 'none',
                markeredgecolor = '0.35',
                markeredgewidth = 0.7,
                markersize = 2.5,
                zorder = 100,
                rasterized = True)
    else:
        line, = ax.plot(data.x, data.y)
        plt.setp(line,
                marker = 'o',
                linestyle = '',
                markerfacecolor = 'none',
                markeredgecolor = '0.35',
                markeredgewidth = 0.7,
                markersize = 2.5,
                zorder = 100,
                rasterized = True)
    if data.has_highlights():
        line, = ax.plot(x = data.highlight_x, y = data.highlight_y)
        plt.setp(line,
                marker = 'o',
                linestyle = '',
                markerfacecolor = 'none',
                markeredgecolor = data.highlight_color,
                markeredgewidth = 0.7,
                markersize = 2.5,
                zorder = 200,
                rasterized = True)
    ax.set_xlim(x_axis_min, x_axis_max)
    ax.set_ylim(y_axis_min, y_axis_max)
    if include_identity_line:
        identity_line, = ax.plot(
                [x_axis_min, x_axis_max],
                [y_axis_min, y_axis_max])
        plt.setp(identity_line,
                color = '0.7',
                linestyle = '-',
                linewidth = 1.0,
                marker = '',
                zorder = 0)
    if include_coverage:
        ax.text(0.02, 0.97,
                "\\normalsize\\noindent$p({0:s} \\in \\textrm{{\\sffamily CI}}) = {1:.3f}$".format(
                        parameter_symbol,
                        proportion_within_ci),
                horizontalalignment = "left",
                verticalalignment = "top",
                transform = ax.transAxes,
                size = 8.0,
                zorder = 300)
    if include_rmse:
        text_y = 0.97
        if include_coverage:
            text_y = 0.87
        ax.text(0.02, text_y,
                "\\normalsize\\noindent RMSE = {0:.2e}".format(
                        rmse),
                horizontalalignment = "left",
                verticalalignment = "top",
                transform = ax.transAxes,
                size = 8.0,
                zorder = 300)
    if x_label is not None:
        ax.set_xlabel(
                x_label,
                fontsize = x_label_size)
    if y_label is not None:
        ax.set_ylabel(
                y_label,
                fontsize = y_label_size)
    if title is not None:
        ax.set_title(plot_title,
                fontsize = title_size)
    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)
    plt.savefig(plot_path, dpi=600)

def generate_model_plot(
        results,
        plot_path,
        number_of_comparisons = 3,
        show_all_models = False,
        plot_title = None,
        include_x_label = True,
        include_y_label = True,
        include_median = True,
        include_cs = True,
        include_prop_correct = True,
        plot_width = 3.5,
        plot_height = 3.0,
        xy_label_size = 16.0,
        title_size = 16.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        lower_annotation_y = 0.02,
        upper_annotation_y = 0.92,
        model_key = "model",
        num_events_key = "num_events"):
    # Only show all models for 3 comparisons
    if show_all_models and (number_of_comparisons != 3):
        show_all_models = False
    cmap = truncate_color_map(plt.cm.binary, 0.0, 0.65, 100)
    model_to_index = {
            "0,0,0": 0,
            "0,0,1": 1,
            "0,1,0": 2,
            "0,1,1": 3,
            "0,1,2": 4,
            }
    index_to_model = {}
    for k, v in model_to_index.items():
        index_to_model[v] = k

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)

    true_map_nevents = []
    true_map_model = []
    true_map_nevents_probs = []
    for i in range(number_of_comparisons):
        true_map_nevents.append([0 for i in range(number_of_comparisons)])
        true_map_nevents_probs.append([[] for i in range(number_of_comparisons)])
    for i in range(5):
        true_map_model.append([0 for i in range(5)])
        # true_map_nevents_probs.append([[] for i in range(5)])
    true_nevents = tuple(int(x) for x in results["true_{num_events}".format(num_events = num_events_key)])
    map_nevents = tuple(int(x) for x in results["map_{num_events}".format(num_events = num_events_key)])
    true_model = tuple(x for x in results["true_{model}".format(model = model_key)])
    map_model = tuple(x for x in results["map_{model}".format(model = model_key)])
    true_nevents_cred_levels = tuple(float(x) for x in results["true_{num_events}_cred_level".format(num_events = num_events_key)])
    true_model_cred_levels = tuple(float(x) for x in results["true_{model}_cred_level".format(model = model_key)])
    assert(len(true_nevents) == len(map_nevents))
    assert(len(true_nevents) == len(true_nevents_cred_levels))
    assert(len(true_nevents) == len(true_model_cred_levels))
    assert(len(true_nevents) == len(true_model))
    assert(len(true_nevents) == len(map_model))

    true_nevents_probs = []
    map_nevents_probs = []
    for i in range(len(true_nevents)):
        true_nevents_probs.append(float(
            results["{num_events}_{n}_p".format(num_events = num_events_key, n = true_nevents[i])][i]))
        map_nevents_probs.append(float(
            results["{num_events}_{n}_p".format(num_events = num_events_key, n = map_nevents[i])][i]))
    assert(len(true_nevents) == len(true_nevents_probs))
    assert(len(true_nevents) == len(map_nevents_probs))

    mean_true_nevents_prob = sum(true_nevents_probs) / len(true_nevents_probs)
    median_true_nevents_prob = pycoevolity.stats.median(true_nevents_probs)

    true_model_probs = tuple(float(x) for x in results["true_{model}_p".format(model = model_key)])
    assert(len(true_nevents) == len(true_model_probs))

    mean_true_model_prob = sum(true_model_probs) / len(true_model_probs)
    median_true_model_prob = pycoevolity.stats.median(true_model_probs)

    nevents_within_95_cred = 0
    model_within_95_cred = 0
    ncorrect = 0
    model_ncorrect = 0
    for i in range(len(true_nevents)):
        true_map_nevents[map_nevents[i] - 1][true_nevents[i] - 1] += 1
        true_map_nevents_probs[map_nevents[i] - 1][true_nevents[i] - 1].append(map_nevents_probs[i])
        if show_all_models:
            true_map_model[model_to_index[map_model[i]]][model_to_index[true_model[i]]] += 1
        if true_nevents_cred_levels[i] <= 0.95:
            nevents_within_95_cred += 1
        if true_model_cred_levels[i] <= 0.95:
            model_within_95_cred += 1
        if true_nevents[i] == map_nevents[i]:
            ncorrect += 1
        if true_model[i] == map_model[i]:
            model_ncorrect += 1
    p_nevents_within_95_cred = nevents_within_95_cred / float(len(true_nevents))
    p_model_within_95_cred = model_within_95_cred / float(len(true_nevents))
    p_correct = ncorrect / float(len(true_nevents))
    p_model_correct = model_ncorrect /  float(len(true_nevents))

    ax = plt.subplot(gs[0, 0])

    if show_all_models:
        ax.imshow(true_map_model,
                origin = 'lower',
                cmap = cmap,
                interpolation = 'none',
                aspect = 'auto'
                )
        for i, row_list in enumerate(true_map_model):
            for j, n in enumerate(row_list):
                ax.text(j, i,
                        str(n),
                        horizontalalignment = "center",
                        verticalalignment = "center",
                        # fontsize = 8,
                        )
    else:
        ax.imshow(true_map_nevents,
                origin = 'lower',
                cmap = cmap,
                interpolation = 'none',
                aspect = 'auto'
                )
        for i, row_list in enumerate(true_map_nevents):
            for j, num_events in enumerate(row_list):
                ax.text(j, i,
                        str(num_events),
                        horizontalalignment = "center",
                        verticalalignment = "center")

    if include_cs:
        if show_all_models:
            ax.text(0.98, lower_annotation_y,
                    "$p(\\mathcal{{T}} \\in \\textrm{{\\sffamily CS}}) = {0:.3f}$".format(
                            p_model_within_95_cred),
                    horizontalalignment = "right",
                    verticalalignment = "bottom",
                    transform = ax.transAxes)
        else:
            ax.text(0.98, lower_annotation_y,
                    "$p(k \\in \\textrm{{\\sffamily CS}}) = {0:.3f}$".format(
                            p_nevents_within_95_cred),
                    horizontalalignment = "right",
                    verticalalignment = "bottom",
                    transform = ax.transAxes)
    if include_prop_correct:
        if show_all_models:
            ax.text(0.02, upper_annotation_y,
                    "$p(\\hat{{\\mathcal{{T}}}} = \\mathcal{{T}}) = {0:.3f}$".format(
                            p_model_correct),
                    horizontalalignment = "left",
                    verticalalignment = "bottom",
                    transform = ax.transAxes)
        else:
            ax.text(0.02, upper_annotation_y,
                    "$p(\\hat{{k}} = k) = {0:.3f}$".format(
                            p_correct),
                    horizontalalignment = "left",
                    verticalalignment = "bottom",
                    transform = ax.transAxes)
    if include_median:
        if show_all_models:
            ax.text(0.98, upper_annotation_y,
                    "$\\widetilde{{p(\\mathcal{{T}}|\\mathbf{{D}})}} = {0:.3f}$".format(
                            median_true_model_prob),
                    horizontalalignment = "right",
                    verticalalignment = "bottom",
                    transform = ax.transAxes)
        else:
            ax.text(0.98, upper_annotation_y,
                    "$\\widetilde{{p(k|\\mathbf{{D}})}} = {0:.3f}$".format(
                            median_true_nevents_prob),
                    horizontalalignment = "right",
                    verticalalignment = "bottom",
                    transform = ax.transAxes)
    if include_x_label:
        if show_all_models:
            ax.set_xlabel("True model ($\\mathcal{{T}}$)",
                    # labelpad = 8.0,
                    fontsize = xy_label_size)
        else:
            ax.set_xlabel("True \\# of events ($k$)",
                    # labelpad = 8.0,
                    fontsize = xy_label_size)
    if include_y_label:
        if show_all_models:
            ax.set_ylabel("MAP model ($\\hat{{\\mathcal{{T}}}}$)",
                    labelpad = 8.0,
                    fontsize = xy_label_size)
        else:
            ax.set_ylabel("MAP \\# of events ($\\hat{{k}}$)",
                    labelpad = 8.0,
                    fontsize = xy_label_size)
    if plot_title:
        ax.set_title(plot_title,
                fontsize = title_size)

    # Make sure ticks correspond only with number of events or model
    if not show_all_models:
        ax.xaxis.set_ticks(range(number_of_comparisons))
        ax.yaxis.set_ticks(range(number_of_comparisons))
    else:
        ax.xaxis.set_ticks(range(5))
        ax.yaxis.set_ticks(range(5))
    xtick_labels = [item for item in ax.get_xticklabels()]
    for i in range(len(xtick_labels)):
        if show_all_models:
            xtick_labels[i].set_text(index_to_model[i])
        else:
            xtick_labels[i].set_text(str(i + 1))
    ytick_labels = [item for item in ax.get_yticklabels()]
    for i in range(len(ytick_labels)):
        if show_all_models:
            ytick_labels[i].set_text(index_to_model[i])
        else:
            ytick_labels[i].set_text(str(i + 1))
    ax.set_xticklabels(xtick_labels)
    ax.set_yticklabels(ytick_labels)

    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)
    plt.savefig(plot_path)

def plot_path_writable(plot_path, force = False):
    if force or (not os.path.exists(plot_path)):
        return True
    sys.stderr.write("Plot path \'{0}\' already exists.\n"
            "  Use \'--force\' option to force overwrite.\n".format(plot_path))
    return False


def main(argv = sys.argv):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('sim_directory',
            metavar = 'PATH-TO-SIM-DIRECTORY',
            type = pycoevolity.argparse_utils.arg_is_dir,
            help = ('Path to directory with simulation files.'))
    parser.add_argument('-b', '--burnin',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
            default = 0,
            help = ('The number of samples to remove from the beginning of '
                    'each log file as burn in. Default: 0.'))
    parser.add_argument('-o', '--output-dir',
            action = 'store',
            type = str,
            default = os.getcwd(),
            help = ('Path to directory in which to write output files. '
                    'Default is the current working directory.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "",
            help = ('A prefix to prepend to all output files.'))
    parser.add_argument('-f', '--force',
            action = 'store_true',
            help = ('Overwrite any existing output files. By default, an error '
                    'is thrown if an output path exists.'))
    parser.add_argument('--no-plot',
            action = 'store_true',
            help = ('Skip plotting; only report summary table.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    sim_dir = args.sim_directory

    true_value_paths = sorted(glob.glob(os.path.join(sim_dir,
            "*simcoevolity-sim-*-true-values.txt")))
    sys.stderr.write("{0} simulations found in \'{1}\'\n".format(
            len(true_value_paths),
            sim_dir))
    
    (results,
            result_keys,
            number_of_comparisons,
            height_keys,
            ancestral_pop_size_keys,
            descendant_pop_size_keys) = parse_simulation_results(
                    true_value_paths,
                    burnin = args.burnin)

    # Write summary table to stdout
    for row in pycoevolity.parsing.dict_line_iter(results, sep = '\t', header = result_keys):
        sys.stdout.write(row)

    if args.no_plot:
        sys.exit(0)

    if not MATPLOTLIB_AVAILABLE:
        sys.stderr.write("ERROR: Unable to import matplotlib for generating plots\n")
        sys.exit(1)


    ###################################################################
    # Plotting
    ###################################################################
    brooks_gelman_1998_recommended_psrf = 1.2

    pad_left = 0.22
    pad_right = 0.94
    pad_bottom = 0.18
    pad_top = 0.965
    plot_width = 4.0
    plot_height = 3.2

    if not os.path.exists(args.output_dir):
        try:
            os.mkdir(args.output_dir)
        except:
            sys.stderr("ERROR: unable to create output directory "
                    "\'{0}\'\n".format(
                        args.output_dir))
            raise
    path_prefix = os.path.join(args.output_dir, args.prefix)

    # Plot event times
    sys.stderr.write("Plotting event times...\n")
    parameter_symbol = "t"
    data = ScatterData.init(
            results = results,
            parameters = height_keys,
            highlight_parameter_prefix = "psrf",
            highlight_threshold = brooks_gelman_1998_recommended_psrf,
            highlight_greater_than = True)
    plot_path = "{prefix}all-comparisons-event-time.pdf".format(
            prefix = path_prefix)
    if plot_path_writable(plot_path, force = args.force):
        generate_scatter_plot(
            data = data,
            plot_path = plot_path,
            parameter_symbol = parameter_symbol,
            title = None,
            title_size = 16.0,
            x_label = "True event time",
            x_label_size = 16.0,
            y_label = "Estimated event time",
            y_label_size = 16.0,
            plot_width = plot_width,
            plot_height = plot_height,
            pad_left = pad_left,
            pad_right = pad_right,
            pad_bottom = pad_bottom,
            pad_top = pad_top,
            force_shared_xy_ranges = True,
            xy_limits = None,
            include_coverage = True,
            include_rmse = True,
            include_identity_line = True,
            include_error_bars = True)
    for parameter in height_keys:
        data = ScatterData.init(
                results = results,
                parameters = [parameter],
                highlight_parameter_prefix = "psrf",
                highlight_threshold = brooks_gelman_1998_recommended_psrf,
                highlight_greater_than = True)
        comparison_label = parameter[12:]
        plot_path = "{prefix}{label}-event-time.pdf".format(
                prefix = path_prefix,
                label = comparison_label)
        if plot_path_writable(plot_path, force = args.force):
            generate_scatter_plot(
                    data = data,
                    plot_path = plot_path,
                    parameter_symbol = parameter_symbol,
                    title = None,
                    title_size = 16.0,
                    x_label = "True event time",
                    x_label_size = 16.0,
                    y_label = "Estimated event time",
                    y_label_size = 16.0,
                    plot_width = plot_width,
                    plot_height = plot_height,
                    pad_left = pad_left,
                    pad_right = pad_right,
                    pad_bottom = pad_bottom,
                    pad_top = pad_top,
                    force_shared_xy_ranges = True,
                    xy_limits = None,
                    include_coverage = True,
                    include_rmse = True,
                    include_identity_line = True,
                    include_error_bars = True)

    # Plot ancestral pop sizes
    sys.stderr.write("Plotting ancestral population sizes...\n")
    parameter_symbol = "N_e"
    data = ScatterData.init(
            results = results,
            parameters = ancestral_pop_size_keys,
            highlight_parameter_prefix = "psrf",
            highlight_threshold = brooks_gelman_1998_recommended_psrf,
            highlight_greater_than = True)
    plot_path = "{prefix}all-comparisons-ancestral-pop-size.pdf".format(
            prefix = path_prefix)
    if plot_path_writable(plot_path, force = args.force):
        generate_scatter_plot(
                data = data,
                plot_path = plot_path,
                parameter_symbol = parameter_symbol,
                title = None,
                title_size = 16.0,
                x_label = "True ancestral $N_e$",
                x_label_size = 16.0,
                y_label = "Estimated ancestral $N_e$",
                y_label_size = 16.0,
                plot_width = plot_width,
                plot_height = plot_height,
                pad_left = pad_left,
                pad_right = pad_right,
                pad_bottom = pad_bottom,
                pad_top = pad_top,
                force_shared_xy_ranges = True,
                xy_limits = None,
                include_coverage = True,
                include_rmse = True,
                include_identity_line = True,
                include_error_bars = True)
    for parameter in ancestral_pop_size_keys:
        data = ScatterData.init(
                results = results,
                parameters = [parameter],
                highlight_parameter_prefix = "psrf",
                highlight_threshold = brooks_gelman_1998_recommended_psrf,
                highlight_greater_than = True)
        comparison_label = parameter[14:]
        plot_path = "{prefix}{label}-ancestral-pop-size.pdf".format(
                prefix = path_prefix,
                label = comparison_label)
        if plot_path_writable(plot_path, force = args.force):
            generate_scatter_plot(
                    data = data,
                    plot_path = plot_path,
                    parameter_symbol = parameter_symbol,
                    title = None,
                    title_size = 16.0,
                    x_label = "True ancestral $N_e$",
                    x_label_size = 16.0,
                    y_label = "Estimated ancestral $N_e$",
                    y_label_size = 16.0,
                    plot_width = plot_width,
                    plot_height = plot_height,
                    pad_left = pad_left,
                    pad_right = pad_right,
                    pad_bottom = pad_bottom,
                    pad_top = pad_top,
                    force_shared_xy_ranges = True,
                    xy_limits = None,
                    include_coverage = True,
                    include_rmse = True,
                    include_identity_line = True,
                    include_error_bars = True)

    # Plot descendant pop sizes
    sys.stderr.write("Plotting descendant population sizes...\n")
    data = ScatterData.init(
            results = results,
            parameters = descendant_pop_size_keys,
            highlight_parameter_prefix = "psrf",
            highlight_threshold = brooks_gelman_1998_recommended_psrf,
            highlight_greater_than = True)
    plot_path = "{prefix}all-comparisons-descendant-pop-size.pdf".format(
            prefix = path_prefix)
    if plot_path_writable(plot_path, force = args.force):
        generate_scatter_plot(
                data = data,
                plot_path = plot_path,
                parameter_symbol = parameter_symbol,
                title = None,
                title_size = 16.0,
                x_label = "True descendant $N_e$",
                x_label_size = 16.0,
                y_label = "Estimated descendant $N_e$",
                y_label_size = 16.0,
                plot_width = plot_width,
                plot_height = plot_height,
                pad_left = pad_left,
                pad_right = pad_right,
                pad_bottom = pad_bottom,
                pad_top = pad_top,
                force_shared_xy_ranges = True,
                xy_limits = None,
                include_coverage = True,
                include_rmse = True,
                include_identity_line = True,
                include_error_bars = True)
    for parameter in descendant_pop_size_keys:
        data = ScatterData.init(
                results = results,
                parameters = [parameter],
                highlight_parameter_prefix = "psrf",
                highlight_threshold = brooks_gelman_1998_recommended_psrf,
                highlight_greater_than = True)
        comparison_label = parameter[9:]
        plot_path = "{prefix}{label}-descendant-pop-size.pdf".format(
                prefix = path_prefix,
                label = comparison_label)
        if plot_path_writable(plot_path, force = args.force):
            generate_scatter_plot(
                    data = data,
                    plot_path = plot_path,
                    parameter_symbol = parameter_symbol,
                    title = None,
                    title_size = 16.0,
                    x_label = "True descendant $N_e$",
                    x_label_size = 16.0,
                    y_label = "Estimated descendant $N_e$",
                    y_label_size = 16.0,
                    plot_width = plot_width,
                    plot_height = plot_height,
                    pad_left = pad_left,
                    pad_right = pad_right,
                    pad_bottom = pad_bottom,
                    pad_top = pad_top,
                    force_shared_xy_ranges = True,
                    xy_limits = None,
                    include_coverage = True,
                    include_rmse = True,
                    include_identity_line = True,
                    include_error_bars = True)

    # Plot model performance
    sys.stderr.write("Plotting event models...\n")
    plot_path = "{prefix}all-comparisons-model.pdf".format(
            prefix = path_prefix)
    if plot_path_writable(plot_path, force = args.force):
        generate_model_plot(
                results = results,
                plot_path = plot_path,
                number_of_comparisons = number_of_comparisons,
                plot_title = None,
                include_x_label = True,
                include_y_label = True,
                include_median = True,
                include_cs = True,
                include_prop_correct = True,
                plot_width = plot_width * 0.90,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.04,
                pad_right = 0.985,
                pad_bottom = pad_bottom,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015)


if __name__ == "__main__":
    main()
