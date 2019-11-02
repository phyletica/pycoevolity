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
        nsamples_per_chain = post_sample.number_of_samples / number_of_runs
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
    return result_keys, results

def get_errors(values, lowers, uppers):
    n = len(values)
    assert(n == len(lowers))
    assert(n == len(uppers))
    return [[values[i] - lowers[i] for i in range(n)],
            [uppers[i] - values[i] for i in range(n)]]

def generate_scatter_plot(
        data,
        plot_path,
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
    
    result_keys, results = parse_simulation_results(true_value_paths,
            burnin = args.burnin)

    # Write summary table to stdout
    for row in pycoevolity.parsing.dict_line_iter(results, sep = '\t', header = result_keys):
        sys.stdout.write(row)

    if args.no_plot:
        sys.exit(0)


    ###################################################################
    # Plotting
    ###################################################################
    import matplotlib as mpl
    
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
    
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    brooks_gelman_1998_recommended_psrf = 1.2

    pad_left = 0.16
    pad_right = 0.94
    pad_bottom = 0.12
    pad_top = 0.965
    plot_width = 2.8
    plot_height = 2.2

    header_keys = tuple(pycoevolity.parsing.parse_header_from_path(true_value_paths[0]))
    height_parameters = [k for k in header_keys if k.startswith("root_height_")]
    ancestral_pop_size_parameters = [k for k in header_keys if k.startswith("pop_size_root_")]
    print height_parameters
    
    


if __name__ == "__main__":
    main()
