#! /usr/bin/env python

import sys
import os
import math
import logging

from pycoevolity import parsing
from pycoevolity import stats

_LOG = logging.getLogger(__name__)


class PosteriorModelSummary(object):
    def __init__(self, nevent_samples, model_samples):
        nevents = tuple(nevent_samples)
        models = tuple(model_samples)
        assert(len(nevents) == len(models))
        self.number_of_samples = len(nevents)
        self.number_of_events_probabilities = stats.get_freqs(nevents)
        self.model_probabilities = stats.get_freqs(models)

    def get_models(self):
        return sorted(self.model_probabilities.items(),
                      reverse = True,
                      key = lambda x: x[1])

    def get_numbers_of_events(self):
        return sorted(self.number_of_events_probabilities.items(),
                      reverse = True,
                      key = lambda x: x[1])

    def get_map_models(self):
        max_p = max(self.model_probabilities.values())
        map_models = []
        for m, p in self.get_models():
            if p < max_p:
                assert(len(map_models) > 0)
                return tuple(map_models)
            map_models.append(m)
        return tuple(map_models)

    def get_map_numbers_of_events(self):
        max_p = max(self.number_of_events_probabilities.values())
        map_numbers_of_events = []
        for m, p in self.get_numbers_of_events():
            if p < max_p:
                assert(len(map_numbers_of_events) > 0)
                return tuple(map_numbers_of_events)
            map_numbers_of_events.append(m)
        return tuple(map_numbers_of_events)

    def get_model_probability(self, model_tuple):
        return self.model_probabilities.get(model_tuple, 0.0)

    def get_number_of_events_probability(self, number_of_events):
        return self.number_of_events_probabilities.get(number_of_events, 0.0)

    def model_in_credibility_set(self, model_tuple, credibility_cut_off = 0.95):
        total_prob = 0.0
        for m, p in self.get_models():
            if m == model_tuple:
                return True
            total_prob += p
            if total_prob > credibility_cut_off:
                return False
        return False

    def number_of_events_in_credibility_set(self, number_of_events,
            credibility_cut_off = 0.95):
        total_prob = 0.0
        for n, p in self.get_numbers_of_events():
            if n == number_of_events:
                return True
            total_prob += p
            if total_prob > credibility_cut_off:
                return False
        return False

    def get_model_credibility_level(self, model_tuple):
        total_prob = 0.0
        for m, p in self.get_models():
            if m == model_tuple:
                return total_prob
            total_prob += p
        assert(math.fabs(1.0 - total_prob) < 1e-6)
        return 1.0

    def get_number_of_events_credibility_level(self, number_of_events):
        total_prob = 0.0
        for n, p in self.get_numbers_of_events():
            if n == number_of_events:
                return total_prob
            total_prob += p
        assert(math.fabs(1.0 - total_prob) < 1e-6)
        return 1.0


class PosteriorSummary(object):
    def __init__(self, paths, burnin = 0):
        self.paths = tuple(paths)
        self.parameter_summaries = {}
        self.model_summary = None
        self.height_index_keys = None
        self.header = None
        self.burnin = burnin
        self.number_of_samples = 0
        self._parse_posterior_paths()
    
    def _parse_posterior_paths(self):
        self.header = tuple(parsing.parse_header_from_path(self.paths[0]))
        d = parsing.get_dict_from_spreadsheets(self.paths, offset = self.burnin)
        n = len(d[self.header[0]])
        for k, v in d.items():
            if (k.startswith('generation') or 
                k.startswith('ln_likelihood') or
                k.startswith('ln_prior') or
                k.startswith('root_height_index')):
                continue
            self.parameter_summaries[k] = stats.get_summary(
                    (float(x) for x in v))
            assert(self.parameter_summaries[k]['n'] == n)
        self.number_of_samples = n
        self.height_index_keys = tuple(h for h in self.header if h.startswith('root_height_index'))
        models = []
        for i in range(n):
            models.append(tuple(int(d[h][i]) for h in self.height_index_keys))
        assert len(models) == self.number_of_samples
        self.model_probabilities = stats.get_freqs(models)
        nevent_samples = (int(x) for x in d['number_of_events'])
        self.model_summary = PosteriorModelSummary(
                nevent_samples,
                models)
        assert(self.model_summary.number_of_samples == self.number_of_samples)

    def get_models(self):
        return self.model_summary.get_models()

    def get_numbers_of_events(self):
        return self.model_summary.get_numbers_of_events()

    def model_in_credibility_set(self, model_tuple, credibility_cut_off = 0.95):
        return self.model_summary.model_in_credibility_set(
                model_tuple,
                credibility_cut_off)

    def number_of_events_in_credibility_set(self, number_of_events,
            credibility_cut_off = 0.95):
        return self.model_summary.number_of_events_in_credibility_set(
                number_of_events,
                credibility_cut_off)

    def get_model_credibility_level(self, model_tuple):
        return self.model_summary.get_model_credibility_level(model_tuple)

    def get_number_of_events_credibility_level(self, number_of_events):
        return self.model_summary.get_number_of_events_credibility_level(number_of_events)

    def get_model_probability(self, model_tuple):
        return self.model_summary.get_model_probability(model_tuple)

    def get_number_of_events_probability(self, number_of_events):
        return self.model_summary.get_number_of_events_probability(
                number_of_events)

    def get_map_models(self):
        return self.model_summary.get_map_models()

    def get_map_numbers_of_events(self):
        return self.model_summary.get_map_numbers_of_events()


class ChainConvergenceSummary(object):
    def __init__(self, paths, burnin_step_size = 100):
        self.paths = list(paths)
        self.posterior_samples = dict(zip(
                self.paths,
                (PosteriorSample(paths = [p], burnin = 0) for p in self.paths)))
        self._vet_posterior_samples()
        self.parameter_keys = []
        self.burnin_step_size = burnin_step_size
        self.chain_lengths = tuple(self.posterior_samples[p].number_of_samples for p in self.paths)
        self.number_of_chains = len(self.paths)
        self.equal_chain_lengths = bool(len(set(self.chain_lengths)) == 1)
        for h in self.posterior_samples[self.paths[0]].header:
            if h.startswith('number_of'):
                continue
            if h in self.posterior_samples[self.paths[0]].parameter_samples:
                self.parameter_keys.append(h)
        self.burnin_values = []
        self.ess_values = {}

    def _vet_posterior_samples(self):
        for i in range(1, len(self.paths)):
            assert (self.posterior_samples[self.paths[0]].header ==
                    self.posterior_samples[self.paths[i]].header), (
                    "Headers of log files {0} and {1} do not match".format(
                            self.paths[0], self.paths[i]))

    def write_summary(self, out = sys.stdout, err = sys.stderr):
        samples_remaining = [x - 1 for x in self.chain_lengths]
        burnin = 0
        out.write("parameter\tburnin\tpsrf\tess_concat\t{0}\n".format(
                "\t".join("ess_" + str(i) for i in range(len(self.paths)))))
        finished = False
        ess_means = {}
        while True:
            for i in range(len(samples_remaining)):
                if ((samples_remaining[i] < self.burnin_step_size) or
                        (samples_remaining[i] < 10)):
                   finished = True 
                samples_remaining[i] -= self.burnin_step_size
            if finished:
                break
            bi = burnin
            if bi < 1:
                bi = 1
            self.burnin_values.append(bi)
            for parameter in self.parameter_keys:
                chains = []
                samples = []
                for path, ps in self.posterior_samples.items():
                    chains.append(ps.parameter_samples[parameter][bi:])
                    samples.extend(ps.parameter_samples[parameter][bi:])
                ess_concat = stats.effective_sample_size(samples)
                if ess_concat == 0.0:
                    continue
                psrf = "-"
                if self.equal_chain_lengths and (self.number_of_chains > 1):
                    psrf = stats.potential_scale_reduction_factor(chains)
                ess = []
                for c in chains:
                    ess.append(stats.effective_sample_size(c))
                out.write("{parameter}\t{burnin}\t{psrf}\t{ess_concat}\t{ess}\n".format(
                        parameter = parameter,
                        burnin = bi,
                        psrf = psrf,
                        ess_concat = ess_concat,
                        ess = "\t".join(str(e) for e in ess)))
                if not parameter in self.ess_values:
                    self.ess_values[parameter] = [ess_concat]
                    ess_means[parameter] = ess_concat
                else:
                    self.ess_values[parameter].append(ess_concat)
                    ess_means[parameter] += ess_concat
            burnin += self.burnin_step_size
        min_ess = float('inf')
        min_ess_key = ""
        for k in ess_means:
            ess_means[k] /= len(self.burnin_values)
            if ess_means[k] < min_ess:
                min_ess = ess_means[k]
                min_ess_key = k
        max_ess = float('-inf')
        max_ess_index = -1
        for i, bi in enumerate(self.burnin_values):
            if self.ess_values[min_ess_key][i] > max_ess:
                max_ess = self.ess_values[min_ess_key][i]
                max_ess_index = i
        err.write("Continuous parameter with lowest mean ESS: {0}\n".format(
                min_ess_key))
        err.write("Burnin value that maximized ESS of {0}: {1} samples (ESS = {2:.2f})\n".format(
                min_ess_key,
                self.burnin_values[max_ess_index],
                max_ess))
        

class PosteriorSample(object):
    def __init__(self, paths, burnin = 0, include_time_in_coal_units = False):
        self.paths = list(paths)
        self.parameter_samples = {}
        self.model_summary = None
        self.height_index_keys = None
        self.height_keys = None
        self.height_labels = None
        self.tip_labels = None
        self.number_of_comparisons = None
        self.header = None
        self.burnin = burnin
        self.number_of_samples = 0
        self._parse_posterior_paths(include_time_in_coal_units)

    def _parse_posterior_paths(self, include_time_in_coal_units = False):
        self.header = tuple(parsing.parse_header_from_path(self.paths[0]))
        d = parsing.get_dict_from_spreadsheets(self.paths, offset = self.burnin)
        n = len(d[self.header[0]])
        for k, v in d.items():
            if (k.startswith('generation') or 
                # k.startswith('ln_likelihood') or
                k.startswith('ln_prior') or
                k.startswith('root_height_index')):
                continue
            if k.startswith('number_of'):
                self.parameter_samples[k] = tuple(int(x) for x in v)
            else:
                self.parameter_samples[k] = tuple(float(x) for x in v)
            assert(len(self.parameter_samples[k]) == n)
        self.number_of_samples = n

        self._parse_height_keys_and_labels()
        self.number_of_comparisons = len(self.height_keys)
        self._parse_tip_labels()
        assert len(self.tip_labels) == self.number_of_comparisons

        models = []
        for i in range(n):
            models.append(tuple(int(d[h][i]) for h in self.height_index_keys))
        assert len(models) == self.number_of_samples
        self.parameter_samples['model'] = tuple(models)
        self.model_summary = PosteriorModelSummary(
                self.parameter_samples['number_of_events'],
                self.parameter_samples['model'])
        assert(self.model_summary.number_of_samples == self.number_of_samples)

        if include_time_in_coal_units:
            for i in range(n):
                for label in self.height_labels:
                    ht_key = "root_height_{0}".format(label)
                    sz_key = "pop_size_{0}".format(label)
                    t = self.parameter_samples[ht_key][i]
                    n = self.parameter_samples[sz_key][i]
                    t_coal = t / (2.0 * n)
                    coal_key = "coal_root_height_{0}".format(label)
                    if coal_key in self.parameter_samples:
                        self.parameter_samples[coal_key].append(t_coal)
                    else:
                        self.parameter_samples[coal_key] = [t_coal]
            assert(len(self.parameter_samples[coal_key]) == self.number_of_samples)

    def _parse_height_keys_and_labels(self):
        ht_index_keys = []
        ht_keys = []
        labels = []
        ht_prefix = "root_height_"
        ht_index_prefix = "root_height_index_"
        for h in self.header:
            if h.startswith(ht_prefix):
                if h.startswith(ht_index_prefix):
                    ht_index_keys.append(h)
                else:
                    ht_keys.append(h)
                    labels.append(h[len(ht_prefix):])
        assert len(ht_keys) == len(ht_index_keys)
        assert len(ht_keys) == len(labels)
        self.height_index_keys = tuple(ht_index_keys)
        self.height_keys = tuple(ht_keys)
        self.height_labels = tuple(labels)

    def _parse_tip_labels(self):
        size_prefix = "pop_size_"
        root_size_prefix = "pop_size_root_"
        labels = []
        current_labels = []
        for h in self.header:
            if h.startswith(size_prefix):
                if h.startswith(root_size_prefix):
                    assert ((len(current_labels) > 0) or (
                            len(current_labels) < 3))
                    assert h[len(root_size_prefix):] == current_labels[0]
                    labels.append(tuple(current_labels))
                    current_labels = []
                else:
                    current_labels.append(h[len(size_prefix):])
        assert len(current_labels) == 0
        self.tip_labels = labels

    def get_models(self):
        return self.model_summary.get_models()

    def get_numbers_of_events(self):
        return self.model_summary.get_numbers_of_events()

    def model_in_credibility_set(self, model_tuple, credibility_cut_off = 0.95):
        return self.model_summary.model_in_credibility_set(
                model_tuple,
                credibility_cut_off)

    def number_of_events_in_credibility_set(self, number_of_events,
            credibility_cut_off = 0.95):
        return self.model_summary.number_of_events_in_credibility_set(
                number_of_events,
                credibility_cut_off)

    def get_model_credibility_level(self, model_tuple):
        return self.model_summary.get_model_credibility_level(model_tuple)

    def get_number_of_events_credibility_level(self, number_of_events):
        return self.model_summary.get_number_of_events_credibility_level(number_of_events)

    def get_model_probability(self, model_tuple):
        return self.model_summary.get_model_probability(model_tuple)

    def get_number_of_events_probability(self, number_of_events):
        return self.model_summary.get_number_of_events_probability(
                number_of_events)

    def get_map_models(self):
        return self.model_summary.get_map_models()

    def get_map_numbers_of_events(self):
        return self.model_summary.get_map_numbers_of_events()

    def get_rank(self, parameter_key, value):
        return stats.rank(self.parameter_samples[parameter_key], value)
    
    def get_heights_2d(self, label_map = {},
            include_model_indices = False):
        labels = []
        heights = []
        map_model = self.get_map_models()[0]
        for i, ht_key in enumerate(self.height_keys):
            tip_label = self.height_labels[i]
            assert ht_key.endswith(tip_label)
            assert tip_label == self.tip_labels[i][0]
            pretty_labels = [label_map.get(l, l) for l in self.tip_labels[i]]
            label = " | ".join(pretty_labels)
            if include_model_indices:
                label = "{0} {1}".format(label, map_model[i])
            hts = list(self.parameter_samples[ht_key])
            labs = [label] * len(hts)
            labels.extend(labs)
            heights.extend(hts)
        return labels, heights

    def get_labels_and_heights(self, label_map = {},
            include_model_indices = False):
        labels = []
        heights = []
        map_model = self.get_map_models()[0]
        for i, ht_key in enumerate(self.height_keys):
            tip_label = self.height_labels[i]
            assert ht_key.endswith(tip_label)
            assert tip_label == self.tip_labels[i][0]
            pretty_labels = [label_map.get(l, l) for l in self.tip_labels[i]]
            label = " | ".join(pretty_labels)
            if include_model_indices:
                label = "{0} {1}".format(label, map_model[i])
            labels.append(label)
            heights.append(self.parameter_samples[ht_key])
        return labels, heights

    def write_height_summary(self, out = sys.stdout):
        out.write("comparison\tmean\tmedian\tci_95_lower\tci_95_upper\thpdi_95_lower\thpdi_95_upper\n")
        for i, ht_key in enumerate(self.height_keys):
            tip_label = self.height_labels[i]
            assert ht_key.endswith(tip_label)
            assert tip_label == self.tip_labels[i][0]
            labels = self.tip_labels[i]
            label = " | ".join(labels)
            summary = stats.get_summary(self.parameter_samples[ht_key])
            out.write("{comparison}\t{mean}\t{median}\t{ci_95_lower}\t{ci_95_upper}\t{hpdi_95_lower}\t{hpdi_95_upper}\n".format(
                    comparison = label,
                    mean = summary["mean"],
                    median = summary["median"],
                    ci_95_lower = summary["qi_95"][0],
                    ci_95_upper = summary["qi_95"][1],
                    hpdi_95_lower = summary["hpdi_95"][0],
                    hpdi_95_upper = summary["hpdi_95"][1]))

    def write_size_summary(self, out = sys.stdout):
        out.write("population\tmean\tmedian\tci_95_lower\tci_95_upper\thpdi_95_lower\thpdi_95_upper\n")
        for tips in self.tip_labels:
            for t in tips:
                size_key = "pop_size_{0}".format(t)
                summary  = stats.get_summary(self.parameter_samples[size_key])
                out.write("{population}\t{mean}\t{median}\t{ci_95_lower}\t{ci_95_upper}\t{hpdi_95_lower}\t{hpdi_95_upper}\n".format(
                        population = t,
                        mean = summary["mean"],
                        median = summary["median"],
                        ci_95_lower = summary["qi_95"][0],
                        ci_95_upper = summary["qi_95"][1],
                        hpdi_95_lower = summary["hpdi_95"][0],
                        hpdi_95_upper = summary["hpdi_95"][1]))
            size_key = "pop_size_root_{0}".format(tips[0])
            summary  = stats.get_summary(self.parameter_samples[size_key])
            out.write("{population}\t{mean}\t{median}\t{ci_95_lower}\t{ci_95_upper}\t{hpdi_95_lower}\t{hpdi_95_upper}\n".format(
                    population = "{0}-root".format(tips[0]),
                    mean = summary["mean"],
                    median = summary["median"],
                    ci_95_lower = summary["qi_95"][0],
                    ci_95_upper = summary["qi_95"][1],
                    hpdi_95_lower = summary["hpdi_95"][0],
                    hpdi_95_upper = summary["hpdi_95"][1]))

    def get_population_sizes_2d(self, label_map = {}):
        labels = []
        sizes = []
        for tips in self.tip_labels:
            for t in tips:
                size_key = "pop_size_{0}".format(t)
                szs = list(self.parameter_samples[size_key])
                pretty_label = label_map.get(t, t)
                labels.extend([pretty_label] * len(szs))
                sizes.extend(szs)
            root_size_key = "pop_size_root_{0}".format(tips[0])
            root_sizes = list(self.parameter_samples[root_size_key])
            root_label = "root-{0}".format(tips[0])
            pretty_root_label = label_map.get(root_label, root_label)
            labels.extend([pretty_root_label] * len(root_sizes))
            sizes.extend(root_sizes)
        return labels, sizes

    def get_relative_population_sizes_2d(self, label_map = {}):
        labels = []
        sizes = []
        for tips in self.tip_labels:
            tip_sizes = []
            for t in tips:
                size_key = "pop_size_{0}".format(t)
                szs = list(self.parameter_samples[size_key])
                tip_sizes.append(szs)
            mean_tip_sizes = []
            if len(tips) == 1:
                mean_tip_sizes = tip_sizes[0]
            else:
                for sample_idx in range(len(tip_sizes[0])):
                    total = 0.0
                    for tip_idx in range(len(tips)):
                        total += tip_sizes[tip_idx][sample_idx]
                    mn = total / float(len(tips))
                    mean_tip_sizes.append(mn)
            root_size_key = "pop_size_root_{0}".format(tips[0])
            root_sizes = list(self.parameter_samples[root_size_key])
            rel_root_sizes = [root_sizes[i] / mean_tip_sizes[i] for i in range(len(root_sizes))]
            root_label = "root-{0}".format(tips[0])
            pretty_root_label = label_map.get(root_label, root_label)
            labels.extend([pretty_root_label] * len(root_sizes))
            sizes.extend(rel_root_sizes)
        return labels, sizes

    def get_labels_and_sizes(self, label_map = {}):
        labels = []
        sizes = []
        for tips in self.tip_labels:
            for t in tips:
                size_key = "pop_size_{0}".format(t)
                szs = list(self.parameter_samples[size_key])
                pretty_label = label_map.get(t, t)
                labels.append(pretty_label)
                sizes.append(szs)
            root_size_key = "pop_size_root_{0}".format(tips[0])
            root_sizes = list(self.parameter_samples[root_size_key])
            root_label = "root-{0}".format(tips[0])
            pretty_root_label = label_map.get(root_label, root_label)
            labels.append(pretty_root_label)
            sizes.append(root_sizes)
        return labels, sizes

    def get_labels_and_relative_sizes(self, label_map = {}):
        labels = []
        sizes = []
        for tips in self.tip_labels:
            tip_sizes = []
            for t in tips:
                size_key = "pop_size_{0}".format(t)
                szs = list(self.parameter_samples[size_key])
                tip_sizes.append(szs)
            mean_tip_sizes = []
            if len(tips) == 1:
                mean_tip_sizes = tip_sizes[0]
            else:
                for sample_idx in range(len(tip_sizes[0])):
                    total = 0.0
                    for tip_idx in range(len(tips)):
                        total += tip_sizes[tip_idx][sample_idx]
                    mn = total / float(len(tips))
                    mean_tip_sizes.append(mn)
            root_size_key = "pop_size_root_{0}".format(tips[0])
            root_sizes = list(self.parameter_samples[root_size_key])
            rel_root_sizes = [root_sizes[i] / mean_tip_sizes[i] for i in range(len(root_sizes))]
            root_label = "root-{0}".format(tips[0])
            pretty_root_label = label_map.get(root_label, root_label)
            labels.append(pretty_root_label)
            sizes.append(rel_root_sizes)
        return labels, sizes


class SumcoevolityNeventsTable(object):
    def __init__(self, sumcoevolity_nevents_table_path):
        self.path = sumcoevolity_nevents_table_path
        self.number_of_elements = None
        self.posterior_probs = None
        self.posterior_probs_annotations = None
        self.prior_probs = None
        self.prior_probs_annotations = None
        self.bayes_factors = None
        self.bayes_factors_annotations = None
        self.no_prior = False
        self._parse_sumcoevolity_table()

    def _parse_sumcoevolity_table(self):
        keys = ["number_of_events",
                "post_prob",
                "cumulative_post_prob",
                "prior_prob",
                "bf"]
        data = parsing.get_dict_from_spreadsheets([self.path],
                sep = "\t",
                offset = 0)
        assert sorted(keys) == sorted(data.keys()), (
                "Unexpected header in sumcoevolity nevents file {0}".format(
                        self.path))
        self.number_of_elements = max(int(x) for x in data["number_of_events"])
        for k in data:
            assert len(data[k]) == self.number_of_elements, (
                    "Unexpected number of rows in column {0} "
                    "of nevents file {1}".format(
                            k, self.path))
        post_probs = [-1] * self.number_of_elements
        post_probs_annotations = [-1] * self.number_of_elements
        pr_probs = [-1] * self.number_of_elements
        pr_probs_annotations = [-1] * self.number_of_elements
        bfs = [-1] * self.number_of_elements
        bfs_annotations = [-1] * self.number_of_elements
        for row_idx in range(self.number_of_elements):
            nevents_idx = int(data["number_of_events"][row_idx]) - 1
            post_p, post_p_annotation = self._parse_entry(data["post_prob"][row_idx])
            if post_p_annotation == "NA":
                raise Exception("Unexpected 'NA' for post_prob in {0}".format(
                        self.path))
            prior_p, prior_p_annotation = self._parse_entry(data["prior_prob"][row_idx])
            if prior_p_annotation == "NA":
                self.no_prior = True
            bf, bf_a = self._parse_entry(data["bf"][row_idx])
            post_probs[nevents_idx] = post_p
            post_probs_annotations[nevents_idx] = post_p_annotation
            pr_probs[nevents_idx] = prior_p
            pr_probs_annotations[nevents_idx] = prior_p_annotation
            bfs[nevents_idx] = bf
            bfs_annotations[nevents_idx] = bf_a
        assert -1 not in post_probs
        assert -1 not in pr_probs
        assert -1 not in bfs
        self.posterior_probs = tuple(post_probs)
        self.posterior_probs_annotations = tuple(post_probs_annotations)
        self.prior_probs = tuple(pr_probs)
        self.prior_probs_annotations = tuple(pr_probs_annotations)
        self.bayes_factors = tuple(bfs)
        self.bayes_factors_annotations = tuple(bfs_annotations)

    def _parse_entry(self, x):
        annotation = None
        value = None
        if x == "NA":
            return None, "NA"
        if x.startswith(">"):
            annotation = ">"
            value = float(x.strip(">"))
        elif x.startswith("<"):
            annotation = "<"
            value = float(x.strip("<"))
        else:
            value = float(x)
        return value, annotation
