#! /usr/bin/env python

import sys
import os
import math
import logging

from sumcoevolity import parsing
from sumcoevolity import stats

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


class PosteriorSample(object):
    def __init__(self, paths, burnin = 0):
        self.paths = list(paths)
        self.parameter_samples = {}
        self.model_summary = None
        self.height_index_keys = None
        self.height_keys = None
        self.comparison_labels = None
        self.number_of_comparisons = None
        self.header = None
        self.burnin = burnin
        self.number_of_samples = 0
        self._parse_posterior_paths()
    
    def get_heights_2d(self, label_map = None,
            include_model_indices = False):
        labels = []
        heights = []
        map_model = self.get_map_models()[0]
        for i, ht_key in enumerate(self.height_keys):
            label = self.comparison_labels[i]
            assert ht_key.endswith(label)
            hts = list(self.parameter_samples[ht_key])
            if label_map:
                label = label_map[label]
            if include_model_indices:
                label = "{0} {1}".format(label, map_model[i])
            labs = [label] * len(hts)
            labels.extend(labs)
            heights.extend(hts)
        return labels, heights

    def _parse_posterior_paths(self):
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
        self.number_of_comparisons = len(ht_keys)
        self.height_index_keys = tuple(ht_index_keys)
        self.height_keys = tuple(ht_keys)
        self.comparison_labels = tuple(labels)

        models = []
        for i in range(n):
            models.append(tuple(int(d[h][i]) for h in self.height_index_keys))
        assert len(models) == self.number_of_samples
        self.parameter_samples['model'] = tuple(models)
        self.model_summary = PosteriorModelSummary(
                self.parameter_samples['number_of_events'],
                self.parameter_samples['model'])
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

    def get_rank(self, parameter_key, value):
        return stats.rank(self.parameter_samples[parameter_key], value)
