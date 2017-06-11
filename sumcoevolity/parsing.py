#! /usr/bin/env python

import sys
import os
import logging
import re

from sumcoevolity.fileio import ReadFile 

_LOG = logging.getLogger(__name__)

HEADER_PATTERN = re.compile(r'^\s*\D.+')
NAN_PATTERN = re.compile(r'^\s*nan\s*$', re.IGNORECASE)

def line_count(stream):
    count = 0
    for line in stream:
        count += 1
    return count

def parse_header(file_stream, sep = '\t', strict = True, seek = True):
    try:
        header_line = next(file_stream)
    except StopIteration:
        if strict:
            raise Exception('did not find header in {0}'.format(file_stream.name))
        else:
            return None
    if not HEADER_PATTERN.match(header_line):
        if strict:
            raise Exception('did not find header in {0}'.format(file_stream.name))
        else:
            return None
    header = header_line.strip().split(sep)
    if seek:
        file_stream.seek(0)
    return header

def parse_header_from_path(path, sep = '\t', strict = True):
    with ReadFile(path) as stream:
        return parse_header(stream, sep = sep, strict = strict, seek = False)

def spreadsheet_iter(spreadsheets, sep = '\t', header = None):
    head_line = False
    if not header:
        head_line = True
        header = parse_header_from_path(spreadsheets[0], sep = sep)
    for sheet_idx, ss in enumerate(spreadsheets):
        with ReadFile(ss) as file_stream:
            if head_line:
                h = next(file_stream).strip().split(sep)
                if header != h:
                    raise Exception('headers do not match')
            for row_idx, row in enumerate(file_stream):
                if row.strip() == '':
                    continue
                r = [el.strip() for el in row.strip().split(sep)]
                if len(r) != len(header):
                    raise Exception('row {0} of spreadsheet {1} has {2} columns, '
                            'header has {3}'.format(row_idx + 1, sheet_idx + 1,
                                    len(r), len(header)))
                yield dict(zip(header, r))

def dict_line_iter(d, sep = '\t', header = None):
    if not header:
        header = sorted(d.keys())
    if sorted(header) != sorted(d.keys()):
        raise ValueError('header does not match dict keys')
    yield '{0}\n'.format(sep.join(header))
    for i in range(len(d[header[0]])):
        yield '{0}\n'.format(sep.join([str(d[h][i]) for h in header]))

def get_dict_from_spreadsheets(spreadsheets, sep = '\t', header = None):
    ss_iter = spreadsheet_iter(spreadsheets, sep = sep, header = header)
    row_dict = next(ss_iter)
    d = dict(zip(row_dict.keys(),
            [[row_dict[k]] for k in row_dict.keys()]))
    for row_dict in ss_iter:
        for k in row_dict.keys():
            d[k].append(row_dict[k])
    return d


class EcoevolityStdOut(object):
    splash_pattern = re.compile(
            r'^\s*estimating\s+evolutionary\s+coevality\s*$',
            re.IGNORECASE)
    patterns = {
# Summary of data from 3 comparisons:
            "number_of_comparisons" :
                re.compile(
                    r'^\s*summary\s+of\s+data\s+from\s+(?P<number_of_comparisons>\d+)\s+comparisons:\s*$',
                    re.IGNORECASE),
            "number_of_sites" :
                re.compile(
                    r'^\s*number\s+of\s+sites:\s+(?P<number_of_sites>\d+)\s*$',
                    re.IGNORECASE),
            "number_of_variable_sites" :
                re.compile(
                    r'^\s*number\s+of\s+variable\s+sites:\s+(?P<number_of_variable_sites>\d+)\s*$',
                    re.IGNORECASE),
            "number_of_patterns" :
                re.compile(
                    r'^\s*number\s+of\s+patterns:\s+(?P<number_of_patterns>\d+)\s*$',
                    re.IGNORECASE),
            "run_time" :
                re.compile(
                    r'^\s*runtime:\s+(?P<run_time>\d+)\s+seconds\.\s*$',
                    re.IGNORECASE),
    }

    def __init__(self, path):
        self._number_of_comparisons = None
        self._numbers_of_sites = tuple()
        self._numbers_of_variable_sites = tuple()
        self._numbers_of_patterns = tuple()
        self._run_time = None
        self._parse_std_out(path)

    def _get_empty_attribute_dict(self):
        return {
                "number_of_comparisons": None,
                "number_of_sites": [],
                "number_of_variable_sites": [],
                "number_of_patterns": [],
                "run_time": None,
                }

    def _parse_std_out(self, path):
        attributes = self._get_empty_attribute_dict()
        with ReadFile(path) as stream:
            for l in stream:
                line = l.strip()
                if self.splash_pattern.match(line):
                    attributes = self._get_empty_attribute_dict()
                for k, pattern in self.patterns.items():
                    m = pattern.match(line)
                    if m:
                        try:
                            attributes[k].append(int(m.group(k)))
                        except AttributeError:
                            attributes[k] = int(m.group(k))
        if attributes["number_of_comparisons"] is None:
            raise Exception(
                    "Unable to parse number of comparisons from {0!r}".format(
                            path))
        self._number_of_comparisons = attributes["number_of_comparisons"]
        if len(attributes["number_of_sites"]) != self._number_of_comparisons:
            raise Exception(
                    "Parsed {0} comparisons, but found number of sites for {1} "
                    "comparisons in {2!r}".format(
                            self._number_of_comparisons,
                            len(attributes["number_of_sites"]),
                            path))
        if len(attributes["number_of_variable_sites"]) != self._number_of_comparisons:
            raise Exception(
                    "Parsed {0} comparisons, but found number of variable "
                    "sites for {1} comparisons in {2!r}".format(
                            self._number_of_comparisons,
                            len(attributes["number_of_variable_sites"]),
                            path))
        if len(attributes["number_of_patterns"]) != self._number_of_comparisons:
            raise Exception(
                    "Parsed {0} comparisons, but found number of patterns for "
                    "{1} comparisons in {2!r}".format(
                            self._number_of_comparisons,
                            len(attributes["number_of_patterns"]),
                            path))
        if attributes["run_time"] is None:
            raise Exception(
                    "Unable to parse run time from {0!r}".format(
                            path))
        self._run_time = attributes["run_time"]
        self._numbers_of_sites = tuple(attributes["number_of_sites"])
        self._numbers_of_variable_sites = tuple(attributes["number_of_variable_sites"])
        self._numbers_of_patterns = tuple(attributes["number_of_patterns"])

    def _get_number_of_comparisons(self):
        return self._number_of_comparisons

    number_of_comparisons = property(_get_number_of_comparisons)

    def _get_run_time(self):
        return self._run_time

    run_time = property(_get_run_time)

    def _get_numbers_of_sites(self):
        return self._numbers_of_sites

    numbers_of_sites = property(_get_numbers_of_sites)

    def _get_numbers_of_variable_sites(self):
        return self._numbers_of_variable_sites

    numbers_of_variable_sites = property(_get_numbers_of_variable_sites)

    def _get_numbers_of_patterns(self):
        return self._numbers_of_patterns

    numbers_of_patterns = property(_get_numbers_of_patterns)

    def get_number_of_sites(self, comparison_index):
        return self._numbers_of_sites[comparison_index]

    def get_min_number_of_sites(self):
        return min(self._numbers_of_sites)

    def get_max_number_of_sites(self):
        return max(self._numbers_of_sites)

    def get_mean_number_of_sites(self):
        return sum(self._numbers_of_sites) / float(len(self._numbers_of_sites))
    
    def get_number_of_variable_sites(self, comparison_index):
        return self._numbers_of_variable_sites[comparison_index]

    def get_min_number_of_variable_sites(self):
        return min(self._numbers_of_variable_sites)

    def get_max_number_of_variable_sites(self):
        return max(self._numbers_of_variable_sites)

    def get_mean_number_of_variable_sites(self):
        return sum(self._numbers_of_variable_sites) / float(len(self._numbers_of_variable_sites))
    
    def get_number_of_patterns(self, comparison_index):
        return self._numbers_of_patterns[comparison_index]

    def get_min_number_of_patterns(self):
        return min(self._numbers_of_patterns)

    def get_max_number_of_patterns(self):
        return max(self._numbers_of_patterns)

    def get_mean_number_of_patterns(self):
        return sum(self._numbers_of_patterns) / float(len(self._numbers_of_patterns))
