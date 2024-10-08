#! /usr/bin/env python

import sys
import os
import logging
import re
import math

from pycoevolity.fileio import ReadFile
from pycoevolity import tempfs

_LOG = logging.getLogger(__name__)

HEADER_PATTERN = re.compile(r'^\s*\D.+')
NAN_PATTERN = re.compile(r'^\s*nan\s*$', re.IGNORECASE)

def line_count(stream):
    count = 0
    for line in stream:
        count += 1
    return count

def data_line_count(path):
    count = 0
    with ReadFile(path) as stream:
        parse_header(stream, strict = True, seek = False)
        for line in stream:
            if line.strip() == '':
                continue
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

def spreadsheet_iter(spreadsheets, sep = '\t', header = None, offset = 0):
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
                if row_idx < offset:
                    continue
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

def get_dict_from_spreadsheets(spreadsheets, sep = '\t', header = None, offset = 0):
    ss_iter = spreadsheet_iter(spreadsheets,
            sep = sep,
            header = header,
            offset = offset)
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


class Loci(object):
    symbol_to_states = {
            'A': ('A',),
            'C': ('C',),
            'G': ('G',),
            'T': ('T',),
            'R': ('A', 'G'),
            'Y': ('C', 'T'),
            'K': ('G', 'T'),
            'M': ('A', 'C'),
            'S': ('C', 'G'),
            'W': ('A', 'T'),
            'V': ('A', 'C', 'G'),
            'H': ('A', 'C', 'T'),
            'D': ('A', 'G', 'T'),
            'B': ('C', 'G', 'T'),
            'N': ('A', 'C', 'G', 'T'),
            '?': tuple(),
            '-': tuple(),
            }

    states_to_symbol = {
            ('A',)               : 'A',
            ('C',)               : 'C',
            ('G',)               : 'G',
            ('T',)               : 'T',
            ('A', 'G')           : 'R',
            ('C', 'T')           : 'Y',
            ('G', 'T')           : 'K',
            ('A', 'C')           : 'M',
            ('C', 'G')           : 'S',
            ('A', 'T')           : 'W',
            ('A', 'C', 'G')      : 'V',
            ('A', 'C', 'T')      : 'H',
            ('A', 'G', 'T')      : 'D',
            ('C', 'G', 'T')      : 'B',
            ('A', 'C', 'G', 'T') : 'N',
            }

    def __init__(self):
        self._labels = set()
        self._population_to_labels = None
        self._populations = None
        self._numbers_of_sites = []
        self._tempfs = tempfs.TempFileSystem()
        self._tmp_locus_paths = []
        self._remove_triallelic_sites = False
        self._convert_to_binary = False
        self._number_of_triallelic_sites = 0
        self._paths = []
        self._sequences_removed = {}
        self._label_suffix = None
        self._label_prefix = None
        self._sample_indices = None
        self._label_change_dict = {}

    def clone(self):
        o = self.__class__()
        o._labels = self._labels
        o._population_to_labels = self._population_to_labels
        o._populations = self._populations
        o._numbers_of_sites = self._numbers_of_sites
        o._tempfs = self._tempfs
        o._tmp_locus_paths = self._tmp_locus_paths
        o._remove_triallelic_sites = self._remove_triallelic_sites
        o._convert_to_binary = self._convert_to_binary
        o._number_of_triallelic_sites = self._number_of_triallelic_sites
        o._paths = self._paths
        o._sequences_removed = self._sequences_removed
        o._label_suffix = self._label_suffix
        o._label_prefix = self._label_prefix
        o._sample_indices = self._sample_indices
        o._label_change_dict = self._label_change_dict
        return o

    @classmethod
    def from_pyrad(cls, loci_path,
            remove_triallelic_sites = False,
            convert_to_binary = False,
            sequence_ids_to_remove = [],
            label_change_map_path = None):
        data = cls()
        data._remove_triallelic_sites = remove_triallelic_sites
        data._convert_to_binary = convert_to_binary
        data._paths = [loci_path]
        if sequence_ids_to_remove:
            data._sequences_removed = dict(zip(
                    sequence_ids_to_remove,
                    (0 for i in range(len(sequence_ids_to_remove)))
                    ))
        if label_change_map_path:
            data._parse_label_change_map(label_change_map_path)
        data._parse_loci_file()
        return data

    @classmethod
    def from_fastas(cls, paths,
            remove_triallelic_sites = False,
            convert_to_binary = False,
            sequence_ids_to_remove = [],
            label_change_map_path = None):
        data = cls()
        data._remove_triallelic_sites = remove_triallelic_sites
        data._convert_to_binary = convert_to_binary
        data._paths = paths
        if sequence_ids_to_remove:
            data._sequences_removed = dict(zip(
                    sequence_ids_to_remove,
                    (0 for i in range(len(sequence_ids_to_remove)))
                    ))
        if label_change_map_path:
            data._parse_label_change_map(label_change_map_path)
        data._parse_fasta_files()
        return data

    @classmethod
    def iter_from_msbayes_config(
        cls,
        msbayes_config,
        remove_triallelic_sites = False,
        convert_to_binary = False,
        recode_ambig_states_as_missing = False,
    ):
        msbayes_cfg = msbayes_config
        if not hasattr(msbayes_cfg, 'sample_table'):
            msbayes_cfg = MsbayesConfig(msbayes_config,
                                        recode_ambig_states_as_missing)
        sample_table = msbayes_cfg.sample_table
        for comp_idx, (comp, loci_data) in enumerate(sample_table.alignments.items()):
            data = cls()
            data._remove_triallelic_sites = remove_triallelic_sites
            data._convert_to_binary = convert_to_binary
            for locus, alignment in loci_data.items():
                sequences = []
                pop_seqs = alignment.sequences
                for i, seqs in enumerate(pop_seqs):
                    pop_suffix = "_comparison{0}pop{1}".format(comp_idx+1, i+1)
                    for name, s in seqs:
                        sequences.append(
                            (f"{name}{pop_suffix}",
                            s)
                        )
                data._process_locus(sequences)
                data._paths.append(alignment.path)
            assert len(data._numbers_of_sites) == len(data._tmp_locus_paths)
            yield comp, data

    def __del__(self):
        self._tempfs.purge()

    def sample_loci(self,
            rng,
            number_of_samples,
            with_replacement = True):
        if not with_replacement:
            # random library will raise ValueError if sample is larger than
            # population
            self._sample_indices = sorted(rng.sample(
                    range(self.number_of_loci),
                    number_of_samples))
            return

        indices = []
        max_index = self.number_of_loci - 1
        for i in range(number_of_samples):
            indices.append(
                    rng.randint(0, max_index))
        self._sample_indices = sorted(indices)

    def get_clone_of_unsampled_loci(self):
        c = self.clone()
        indices = []
        for i in range(len(self._numbers_of_sites)):
            if i not in self._sample_indices:
                indices.append(i)
        c._sample_indices = indices
        return c

    def split_loci(self, rng,
            auto_annotate_labels = True):
        self.clear_loci_samples()
        half_n_loci = int(math.ceil(self.number_of_loci / 2.0))
        self.sample_loci(
                rng = rng,
                number_of_samples = half_n_loci,
                with_replacement = False)
        c = self.get_clone_of_unsampled_loci()
        if auto_annotate_labels:
            if self._label_prefix is None:
                self._label_prefix = "set1"
            else:
                self._label_prefix += "set1"
            if self._label_suffix is None:
                self._label_suffix = "set1"
            else:
                self._label_suffix += "set1"
            if c._label_prefix is None:
                c._label_prefix = "set2"
            else:
                c._label_prefix += "set2"
            if c._label_suffix is None:
                c._label_suffix = "set2"
            else:
                c._label_suffix += "set2"
        return c

    def clear_loci_samples(self):
        self._sample_indices = None

    def group_labels_by_population(self,
            population_name_delimiter = "-",
            population_name_is_prefix = True):
        self._population_to_labels = {}
        for label in self._labels:
            label_parts = label.split(population_name_delimiter)
            if population_name_is_prefix:
                pop = label_parts[0]
            else:
                pop = label_parts[-1]
            if pop in self._population_to_labels:
                self._population_to_labels[pop].add(label)
            else:
                self._population_to_labels[pop] = set([label])
        self._populations = sorted(self._population_to_labels.keys())

    def _process_locus(self, sequences):
        if self._sequences_removed:
            seqs = []
            for seq_label, seq_chars in sequences:
                if seq_label in self._sequences_removed:
                    self._sequences_removed[seq_label] += 1
                else:
                    seqs.append((seq_label, seq_chars))
            sequences = seqs
        if self.convert_to_binary or self._remove_triallelic_sites:
            sequences = [[label, list(s)] for label, s in sequences]
        nsites = len(sequences[0][1])
        residues = [set() for i in range(nsites)]
        for site_idx in range(nsites):
            for seq_idx in range(len(sequences)):
                symbol = sequences[seq_idx][1][site_idx]
                r = self.symbol_to_states[symbol]
                if len(r) > 2:
                    raise Exception(
                            "Polymorphic site with more than two states: {0}".format(
                                    symbol))
                residues[site_idx].update(r)
        triallelic_site_indices = []
        for i, r in enumerate(residues):
            if len(r) > 2:
                triallelic_site_indices.append(i)
        self._number_of_triallelic_sites += len(triallelic_site_indices)
        if self._remove_triallelic_sites:
            for site_idx in sorted(triallelic_site_indices, reverse = True):
                for seq_idx in range(len(sequences)):
                    sequences[seq_idx][1].pop(site_idx)

        ######################################################################
        # Ecoevolity has option to automatically re-code triallelic sites, so
        # not doing it here.
        #
        # else:
        #     for site_idx in triallelic_site_indices:
        #         seq_idx = 0
        #         state0 = None
        #         state1 = None
        #         while state0 is None:
        #             states = self.symbol_to_states[sequences[seq_idx][1][site_idx]]
        #             if len(states) > 0:
        #                 state0 = states[0]
        #             if len(states) > 1:
        #                 state1 = states[1]
        #             seq_idx += 1
        #         while state1 is None:
        #             for s in self.symbol_to_states[sequences[seq_idx][1][site_idx]]:
        #                 if s != state0:
        #                     state1 = s
        #                     break
        #             seq_idx +=1
        #         possible_states = (state0, state1)
        #         for seq_idx in range(len(sequences)):
        #             replace_symbol = False
        #             states = list(self.symbol_to_states[sequences[seq_idx][1][site_idx]])
        #             if not states:
        #                 continue
        #             for state_idx in range(len(states)):
        #                 if not states[state_idx] in possible_states:
        #                     states[state_idx] = state1
        #                     replace_symbol = True
        #             if replace_symbol:
        #                 new_symbol = self.states_to_symbol[tuple(sorted(set(states)))]
        #                 sequences[seq_idx][1][site_idx] = new_symbol
        ######################################################################

        if self.convert_to_binary:
            for site_idx in range(nsites):
                seq_idx = 0
                state0 = None
                state1 = None
                while (seq_idx < len(sequences)) and (state0 is None):
                    states = self.symbol_to_states[sequences[seq_idx][1][site_idx]]
                    if not states:
                        seq_idx += 1
                        continue
                    if len(states) > 0:
                        state0 = states[0]
                    if len(states) > 1:
                        state1 = states[1]
                    seq_idx += 1
                while (seq_idx < len(sequences)) and (state1 is None):
                    for s in self.symbol_to_states[sequences[seq_idx][1][site_idx]]:
                        if s != state0:
                            state1 = s
                            break
                    seq_idx += 1
                if state0 is None:
                    continue
                possible_states = (state0, state1)
                if state1 is None:
                    possible_states = (state0,)
                for seq_idx in range(len(sequences)):
                    states = list(self.symbol_to_states[sequences[seq_idx][1][site_idx]])
                    if not states:
                        continue
                    for state_idx in range(len(states)):
                        if not states[state_idx] in possible_states:
                            states[state_idx] = state1
                    new_symbol = 0
                    if (len(states) < 2):
                        if states[0] == state1:
                            new_symbol = 2
                    else:
                        for s in states:
                            if s != state0:
                                new_symbol += 1
                    sequences[seq_idx][1][site_idx] = str(new_symbol)

        recorded_length = len(sequences[0][1])
        self._numbers_of_sites.append(recorded_length)
        tmp_path = self._tempfs.get_file_path()
        self._tmp_locus_paths.append(tmp_path)
        expected_length = nsites
        if self._remove_triallelic_sites:
            expected_length -= len(triallelic_site_indices) 
        with open(tmp_path, "w") as out:
            for label, seq in sequences:
                assert len(seq) == expected_length, (
                        "Seq length of {0}; expecting {1}".format(
                                len(seq),
                                expected_length))
                assert recorded_length == expected_length
                l = self._label_change_dict.get(label, label)
                self._labels.add(l)
                out.write("{0}\t{1}\n".format(l, "".join(seq)))

    def _parse_loci_file(self):
        with ReadFile(self._paths[0]) as stream:
            seqs = []
            for i, line in enumerate(stream):
                if line.startswith("//"):
                    self._process_locus(seqs)
                    seqs = []
                    continue
                try:
                    label, seq = line.strip().split()
                except ValueError as e:
                    sys.stderr.write("ERROR: Problem parsing line {0} of {1}:\n{2}".format(
                            i + 1, self._paths[0], line))
                    raise
                if label in self._sequences_removed:
                    self._sequences_removed[label] += 1
                    continue
                seq = seq.replace("N", "?")
                seqs.append([label, [c for c in seq]])
        if seqs:
            self._process_locus(seqs)
        assert len(self._numbers_of_sites) == len(self._tmp_locus_paths)

    def _get_label_suffix(self):
        return self._label_suffix
    def _set_label_suffix(self, s):
        self._label_suffix = str(s)

    label_suffix = property(_get_label_suffix, _set_label_suffix)

    def _get_label_prefix(self):
        return self._label_prefix
    def _set_label_prefix(self, s):
        self._label_prefix = str(s)

    label_prefix = property(_get_label_prefix, _set_label_prefix)

    def _get_path(self):
        return self._paths[0]

    path = property(_get_path)

    def _get_total_number_of_sites(self):
        if not self._sample_indices:
            return sum(self._numbers_of_sites)
        else:
            return sum(ns for i, ns in enumerate(self._numbers_of_sites) if i in self._sample_indices)

    number_of_sites = property(_get_total_number_of_sites)

    def _get_number_of_loci(self):
        if not self._sample_indices:
            return len(self._numbers_of_sites)
        else:
            return len(self._sample_indices)

    number_of_loci = property(_get_number_of_loci)

    def _get_number_of_taxa(self):
        return len(self._labels)

    number_of_taxa = property(_get_number_of_taxa)

    def _get_labels(self):
        return sorted(self._labels)

    labels = property(_get_labels)

    def _get_convert_to_binary(self):
        return self._convert_to_binary
    convert_to_binary = property(_get_convert_to_binary)

    def _get_number_of_triallelic_sites_found(self):
        return self._number_of_triallelic_sites

    number_of_triallelic_sites_found = property(_get_number_of_triallelic_sites_found)

    def _get_number_of_triallelic_sites_removed(self):
        if self._remove_triallelic_sites:
            return self._number_of_triallelic_sites
        return 0

    number_of_triallelic_sites_removed = property(_get_number_of_triallelic_sites_removed)

    def _get_removed_sequence_counts(self):
        return self._sequences_removed

    removed_sequence_counts = property(_get_removed_sequence_counts)

    def get_nexus_taxa_block(self):
        if self._label_suffix or self._label_prefix:
            p = ""
            s = ""
            if self._label_prefix:
                p = self._label_prefix
            if self._label_suffix:
                s = self._label_suffix
            b = ("BEGIN TAXA;\n"
                 "    DIMENSIONS NTAX={ntax};\n"
                 "    TAXLABELS\n"
                 "        {labels}\n"
                 "    ;\n"
                 "END;".format(
                     ntax = self.number_of_taxa,
                     labels = "\n        ".join(p + l + s for l in self.labels)))
            return b
        b = ("BEGIN TAXA;\n"
             "    DIMENSIONS NTAX={ntax};\n"
             "    TAXLABELS\n"
             "        {labels}\n"
             "    ;\n"
             "END;".format(
                 ntax = self.number_of_taxa,
                 labels = "\n        ".join(self.labels)))
        return b

    def get_phylip_header(self):
        return "{ntax} {nchar}".format(
                 ntax = self.number_of_taxa,
                 nchar = self.number_of_sites)

    def get_nexus_characters_block_preamble(self):
        if self.convert_to_binary:
            s = ("BEGIN CHARACTERS;\n"
                 "    DIMENSIONS NCHAR={nchar};\n"
                 "    FORMAT DATATYPE=STANDARD SYMBOLS=\"012\" MISSING=? GAP=- INTERLEAVE=YES;\n"
                 "    MATRIX".format(
                     nchar = self.number_of_sites))
        else:
            s = ("BEGIN CHARACTERS;\n"
                 "    DIMENSIONS NCHAR={nchar};\n"
                 "    FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;\n"
                 "    MATRIX".format(
                     nchar = self.number_of_sites))
        return s

    def write_nexus_charset_block(self, stream):
        nsites = self._numbers_of_sites
        if self._sample_indices:
            nsites = [n for i, n in enumerate(self._numbers_of_sites) if i in self._sample_indices]
        stream.write("BEGIN SETS;\n")
        total_n_sites = 0
        for i, n in enumerate(nsites):
            stream.write("    CHARSET locus{locus_id}={start}-{end};\n".format(
                    locus_id = i + 1,
                    start = total_n_sites + 1,
                    end = total_n_sites + n))
            total_n_sites += n
        stream.write("END;\n")

    def _parse_tmp_locus_file(self, path):
        seqs = {}
        with ReadFile(path) as stream:
            for line in stream:
                label, seq = line.strip().split()
                seqs[label] = seq
        return seqs

    def _parse_tmp_locus_file_as_list(self, path):
        seqs = []
        with ReadFile(path) as stream:
            for line in stream:
                label, seq = line.strip().split()
                seqs.append((label, seq))
        return seqs

    def get_label_buffer_size(self):
        mx = -1
        for l in self._labels:
            if len(l) > mx:
                mx = len(l)
        if self._label_prefix:
            mx += len(self._label_prefix)
        if self._label_suffix:
            mx += len(self._label_suffix)
        return mx + 4

    def write_interleaved_sequences(self, stream = None, indent = "",
            use_names_at_interleaves = True):
        prefix = ""
        if self._label_prefix:
            prefix = self.label_prefix
        suffix = ""
        if self._label_suffix:
            suffix = self.label_suffix
        if stream is None:
            stream = sys.stdout
        nsites = 0
        label_buffer = self.get_label_buffer_size()
        labels = self._get_labels()
        paths = self._tmp_locus_paths
        if self._sample_indices:
            paths = [p for i, p in enumerate(self._tmp_locus_paths) if i in self._sample_indices]
        for i, tmp_path in enumerate(paths):
            if i > 0:
                stream.write("\n")
            seqs = self._parse_tmp_locus_file(tmp_path)
            seq_length = len(list(seqs.values())[0])
            nsites += seq_length
            for l in labels:
                s = seqs.get(l, "?" * seq_length)
                if (not use_names_at_interleaves) and (i > 0):
                    stream.write("{sequence}\n".format(sequence = s))
                else:
                    stream.write("{indent}{label:<{fill}}{sequence}\n".format(
                            indent = indent,
                            label = prefix + l + suffix,
                            fill = label_buffer,
                            sequence = s))
        assert nsites == self.number_of_sites

    def write_fasta_files(self,
            directory,
            pair_label = "0",
            population_name_delimiter = "-",
            population_name_is_prefix = True,
            write_sample_table = False,
            sample_table_stream = None,
            sample_table_directory = None):
        if not self._population_to_labels:
            self.group_labels_by_population(
                    population_name_delimiter = population_name_delimiter,
                    population_name_is_prefix = population_name_is_prefix)
        if sample_table_stream is None:
            sample_table_stream = sys.stdout
        if not os.path.isdir(directory):
            raise Exception("{0!r} is not a valid directory".format(directory))
        locus_number_buffer = len(str(self.number_of_loci))
        prefix = ""
        if self._label_prefix:
            prefix = self.label_prefix
        suffix = ""
        if self._label_suffix:
            suffix = self.label_suffix
        paths = self._tmp_locus_paths
        if self._sample_indices:
            paths = [p for i, p in enumerate(self._tmp_locus_paths) if i in self._sample_indices]
        for i, tmp_path in enumerate(paths):
            locus_number = "{locus_num:0{buffer_size}d}".format(
                    locus_num = i,
                    buffer_size = locus_number_buffer)
            path = os.path.join(directory,
                    "pair-{0}-locus-{1}.fasta".format(pair_label, locus_number))
            path_for_table = path
            if sample_table_directory:
                path_for_table = os.path.relpath(path, start = sample_table_directory)
            if os.path.exists(path):
                raise Exception("The path {0!r} already exists. "
                        "Please designate a different directory.".format(path))
            seqs = self._parse_tmp_locus_file(tmp_path)
            locus_length = None
            pop_sample_sizes = [0 for p in self._populations]
            with open(path, "w") as out:
                for pop_index, pop in enumerate(self._populations):
                    for l in self._population_to_labels[pop]:
                        if l in seqs:
                            s = seqs.pop(l)
                            if locus_length is None:
                                locus_length = len(s)
                            else:
                                assert locus_length == len(s)
                            out.write(">{label}\n{sequence}\n".format(
                                    label = prefix + l + suffix,
                                    sequence = s))
                            pop_sample_sizes[pop_index] += 1
                if len(seqs) > 0:
                    raise Exception("Unexpected sequence labels in temp locus file")
            if write_sample_table:
                sample_size_str = "\t".join(str(n) for n in pop_sample_sizes)
                sample_table_stream.write("pair-{pair_label}\tlocus-{locus_number}\t1.0\t1.0\t{sample_sizes}\t1.0\t{locus_length}\t0.25\t0.25\t0.25\t{path}\n".format(
                        pair_label = pair_label,
                        locus_number = locus_number,
                        sample_sizes = sample_size_str,
                        locus_length = locus_length,
                        path = path_for_table))

    def write_nexus(self, stream = None,
            include_charset_block = False):
        if stream is None:
            stream = sys.stdout
        stream.write("#NEXUS\n\n")
        stream.write("{0}\n\n".format(self.get_nexus_taxa_block()))
        stream.write("{0}\n".format(
                self.get_nexus_characters_block_preamble()))
        self.write_interleaved_sequences(stream, indent = " " * 8)
        stream.write("    ;\nEND;\n")
        if include_charset_block:
            stream.write("\n")
            self.write_nexus_charset_block(stream)

    def write_phylip(self, stream = None):
        if stream is None:
            stream = sys.stdout
        stream.write("{0}\n".format(self.get_phylip_header()))
        self.write_interleaved_sequences(stream,
                indent = "",
                use_names_at_interleaves = False)

    def write_union_fasta_files(self,
            directory):
        if not os.path.isdir(directory):
            raise Exception("{0!r} is not a valid directory".format(directory))
        locus_number_buffer = len(str(self.number_of_loci))
        prefix = ""
        if self._label_prefix:
            prefix = self.label_prefix
        suffix = ""
        if self._label_suffix:
            suffix = self.label_suffix
        labels = self._get_labels()
        path_indices = range(len(self._numbers_of_sites))
        if self._sample_indices:
            path_indices = self._sample_indices
        for counter, locus_idx in enumerate(path_indices):
            tmp_path = self._tmp_locus_paths[locus_idx]
            if len(self._paths) == len(self._numbers_of_sites):
                path = os.path.join(directory,
                        os.path.basename(self._paths[locus_idx]))
            else:
                locus_number = "{locus_num:0{buffer_size}d}".format(
                        locus_num = counter,
                        buffer_size = locus_number_buffer)
                path = os.path.join(directory,
                        "locus-{0}.fasta".format(locus_number))
            if os.path.exists(path):
                raise Exception("The path {0!r} already exists. "
                        "Please designate a different directory.".format(path))
            seqs = self._parse_tmp_locus_file(tmp_path)
            seq_length = len(list(seqs.values())[0])
            with open(path, "w") as out:
                for label in labels:
                    seq = seqs.get(label, "?" * seq_length)
                    out.write(">{label}\n{sequence}\n".format(
                            label = prefix + label + suffix,
                            sequence = seq))

    def _parse_fasta_files(self):
        for path in self._paths:
            seqs = self.parse_fasta_file(path)
            self._process_locus(seqs)
        assert len(self._numbers_of_sites) == len(self._tmp_locus_paths)

    def _parse_label_change_map(self, path):
        d = {}
        with ReadFile(path) as stream:
            for line_idx, line in enumerate(stream):
                labels = line.strip().split(",")
                if len(labels) != 2:
                    msg = ("Line {ln_idx} of {path} has {ncols} columns. "
                            "It should have 2.".format(
                                ln_idx = line_idx,
                                path = path,
                                ncols = len(labels)))
                    raise Exception(msg)
                current_label = labels[0].strip()
                new_label = labels[1].strip()
                if current_label in d:
                    msg = ("Current seq label '{label}' found more than once "
                            "in {path}".format(
                                label = current_label,
                                path = path))
                    raise Exception(msg)
                if new_label in d.values():
                    msg = ("Replacement seq label '{label}' found more than once "
                            "in {path}".format(
                                label = new_label,
                                path = path))
                    raise Exception(msg)
                d[current_label] = new_label
        self._label_change_dict = d

    @classmethod
    def parse_fasta_file(cls, path):
        seqs = []
        current_seq = []
        label = None
        with ReadFile(path) as stream:
            for i, line in enumerate(stream):
                stripped_line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        seqs.append((label, tuple(current_seq)))
                    label = stripped_line[1:]
                    current_seq = []
                else:
                    if stripped_line:
                        current_seq.extend([c for c in stripped_line])
        if current_seq:
            seqs.append((label, tuple(current_seq)))
        return seqs


class MsbayesConfig(object):

    def __init__(self, cfg_path, recode_ambig_states_as_missing = False):
        if not self.is_config(cfg_path):
            raise Exception(
                'The following file does not appear to be a valid '
                f'msbayes config: {cfg_path}')
        self.sample_table = None
        self.prior_settings = {}
        self._parse_config(cfg_path, recode_ambig_states_as_missing)

    @classmethod
    def is_config(cls, cfg_path, max_lines_to_sample_table = 10000):
        with ReadFile(cfg_path) as cfg_stream:
            for i, line in enumerate(cfg_stream):
                if i >= max_lines_to_sample_table:
                    return False
                if MsbayesSampleTable.begin_pattern.match(line.strip()):
                    return True
        return False

    def _parse_config(self, cfg_path, recode_ambig_states_as_missing = False):
        self.sample_table = MsbayesSampleTable(
            cfg_path,
            recode_ambig_states_as_missing,
        )
        self.prior_settings = self._parse_preamble(cfg_path)

    def _get_taxa(self):
        if not self.sample_table:
            return None
        if not self.sample_table.taxa:
            return None
        return self.sample_table.taxa

    taxa = property(_get_taxa)

    def _get_npairs(self):
        if self.taxa is None:
            return 0
        return len(self.taxa)

    npairs = property(_get_npairs)

    def _parse_preamble(self, cfg_path):
        d = {}
        with ReadFile(cfg_path) as cfg_stream:
            for i, line in enumerate(cfg_stream):
                if MsbayesSampleTable.begin_pattern.match(line.strip()):
                    break
                line = line.strip()
                if (line == '') or (line.startswith('#')):
                    continue
                line = line.replace("=", " ")
                key_val = [x.strip() for x in line.split()]
                if len(key_val) != 2:
                    _LOG.error(
                        'msbayes preamble item on line {0} is invalid:\n{1}\n'.format(
                            i+1, l))
                d[key_val[0].lower()] = key_val[1]
        return d

    def using_dpp_model(self):
        dpp_keys = (
            'concentrationshape',
            'concentrationscale',
            'taushape',
            'tauscale',
            'thetashape',
            'thetascale',
            'ancestralthetashape',
            'ancestralthetascale',
            'thetaparameters',
        )
        for k in dpp_keys:
            if k in self.prior_settings:
                return True
        return False

    def using_time_in_subs_per_site(self):
        value = int(self.prior_settings.get('timeinsubspersite', 0))
        return bool(value)

    def write_ecoevolity_model_settings(self, out_stream = None):
        if not out_stream:
            out_stream = sys.stdout
        if self.using_dpp_model():
            if (('concentrationshape' in self.prior_settings) and (
                'concentrationscale' in self.prior_settings)):
                shape = self.prior_settings['concentrationshape']
                scale = self.prior_settings['concentrationscale']
                out_stream.write('event_model_prior:\n')
                out_stream.write('    dirichlet_process:\n')
                out_stream.write('        parameters:\n')
                out_stream.write('            concentration:\n')
                out_stream.write('                estimate: true\n')
                out_stream.write('                prior:\n')
                out_stream.write('                    gamma_distribution:\n')
                out_stream.write(f'                        shape: {shape}\n')
                out_stream.write(f'                        scale: {scale}\n')
                out_stream.write('\n')

            if (('taushape' in self.prior_settings) and (
                'tauscale' in self.prior_settings)):
                tau_shape = self.prior_settings['taushape']
                tau_scale = self.prior_settings['tauscale']

                if self.using_time_in_subs_per_site():
                    out_stream.write('event_time_prior:\n')
                    out_stream.write('    gamma_distribution:\n')
                    out_stream.write(f'        shape: {tau_shape}\n')
                    out_stream.write(f'        scale: {tau_scale}\n')
                    out_stream.write('\n')

                elif (('thetashape' in self.prior_settings) and (
                    'thetascale' in self.prior_settings)):
                    theta_shape = float(self.prior_settings['thetashape'])
                    theta_scale = float(self.prior_settings['thetascale'])
                    theta_mean = theta_shape * theta_scale
                    tau_scale = float(tau_scale) * theta_mean
                    out_stream.write('event_time_prior:\n')
                    out_stream.write('    gamma_distribution:\n')
                    out_stream.write(f'        shape: {tau_shape}\n')
                    out_stream.write(f'        scale: {tau_scale}\n')
                    out_stream.write('\n')

        ne_shape = '5.0'
        ne_scale = '0.0002'
        if (('thetashape' in self.prior_settings) and (
            'thetascale' in self.prior_settings)):
            ne_shape = self.prior_settings['thetashape']
            # msbayes pop size prior is on 4*Ne*mu, whereas ecoevolity is on
            # Ne*mu
            ne_scale = float(self.prior_settings['thetascale']) / 4.0
        out_stream.write('global_comparison_settings:\n')
        out_stream.write('    ploidy: 1\n')
        out_stream.write('    genotypes_are_diploid: false\n')
        out_stream.write('    markers_are_dominant: false\n')
        out_stream.write('    constant_sites_removed: false\n')
        out_stream.write('    population_name_delimiter: \" \"\n')
        out_stream.write('    population_name_is_prefix: false\n')
        out_stream.write('    equal_population_sizes: false\n')
        out_stream.write('    parameters:\n')
        out_stream.write('        freq_1:\n')
        out_stream.write('            value: 0.5\n')
        out_stream.write('            estimate: false\n')
        out_stream.write('        mutation_rate:\n')
        out_stream.write('            value: 1.0\n')
        out_stream.write('            estimate: false\n')
        out_stream.write('        root_relative_population_size:\n')
        out_stream.write('            value: 1.0\n')
        out_stream.write('            estimate: true\n')
        out_stream.write('            prior:\n')
        out_stream.write('                gamma_distribution:\n')
        out_stream.write('                    shape: 100.0\n')
        out_stream.write('                    scale: 0.01\n')
        out_stream.write('        population_size:\n')
        out_stream.write('            estimate: true\n')
        out_stream.write('            prior:\n')
        out_stream.write('                gamma_distribution:\n')
        out_stream.write(f'                    shape: {ne_shape}\n')
        out_stream.write(f'                    scale: {ne_scale}\n')


class MsbayesSampleTable(object):
    begin_pattern = re.compile(r'^begin\s*sample_tbl$', re.IGNORECASE)
    end_pattern = re.compile(r'^end\s*sample_tbl$', re.IGNORECASE)

    def __init__(self, config_path, recode_ambig_states_as_missing = False):
        self.alignments = None
        self.ordering = []
        self._parse_table(config_path, recode_ambig_states_as_missing)

    def _parse_table(self, config_path, recode_ambig_states_as_missing = False):
        self.alignments = {}
        table_started = False
        table_finished = False
        row_num = 0
        with ReadFile(config_path) as config_stream:
            for i, l in enumerate(config_stream):
                line = l.strip()
                if self.end_pattern.match(line):
                    if not table_started:
                        raise Exception(
                            'hit end of sample table before beginning')
                    if len(self.alignments) < 1:
                        raise Exception(
                            'no rows found in sample table')
                    table_finished = True
                    break
                if self.begin_pattern.match(line):
                    table_started = True
                    continue
                if not table_started:
                    continue
                if (line == '') or (line.startswith('#')):
                    continue
                row_num += 1
                al = MsbayesAlignment(
                    line,
                    config_path,
                    recode_ambig_states_as_missing,
                )
                if not al.taxon_name in self.alignments:
                    self.alignments[al.taxon_name] = {}
                    self.alignments[al.taxon_name][al.locus_name] = al
                    self.ordering.append((al.taxon_name, al.locus_name))
                    continue
                if al.locus_name in self.alignments[al.taxon_name]:
                    raise Exception('locus {0} found twice '
                            'for taxon {1} at row {2} of sample '
                            'table'.format(al.locus_name, al.taxon_name,
                                    row_num))
                self.alignments[al.taxon_name][al.locus_name] = al
                self.ordering.append((al.taxon_name, al.locus_name))
            if not table_started:
                raise Exception('no sample table found')
            if not table_finished:
                raise Exception('no end of table found')

    def _get_taxa(self):
        return tuple(self.alignments.keys())

    taxa = property(_get_taxa)

    def _get_loci(self):
        l = []
        for t, d in self.alignments.items():
            for locus in d.keys():
                if not locus in l:
                    l.append(locus)
        return l

    loci = property(_get_loci)

    def _get_number_of_taxa(self):
        return len(self.taxa)

    npairs = property(_get_number_of_taxa)

    def get_sample_table_string(self):
        return '\n'.join(('BEGIN SAMPLE_TBL',
                '\n'.join((str(self.alignments[t][l]) for t, l in self.ordering)),
                'END SAMPLE_TBL'))

    def __str__(self):
        return self.get_sample_table_string()

    def equals(self, other):
        if not isinstance(other, MsbayesSampleTable):
            return False
        if len(self.alignments) != len(other.alignments):
            return False
        for i1, i2 in zip(self.alignments.items(), other.alignments.items()):
            if i1[0] != i2[0]:
                return False
            if len(i1[1]) != len(i2[1]):
                return False
            for t1, t2 in zip(i1[1].items(), i2[1].items()):
                if t1[0] != t2[0]:
                    return False
                if not t1[1].equals(t2[1]):
                    return False
        if self.ordering != self.ordering:
            return False
        return True


class MsbayesAlignment(object):

    def __init__(
        self,
        sample_table_row = None,
        config_path = None,
        recode_ambig_states_as_missing = False
    ):
        self.taxon_name = None
        self.locus_name = None
        self._ploidy_multiplier = None
        self._mutation_rate_multiplier = None
        self._number_of_gene_copies = None
        self._kappa = None
        self._length = None
        self._base_frequencies = None
        self.path = None
        self.sequences = None
        if sample_table_row != None:
            self._parse_sample_table_row(
                sample_table_row,
                config_path,
                recode_ambig_states_as_missing)

    def _parse_sample_table_row(
        self,
        sample_table_row,
        config_path,
        recode_ambig_states_as_missing = False
    ):
        row = [e.strip() for e in sample_table_row.strip().split()]
        if len(row) != 12:
            raise Exception(
                'sample table row has {0} columns'.format(len(row)))
        self.taxon_name, self.locus_name = row[0:2]
        self._ploidy_multiplier, self._mutation_rate_multiplier = row[2:4]
        self._number_of_gene_copies = tuple(row[4:6])
        self._kappa, self._length = row[6:8]
        self._base_frequencies = tuple(row[8:11])
        self.path = row[11]
        if config_path:
            self.path = os.path.join(
                os.path.dirname(config_path),
                self.path,
            )
        self._parse_fasta_file(self.path, recode_ambig_states_as_missing)

    def _get_ploidy_multiplier(self):
        return float(self._ploidy_multiplier)

    ploidy_multiplier = property(_get_ploidy_multiplier)

    def _get_mutation_rate_multiplier(self):
        return float(self._mutation_rate_multiplier)

    mutation_rate_multiplier = property(_get_mutation_rate_multiplier)

    def _get_number_of_gene_copies(self):
        return tuple(int(x) for x in self._number_of_gene_copies)

    number_of_gene_copies = property(_get_number_of_gene_copies)

    def _get_kappa(self):
        return float(self._kappa)

    kappa = property(_get_kappa)

    def _get_length(self):
        return int(self._length)

    length = property(_get_length)

    def _get_base_frequencies(self):
        return tuple(float(x) for x in self._base_frequencies)

    base_frequencies = property(_get_base_frequencies)

    def get_sample_table_row_elements(self):
        return(self.taxon_name,
                self.locus_name,
                self.ploidy_multiplier,
                self.mutation_rate_multiplier,
                self.number_of_gene_copies[0],
                self.number_of_gene_copies[1],
                self.kappa,
                self.length,
                self.base_frequencies[0],
                self.base_frequencies[1],
                self.base_frequencies[2],
                self.path)

    def get_sample_table_row_element_strings(self):
        return(self.taxon_name,
                self.locus_name,
                self._ploidy_multiplier,
                self._mutation_rate_multiplier,
                self._number_of_gene_copies[0],
                self._number_of_gene_copies[1],
                self._kappa,
                self._length,
                self._base_frequencies[0],
                self._base_frequencies[1],
                self._base_frequencies[2],
                self.path)

    def get_sample_table_row_string(self):
        return '\t'.join([str(e) for e in self.get_sample_table_row_element_strings()])

    def equals(self, other):
        return self.__dict__ == other.__dict__

    def __str__(self):
        return self.get_sample_table_row_string()

    def _parse_fasta_file(self, path, recode_ambig_states_as_missing = False):
        seqs = Loci.parse_fasta_file(path)
        if not len(seqs) == sum(self.number_of_gene_copies):
            raise Exception(
                "Expected {0} seqs in fasta file {1}, but only found {2}".format(
                    sum(self.number_of_gene_copies),
                    self.path,
                    len(seqs),
                )
            )
        if recode_ambig_states_as_missing:
            symbol_map = {}
            for symbol, state in Loci.symbol_to_states.items():
                if len(state) > 1:
                    symbol_map[symbol] = '?'
                else:
                    symbol_map[symbol] = symbol
            new_seqs = []
            for label, s in seqs:
                new_s = tuple(symbol_map[c] for c in s)
                new_seqs.append( (label, new_s) )
            seqs = new_seqs
        n_pop_1_seqs = self.number_of_gene_copies[0]
        self.sequences = tuple(seqs[0:n_pop_1_seqs]), tuple(seqs[n_pop_1_seqs:])
        assert len(self.sequences) == 2
        assert len(self.sequences[0]) == self.number_of_gene_copies[0]
        assert len(self.sequences[1]) == self.number_of_gene_copies[1]
