#!/usr/bin/env python

import os
import sys
import json

def parse_ipyrad_loci_file(path, treat_n_as_missing = True):
    loci = []
    with open(path, "r") as in_stream:
        locus = {}
        labels = set()
        for line in in_stream:
            l = line.strip()
            if l.startswith("//") and l.endswith("|"):
                assert len(locus) > 0
                loci.append(locus)
                locus = {}
                continue
            label, seq = l.split()
            assert label not in locus
            s = seq.strip()
            if treat_n_as_missing:
                s = s.replace("N", "?").replace("n", "?")
            locus[label.strip()] = s
            labels.add(label.strip())
        assert l.startswith("//")
    return loci, sorted(labels)

def parse_ipyrad_json_file(path):
    with open(path, "r") as in_stream:
        data = json.load(in_stream)
    return data

def parse_ipyrad_json_sample_stats(path, stat_key = "reads_passed_filter"):
    d = parse_ipyrad_json_file(path)
    labels = d["assembly"]["samples"]
    stats = {}
    for label in labels:
        s = d["samples"][label]["stats"][stat_key]
        stats[label] = s
    return stats

def get_overlap_and_div(seq1, seq2, missing = ("?", "N", "n", "-")):
    assert len(seq1) == len(seq2)
    seq_length = len(seq1)
    overlap_count = 0
    diff_count = 0
    for site_idx in range(seq_length):
        if (seq1[site_idx] not in missing) and (seq2[site_idx] not in missing):
            overlap_count += 1
            if seq1[site_idx].upper() != seq2[site_idx].upper():
                diff_count += 1
    overlap = overlap_count / float(seq_length)
    div = None
    if overlap_count > 0:
        div = diff_count / float(overlap_count)
    return overlap, div

def get_pairwise_overlaps_and_divs(labels, loci):
    assert len(set(labels)) == len(labels)
    pair_overlaps_and_divs = {}
    for i in range(len(labels) - 1):
        for j in range(i + 1, len(labels)):
            pair = labels[i], labels[j]
            assert pair not in pair_overlaps_and_divs
            pair_overlaps_and_divs[pair] = []
            for locus in loci:
                overlap = 0.0
                div = None
                if (pair[0] in locus) and (pair[1] in locus):
                    overlap, div = get_overlap_and_div(locus[pair[0]], locus[pair[1]])
                pair_overlaps_and_divs[pair].append((overlap, div))
    return pair_overlaps_and_divs

def get_pairwise_overlap_count_and_mean_div(
    labels,
    loci,
    min_overlap,
):
    assert(min_overlap > 0.0)
    pair_overlaps_and_divs = get_pairwise_overlaps_and_divs(labels, loci)
    overlap_counts_and_mean_divs = {}
    for pair in pair_overlaps_and_divs.keys():
        overlap_div_tups = pair_overlaps_and_divs[pair]
        overlap_count = 0
        div_sum = 0.0
        for overlap, div in overlap_div_tups:
            if overlap >= min_overlap:
                overlap_count += 1
                div_sum += div
        div_mean = None
        if overlap_count > 0:
            div_mean = div_sum / overlap_count
        overlap_counts_and_mean_divs[pair] = overlap_count, div_mean
    return overlap_counts_and_mean_divs

def get_per_sample_mean_overlap(labels, pairwise_overlap_counts_and_mean_divs):
    sample_overlap_sums = {l : 0 for l in labels}
    for pair_tup in pairwise_overlap_counts_and_mean_divs.keys():
        overlap_count, div_mean = pairwise_overlap_counts_and_mean_divs[pair_tup]
        sample_overlap_sums[pair_tup[0]] += overlap_count
        sample_overlap_sums[pair_tup[1]] += overlap_count
    sample_overlap_means = {}
    for label in sample_overlap_sums.keys():
        olap_sum = sample_overlap_sums[label]
        sample_overlap_means[label] = olap_sum / float(len(labels) - 1)
    return sample_overlap_means

def get_expected_pairwise_overlap_div_table(
    ipyrad_loci_path,
    min_overlap = 0.5,
    ipyrad_json_path = None,
    stat_key = "reads_passed_filter",
):
    loci, labels = parse_ipyrad_loci_file(ipyrad_loci_path)
    read_cov = None
    if ipyrad_json_path:
        read_cov = parse_ipyrad_json_sample_stats(ipyrad_json_path, stat_key)
    overlap_counts_and_mean_divs = get_pairwise_overlap_count_and_mean_div(labels, loci, min_overlap)
    lines = []
    header = "sample_1\tsample_2\tnum_loci_shared\tmean_divergence\n"
    if read_cov:
        header = "sample_1\tsample_2\tnum_loci_shared\tmean_divergence\tsample_1_num_reads_passed_filter\tsample_2_num_reads_passed_filter\n"
    lines.append(header)
    for pair in overlap_counts_and_mean_divs.keys():
        overlap, div = overlap_counts_and_mean_divs[pair]
        line = f"{pair[0]}\t{pair[1]}\t{overlap}\t{div}\n"
        if read_cov:
            line = f"{pair[0]}\t{pair[1]}\t{overlap}\t{div}\t{read_cov[pair[0]]}\t{read_cov[pair[1]]}\n"
        lines.append(line)
    return lines

def get_expected_per_sample_overlap_table(
    ipyrad_loci_path,
    min_overlap = 0.5,
    ipyrad_json_path = None,
    stat_key = "reads_passed_filter",
):
    loci, labels = parse_ipyrad_loci_file(ipyrad_loci_path)
    read_cov = None
    if ipyrad_json_path:
        read_cov = parse_ipyrad_json_sample_stats(ipyrad_json_path, stat_key)
    overlap_counts_and_mean_divs = get_pairwise_overlap_count_and_mean_div(labels, loci, min_overlap)
    sample_overlap_means = get_per_sample_mean_overlap(labels, overlap_counts_and_mean_divs)
    lines = []
    header = "sample\tmean_num_loci_shared\n"
    if read_cov:
        header = "sample\tmean_num_loci_shared\tnum_reads_passed_filter\n"
    lines.append(header)
    for label in sample_overlap_means.keys():
        mean_overlap = sample_overlap_means[label]
        line = f"{label}\t{mean_overlap}\n"
        if read_cov:
            line = f"{label}\t{mean_overlap}\t{read_cov[label]}\n"
        lines.append(line)
    return lines

def get_missing_proportion(seq, missing = ("?", "-", "N", "n")):
    assert seq
    seq_missing = [c for c in seq if c in missing]
    return len(seq_missing) / float(len(seq))

def get_expected_missing_data_proportions(
    ipyrad_loci_path,
    missing = ("?", "-", "N", "n"),
):
    loci, labels = parse_ipyrad_loci_file(ipyrad_loci_path)
    missing_vecs = []
    for locus in loci:
        missing_props = []
        for label in labels:
            p_missing = 1.0
            if label in locus:
                seq = locus[label]
                p_missing = get_missing_proportion(seq, missing)
            missing_props.append(p_missing)
        assert len(missing_props) == len(labels)
        missing_vecs.append(missing_props)
    assert len(missing_vecs) == len(loci)
    return missing_vecs

