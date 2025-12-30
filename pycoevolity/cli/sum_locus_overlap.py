#! /usr/bin/env python

"""
CLI program for summarizing sample overlap across loci.
"""

import os
import sys
import random
import argparse

import pycoevolity


def main(argv = sys.argv):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_paths',
            metavar = 'INPUT-FILE-PATH',
            nargs = '+',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Path to ipyrad locus file or paths to fasta files. '
                    'If one path is provided, it is assumed to by an ipyrad '
                    'loci file; if more then one path is provided, they are '
                    'assumed to be fasta files.'))
    parser.add_argument('-j', '--json-path',
            metavar = 'IPYRAD_JSON_PATH',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Path to ipyrad JSON file. If provided, read depth will '
                    'be summarized.'))
    parser.add_argument('-o', '--min-overlap',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_float,
            default = 0.5,
            help = ('The minimum proportion of sites of the aligned sequence '
                    'length two samples need to both have non-missing '
                    'character states.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "",
            help = ('A prefix to prepend to all output files.'))
    parser.add_argument('-f', '--force',
            action = 'store_true',
            help = ('Overwrite any existing output files. By default, an error '
                    'is thrown if an output path exists.'))
    parser.add_argument('-b', '--batch-size',
            action = 'store',
            default = 100,
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            help = ('The number of loci to process in memory.'))
    parser.add_argument('-n', '--num-procs',
            action = 'store',
            default = 1,
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            help = ('The number of processes to use.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    if args.min_overlap > 1.0:
        sys.stderr.write(
            "ERROR: min-overlap must be greater than 0 and less than or equal to 1"
        )
        sys.exit(1)

    prefix = args.prefix
    if len(prefix.split(os.path.sep)) < 2:
        prefix = os.path.join(os.curdir, prefix)

    pairwise_path = prefix + "pycoevolity-pairwise-overlap-and-div.tsv"
    per_sample_path = prefix + "pycoevolity-sample-mean-overlap.tsv"
    output_dir = os.path.dirname(pairwise_path)
    if not output_dir:
        output_dir = os.curdir
    if not args.force:
        if os.path.exists(pairwise_path):
            raise Exception(
                    "\nERROR: File {0!r} already exists.\n"
                    "Use \'-p/--prefix\' option to specify a different prefix,\n"
                    "or the \'-f/--force\' option to overwrite existing "
                    "files.".format(pairwise_path))
        if os.path.exists(per_sample_path):
            raise Exception(
                    "\nERROR: File {0!r} already exists.\n"
                    "Use \'-p/--prefix\' option to specify a different prefix,\n"
                    "or the \'-f/--force\' option to overwrite existing "
                    "files.".format(per_sample_path))

    if len(args.input_paths) == 1:
        try:
            data = pycoevolity.parsing.Loci.from_pyrad(args.input_paths[0],
                    remove_triallelic_sites = False,
                    convert_to_binary = False,
                    sequence_ids_to_remove = None,
                    label_change_map_path = None,
                    treat_n_as_missing = True)
        except:
            sys.stderr.write("ERROR: there was a problem parsing the loci file")
            raise
    else:
        try:
            data = pycoevolity.parsing.Loci.from_fastas(args.input_paths,
                    remove_triallelic_sites = False,
                    convert_to_binary = False,
                    sequence_ids_to_remove = None,
                    label_change_map_path = None,
                    treat_n_as_missing = True)
        except:
            sys.stderr.write("ERROR: there was a problem parsing the fasta files")
            raise

    sys.stderr.write("Command: {0}\n".format(" ".join(argv)))
    sys.stderr.write("\tNumber of taxa:  {0}\n".format(data.number_of_taxa))
    sys.stderr.write("\tNumber of loci:  {0}\n".format(data.number_of_loci))
    sys.stderr.write("\tNumber of sites: {0}\n".format(data.number_of_sites))

    (
        pairwise_summary,
        sample_nshared_mean,
    ) = data.get_pairwise_sample_overlap_and_div(
        
        min_overlap = args.min_overlap,
        num_processes = args.num_procs,
        num_loci_per_batch = args.batch_size,
        missing_symbols = ("?", "-", "N", "n"),
    )

    read_data = None
    if args.json_path:
        read_data = pycoevolity.parsing.parse_sample_stat_from_ipyrad_json(
            args.json_path,
            stat_key = "reads_passed_filter")

    with open(pairwise_path, "w") as out:
        if read_data:
            out.write("sample_1\tsample_2\tnum_loci_shared\tmean_divergence\tsample_1_num_reads_passed_filter\tsample_2_num_reads_passed_filter\n")
        else:
            out.write("sample_1\tsample_2\tnum_loci_shared\tmean_divergence\n")
        for sample_pair, shared_div_tup in pairwise_summary.items():
            assert len(sample_pair) == 2
            assert len(shared_div_tup) == 2
            if read_data:
                out.write(f"{sample_pair[0]}\t{sample_pair[1]}\t{shared_div_tup[0]}\t{shared_div_tup[1]}\t{read_data[sample_pair[0]]}\t{read_data[sample_pair[1]]}\n")

            else:
                out.write(f"{sample_pair[0]}\t{sample_pair[1]}\t{shared_div_tup[0]}\t{shared_div_tup[1]}\n")

    with open(per_sample_path, "w") as out:
        if read_data:
            out.write("sample\tmean_num_loci_shared\tnum_reads_passed_filter\n")
        else:
            out.write("sample\tmean_num_loci_shared\n")
        for sample, mean_num_shared in sample_nshared_mean.items():
            if read_data:
                out.write(f"{sample}\t{mean_num_shared}\t{read_data[sample]}\n")
            else:
                out.write(f"{sample}\t{mean_num_shared}\n")

    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib import gridspec
    except ImportError as e:
        sys.stderr.write('Could not import matplotlib; skipping plotting')
        sys.exit(1)

    have_scipy = False
    try:
        import scipy.stats
        have_scipy = True
    except:
        pass

    # Use TrueType (42) fonts rather than Type 3 fonts
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    plt.style.use('tableau-colorblind10')

    keys = sorted(pairwise_summary.keys())
    # Need condition to avoid divergences of `None`
    num_loci_shared = [pairwise_summary[k][0] for k in keys if pairwise_summary[k][0] > 0]
    div = [pairwise_summary[k][1] for k in keys if pairwise_summary[k][0] > 0]

    pearson_corr = None
    pearson_pval = None
    spearman_corr = None
    spearman_pval = None
    if have_scipy:
        pearson_corr, pearson_pval = scipy.stats.pearsonr(div, num_loci_shared)
        spearman_corr, spearman_pval = scipy.stats.spearmanr(div, num_loci_shared)

    plot_width = 4.0
    plot_height = 3.5
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    line, = ax.plot(div, num_loci_shared)
    plt.setp(line,
            marker = 'o',
            linestyle = '',
            markerfacecolor = 'none',
            markeredgewidth = 1.0,
            markersize = 3.5,
            rasterized = False)
    ax.set_xlabel("Pairwise mean divergence")
    ax.set_ylabel("Pairwise number of shared loci")

    if have_scipy:
        ax.text(0.99, 0.99,
                f"Pearson's corr = {pearson_corr:.2g}",
                horizontalalignment = "right",
                verticalalignment = "top",
                transform = ax.transAxes,
                zorder = 1000,
                fontsize = 'x-small',
               )
        ax.text(0.99, 0.94,
                f"Spearman's corr = {spearman_corr:.2g}",
                horizontalalignment = "right",
                verticalalignment = "top",
                transform = ax.transAxes,
                zorder = 1000,
                fontsize = 'x-small',
               )

    fig.tight_layout()
    pdf_path = prefix + "pycoevolity-pairwise-mean-div-vs-num-loci-shared.pdf"
    plt.savefig(pdf_path)

    if read_data:
        num_loci_shared = [pairwise_summary[k][0] for k in keys]
        min_num_reads = [min((read_data[k[0]], read_data[k[1]])) for k in keys]
        # mean_num_reads = [sum((read_data[k[0]], read_data[k[1]])) / 2.0 for k in keys]

        if have_scipy:
            pearson_corr, pearson_pval = scipy.stats.pearsonr(min_num_reads, num_loci_shared)
            spearman_corr, spearman_pval = scipy.stats.spearmanr(min_num_reads, num_loci_shared)

        plt.close('all')
        fig = plt.figure(figsize = (plot_width, plot_height))
        gs = gridspec.GridSpec(1, 1,
                wspace = 0.0,
                hspace = 0.0)
        ax = plt.subplot(gs[0, 0])
        line, = ax.plot(min_num_reads, num_loci_shared)
        plt.setp(line,
                marker = 'o',
                linestyle = '',
                markerfacecolor = 'none',
                markeredgewidth = 1.0,
                markersize = 3.5,
                rasterized = False)
        ax.set_xlabel("Pairwise minimum number of reads")
        ax.set_ylabel("Pairwise number of shared loci")

        if have_scipy:
            ax.text(0.99, 0.99,
                    f"Pearson's corr = {pearson_corr:.2g}",
                    horizontalalignment = "right",
                    verticalalignment = "top",
                    transform = ax.transAxes,
                    zorder = 1000,
                    fontsize = 'x-small',
                   )
            ax.text(0.99, 0.94,
                    f"Spearman's corr = {spearman_corr:.2g}",
                    horizontalalignment = "right",
                    verticalalignment = "top",
                    transform = ax.transAxes,
                    zorder = 1000,
                    fontsize = 'x-small',
                   )

        fig.tight_layout()
        pdf_path = prefix + "pycoevolity-pairwise-min-reads-vs-num-loci-shared.pdf"
        plt.savefig(pdf_path)

if __name__ == "__main__":
    main()
