Version 0.3.1
=============

Bug fixes
---------

-   Lower case nucleotide character codes were causing errors. Fixing that bug
    to accommodate lower case character codes.

Changes
-------

-   Adding option to treat "N" characters as missing data.


Version 0.3.0
=============

Changes
-------

-   Adding pyco-multi-sumevents CLI tool, which allows for plotting multiple
    nevents probability results for easy comparison.

-   Fixing issues with pyco-sumevents CLI where command line args were being
    ignored.


Version 0.2.13
==============

Changes
-------

-   Updating SumcoevolityNeventsTable to be more flexible when parsing output
    file from sumceovolity.


Version 0.2.12
==============

Changes
-------

-   Updating sumevents CLI tool. Adding '--legend-in-plot' and
    '--full-prob-axis' command line options.  The former puts the
    prior/posterior legend inside the plot box and the latter forces the full
    probability (y) axis of 0--1.


Version 0.2.11
==============

Changes
-------
-   Adding command-line tool ``msb2nex`` (and supporting code in ``parsing.py``
    module) for converting msbayes-formatted configs and data to
    ecoevolity-formatted configs and data.


Version 0.2.10
==============

Changes
-------
-   Removing matplotlib dependencies on Latex.


Version 0.2.9
=============

Changes
-------
-   Adding option to command-line tools ``loci2nex`` and ``loci2phy`` to
    randomly subsample loci.


Version 0.2.8
=============

Bug fixes
---------
-   Fixing issue with pyco-sumevents where it would crash if any Bayes factors
    were "NA".


Version 0.2.7
=============

Changes
-------
-   Adding option to command-line tools ``loci2nex`` and ``loci2phy`` to
    convert names of sequence label. A CSV file can now be provided with the
    current seq labels in the first colum and the desired labels in the second.


Version 0.2.6
=============

Changes
-------
-   Adding pyco-dummy-data command-line tool for generating dummy biallelic
    data in YAML format. Such files are often useful as templates for running
    simulations with simcoevolity.


Version 0.2.5
=============

Changes
-------
-   Updating pyco-sumsims to allow mutlple summary-table files to plot from.
-   Updating pyco-sumevents to plot via matplotlib by default; can still choose
    to use R.
-   Updating posterior classes to summarize model distances (set-partitions
    distances); this uses the Hungarian (Munkres) algorithm via the Munkres
    package.
-   Updating logic in setup script to be more robust/portable.


Version 0.2.4
=============

Changes
-------
-   Fixing highlighting bug in plots of sumsims tool (highlighting estimates
    with poor ESS values, rather than good).


Version 0.2.3
=============

Changes
-------
-   Updating sumsims tool to make plotting more robust across matplotlib
    environments.
-   Adding option to plot directly from previous summary table produced by
    sumsims.


Version 0.2.2
=============

Changes
-------
-   Updating sumsims tool to avoid LaTeX dependencies.


Version 0.2.1
=============

Changes
-------
-   Updating sumsims tool to specify a backend for matplotlib (mpl.use("agg")).
    This should avoid the invalid DISPLAY environmental variable RuntimeError
    on some systems.


Version 0.2.0
=============

Changes
-------
-   Adding sumsims CLI tool for summarizing and plotting results of ecoevolity
    analyses on datasets simulated with simcoevolity.

-   Adding formatting options to CLI plotting tools for greater control over
    output plots.

-   Adding option to output violin plots for sumsizes and sumtimes tools.

-   Adding option to plotting tools to ignore (i.e., do not plot) specified
    comparisons.

-   Adding option to PosteriorSample class to store event times in coalescent
    units.

-   Posterior classes updated to keep track of models and number of events
    separately for divergence and demographic comparisons.

-   Adding option to randomly split loci from one dataset into two datasets.
    This option is offered by CLI tools for converting pyrad loci file or fasta
    files to a nexus or phylip alignment.


Version 0.1.2
=============

Changes
-------

-   Adding CLI option to specify colors for plotting the comparison of
    divergence-time densities.

-   Adding '--charsets' option to include include a 'sets' block that defines
    locus boundaries in output nexus files.

-   Adding CLI tools for converting fasta files (rather than just .loci files
    from iPyrad). There is now a fastas2alignment tool, which converts a bunch
    of fastas into a concatenated nexus or phylip file, and a fastas2 fastas
    tool which 'normalizes' a bunch of fasta files such that they all have the
    same set of sequences.


Version 0.1.1
=============

Changes
-------
-   Adding CLI tool for converting iPyrad loci files into input files for
    PyMsBayes.


Version 0.1.0
=============

-   Initial release.
