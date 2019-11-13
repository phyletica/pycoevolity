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
