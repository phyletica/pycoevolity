Table of Contents
=================

 -  [Pycoevolity](#pycoevolity)
 -  [Installation](#installation)
 -  [Usage](#usage)
 -  [Acknowledgments](#acknowledgments)
 -  [License](#license)

Pycoevolity
============

A package for summarizing the output of Ecoevolity.

Installation
============

Prerequisites
-------------

For some of the plots that pycoevolity produces, there is the option of using
the
[R](https://www.r-project.org/)
packages
[ggplot2](http://ggplot2.tidyverse.org/)
and
[ggridges](https://github.com/clauswilke/ggridges).
So, if you want plotting by pycoevolity to be fully functional,
and you don't already have
[R](https://www.r-project.org/)
installed, you'll need to install it.
Once
[R](https://www.r-project.org/)
is in place, you can install the necessary packages from the
[R](https://www.r-project.org/)
prompt using:

    install.packages(c("ggplot2", "ggridges"))

However, there are python (matplotlib) options for plotting, so the R packages
are not required.

Installing pycoevolity with pip
-------------------------------

If you have the Python package manager [pip](https://pypi.org/project/pip/)
installed, you should be able to install pycoevolity using:

    pip install git+https://github.com/phyletica/pycoevolity.git@master


Manual install
--------------

If you don't have [pip](https://pypi.org/project/pip/), I recommend you install
it and use the command above to install pycoevolity.
If that is not an option, you can download and install manually:

    git clone https://github.com/phyletica/pycoevolity.git
    cd pycoevolity
    pip install .

If the last command fails with a "permission denied" error, you can install
within your home directory using:

    pip install . --user

If you use the last command, pycoevolity will be installed in a directory
within your home directory.
Specifically, the command-line tools will probably be installed to a
`.local/bin` directory within your home folder.
You will probably have to add this directory to your PATH for your shell to
find and run the pycoevolity tools.


Example usage
=============

    pyco-sumtimes --help

Acknowledgments
================

This work was supported by funding provided to Jamie Oaks from the National
Science Foundation (grant numbers DBI 1308885 and DEB 1656004).

License
=======

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

See "LICENSE.txt" for full terms and conditions of usage.
