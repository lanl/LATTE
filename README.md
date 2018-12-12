[![Build Status](https://travis-ci.org/lanl/LATTE.svg?branch=master)](https://travis-ci.org/lanl/LATTE)
[![codecov](https://codecov.io/gh/lanl/LATTE/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl/LATTE)

# LATTE

Open source density functional tight binding molecular dynamics.

#LA-CC-10-004

This material was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, this program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2.0 of the License. Accordingly, this program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


# Compiling LATTE (Makefile):

Enter in the LATTE directory (`cd ~/LATTE`) and
modify the `makefile.CHOICES` file according to your operative system,
libraries, compiler, etc. There are several examples of
`makefile.CHOICES` files inside the makefiles directory
(`~/LATTE/makefiles`) that could be used as a template to replace the
`makefile.CHOICES` that is located in the main directory.

To build the code just type make inside the main directory as follows:

          $ cd; cd ~/LATTE
          $ make; make test


The latter should build the code and test it with some examples that are
located in `~/LATTE/tests`.

# Compiling LATTE (CMake):

To build LATTE with CMake:

```
$ cd; cd ~/LATTE
$ mkdir build
$ cd build
$ cmake ../cmake -D<OPTIONS>=<VALUE>
$ make
```

Some useful options are `OPENMP` (`ON` or `OFF`), `DO_MPI` and `PROGRESS`.
For a full list of options use a CMake GUI, e.g. `ccache` or `cache-gui`.


# Authors

Nicolas Bock (T-1)
Marc Cawkwell (T-1)
Josh D. Coe (T-1)
Aditi Krishnapriyan (T-1 & Stanford)
Matthew P. Kroonblawd (T-1 & Mizzou)
Adam Lang (T-1)
Chang Liu (T-1 & NCSU)
Enrique Martinez Saez (MST-8)
Susan M. Mniszewski (CCS-3)
Christian F. A. Negre (T-1)
Anders M. N. Niklasson (T-1)
Edward Sanville (T-1)
Mitchell A. Wood (T-1 &  Purdue)
Ping Yang (T-1)

Los Alamos National Laboratory


# Citing

To cite the code, please proceed as follows:

[![DOI](https://zenodo.org/badge/75976231.svg)](https://zenodo.org/badge/latestdoi/75976231)

with the following `bibtex` citation:

    @misc{LATTE,
      title = {LATTE},
      url = {https://github.com/lanl/LATTE},
      author = {N. Bock and M. J. Cawkwell and J. D. Coe and A. Krishnapriyan and M. P. Kroonblawd and A. Lang and and C. Liu and E. Martinez Saez and S. M. Mniszewski and C. F. A. Negre and A. M. N. Niklasson and E. Sanville and M. A. Wood and P. Yang},
      year = {2008}
    }

# Contacts

Please subscribe to our [mailing list](https://lanl.github.io/LATTE/) to
request information or send us feedback.
