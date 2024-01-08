Rosetta Biomolecular Modeling Library
=====================================

The Rosetta software suite includes algorithms for computational modeling and analysis of protein structures. It has enabled notable scientific advances in computational biology, including de novo protein design, enzyme design, ligand docking, and structure prediction of biological macromolecules and macromolecular complexes.

Rosetta is maintained by the RosettaCommons, a collaboration of 50+ academic research groups, who have been developing Rosetta for over 20 years. See <https://ww.rosettacommons.org> for more information about Rosetta and the RosettaCommons.

Rosetta Code
============

While the Rosetta source code is published on GitHub, it is not "Open Source" (according to the OSI definition). Most notably, use for commercial purposes requires purchase of a separate license. See LICENSE.md for further information. 

Getting Started Using Rosetta
============================

Start here: https://www.rosettacommons.org/docs/latest/getting_started/Getting-Started

Compilation quick-start (see <https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation> for more details):

``` sh
$ cd Rosetta/main/source
$ ./scons.py -j<NumOfJobs> mode=release bin
```

Questions about how to use Rosetta are best directed to the RosettaCommons forums <https://www.rosettacommons.org/forum>

PyRosetta
=========

PyRosetta are Python bindings to the Rosetta library. These can be built from the Rosetta source code. 

See <https://www.pyrosetta.org/> for more information about PyRosetta.

Developing Rosetta
==================

We welcome contributions to improve Rosetta. We use a fork-and-PR system for contribution. 
To contribute to Rosetta, please fork the Rosetta repo(s) under your own Github user space. 
You can then develop your additions in your own space. Once you're ready to contribute it back, open a PR agaist the main Rosetta repos.
You will need to sign the Rosetta Contributor License Agreement before your contribution can be accepted.

See CONTRIBUTING.md for more details.

Rosetta Code Organization
=========================

Due to its size, Rosetta uses git submodules to help in organization.

The main repository (RosettaCommons/rosetta) contains the Rosetta source code, database, unit test and integration tests
* rosetta/source/src -- The Rosetta source
* rosetta/database/ -- The Rosetta database (used during runtime)
* rosetta/source/test/ -- The compiled unit tests
* rosetta/test/integration/ -- The integration tests
* rosetta/source/bin/ -- The location of the (symlinks to) the Rosetta executables
* rosetta/source/build/ -- The location of the built libraries

Additional information is located in submodules:
* rosetta/documentation/ -- https://github.com/RosettaCommons/documentation -- Source for the online documentation
* rosetta/demos/ -- https://github.com/RosettaCommons/demos -- Various demos on using Rosetta
* rosetta/tools/ -- https://github.com/RosettaCommons/tools -- Additional helper scripts and protocols
* rosetta/rosetta_scripts_scripts -- https://github.com/RosettaCommons/rosetta_scripts_scripts -- Example XML scripts for use with RosettaScripts
* rosetta/pyrosetta_scripts -- https://github.com/RosettaCommons/pyrosetta_scripts -- Example PyRosetta scripts
* rosetta/PyRosetta.notebooks -- https://github.com/RosettaCommons/PyRosetta.notebooks -- Example Jupyter notebooks using PyRosetta.
* rosetta/source/external/... -- Various 'venderized' external dependencies

