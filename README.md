Don't Fork Rosetta main
=======================

In our current development workflow, all work should be done in branches. Do not fork this repository.  Read [our workflow Documentation](https://wiki.rosettacommons.org/index.php/GithubWorkflow) to learn how to work with the Rosetta repositories.

We avoid forking, using a branching-only method, to ensure that all the Rosetta code stays accessible.  We would lose access to unmerged forked branches as developers leave the community.  It also ensures we have a single backup of all pushed branches in the event that GitHub were to disappear.

Rosetta main
============

Rosetta/main (the repository you are looking at) contains the Rosetta source code, database, unit tests and integration tests. The source code is located in source/src can be compiled with SCons using the following commands:

``` sh
$ cd Rosetta/main/source
$ ./scons.py -j<NumOfJobs> mode=[debug/release] [bin]
```
