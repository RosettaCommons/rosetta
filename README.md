Don't Fork Rosetta main
=======================

In our current development workflow, all work should be done in branches. Do not fork this repository.  Read [Our workflow Documentation](https://wiki.rosettacommons.org/index.php/GithubWorkflow) to learn how to work with the Rosetta repositories. Rosetta avoids PUBLIC forks because then we'll lose control of the Rosetta IP and/or Rosetta becomes de facto free and we can no longer afford RosettaCON.  Rosetta avoids PRIVATE forks inside RosettaCommons because that costs money to no purpose.  PRIVATE forks to approved outside groups are an exception.

Rosetta main
============

Rosetta/main contains the Rosetta source code, database, unit tests and integration tests. The source code is located in source/src can be compiled with SCons using the following commands:

``` sh
$ cd Rosetta/main/source
$ ./scons.py -j<NumOfJobs> mode=[debug/release] [bin]
```
