# integration test suite
These tests are mainly focused on running Rosetta and making sure the results do not change unexpectedly.

-----
### (plain)
Run Rosetta's integration tests in release mode and compare their output with a previous revision.

-----
### debug
Run Rosetta's integration tests in debug mode and compare their output with a previous revision.

-----
### release_debug
Run Rosetta's integration tests in release_debug mode and compare their output with a previous revision.

-----
### release_debug_with_symbols
Run Rosetta's integration tests in release_debug mode and compare their output with a previous revision.
Has more informative diagnostic output than the plain release_debug mode.

-----
### mpi
Run Rosetta's MPI integration tests (command.mpi) with a release_debug mode MPI-enabled build and compare their output with a previous revision.

-----
### pytorch
Run Rosetta's Pytorch integration tests (command.pytorch) with a debug mode Pytorch-enabled build and compare their output with a previous revision.

-----
### tensorflow
Run Rosetta's Tensorflow integration tests (command.tensorflow) with a debug mode Tensorflow-enabled build and compare their output with a previous revision.

-----
### thread
Run Rosetta's threading integration tests (command.thread) with a debug mode threading-enabled build and compare their output with a previous revision.

-----
### bcl
Run Rosetta's BCL integration tests (command.bcl) with a debug mode threading-enabled build and compare their output with a previous revision.
