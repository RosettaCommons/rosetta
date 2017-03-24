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
Run Rosetta's MPI integration tests with a release mode MPI-enabled build and compare their output with a previous revision.
