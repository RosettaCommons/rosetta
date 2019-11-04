Integration test for invoking the multithreaded packer setup in the fixbb application.

Note that this test only runs for the multithreaded build of Rosetta.  It tests three things:

- Packing when as many threads are requested as are available.
- Packing when more threads are requested than are available.
- Packing when only one thread (the calling thread) is available.

The test also checks that the output pdb files are identical in each case.  We're using the Trp cage miniprotein as our
design test case.

Author: Vikram K. Mulligan, Ph.D., (vmulligan@flatironinstitute.org), Center for Computational Biology, Flatiron Insitute. 

