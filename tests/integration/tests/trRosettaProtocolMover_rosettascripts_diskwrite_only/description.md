# trRosettaProtocolMover\_rosettascripts\_diskwrite\_only

## Author

Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org).

## Description

This test ensures that the trRosettaProtocolMover can be used within RosettaScripts.  This test uses the mover to write constraints to disk, skipping the structure prediction phase.  For the non-Tensorflow build of Rosetta, this test only tests whether the -info flag works.
