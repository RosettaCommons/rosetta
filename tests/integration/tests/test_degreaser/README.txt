Two integration tests for the Degreaser, AKA the SecretionOptimizationMover.

First integration test is in 'symmetric mode,' which takes in a symdef file and passes it to SetupForSymmetryMover
Second integration test is in 'asymmetric mode,' and is symmetry unaware, but is removed to have this test take less than 30 seconds to run. enable it if you're curious to see if for some reason the asymmetric mode would break ?

Both should be 'deterministic' in that they should converge on the positions to mutate but not necessarily on amino acid identity. Small differences in repacking could also account for changes.

Author: John Wang, Neil King lab, IPD, UW Seattle

