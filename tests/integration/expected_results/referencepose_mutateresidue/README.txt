Integration test for the MutateResidue mover working with ReferencePose objects stored by the StorePoseSnapshot mover.  Changes to this test could mean:
-- The reference pose machinery in core::pose::reference_pose isn't working properly or has changed (though this is covered by unit tests).
-- The links between the core::pose::Pose class, the core::pose::reference_pose::ReferencePoseSet class, and/or the core::pose::reference_pose::ReferencePose class have changed or broken (though this is covered by unit tests).
-- The StorePoseSnapshot mover is not working properly.
-- The MutateResidue mover is not working properly, or is not parsing reference pose strings properly.
-- Code in core/pose/sequence.cc that parses residue strings (Rosetta numbering, PDB numbering, or reference pose numbering) is not parsing these strings properly (though this is covered by unit tests).

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

