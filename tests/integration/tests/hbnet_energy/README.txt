Tests of the hbnet energy term with design.   This energy term is intended
to approximate what Scott's HBNet protocol does: find side-chain rotamers
that yield good hydrogen bond networks.  Unlike the HBNet protocol, the
hbnet energy is a non-pairwise-decomposible score term that works with the
packer to guide it to solutions with hydrogen bond networks.  This allows
simultaneous optimization of hydrogen bond networks and all of the other
things (hydrophobic interactions, Van der Waals packing, etc.) that the
packer optimizes.

This version tests the asymmetric case.

Author Vikram K. Mulligan (vmullig@uw.edu)
