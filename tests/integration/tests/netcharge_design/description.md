#Netcharge integration test

Tests of the netcharge energy term.  This term is used to require that a
pose or region have a particular net charge.  It is a non-pairwise-decomposible,
yet packer-compatible, score term that provides a nonlinearly-ramping penalty for
deviation from a desired net charge.

This version tests the asymmetric case.

Expected output: a pose in which the surface has a predminantly positive charges at one
end and predominantly negative charges at the other end, with the two fading into one
another in the middle.  The net charge should be very close to zero.

Author Vikram K. Mulligan (vmullig@uw.edu)
