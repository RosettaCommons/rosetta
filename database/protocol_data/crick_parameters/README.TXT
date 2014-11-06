This directory contains parameter files for the MakeBundle and MakeBundleHelix movers.  Parameters in these files define the minor helix in a helical bundle.

The .crick_parameters files have the following format:
--ASCII text.
--Any line starting with a pound sign (#) is ignored, and may be used to comment the file.
--Lines that are read start with r1, omega1, z1, delta_omega1, or delta_z1.
--A single omega1 line and a single z1 line define values of omega1 and z1, respectively.  (More than one such line makes no sense.)
--Lines starting with r1, delta_omega1, and delta_z1 define these values.  One of each line is included for each mainchain atom, in the order of mainchain atoms.
--Lines starting with r1, delta_omega1, and delta_z1 are followed by a tab, then a single floating-point number.

Example:

This is my file.  It is for a polymer with four mainchain atoms.
All values are completely made up.  This is a comment line.
r1 1.252525
r1 2.004426
r1 1.833133
r1 1.250311
omega1 1.431353
z1 1.809009
delta_omega1 -0.3031513
delta_omega1 0.0000000
delta_omega1 0.5698985
delta_omega1 1.0021200
delta_z1 -0.8315133
delta_z1 0.0000000
delta_z1 -1.315315
delta_z1 0.2123560
