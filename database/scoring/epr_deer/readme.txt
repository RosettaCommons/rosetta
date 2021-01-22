DEFAULT SPIN LABELS

These coordinates describe the electron position on a dynamic MTSSL side chain in a coarse-grained manner.
The XYZ coordinates are each assigned a weight (second member of the pair) to characterize a given electron's
propensity to occupy that given area. This works as well as existing methods without relying on more expensive
fullatom representations. Clashes are evaluated using a dummy atom (see below). Custom positions may also be
provided. These are introduced using a homogeneous transform object - see corresponding cc file.
The larger vector is more accurate but about 4-fold slower. For de novo folding, users are encouraged to use
the spin label DEFAULT_FAST (the first one), whereas for things like comparative modeling DEFAULT (the second)
is more appropriate.