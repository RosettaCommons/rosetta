This integration test runs a packing and minimization trajectory with the mhc_epitope scoreterm
turned on.  Globally, the scoreterm uses the propred matrix to prevent the epitope counts from
increasing beyond native.  An MHC constraint is also applied, which uses an external NetMHCII
database file (generated specifically for this hotspot) to specifically reduce immunogenicity at
that hotspot.  Catalytic residue restrictions are applied, and design is allowed only around the
hotspot.

Author: Brahm J. Yachnin, Ph.D., Khare laboratory, Rutgers University (brahm.yachnin@rutgers.edu).
Chris Bailey-Kellogg, Ph.D., Dartmouth College (cbk@cs.dartmouth.edu)

The following manuscript is relevant to this test:
Yachnin BJ, Mulligan VK, Khare SD, and Bailey-Kellogg C.  (2021).  MHCEpitopeEnergy, a flexible Rosetta-based
biotherapeutic deimmunization platform.  J Chem Inf Model 61(5):2368-2382.  doi: 10.1021/acs.jcim.1c00056.