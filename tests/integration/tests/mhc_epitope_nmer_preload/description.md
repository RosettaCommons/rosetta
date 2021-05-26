This integration test runs a packing and minimization trajectory with the mhc_epitope scoreterm
turned on.  An MHC constraint is applied, which uses a pre-loaded NetMHCII CSV database file
(generated specifically for this hotspot) to specifically reduce immunogenicity at
that hotspot.  Catalytic residue restrictions are applied, and design is allowed only around the
hotspot.

After execution, the final scoring is done with nmer-based scores that cover the entire pose in the OUTPUT
tag.  This is not used during packing as it is very slow.  If faster SVMs are added in the future, it would
be worthwhile to add this to the packing trajectory as well.

Author: Brahm J. Yachnin, Ph.D., Khare laboratory, Rutgers University (brahm.yachnin@rutgers.edu).
Chris Bailey-Kellogg, Ph.D., Dartmouth College (cbk@cs.dartmouth.edu)

The following manuscript is relevant to this test:
Yachnin BJ, Mulligan VK, Khare SD, and Bailey-Kellogg C.  (2021).  MHCEpitopeEnergy, a flexible Rosetta-based
biotherapeutic deimmunization platform.  J Chem Inf Model 61(5):2368-2382.  doi: 10.1021/acs.jcim.1c00056.