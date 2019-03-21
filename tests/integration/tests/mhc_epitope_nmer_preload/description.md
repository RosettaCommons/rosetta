This integration test runs a packing and minimization trajectory with the mhc_epitope scoreterm
turned on.  Globally, the scoreterm uses the "nmer" SVM method (SEE NOTE) to prevent the epitope counts from
increasing beyond native.  An MHC constraint is also applied, which uses a pre-loaded NetMHCII CSV
database file (generated specifically for this hotspot) to specifically reduce immunogenicity at
that hotspot.  Catalytic residue restrictions are applied, and design is allowed only around the
hotspot.

NOTE: Currently, this is using the Propred instead of nmer, because nmer has not yet been implemented.
Once nmer is working, we will switch this over to use that predictor.

Author: Brahm J. Yachnin, Ph.D., Khare laboratory, Rutgers University (brahm.yachnin@rutgers.edu).
Chris Bailey-Kellogg, Ph.D., Dartmouth College (cbk@cs.dartmouth.edu)
