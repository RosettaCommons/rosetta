// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_coarse_Coarse.hh
#define INCLUDED_core_coarse_Coarse.hh

/*
different design ideas for Coarse/Centroid residues:


A LiteResidue could also contain the information about its resolution
level (int / enum)
0 = full-atom
2 = custom1
3 = custom2

last = centroid
...
but otherwise just contains the atoms/beads which it has in its current resolution.
The aa-names would stay untouched

a coarsifier class would be used to translate between different levels
of resolution

 to lower resolution -> directly compute new bead centers
 to higher resolution --> do repacking

 are resolution always in such a hirachie?
 should one do rotamer-trials or full-repack for going up in the hirachie?

the score funtion might query the resolution level to determine how to
score things.

how would the interface/class hirachy of the coarsifier look like ?


*****************
only one level of coars-graining is supported at any given time:
--> residues do not need to keep track of their coarse graining level

to switch the level
the a new residue-set (the platonic-types) is created (or already there)
old residues are replaced by new ones,
the respective parts of scoring functions get replaced and a new rotamerlibrary
is plugged in.

The coarse-grainer must be able to do the following things:
take a full-atom (or high-level) rotamerlibrary and produce a coarse-one
take full-atom residues and make coarse ones.


evtl. new classes have to be derived from rotamer-library and score-functions
to get the new logic in there if necessary.
for instance, sequence dependence of vdw energy or
angle dofs in rotamer packing.

for starters there will be fixed angles...




 *************
new residue types
*********
the full-atom lite-residue carries also the
other resolutions with it. --> this might give us problems
with atom-trees etc.


****************
The translation-process:

A Meta-Translator will contain the rules how to translate
One Meta-translator contains one rule set.
It is matched to the residue type by aa(). Multiple translators might
be matched to one aa(). In this case there should be a default Translator which is overloaded by
specializations that apply to special modifications.

The Meta-Translator will be used to create Translator objects.
There will be exactly one Translator per Residue-Type in the ResidueSet.
It contains a mapping bead index <--> atom_indices in Residue-Type

The container of Translators should be deleted after use. (since the residue-types might change, etc etc)
or we need mechanism to notify if residue-types get atoms added.

A coarse-atom-iterator as in pose_coarse.cc would use the mapping encoded by the Translators...

*******************



*/

#endif // INCLUDED_core_coarse_Coarse.hh
