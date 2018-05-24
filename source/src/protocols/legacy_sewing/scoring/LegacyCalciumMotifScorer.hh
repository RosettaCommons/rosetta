// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LegacyCalciumMotifScorer.hh
///
/// @brief Favors interactions between all segments of pdb substructure and all other substructures
/// @author Sharon Guffy

#ifndef INCLUDED_protocols_legacy_sewing_scoring_LegacyCalciumMotifScorer_hh
#define INCLUDED_protocols_legacy_sewing_scoring_LegacyCalciumMotifScorer_hh

//Unit headers
#include <protocols/legacy_sewing/scoring/LegacyCalciumMotifScorer.fwd.hh>
#include <protocols/legacy_sewing/scoring/LegacyMotifScorer.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>

//Core headers
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/chemical/ResidueTypeSet.hh>

//Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace legacy_sewing {
namespace scoring {

class LegacyCalciumMotifScorer : public LegacyMotifScorer {

public:

	///@brief default construct
	LegacyCalciumMotifScorer();

	virtual ~LegacyCalciumMotifScorer(){}

	virtual
	core::Real
	score(
		legacy_sewing::AssemblyCOP assembly
	);
	//Does the assembly know what the first node was?



	//Don't know what these methods do; they might not be applicable
	core::Real
	full_motif_score(
		legacy_sewing::AssemblyCOP assembly
	);

	/* core::Real
	norm_motif_score(
	legacy_sewing::AssemblyCOP assembly
	);*/ //I don't think this ever even gets called (or defined???) in the regular InterModelMotifScore code

private:
	//If not, need to store it as private data & make a new constructor that stores its ID
	//check assembly->segments()
	//Does the AppendAssemblyMover always give the supplied motif the same ID?
};


} //scoring namespace
} //legacy_sewing namespace
} //protocols namespace


#endif
