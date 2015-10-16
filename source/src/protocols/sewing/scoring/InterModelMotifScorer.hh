// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file InterModelMotifScorer.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_scoring_InterModelMotifScorer_hh
#define INCLUDED_protocols_sewing_scoring_InterModelMotifScorer_hh

//Unit headers
#include <protocols/sewing/scoring/InterModelMotifScorer.fwd.hh>
#include <protocols/sewing/scoring/MotifScorer.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.hh>

//Core headers
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/chemical/ResidueTypeSet.hh>

//Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

class InterModelMotifScorer : public MotifScorer {

public:

	///@brief default construct
	InterModelMotifScorer();

	virtual ~InterModelMotifScorer(){}

	virtual
	core::Real
	score(
		AssemblyCOP assembly
	);

	core::Real
	full_motif_score(
		AssemblyCOP assembly
	);

	core::Real
	norm_motif_score(
		AssemblyCOP assembly
	);

private:

};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
