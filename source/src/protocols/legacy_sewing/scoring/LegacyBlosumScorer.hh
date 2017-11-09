// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyBlosumScorer.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_legacy_sewing_scoring_LegacyBlosumScorer_hh
#define INCLUDED_protocols_legacy_sewing_scoring_LegacyBlosumScorer_hh

//Unit headers
#include <protocols/legacy_sewing/scoring/LegacyBlosumScorer.fwd.hh>
#include <protocols/legacy_sewing/scoring/LegacyAssemblyScorer.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>

//Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace legacy_sewing  {
namespace scoring {

class LegacyBlosumScorer : public LegacyAssemblyScorer {

public:

	///@brief default construct
	LegacyBlosumScorer();

	virtual ~LegacyBlosumScorer(){}

	virtual
	core::Real
	score( AssemblyCOP assembly ) = 0;

	// bool
	// check_blosum(
	//  AtomMap const & atom_map,
	//  std::set<SewSegment> const & reference_segments,
	//  Model const & mobile_model
	// );

private:

};


} //scoring namespace
} //legacy_sewing namespace
} //protocols namespace

#endif
