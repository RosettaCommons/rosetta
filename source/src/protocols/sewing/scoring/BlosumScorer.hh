// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file BlosumScorer.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_scoring_BlosumScorer_hh
#define INCLUDED_protocols_sewing_scoring_BlosumScorer_hh

//Unit headers
#include <protocols/sewing/scoring/BlosumScorer.fwd.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.hh>

//Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

class BlosumScorer : public AssemblyScorer {

public:

	///@brief default construct
	BlosumScorer();

	virtual ~BlosumScorer(){}

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
} //sewing namespace
} //protocols namespace

#endif
