// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file tags that go with MCA status.
/// @brief
/// @author Rhiju Das

// Unit Headers
#include <protocols/moves/MonteCarloStatus.hh>
#include <map>

namespace protocols {
namespace moves {

	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	std::string
	to_string( MCA const & mc_accepted ){

		static bool init( false );
		static std::map< MCA, std::string> mc_accepted_name;

		if ( !init ){
			mc_accepted_name[ MCA_accepted_score_beat_low ] = "accepted score beat low";
			mc_accepted_name[ MCA_accepted_score_beat_last ] = "accepted score beat last";
			mc_accepted_name[ MCA_accepted_thermally ] = "accepted thermally";
			mc_accepted_name[ MCA_rejected ] = "rejected";
			init = true;
		}

		return mc_accepted_name[ mc_accepted ];
	}


} // namespace moves
} // namespace core
