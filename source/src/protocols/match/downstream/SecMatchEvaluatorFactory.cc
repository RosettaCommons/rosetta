// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @filesrc/protocols/match/downstream/SecMatchEvaluatorFactory.cc
/// @brief
/// @author Kui Chan (kuichan@uw.edu), oct 09


// Unit headers
#include <protocols/match/downstream/SecMatchEvaluatorFactory.hh>
#include <protocols/match/downstream/GeometrySecMatchRPE.hh>
#include <protocols/match/downstream/ScoringSecMatchRPE.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


//#include <basic/Tracer.hh>

namespace protocols {
namespace match {
namespace downstream {

//static THREAD_LOCAL basic::Tracer TR( "protocols.match.downstream.SecMatchEvaluatorFactory" );


SecMatchResiduePairEvaluatorOP
SecMatchEvaluatorFactory::create_SecMatchResiduePairEvaluatorOP(
	protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo const & mcfi,
	utility::vector1< core::Size > const & downstream_inds,
	utility::vector1< core::Size > const & upstream_inds,
	std::string SecMatchStr,
	core::pose::Pose const & upstream_pose
){
	Size found = SecMatchStr.find("SCORING_SECMATCH");
	if ( found!=std::string::npos ) {
		ScoringSecMatchRPEOP ssmOP( new ScoringSecMatchRPE ( SecMatchStr, upstream_pose ) );
		return ssmOP;
	}

	GeometrySecMatchRPEOP gsmOP( new GeometrySecMatchRPE ( mcfi, downstream_inds, upstream_inds ) );
	return gsmOP;
}


} // namespace downstream
} // namespace scoring
} // namespace core


