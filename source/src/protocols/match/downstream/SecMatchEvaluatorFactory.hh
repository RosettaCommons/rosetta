// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/match/downstream/SecMatchEvaluatorFactory.hh
/// @brief
/// @author Kui Chan (kuichan@uw.edu), oct 09


#ifndef INCLUDED_protocols_match_downstream_SecMatchEvaluatorFactory_hh
#define INCLUDED_protocols_match_downstream_SecMatchEvaluatorFactory_hh

// Unit headers
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.fwd.hh>
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <core/types.hh>

// C++ headers
//#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

/// @brief a collection of functions making a single score_function
class SecMatchEvaluatorFactory
{
public:
	typedef core::Size Size;
	typedef core::Real Real;

public:

	static
	SecMatchResiduePairEvaluatorOP
	create_SecMatchResiduePairEvaluatorOP(
		protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo const & mcfi,
		utility::vector1< core::Size > const & downstream_inds,
		utility::vector1< core::Size > const & upstream_inds,
		std::string SecMatchStr,
		core::pose::Pose const & upstream_pose
	);


};


} // namespace downstream
} // namespace match
} // namespace protocols

#endif
