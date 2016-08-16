// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_dna_setup_hh
#define INCLUDED_core_scoring_dna_setup_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


//#include <map>

namespace core {
namespace scoring {
namespace dna {

void
set_base_partner( pose::Pose & pose );


void
find_basepairs(
	pose::Pose const & pose,
	//utility::vector1< std::pair< int, int > > & pairs,
	utility::vector1< Size > & partner
);


} // namespace dna
} // scoring
} // core

#endif
