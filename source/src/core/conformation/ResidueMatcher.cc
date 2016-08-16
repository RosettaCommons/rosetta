// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/Residue.fwd.hh
/// @author Phil Bradley


// Project headers
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/Residue.hh>


// Utility headers

// C++ headers


namespace core {
namespace conformation {

/// @details Auto-generated virtual destructor
ResidueMatcher::~ResidueMatcher() {}


bool
WatsonCrickResidueMatcher::operator()( Residue const & rsd1, Residue const & rsd2 ) const
{
	using namespace chemical;
	switch ( rsd1.aa() ) {
	case na_thy :
		return ( rsd2.aa() == na_ade );
	case na_ade :
		return ( rsd2.aa() == na_thy );
	case na_cyt :
		return ( rsd2.aa() == na_gua );
	case na_gua :
		return ( rsd2.aa() == na_cyt );
	default :
		return false;
	}
}

bool
ExactResidueMatcher::operator()( Residue const & rsd1, Residue const & rsd2 ) const
{
	return ( rsd1.name3() == rsd2.name3() );
}

} // namespace conformation
} // namespace core


