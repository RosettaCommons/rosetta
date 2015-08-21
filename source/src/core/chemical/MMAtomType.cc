// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/MMAtomType.cc
/// @brief  Molecular mechanics atom type class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Project headers
#include <core/chemical/MMAtomType.hh>

// Utility headers
#include <utility/exit.hh>

namespace core {
namespace chemical {

void
MMAtomType::set_parameter( std::string const & param, Real const setting )
{
	if ( param == "LJ_RADIUS" ) {
		lj_radius_ = setting;
	} else if ( param == "LJ_WDEPTH" ) {
		lj_wdepth_ = setting;
	} else if ( param == "LJ_3B_RADIUS" ) {
		lj_three_bond_radius_ = setting;
	} else if ( param == "LJ_3B_WDEPTH" ) {
		lj_three_bond_wdepth_ = setting;
	} else {
		utility_exit_with_message( "unrecognized atomtype parameter "+param );
	}
}

} // chemical
} // core

