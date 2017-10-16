// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pockets/util.cc
/// @brief protocols::pockets::util functions
/// @author David Johnson
/// @author Ragul Gowthaman

#include <protocols/pockets/CCluster.hh>
#include <protocols/pockets/PCluster.hh>

// Core Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>

#include <algorithm>


namespace protocols {
namespace pockets {

using namespace core;

void min_values_versus_cluster( Size & x, Size & y, Size & z, PCluster const & c2 ) {
	x = std::min( x, c2.get_minX() );
	y = std::min( y, c2.get_minY() );
	z = std::min( z, c2.get_minZ() );
}
void max_values_versus_cluster( Size & x, Size & y, Size & z, PCluster const & c2 ) {
	x = std::max( x, c2.get_maxX() );
	y = std::max( y, c2.get_maxY() );
	z = std::max( z, c2.get_maxZ() );
}

void min_values_versus_cluster( Size & x, Size & y, Size & z, CCluster const & c2 ) {
	x = std::min( x, c2.get_minX() );
	y = std::min( y, c2.get_minY() );
	z = std::min( z, c2.get_minZ() );
}
void max_values_versus_cluster( Size & x, Size & y, Size & z, CCluster const & c2 ) {
	x = std::max( x, c2.get_maxX() );
	y = std::max( y, c2.get_maxY() );
	z = std::max( z, c2.get_maxZ() );
}

// AMW TODO: refactor so as not to rely on options system reads on-the-fly
Size
counting_atoms_in_residue( core::conformation::Residue const & rsd ) {
	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::include_hydrogens ]() ) {
		return rsd.natoms();
	} else {
		return rsd.nheavyatoms();
	}
}


} // pockets
} // protocols

