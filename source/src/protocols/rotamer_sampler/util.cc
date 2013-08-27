// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/util.cc
/// @brief Useful functions for rotamer generator.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/util.hh>

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

using namespace core;
static basic::Tracer TR( "protocols.rotamer_sampler.util" );

namespace protocols {
namespace rotamer_sampler {

void add_values_from_center(
	utility::vector1<core::Real> & torsions,
	Real const center,
	Real const max_range,
	Real const bin_size
) {
	Real const epsilon = 0.0001; //A small number to avoid numerical error for 'Real'.
	for ( Real offset = 0; offset <= max_range + epsilon; offset += bin_size ) {
		torsions.push_back( center + offset );
		if ( offset != 0 ) torsions.push_back( center - offset );
	}
	std::sort( torsions.begin(), torsions.end() );
}

}
}
