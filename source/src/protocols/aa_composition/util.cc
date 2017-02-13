// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/aa_composition/util.cc
/// @brief  Utility functions for aa_composition-related protocols.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/aa_composition/util.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>

// Auto Headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.aa_composition.util" );

namespace protocols {
namespace aa_composition {

/// @brief Given a pose, run DSSP and populate the helices list.
/// @details Ignores helices shorter than min_length.
void
find_helices_over_length(
	core::pose::Pose const &pose,
	utility::vector1< std::pair < core::Size, core::Size > > &helices,
	core::Size const min_length
) {
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	std::string const secstruct( dssp.get_dssp_secstruct() );

	TR.Debug << "Secondary structure: " << secstruct << std::endl;

	bool inhelix(false);
	core::Size curstart(0);
	for ( core::Size i(0), imax(secstruct.length()); i<imax; ++i ) {
		if ( !inhelix ) {
			if ( secstruct[i] == 'H' ) {
				inhelix = true;
				curstart = i+1;
			}
		} else { //if inhelix
			if ( secstruct[i] != 'H' ) {
				inhelix = false;
				if ( ( i - curstart + 1 ) >= min_length ) {
					helices.push_back( std::pair< core::Size, core::Size >( curstart, i ) );
				}
			}
		}
	}
}


} //namespace aa_composition
} //namespace protocols
