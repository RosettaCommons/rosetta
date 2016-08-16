// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/Loops.tmpl.hh
/// @brief
/// @author Brian D. Weitzner

#ifndef INCLUDED_protocols_loops_Loops_TMPL_HH
#define INCLUDED_protocols_loops_Loops_TMPL_HH

// Unit header
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace loops {

/// @brief set each loop-residue in the vector to val.
/// input vector of nres length ( if shorter last residues of loop are ignored )
template< class T >
void Loops::transfer_to_residue_vector( utility::vector1< T > & vector, T val ) const {
	core::Size nres = vector.size();
	for ( const_iterator it = begin(); it  != end(); ++it ) {
		if ( it->start() <= nres ) {
			for ( core::Size pos = it->start(); pos <= std::min( it->stop(), nres ); pos++ ) {
				vector[ pos ] = val;
			}
		}
	}
}

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_Loops_TMPL_HH
