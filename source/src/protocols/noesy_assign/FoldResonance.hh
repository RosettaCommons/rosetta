// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_FoldResonance_hh
#define INCLUDED_protocols_noesy_assign_FoldResonance_hh


// Unit Headers
//#include <protocols/noesy_assign/FoldResonance.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/assert.hh>

// C++ headers
#include <string>
#include <cmath>
#include <iostream>

namespace protocols {
namespace noesy_assign {

class FoldResonance {
public:
	FoldResonance() : start_( 0 ), window_( 0 ) {}

	void set_window( core::Real start, core::Real end ) {
		start_ = start;
		window_ = end-start;
	}

	bool is_folded_down( core::Real freq ) const {
		return window_ > 0.1 && freq < start();
	}

	bool is_folded_up( core::Real freq ) const {
		return window_ > 0.1 && freq > end();
	}

	bool is_folded( core::Real freq ) const {
		return is_folded_down(freq) || is_folded_up( freq );
	}

	core::Real start() const { return start_; }
	core::Real end() const { return start_+window_; }

	core::Real operator() ( core::Real freq ) const {
		if ( !is_folded( freq ) ) return freq;
		//e return std::fmod( freq-start_, window_ ) + start_;   //modulus does not necessarily work right with negative values
		while ( is_folded_down( freq ) ) { freq+=window_; }
		while ( is_folded_up( freq ) ) { freq-=window_; }
		debug_assert( !is_folded( freq ) );
		return freq;
	}

	bool is_folded() const {
		return window_ > 0; //< BOGUS_SW-1;
	}

	void show( std::ostream& os ) const {
		if ( is_folded() ) {
			os << "FOLDED with " << window_ << " " << start() << " " << end() << std::endl;
		} else {
			os << "UNFOLDED" << std::endl;
		}
	}

private:
	core::Real start_;
	core::Real window_;
};

}
}

#endif
