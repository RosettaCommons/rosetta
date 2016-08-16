// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/LAMBEGO_IO.hh
/// @brief reader for a probability distribution over A,B,E,G,O torsion bins.
/// @author James Thompson

#ifndef INCLUDED_protocols_frag_picker_LAMBEGO_IO_hh
#define INCLUDED_protocols_frag_picker_LAMBEGO_IO_hh

// Package headers
#include <core/types.hh>
#include <iostream>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace frag_picker {

class LAMBEGO_IO {
public:

	LAMBEGO_IO();

	void read( std::istream & input );
	void write( std::ostream & output );

	utility::vector1< core::Real > const & prof_row( core::Size const idx ) const;

	utility::vector1< utility::vector1< core::Real > > const & matrix() const;

	core::Size nrows() const;

private:
	std::string sequence_;
	utility::vector1< utility::vector1< core::Real > > probs_;
	utility::vector1< char > bin_names_;
}; // LAMBEGO_IO

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_LAMBEGO_IO_HH */
