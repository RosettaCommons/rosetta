// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/LAMBEGO_IO.cc
/// @brief
/// @author James Thompson

// unit headers
#include <protocols/frag_picker/LAMBEGO_IO.hh>

// project headers
#include <core/types.hh>
#include <basic/Tracer.hh>

// utility headers
#include <utility/exit.hh>

#include <iostream>
#include <string>

namespace protocols {
namespace frag_picker {

using namespace core;

static THREAD_LOCAL basic::Tracer tracer( "protocols.frag_picker.LAMBEGO_IO" );

LAMBEGO_IO::LAMBEGO_IO() {
	bin_names_.push_back( 'L' );
	bin_names_.push_back( 'A' );
	bin_names_.push_back( 'M' );
	bin_names_.push_back( 'B' );
	bin_names_.push_back( 'E' );
	bin_names_.push_back( 'G' );
	bin_names_.push_back( 'O' );
}

void LAMBEGO_IO::read(
	std::istream & input
) {
	using core::Real;
	using std::string;
	using utility::vector1;
	using std::istringstream;

	string line;
	getline( input, line );
	tracer.Debug << "read header " << line << std::endl;

	probs_.clear();
	sequence_.clear();
	while ( getline( input, line ) ) {
		if ( line.substr(0,1) == "#" ) continue;
		istringstream line_stream( line );
		tracer.Debug << "line " << line << std::endl;

		char aa;
		Size resi;
		vector1< Real > per_residue_probs( bin_names_.size(), 0.0 );
		line_stream >> resi;
		line_stream >> aa;
		sequence_ += aa;
		for ( Size ii = 1, n_bins = bin_names_.size(); ii <= n_bins; ++ii ) {
			line_stream >> per_residue_probs[ii];
		}
		if ( line_stream.fail() ) {
			utility_exit_with_message( "Error reading in LAMBEGO_IO::read()!" );
		}
		runtime_assert( per_residue_probs.size() == bin_names_.size() );
		probs_.push_back( per_residue_probs );
	}
}

void LAMBEGO_IO::write( std::ostream & output ) {
	using core::Real;
	using utility::vector1;

	output << "resi aa ";
	for ( vector1< char >::const_iterator it = bin_names_.begin(),
			end = bin_names_.end(); it != end; ++it
			) {
		output << *it;
		if ( it + 1 != end ) output << ' ';
	}
	output << std::endl;

	Size seq_idx(1);
	for ( vector1< vector1< Real > >::const_iterator row_it = probs_.begin(),
			row_end = probs_.end(); row_it != row_end; ++row_it
			) {
		output << seq_idx << ' ' << sequence_[seq_idx-1] << ' ';
		for ( vector1< Real >::const_iterator it = row_it->begin(),
				end = row_it->end(); it != end; ++it
				) {
			output << *it;
			if ( it + 1 != end ) output << ' ';
		}
		output << std::endl;
		seq_idx++;
	}
}

utility::vector1< core::Real > const &
LAMBEGO_IO::prof_row( Size const idx ) const {
	runtime_assert( idx <= probs_.size() );
	return probs_[ idx ];
}

utility::vector1< utility::vector1< core::Real > > const & LAMBEGO_IO::matrix() const {
	return probs_;
}

Size LAMBEGO_IO::nrows() const {
	return probs_.size();
}

} // frag_picker
} // protocols
