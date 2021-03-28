// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaMultipleSequenceAlignment.cc
/// @brief A class to store multiple sequence alignment data for input into trRosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef USE_TENSORFLOW

// Project headers:
#include <protocols/trRosetta/trRosettaMultipleSequenceAlignment.hh>

// Core headers:
#include <core/types.hh>

// Basic headers:
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.trRosetta.trRosettaMultipleSequenceAlignment" );


namespace protocols {
namespace trRosetta {

/// @brief Options constructor: initialize this object from the contents
/// of a multiple sequence alignment.
trRosettaMultipleSequenceAlignment::trRosettaMultipleSequenceAlignment(
	std::string const & file_contents
) {
	initialize_from_file_contents( file_contents );
}

/// @brief Destructor.
trRosettaMultipleSequenceAlignment::~trRosettaMultipleSequenceAlignment() = default;

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
trRosettaMultipleSequenceAlignmentOP
trRosettaMultipleSequenceAlignment::clone() const {
	return utility::pointer::make_shared< trRosettaMultipleSequenceAlignment >( *this );
}

////////////////////////////////
// PUBLIC MEMBER FUNCTIONS
////////////////////////////////

/// @brief Return the contents as input to a trRosetta model.
utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< int32_t > > 
trRosettaMultipleSequenceAlignment::construct_input_tensors() const {
	debug_assert( input_tensors_.size() == 3 );
	return input_tensors_;
}

////////////////////////////////
// PRIVATE MEMBER FUNCTIONS
////////////////////////////////

/// @brief Initialize this object from the contents of a multiple sequence alignment file.
void
trRosettaMultipleSequenceAlignment::initialize_from_file_contents(
	std::string const & file_contents
) {
	using namespace basic::tensorflow_manager;

	debug_assert(input_tensors_.size() == 0); //Should not already be initialized.

	utility::vector1< std::string > const lines( utility::string_split(file_contents, '\n') );
	utility::vector1< std::string > actual_sequences;
	actual_sequences.reserve( lines.size() );

	std::string const allowed_chars( "ARNDCQEGHILKMFPSTWYV-" ); //Looks like it was alphabetical order of three-letter codes in the original Python script.

	// First, find the minimum sequence length and count sequences.
	core::Size minlen(0);
	core::Size seqcount(0);
	for( std::string linestripped /*Deliberately copied*/ : lines) {
		if(linestripped.empty()) continue;
		remove_lcase( utility::rstrip_whitespace( linestripped ) ); //Strip only trailing whitespace and remove lowercase characters.
		if( linestripped.empty() || linestripped[0] == '>' || linestripped[0] == '#' ) {
			continue; //Skip empty lines or annotation lines.
		}
		if( minlen == 0 || linestripped.size() < minlen ) {
			minlen = linestripped.size();
		}
		++seqcount;
	}

	// Second, allocate storage space:
	input_tensors_.emplace_back( RosettaTensorflowTensorContainer< int32_t >(utility::vector1<int64_t>{}, 0) );
	input_tensors_.emplace_back( RosettaTensorflowTensorContainer< int32_t >(utility::vector1<int64_t>{}, 0) );
	input_tensors_.emplace_back( RosettaTensorflowTensorContainer< int32_t >(utility::vector1<int64_t>{ static_cast<int64_t>(seqcount), static_cast<int64_t>(minlen) }, 0) );

	// Third, store data:
	input_tensors_[1](1) = seqcount;
	input_tensors_[2](1) = minlen;
	core::Size seqcount2(0); 
	for( std::string linestripped /*Deliberately copied*/ : lines) {
		if(linestripped.empty()) continue;
		remove_lcase( utility::rstrip_whitespace( linestripped ) ); //Strip only trailing whitespace and remove lowercase characters.
		if( linestripped.empty() || linestripped[0] == '>' ) {
			continue; //Skip empty lines or annotation lines.
		}

		++seqcount2;
		for( core::Size i(0); i<minlen; ++i ) {
			core::Size char_as_int( allowed_chars.find( linestripped[i] ) );
			if( char_as_int > 20 ) {
				char_as_int = 20; //Treat all unknown characters as gap characters.
			}
			input_tensors_[3](seqcount2, i+1) = char_as_int;
		}
	}
	debug_assert(seqcount2 == seqcount); //Should be true.
}

/// @brief Remove the lowercase characters (i.e. delete them and ligate the flanking regions) from a string.
/// @details Replaces the input.
std::string &
trRosettaMultipleSequenceAlignment::remove_lcase(
	std::string & s
) const {
	if(s.empty()) return s;

	std::string o;
	for( core::Size i(0), imax(s.size()); i<imax; ++i ) {
		if( !std::islower(s[i]) ) {
			o += s[i];
		}
	}
	s = o;
	return s;	
}

} //trRosetta
} //protocols

#endif //USE_TENSORFLOW
