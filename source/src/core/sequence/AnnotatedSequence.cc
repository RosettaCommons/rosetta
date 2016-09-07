// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sequence/AnnotatedSequence.cc
/// @brief
/// @author Oliver Lange

// C/C++ headers
#include <string>
#include <map>

#include <core/chemical/AA.hh>
// unit headers

#include <core/sequence/AnnotatedSequence.hh>

#include <utility>
#include <utility/exit.hh>

namespace core {
namespace sequence {

AnnotatedSequence::AnnotatedSequence()
: map_is_clean_( true )
{}

AnnotatedSequence::AnnotatedSequence( std::string  str_in )
: std::string(std::move( str_in )),
	map_is_clean_( false )
{}

AnnotatedSequence::AnnotatedSequence( AnnotatedSequence const& ) = default;

AnnotatedSequence& AnnotatedSequence::operator=( AnnotatedSequence const& other ) {
	if ( this==&other ) return *this;
	std::string::operator=( other );
	map_is_clean_ = other.map_is_clean_;
	pos_map_ = other.pos_map_;
	one_letter_sequence_ = other.one_letter_sequence_;
	length_ = other.length_;
	return *this;
}

void AnnotatedSequence::operator=( std::string const& other ) {
	std::string::operator=( other );
	map_is_clean_ = false;
}

//@brief sequence position is patched
bool AnnotatedSequence::is_patched( core::Size seqpos ) const {
	if ( !map_is_clean_ ) calculate_map();
	runtime_assert( seqpos <= pos_map_.size() );
	Size apos( pos_map_[ seqpos ] );
	std::string const& annotated_seq = *this;
	if ( apos < annotated_seq.length()-1 ) {
		return annotated_seq[ apos+1 ]=='[';
	}
	return false;
}

//@brief return the string within the [ ] or ""
std::string AnnotatedSequence::patch_str( core::Size seqpos ) const {
	if ( !map_is_clean_ ) calculate_map();
	runtime_assert( seqpos <= pos_map_.size() );
	Size apos( pos_map_[ seqpos ] );
	std::string const& annotated_seq = *this;
	std::string patch_txt;
	if ( apos < annotated_seq.length()-2 ) {
		if ( annotated_seq[ apos+1 ] != '[' ) return "";
		for ( Size i = apos+2; i < annotated_seq.length() && annotated_seq[ i ]!= ']'; ++i ) {
			patch_txt = patch_txt + annotated_seq[ i ];
		}
	}
	return patch_txt;
}

char AnnotatedSequence::one_letter( core::Size seqpos ) const {
	if ( !map_is_clean_ ) calculate_map();
	runtime_assert( seqpos <= pos_map_.size() );
	Size apos( pos_map_[ seqpos ] );
	std::string const& annotated_seq = *this;
	return annotated_seq[ apos ];
}

core::chemical::AA AnnotatedSequence::aa( core::Size seqpos ) const {
	return core::chemical::aa_from_oneletter_code( one_letter( seqpos ) );
}


std::string AnnotatedSequence::one_letter_sequence() const {
	if ( !map_is_clean_ ) calculate_map();
	return one_letter_sequence_;
}

core::Size AnnotatedSequence::length() const {
	if ( !map_is_clean_ ) calculate_map();
	return one_letter_sequence_.size();
}

void AnnotatedSequence::calculate_map() const {
	if ( map_is_clean_ ) return;
	map_is_clean_ = true;

	std::string const& annotated_seq = *this;
	std::string sequence("");

	bool in_bracket = false;
	pos_map_.clear();
	pos_map_.reserve( annotated_seq.length() ); //a little more ... but avoiding resizing...
	for ( Size i = 0, ie = annotated_seq.length(); i < ie; ++i ) {
		char c = annotated_seq[i];
		if ( c == '[' ) {
			in_bracket = true;
			continue;
		} else if ( c == ']' ) {
			in_bracket = false;
			continue;
		} else {
			if ( in_bracket ) {
				continue;
			} else {
				pos_map_.push_back( i );
				sequence = sequence + c;
			}
		}
	}
	one_letter_sequence_ = sequence;
	runtime_assert( !in_bracket );
}


}
}
