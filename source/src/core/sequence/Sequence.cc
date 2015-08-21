// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Sequence.cc
/// @brief method definitions for Sequence class
/// @author James Thompson

// Unit headers
#include <core/sequence/Sequence.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility/exit.hh>

#include <iostream>
#include <string>

#include <utility/vector1.fwd.hh>
#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

using core::Size;
using std::string;
using utility::vector1;

Sequence::Sequence( core::pose::Pose const & pose ):
	start_(1),
	gap_char_('-')
{
	if ( ! pose.pdb_info() ) {
		id_ = pose.pdb_info()->name();
	} else {
		id_ = "unknown";
	}
	sequence( pose.sequence() );
}

Sequence::~Sequence() {}

SequenceOP Sequence::clone() const {
	SequenceOP new_seq_op( new Sequence( *this ) );
	return new_seq_op;
}

/// @brief initializes this sequence object from a file.
void Sequence::read_from_file( utility::file::FileName const & /*fn*/ ) {
	utility_exit_with_message(
		"Error: class doesn't define method read_from_file!"
	);
}


void Sequence::sequence( string sequence ) {
	seq_ = sequence;
}

void Sequence::id( std::string new_id ) {
	id_ = new_id;
}

void Sequence::start( core::Size new_start ) {
	start_ = new_start;
}

void Sequence::gap_char( char new_gap_char ) {
	gap_char_ = new_gap_char;
}

core::Size Sequence::start() const {
	return start_;
}

std::string Sequence::id() const {
	return id_;
}

char Sequence::gap_char() const {
	return gap_char_;
}

Size Sequence::length() const {
	return seq_.length();
}

Size Sequence::ungapped_length() const {
	return ungapped_sequence().length();
} // ungapped_length

string Sequence::ungapped_sequence() const {
	string ungapped("");
	for ( core::Size i = 1; i <= length(); ++i ) {
		if ( !is_gap(i) ) ungapped += (*this)[i];
	}

	return ungapped;
}

string Sequence::sequence() const {
	return seq_;
}

char Sequence::operator[]( core::Size pos ) const {
	runtime_assert( pos > 0 );
	runtime_assert( pos <= seq_.size() );
	return seq_.at( pos - 1 ); // strings are indexed by zero
}

char Sequence::at( core::Size pos ) const {
	//runtime_assert( pos > 0 );
	//runtime_assert( pos <= seq_.size() );
	//return seq_.at( pos - 1 ); // strings are indexed by zero
	return( (*this)[pos] );
}

void Sequence::insert_char( Size pos, char new_char ) {
	runtime_assert( pos <= length() + 1 );

	std::string new_seq( "" );
	for ( Size i = 0; i <= length() + 1; ++i ) {
		if ( i == pos )                new_seq += new_char;
		if ( i >= 1 && i <= length() ) new_seq += (*this)[i];
	}

	sequence( new_seq );
} // insert_char

void Sequence::delete_position( core::Size pos ) {

	std::string new_seq( "" );
	for ( core::Size i = 1; i <= length(); ++i ) {
		if ( i != pos ) new_seq += (*this)[i];
	}

	sequence( new_seq );
}

void Sequence::insert_gap( core::Size pos ) {
	insert_char( pos, gap_char() );
}

void Sequence::append_char( char new_char ) {
	seq_ += new_char;
} // insert_char

void Sequence::append_gap() {
	seq_ += gap_char();
}

bool Sequence::is_gap( core::Size pos ) const {
	// define anything outside the length of the sequence as a gap.
	if ( pos < 1 || pos > length() ) return true;

	return ( (*this)[pos] == gap_char_ );
}

core::Size Sequence::resnum( core::Size idx ) const {
	runtime_assert( idx <= length() );
	runtime_assert( idx > 0 );

	if ( is_gap( idx ) ) return 0;

	core::Size num( start() );
	for ( core::Size i = 1; i < idx; ++i ) {
		if ( !is_gap(i) ) ++num;
	}

	return num;
}

void Sequence::read_data( std::istream & in ) {
	std::string seq, name;
	core::Size begin;

	in >> name >> begin >> seq;

	id      ( name  );
	start   ( begin );
	sequence( seq   );
}

std::ostream & operator<<(
	std::ostream & out,
	const Sequence & seq
) {
	out << seq.to_string();
	// out << std::endl;
	return out;
} // operator <<

std::istream & operator>> (
	std::istream & in,  Sequence & seq
) {
	seq.read_data( in );
	return in;
}

std::string Sequence::to_string() const {
	std::string retval("");
	Size const id_width( 20 );
	Size const start_width( 8 );

	retval += ObjexxFCL::format::A( id_width, id() );
	retval += ObjexxFCL::format::I( start_width, start() );
	retval += ' ' + sequence();
	return retval;
}

std::string Sequence::type() const {
	return "sequence";
}

} // sequence
} // core
