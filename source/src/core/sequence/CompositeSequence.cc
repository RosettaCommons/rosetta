// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CompositeSequence.cc
/// @author James Thompson

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/CompositeSequence.hh>

#include <utility/pointer/owning_ptr.hh>


#include <string>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>


namespace core {
namespace sequence {

static basic::Tracer tr( "core.sequence.CompositeSequence" );

std::ostream & operator<<(
	std::ostream & out, const CompositeSequence & p
) {
	out << p.to_string();
	return out;
}

core::Size CompositeSequence::n_seqs() const {
	return seqs_.size();
}

utility::vector1< SequenceOP > CompositeSequence::seqs() const {
	return seqs_;
}

void CompositeSequence::add_sequence( SequenceOP seq ) {
	if ( seqs_.size() > 0 ) runtime_assert( length() == seq->length() );
	id      (seq->id()      );
	gap_char(seq->gap_char());
	sequence(seq->sequence());
	seqs_.push_back(seq);
}

Size CompositeSequence::length() const {
	if ( seqs_.size() == 0 ) return 0;
	return seqs_[1]->length();
}

void CompositeSequence::delete_position( core::Size pos ) {
	for ( Size ii = 1; ii <= n_seqs(); ++ii ) {
		seqs_[ii]->delete_position(pos);
	}
	Sequence::delete_position( pos );
}

void CompositeSequence::insert_char( core::Size pos, char new_char ) {
	Sequence::insert_char( pos, new_char );
	for ( Size ii = 1; ii <= n_seqs(); ++ii ) {
		seqs_[ii]->insert_char(pos,new_char);
	}
}

SequenceOP CompositeSequence::seq( core::Size idx ) const {
	//runtime_assert( idx <= seqs_.size() );
	return seqs_[idx];
}

/// @brief dtor
CompositeSequence::~CompositeSequence() = default;

/// @brief Returns an owning pointer to a new CompositeSequence object,
/// with data that is a deep copy of the information in this object.
SequenceOP CompositeSequence::clone() const {
	SequenceOP new_seq_op( new CompositeSequence( *this ) );
	return new_seq_op;
}

std::string CompositeSequence::to_string() const {
	using ObjexxFCL::string_of;
	std::string retval("");
	retval = Sequence::to_string();
	//retval += "*** IN CompositeSequence::to_string() ***\n";
	//retval += "*** print out " + string_of(n_seqs()) + " sequences\n";
	//retval += "*** base sequence " + sequence() + "\n";
	//for ( Size ii = 1; ii <= n_seqs(); ++ii ) {
	// retval += seqs_[ii]->to_string() + "\n";
	//}
	return retval;
	//retval += "*** DONE WITH CompositeSequence::to_string() ***\n";
}

std::string CompositeSequence::type() const {
	return "composite_sequence";
}

} // sequence
} // core
