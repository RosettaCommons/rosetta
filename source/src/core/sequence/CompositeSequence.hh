// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CompositeSequence.hh
/// @brief class definition for a sequence profile that represents each
/// position in a sequence as a probability distribution over the allowed amino
/// acids at that position.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_CompositeSequence_hh
#define INCLUDED_core_sequence_CompositeSequence_hh

#include <core/sequence/CompositeSequence.fwd.hh>

#include <core/types.hh>
#include <core/sequence/Sequence.hh>

#include <utility/file/FileName.fwd.hh>

#include <utility/vector1.hh>

//Auto Headers
namespace core {
namespace sequence {

class CompositeSequence : public Sequence {
	typedef std::string string;

public:

	/// @brief ctors
	CompositeSequence() {}

	/// @brief copy ctor
	CompositeSequence( CompositeSequence const & src )
	: Sequence()
	{
		*this = src;
	}

	core::Size n_seqs() const;

	utility::vector1< SequenceOP > seqs() const;

	void add_sequence( SequenceOP seq );

	Size length() const override;

	void delete_position( core::Size pos ) override;

	void insert_char( core::Size pos, char new_char ) override;

	SequenceOP seq( core::Size idx ) const;

	/// @brief assignment operator.
	CompositeSequence & operator = ( CompositeSequence const & rhs ) {
		if ( this == &rhs ) return *this;

		seqs_.clear();
		for ( core::Size ii = 1; ii <= rhs.seqs().size(); ++ii ) {
			//seqs_.push_back( rhs.seq(ii)->clone() );
			add_sequence( rhs.seq(ii)->clone() );
		}

		id      ( rhs.id() );
		start   ( rhs.start() );
		gap_char( rhs.gap_char() );
		sequence( rhs.sequence() );
		return *this;
	}

	/// @brief dtor
	~CompositeSequence() override;

	/// @brief Returns an owning pointer to a new CompositeSequence object,
	/// with data that is a deep copy of the information in this object.
	SequenceOP clone() const override;

	std::string to_string() const override;

	std::string type() const override;

private:
	utility::vector1< SequenceOP > seqs_;
}; // class CompositeSequence

} // sequence
} // core

#endif
