// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ChemicalShiftSequence.hh
/// @brief class definition for a sequence profile that represents each
/// position in a sequence as a probability distribution over the allowed amino
/// acids at that position.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_ChemicalShiftSequence_hh
#define INCLUDED_core_sequence_ChemicalShiftSequence_hh

#include <core/sequence/ChemicalShiftSequence.fwd.hh>

#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>

#include <utility/file/FileName.fwd.hh>

#include <utility/vector1.hh>

namespace core {
namespace sequence {

class ChemicalShiftSequence : public SequenceProfile {
	typedef std::string string;
	typedef utility::file::FileName FileName;

public:
	/// @brief ctors
	ChemicalShiftSequence() {}

	ChemicalShiftSequence(
		FileName const & fn
	) {
		read_from_file( fn );
	}

	ChemicalShiftSequence(
		utility::vector1< utility::vector1< core::Real > > cs_prof,
		std::string const & seq,
		std::string const & ident,
		Size start_idx = 1
	) {
		profile ( cs_prof );
		sequence( seq );
		id      ( ident );
		start   ( start_idx );

		debug_assert( profile().size() == length() );
	}

	/// @brief copy ctor
	ChemicalShiftSequence( ChemicalShiftSequence const & src ):
		SequenceProfile()
	{
		*this = src;
	}

	/// @brief assignment operator.
	ChemicalShiftSequence & operator = ( ChemicalShiftSequence const & rhs ) {
		if ( this == &rhs ) return *this;

		id      ( rhs.id() );
		start   ( rhs.start() );
		gap_char( rhs.gap_char() );
		sequence( rhs.sequence() );

		profile ( rhs.profile() );
		alphabet( rhs.alphabet() );

		return *this;
	}

	/// @brief dtor
	virtual ~ChemicalShiftSequence() {}

	/// @brief Returns an owning pointer to a new ChemicalShiftSequence object,
	/// with data that is a deep copy of the information in this object.
	virtual SequenceOP clone() const {
		SequenceOP new_seq_op( new ChemicalShiftSequence( *this ) );
		return new_seq_op;
	}

	/// @brief Read an profile matrix from the given filename using the
	/// 2nd_inCS.tab format.
	void read_from_file( FileName const & fn );

	/// @brief Return the alphabet used by this sequence profile. This is an
	/// N-dimensional vector1 where N is the width of the profile, and the ith
	/// entry of any row in the profile represents the probability of ith
	/// character at that row in the sequence.
	utility::vector1< std::string > alphabet() const {
		return alphabet_;
	}

	void alphabet( utility::vector1< std::string > new_alphabet ) {
		alphabet_ = new_alphabet;
	}

	/// @brief Print this ChemicalShiftSequence object to the given std::ostream.
	friend std::ostream & operator<<(
		std::ostream & out, const ChemicalShiftSequence & p
	);

	Real sigma(
		Size /*pos*/, Size idx
	);

private:
	void check_internals_() const;
	utility::vector1< std::string > alphabet_;
}; // class ChemicalShiftSequence

} // sequence
} // core

#endif
