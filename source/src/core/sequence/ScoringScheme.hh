// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ScoringScheme.hh
/// @brief abstract base class for representing scoring schemes for alignments.
/// @details ScoringScheme objects contain two scoring rules:
/// - a rule for comparing any character from a sequence to a gap character (usually
/// liste as gap insertion and gap extension, or d and e)
/// - a rule for scoring two elements from a sequence (simply S)
/// Generally the rule for scoring gaps is composed of a gap insertion and a gap
/// extension parameter. It's important to note that alignments derived using
/// the ScoringScheme and Aligner objects are only guaranteed to be optimal if
/// every element of S is bigger than -2 * e.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_ScoringScheme_hh
#define INCLUDED_core_sequence_ScoringScheme_hh

// Unit headers
#include <core/sequence/ScoringScheme.fwd.hh>

#include <core/types.hh>
#include <core/sequence/Sequence.fwd.hh>

#include <utility/io/izstream.fwd.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>

namespace core {
namespace sequence {

class ScoringScheme : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor
	ScoringScheme();

	/// @brief dtor
	virtual ~ScoringScheme();

	/// @brief clone method.
	virtual ScoringSchemeOP clone() const = 0;

	/// @brief Initialize from a file.
	virtual void read_from_file( utility::file::FileName const & /*fn*/ );

	virtual void read_data( utility::io::izstream & /*input*/ );

	/// @brief Gets the gap opening penalty.
	virtual Real gap_open() const;

	/// @brief Gets the gap extension penalty.
	virtual Real gap_extend() const;

	/// @brief Sets the gap opening penalty.
	void gap_open( Real const gap_open );

	/// @brief Sets the gap extension penalty.
	void gap_extend( Real const gap_extend );

	/// @brief getters for type, which is a unique string name for this object.
	std::string type() const;

	/// @brief getters for type, which is a unique string name for this object.
	void type( std::string new_type );

	virtual Real score(
		SequenceOP seq1,
		SequenceOP seq2,
		Size pos1,
		Size pos2
	) = 0;

	/// @brief Utility method for producing useful error messages and exiting
	/// from program. Declared const which is funny, because exiting the program
	/// certainly changes the state of this object! This might be replaced with
	/// exception handling if we ever start using those.
	void unimplemented_method_error( std::string const & method_name ) const;

	bool is_good(
		Real const & num
	);

private:
	Real gap_open_;
	Real gap_extend_;
	std::string type_;
}; // class ScoringScheme

} // sequence
} // core

#endif
