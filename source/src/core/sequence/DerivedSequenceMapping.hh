// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_sequence_DerivedSequenceMapping_hh
#define INCLUDED_core_sequence_DerivedSequenceMapping_hh

// Unit headers
#include <core/id/SequenceMapping.hh>
#include <core/id/SequenceMapping.fwd.hh>

// Project headers
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>

namespace core {
namespace sequence {

class DerivedSequenceMapping : public core::id::SequenceMapping {

public:
// constructors, destructors and assignment operator
	/// @brief ctor
	DerivedSequenceMapping() :
		core::id::SequenceMapping(),
		seq1_(""),
		seq2_(""),
		start_seq2_(0)
	{}

	/// @brief ctor
	DerivedSequenceMapping( Size s1, Size s2 ) :
		core::id::SequenceMapping(s1,s2)
	{}

	/// @brief dtor
	~DerivedSequenceMapping();

	/// @brief copy constructor
	DerivedSequenceMapping( DerivedSequenceMapping const & src );

	DerivedSequenceMapping &
	operator=( DerivedSequenceMapping const & src );

public:

	std::string & seq1() {
		return seq1_;
	}

	std::string & seq2() {
		return seq2_;
	}

	std::string const & seq1() const {
		return seq1_;
	}

	std::string const & seq2() const {
		return seq2_;
	}

	Size const & start_seq2() const {
		return start_seq2_;
	}

	void seq1( std::string const & s ) {
		seq1_ = s;
	}

	void seq2( std::string const & s ) {
		seq2_ = s;
	}

	void start_seq2( Size s ) {
		start_seq2_ = s;
	}


private:
	std::string seq1_; //target=query sequence ... might be empty if not available
	std::string seq2_; //template sequence
 	Size start_seq2_; // which is the first seq_pos of seq2_
}; // class SequenceMapping

} // sequence
} // core

#endif
