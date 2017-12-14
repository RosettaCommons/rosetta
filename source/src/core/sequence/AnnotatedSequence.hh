// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sequence/AnnotatedSequence.hh
/// @brief
/// @author Oliver Lange

#ifndef INCLUDED_core_sequence_AnnotatedSequence_hh
#define INCLUDED_core_sequence_AnnotatedSequence_hh

#include <utility/vector1.hh>

// C/C++ headers
#include <string>
#include <map>

// Project headers
#include <core/types.hh>

#include <core/chemical/AA.hh>


namespace core {
namespace sequence {

class AnnotatedSequence : public std::string {
public:

	AnnotatedSequence();
	AnnotatedSequence( std::string const & );
	AnnotatedSequence( AnnotatedSequence const& other );

	void operator=( std::string const& );
	AnnotatedSequence& operator=( AnnotatedSequence const& other );

	//@brief sequence position is patched
	bool is_patched( core::Size seqpos ) const;

	//@brief return the string within the [ ] or ""
	std::string patch_str( core::Size seqpos ) const;

	char one_letter( core::Size seqpos ) const;

	core::chemical::AA aa( core::Size seqpos ) const;

	std::string one_letter_sequence() const;

	core::Size length() const;

private:

	void calculate_map() const;

	typedef utility::vector1< core::Size > PositionMap;

	mutable PositionMap pos_map_;
	mutable bool map_is_clean_;
	mutable std::string one_letter_sequence_;
	mutable core::Size length_;
};

}
}

#endif
