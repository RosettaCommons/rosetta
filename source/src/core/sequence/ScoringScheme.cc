// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ScoringScheme.cc
/// @brief abstract base class for representing scoring schemes for alignments.

// Unit headers
#include <core/sequence/ScoringScheme.fwd.hh>

#include <core/types.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/ScoringScheme.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>
#include <complex>

namespace core {
namespace sequence {

/// @brief ctor
ScoringScheme::ScoringScheme() = default;

/// @brief dtor
ScoringScheme::~ScoringScheme() = default;

/// @brief Initialize from a file.
void ScoringScheme::read_from_file(
	utility::file::FileName const & /*fn*/
) {
	unimplemented_method_error( "read_from_file" );
}

void ScoringScheme::read_data( utility::io::izstream & /*input*/ ) {
	unimplemented_method_error( "read_data" );
}

/// @brief Gets the gap opening penalty.
Real ScoringScheme::gap_open() const {
	return gap_open_;
}

/// @brief Gets the gap extension penalty.
Real ScoringScheme::gap_extend() const {
	return gap_extend_;
}

/// @brief Sets the gap opening penalty.
void ScoringScheme::gap_open( Real const gap_open ) {
	gap_open_ = gap_open;
}

/// @brief Sets the gap extension penalty.
void ScoringScheme::gap_extend( Real const gap_extend ) {
	gap_extend_ = gap_extend;
}

/// @brief getters for type, which is a unique string name for this object.
std::string ScoringScheme::type() const {
	return type_;
}

/// @brief getters for type, which is a unique string name for this object.
void ScoringScheme::type( std::string new_type ) {
	type_ = new_type;
}

/// @brief Utility method for producing useful error messages and exiting
/// from program. Declared const which is funny, because exiting the program
/// certainly changes the state of this object! This might be replaced with
/// exception handling if we ever start using those.
void ScoringScheme::unimplemented_method_error(
	std::string const & method_name
) const {
	utility_exit_with_message(
		"Called ScoringScheme::" + method_name + " method from derived class " +
		type() + "," + "ended up in ScoringScheme::" + method_name + "\n"
	);
}

bool ScoringScheme::is_good(
	Real const & num
) {
	static Real const TOL(1e-5);
	using std::abs;
	return (
		std::abs( num - 9999.000 ) > TOL && std::abs( num - 0 ) > TOL
	);
}

} // sequence
} // core
