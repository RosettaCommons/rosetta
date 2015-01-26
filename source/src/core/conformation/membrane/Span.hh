// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/conformation/membrane/Span.fwd.hh
///
/// @brief		Object for describing start and end of a transmembrane span
/// @details	The Span object stores 2 SSizes - a stard and end position of a transmembrane span.
///				Should be kept in a vector of spans toward describing the total spanning topology of a
///				membrane protein.
///				Last Modified: 7/23/14
///
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_Span_hh
#define INCLUDED_core_conformation_membrane_Span_hh

// Unit headers
#include <core/conformation/membrane/Span.fwd.hh>

// Package Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace conformation {
namespace membrane {

using namespace core;

class Span : public utility::pointer::ReferenceCount {

public: // constructors

	/// @brief	Default Constructor
	/// @details Construct a default span object representing a span from 1-1
	/// this constructor should eventually be made private because it doesn't build a real thing
	Span();

	/// @brief Custom Constructor - Construct new span
	/// @details Constructor from start and end
	Span( Size start, Size end );

	/// @brief Copy Consturctor
	/// @details Make a deep copy of this object
    Span( Span const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this object
    Span &
    operator=( Span const & src );

	/// @brief	Destructor
	~Span();

public: // getters

    /// @brief Get start position
	/// @details Get the Starting Position of a transmembrane span
    Size start() const;

    /// @brief Get end position
	/// @details Get the end position of a transmembrane span
    Size end() const;

	
	/// @brief Get residue closest to center
	Size center() const;
	
	/// @brief Shift by offset
	/// @details Shift the transmembrane span by a user-provided offset
	void shift( Size offset );

	/// @brief Show This Span
	/// @details Show the information in this span. TODO: Should override base method
	void show();

	/// @brief Check that this Span is Valid
	/// @details Check that this span describes a consecutive transmembrane span
	/// of nonzero length.
	bool is_valid();

	// TODO: get rid of this guy
	void not_valid();

private: // data

	// Specify start/end position of a transmembrane span
    Size start_;
    Size end_;

};

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_Span_hh
