// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/conformation/membrane/Span.fwd.hh
///
/// @brief  Object for describing start and end of a transmembrane span
/// @details The Span object stores 2 SSizes - a stard and end position of a transmembrane span.
///    Should be kept in a vector of spans toward describing the total spanning topology of a
///    membrane protein.
///    Last Modified: 7/23/14
///
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <core/conformation/membrane/Span.hh>

// Package Headers
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

static basic::Tracer TR( "core.conformation.membrane.Span" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace membrane {

////////////////////
/// Constructors ///
////////////////////

/// @brief Default Constructor
/// @details Construct a default span object representing a span from 1-1
/// this constructor should eventually be made private because it doesn't build a real thing
Span::Span() :
	utility::pointer::ReferenceCount(),
	start_(1),
	end_(1),
	orientation_( in )
{}

/// @brief Custom Constructor - Construct new span
/// @details Constructor from start and end
Span::Span( core::Size start, core::Size end, Orientation orient ) :
	utility::pointer::ReferenceCount(),
	start_( start ),
	end_( end ),
	orientation_( orient )
{}

/// @brief Copy Consturctor
/// @details Make a deep copy of this object
Span::Span( Span const & src ) :
	utility::pointer::ReferenceCount(),
	start_( src.start_ ),
	end_( src.end_ ),
	orientation_( src.orientation_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this object
Span &
Span::operator=( Span const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Make a deep copy of everything
	this->start_ = src.start_;
	this->end_ = src.end_;
	this->orientation_ = src.orientation_;

	return *this;
}

/// @brief Destructor
Span::~Span(){}

//////////////////
/// Accessors  ///
//////////////////

/// @brief Get start position
/// @details Get the Starting Position of a transmembrane span
core::Size
Span::start() const {
	return start_;
}

/// @brief Get end position
/// @details Get the end position of a transmembrane span
core::Size
Span::end() const {
	return end_;
}

/// @brief Get the orientation of the span
Orientation
Span::orientation() const {
	return orientation_;
}

/// @brief Set the orientation of the span
void Span::orientation( Orientation orient_in ) {
	orientation_ = orient_in;
}

/// @brief get residue closest to center
core::Size Span::center() const {
	return ( ( start_ + end_ ) / 2 );
}

/// @brief Shift by offset
/// @details Shift the transmembrane span by a user-provided offset
void Span::shift( core::Size offset ) {
	start_ += offset;
	end_ += offset;
	not_valid();
}

/// @brief Show This Span
/// @details Show the information in this span. TODO: Should override base method
void
Span::show() const {
	TR << "Span: start: " << start_ << ", end: " << end_ << std::endl;
}

/// @brief Check that this Span is Valid
/// @details Check that this span describes a consecutive transmembrane span
/// of nonzero length.
bool Span::is_valid() const {

	bool valid( false );

	// span start > span end
	if ( start_ > end_ ) {
		TR << "Span is invalid: start > end!" << std::endl;
		return false;
	}

	// start or end are zero
	if ( start_ == 0 || end_ == 0 ) {
		TR << "Span is invalid: either start or end = 0!" << std::endl;
		return false;
	}

	// define length
	core::Size length = end_ - start_ + 1;

	// short span
	if ( length <= 5 ) {
		TR.Warning << "SHORT SPAN: SPAN IS ONLY " << length << " RESIDUES LONG!!!" << std::endl;
	}

	if ( length >= 30 ) {
		TR.Warning << "LONG SPAN: SPAN IS " << length << " RESIDUES LONG!!!" << std::endl;
	}

	valid = true;
	return valid;

} // is valid

void Span::not_valid() const {

	if ( ! is_valid() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Span is invalid!" );
	}

} // not valid

} // membrane
} // conformation
} // core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::membrane::Span::save( Archive & arc ) const {
	arc( CEREAL_NVP( start_ ) ); // Size
	arc( CEREAL_NVP( end_ ) ); // Size
	arc( CEREAL_NVP( orientation_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::membrane::Span::load( Archive & arc ) {
	arc( start_ ); // Size
	arc( end_ ); // Size
	arc( orientation_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::membrane::Span );
CEREAL_REGISTER_TYPE( core::conformation::membrane::Span )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_membrane_Span )
#endif // SERIALIZATION
