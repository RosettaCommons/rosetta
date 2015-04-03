// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/gen/SecStructGen.cc
/// @brief  default constant length fragment VallExtentGenerator
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/gen/SecStructGen.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace gen {


/// @brief default constructor
SecStructGen::SecStructGen() :
	Super()
{}


/// @brief fragment size constructor
/// @param[in] ss the required secondary structure string of the fragment
SecStructGen::SecStructGen( String const & ss ) :
	Super(),
	ss_( ss )
{}


/// @brief copy constructor
SecStructGen::SecStructGen( SecStructGen const & rval ) :
	Super( rval ),
	ss_( rval.ss_ )
{}


/// @brief default destructor
SecStructGen::~SecStructGen()
{}


/// @brief copy assignment
SecStructGen & SecStructGen::operator =( SecStructGen const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
		ss_ = rval.ss_;
	}
	return *this;
}


/// @brief clone this object
VallFragmentGenOP SecStructGen::clone() const {
	return VallFragmentGenOP( new SecStructGen( *this ) );
}


/// @brief return the desired fragment extent w/ length equal to the
///  secondary structure string
/// @return Valid (true) extent if the extent has exactly the required
///  secondary structure string and the end of the extent does not go past
///  section_end.  Invalid (false) extent otherwise.
/// @remarks we assume VallResidueIterator is a type of RandomAccessIterator, such as
///  those used in std::vector
SecStructGen::Extent SecStructGen::operator ()( VallResidueIterator extent_begin, VallResidueIterator section_end ) const {
	Extent extent;
	extent.begin = extent_begin;
	extent.end = extent_begin + ss_.length();
	extent.valid = ( extent.end <= section_end );

	if ( extent.valid ) {

		// check if secondary structure string matches
		Size str_idx = 0;
		for ( VallResidueIterator i = extent.begin; i != extent.end; ++i, ++str_idx ) {
		debug_assert( str_idx != ss_.length() );

			// if secondary structure string doesn't match, set invalid
			// Extent and break out of loop
			if ( i->ss() != ss_.at( str_idx ) ) {
				extent.valid = false;
				break;
			}
		}

	} // if extent.valid

	return extent;
}


} // namespace gen
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core
