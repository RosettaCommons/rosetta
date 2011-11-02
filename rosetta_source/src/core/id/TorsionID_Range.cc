// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/TorsionID_Range.cc
/// @brief  Conformation torsion identifier class
/// @author Colin Smith


// Unit headers
#include <core/id/TorsionID_Range.hh>

// C++ headers
// AUTO-REMOVED #include <ostream>


namespace core {
namespace id {


/// @brief stream << TorsionID_Range
std::ostream &
operator <<( std::ostream & os, TorsionID_Range const & a )
{
	os << a.torsion_id() << " min= " << a.min() << " max= " << a.max() << ' ';
	return os;
}


} // namespace id
} // namespace core
