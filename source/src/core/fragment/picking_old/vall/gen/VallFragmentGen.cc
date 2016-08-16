// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/gen/VallFragmentGen.cc
/// @brief  base class Vall ExtentGenerator
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/gen/VallFragmentGen.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace gen {


/// @brief default constructor
VallFragmentGen::VallFragmentGen() :
	Super()
{}


/// @brief copy constructor
VallFragmentGen::VallFragmentGen( VallFragmentGen const & rval ) :
	Super( rval )
{}


/// @brief default destructor
VallFragmentGen::~VallFragmentGen()
{}


/// @brief copy assignment
VallFragmentGen & VallFragmentGen::operator =( VallFragmentGen const & rval ) {
	if ( this != &rval ) {}
	return *this;
}


} // namespace gen
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core
