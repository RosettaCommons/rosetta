// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/eval/VallFragmentEval.cc
/// @brief  base class for Vall ExtentEvaluator
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.hh>

// package headers
#include <core/fragment/picking_old/vall/scores/VallFragmentScore.hh>

#include <core/fragment/picking_old/concepts/Extent.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


/// @brief default constructor
VallFragmentEval::VallFragmentEval() :
	Super()
{}


/// @brief default copy constructor
VallFragmentEval::VallFragmentEval( VallFragmentEval const & rval ) :
	Super( rval )
{}


/// @brief default destructor
VallFragmentEval::~VallFragmentEval()
{}


/// @brief copy assignment
VallFragmentEval & VallFragmentEval::operator =( VallFragmentEval const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
	}
	return *this;
}


/// @brief called by VallLibrarian: for a fragment extent, evaluate and store
///  results in a VallFragmentScore
/// @return true if score should be stored, false otherwise
bool VallFragmentEval::operator ()(
	Extent const & extent,
	VallFragmentScore & fs
)
{
	// store part of the extent data
	fs.extent_begin = extent.begin;
	fs.extent_end = extent.end;

	// evaluate the extent
	return eval_impl( extent, fs );
}


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core
