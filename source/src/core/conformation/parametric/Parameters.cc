// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/conformation/parametric/Parameters.cc
/// @brief  A class for holding sets of parameters for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <core/conformation/parametric/Parameters.hh>

// Package headers

// Project headers

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/assert.hh>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>


namespace core {
namespace conformation {
namespace parametric {

static thread_local basic::Tracer TR( "core.conformation.parametric.Parameters" );

/// @brief Constructor.
///
Parameters::Parameters() :
	residue_list_()
{
}

Parameters::Parameters( Parameters const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< Parameters >()
{
	residue_list_.clear();
	if ( src.residue_list_.size()>0 ) {
		for ( core::Size i=1, imax=src.residue_list_.size(); i<=imax; ++i ) {
			residue_list_.push_back( src.residue_list_[i]->clone() ); //This copies the residue that was being pointed at.
			//Note that when copying a Conformation, I need to add logic that will ensure that the Parameters objects that result have owning pointers to the residues in the Conformation,
			//rather than to residues that only exist in the Parameters object.
		}
	}
}

Parameters::~Parameters() {}


/// @brief make a copy of this residue( allocate actual memory for it )
///
ParametersOP
Parameters::clone() const
{
	return ParametersOP( new Parameters( *this ) );
}

} // namespace parametric
} // namespace conformation
} // namespace core

