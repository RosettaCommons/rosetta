// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMLJLibary.cc
/// @brief  Molecular mechanics lj library
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/mm/MMLJLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/MMAtomType.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers

// C++ headers
#include <string>
#include <map>

//Auto Headers
namespace core {
namespace scoring {
namespace mm {

/// @details Auto-generated virtual destructor
MMLJLibrary::~MMLJLibrary() {}

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "core.mm.MMLJLibrary" );

/// @details Constructs a MMLJLibrary instance from a filename string and constant access pointer to an MMAtomTypeSet
MMLJLibrary::MMLJLibrary(
	core::chemical::MMAtomTypeSetCOP mm_atom_set
):
	// hard coding these here for now
	//  nblist_dis2_cutoff_XX_( 60.0 ),
	//  nblist_dis2_cutoff_XH_( 60.0 ),
	//  nblist_dis2_cutoff_HH_( 60.0 )
	// TESTING THESE
	nblist_dis2_cutoff_XX_( 72.25 ),
	nblist_dis2_cutoff_XH_(  3.19 ),
	nblist_dis2_cutoff_HH_( 19.36 )
{
	// set the MM atom type set
	mm_atom_set_ = core::chemical::MMAtomTypeSetCAP( mm_atom_set );
	core::chemical::MMAtomTypeSetCOP mm_atom_set_op( mm_atom_set_ );

	// add the lj params
	for ( Size i = 1; i <= mm_atom_set_op->n_atomtypes(); ++i ) {

		// get lj radius and well depth
		Real radius = (*mm_atom_set_op)[ i ].lj_radius();
		Real wdepth = (*mm_atom_set_op)[ i ].lj_wdepth();
		Real radius_3b = (*mm_atom_set_op)[ i ].lj_three_bond_radius();
		Real wdepth_3b = (*mm_atom_set_op)[ i ].lj_three_bond_wdepth();

		// add to correct library
		mm_lj_library_.push_back( mm_lj_param_set( radius, wdepth ) );
		mm_lj_three_bond_library_.push_back( mm_lj_param_set( radius_3b, wdepth_3b ) );
	}

	// print number lj params added
	TR << "MM lj sets added: " << mm_lj_library_.size() << std::endl;
}

} // namespace mm
} // namespace scoring
} // namespace core
