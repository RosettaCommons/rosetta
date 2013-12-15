// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/pose/ncbb/util.cc
/// @brief   Utility function definitions for poses with noncanonical backbones
/// @author  kdrew

// Unit headers
#include <core/pose/ncbb/util.hh>
#include <core/pose/Pose.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>

// Utility headers
//#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>


// Construct tracer.
static basic::Tracer TR("core.pose.ncbb.util");


namespace core {
namespace pose {
namespace ncbb {

using namespace std;
using namespace core;

// initializes oops in poses
/// @details loops through residues searching for oop variants and sets up constraints
utility::vector1< core::Size >
initialize_oops(
	Pose & pose
) {
	utility::vector1< core::Size > oop_seq_positions;
	for ( Size i=1; i<= pose.total_residue(); ++i )
	{
		if( pose.residue(i).has_variant_type(chemical::OOP_PRE) == 1 )
		{
			oop_seq_positions.push_back( i );
			core::pose::ncbb::add_oop_constraint(pose, i);
		}
	}
	return oop_seq_positions;
}

// Add constraints to keep oligooxopiperazine (oop) ring closed, default values (distance = 1.5, std = 0.05)
/// @details Overloaded function which defines default values, calls more general add_oop_constraint function
void add_oop_constraint( core::pose::Pose & pose, core::Size oop_seq_position )
{
	add_oop_constraint( pose, oop_seq_position, 1.5, 0.05 );
}
// Add constraints to keep oligooxopiperazine (oop) ring closed
/// @details General function to add atom pair constraints to CYP and CZP for given distance \n
/// Other constraints are setup to keep virtual atoms VYP fixed near CYP and VZP fixed near CZP
void add_oop_constraint( core::pose::Pose & pose, core::Size oop_seq_position, core::Real distance, core::Real std )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	runtime_assert_msg ( pose.residue(oop_seq_position).has_variant_type(chemical::OOP_PRE) == 1, "residue must have OOP_PRE variant type" );

	//kdrew: add constraint
	core::scoring::func::HarmonicFuncOP harm_func  (new core::scoring::func::HarmonicFunc( distance, std ) );

	//kdrew: constrain: CYP VYP  and CZP VZP, hard coded to have zero distance
	core::scoring::func::HarmonicFuncOP virtual_atom_overlap_harm_func  (new core::scoring::func::HarmonicFunc( 0.0, 0.1 ) );

	AtomID aidCYP( pose.residue( oop_seq_position ).atom_index("CYP"), oop_seq_position );
	AtomID aidCZP( pose.residue( oop_seq_position+1 ).atom_index("CZP"), oop_seq_position+1 );
	AtomID aidVYP( pose.residue( oop_seq_position+1 ).atom_index("VYP"), oop_seq_position+1 );
	AtomID aidVZP( pose.residue( oop_seq_position ).atom_index("VZP"), oop_seq_position );

	AtomPairConstraintCOP CYP_CZP_atompair = new AtomPairConstraint( aidCYP, aidCZP, harm_func );
	//kdrew: setup virtual atoms constraints
	AtomPairConstraintCOP CYP_VYP_atompair = new AtomPairConstraint( aidCYP, aidVYP, virtual_atom_overlap_harm_func );
	AtomPairConstraintCOP CZP_VZP_atompair = new AtomPairConstraint( aidCZP, aidVZP, virtual_atom_overlap_harm_func );

	//kdrew: remove old constraints that are identical
	core::scoring::constraints::ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();
	for(Size i = 1; i <= cs.size(); i++)
	{
		Constraint const & other_cst = *cs[i];
		if( !dynamic_cast< AtomPairConstraint const * > ( &other_cst ) )
		{ continue; }

		AtomPairConstraint const & constraint_i( static_cast< AtomPairConstraint const & > (other_cst) );

		if ( (constraint_i.atom(1) == CYP_CZP_atompair->atom(1) && constraint_i.atom(2) == CYP_CZP_atompair->atom(2))
			|| (constraint_i.atom(1) == CYP_VYP_atompair->atom(1) && constraint_i.atom(2) == CYP_VYP_atompair->atom(2))
			|| (constraint_i.atom(1) == CZP_VZP_atompair->atom(1) && constraint_i.atom(2) == CZP_VZP_atompair->atom(2))  )
		{
			pose.remove_constraint( cs[i], true );
			TR << "found and removed atom pair constraint from oop at residue: " << oop_seq_position << std::endl;
		}
	}

	pose.add_constraint( CYP_CZP_atompair );
	pose.add_constraint( CYP_VYP_atompair );
	pose.add_constraint( CZP_VZP_atompair );

	TR << "added atom pair constraints to oop at residue: " << oop_seq_position << " with distance: " << distance << " and std: "<< std << std::endl;

}


}  // namespace ncbb
}  // namespace pose
}  // namespace core
