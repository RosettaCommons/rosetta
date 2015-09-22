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
/// @author  Andy Watkins

// Unit headers
#include <core/pose/ncbb/util.hh>
#include <core/pose/Pose.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
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
static THREAD_LOCAL basic::Tracer TR( "core.pose.ncbb.util" );


namespace core {
namespace pose {
namespace ncbb {

using namespace std;
using namespace core;

utility::vector1< core::Size >
initialize_ncbbs (
	Pose & pose
) {
	utility::vector1< core::Size > ncbb_seq_positions;
	core::Size hbs = 0;

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		// now we return both PRE and POST locations...
		if ( pose.residue(i).has_variant_type(chemical::OOP_PRE) == 1 ) {
			ncbb_seq_positions.push_back( i );
			core::pose::ncbb::add_oop_constraint(pose, i);
		}
		if ( pose.residue(i).has_variant_type(chemical::OOP_POST)== 1 ) {
			ncbb_seq_positions.push_back( i );
		}
		if ( pose.residue(i).has_variant_type(chemical::TRIAZOLAMERC) == 1 ) {
			ncbb_seq_positions.push_back( i );
			core::pose::ncbb::add_triazole_constraint(pose, i);
		}
		if ( pose.residue(i).has_variant_type(chemical::TRIAZOLAMERN)== 1 ) {
			ncbb_seq_positions.push_back( i );
		}
		if ( pose.residue(i).has_variant_type(chemical::HBS_PRE) == 1 ) {
			ncbb_seq_positions.push_back( i );
			hbs = 1;
			core::pose::ncbb::add_hbs_constraint(pose, i);
		}
		if ( pose.residue(i).has_variant_type(chemical::HBS_POST) == 1 ) {
			ncbb_seq_positions.push_back( i );
			hbs = 0;
		}

		if ( hbs == 1 ) { // we're inside the hbs macrocycle; even though these residues are unpatched
			ncbb_seq_positions.push_back( i );
		}
	}

	return ncbb_seq_positions;
}

utility::vector1< core::Size >
initialize_hbs (
	Pose & pose
) {
	utility::vector1< core::Size > hbs_seq_positions;
	core::Size hbs = 0;

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).has_variant_type(chemical::HBS_PRE) == 1 ) {
			hbs_seq_positions.push_back( i );
			hbs = 1;
		}
		if ( pose.residue(i).has_variant_type(chemical::HBS_POST)== 1 ) {
			hbs_seq_positions.push_back( i );
			hbs = 0;
		}

		if ( hbs == 1 ) { // we're inside the hbs macrocycle; even thoough these residues are unpatched
			hbs_seq_positions.push_back( i );
		}
	}
	return hbs_seq_positions;
}

void add_a3b_hbs_constraint( core::pose::Pose & pose, core::Size a3b_hbs_pre_position )
{
	add_a3b_hbs_constraint( pose, a3b_hbs_pre_position, 1.479871, 0.05 );
}

void add_a3b_hbs_constraint( core::pose::Pose & pose, core::Size hbs_pre_position, core::Real distance, core::Real std )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;

	HarmonicFuncOP harm_func  (new HarmonicFunc( distance, std ) );
	HarmonicFuncOP harm_func_0  (new HarmonicFunc( 0, std ) );
	CircularHarmonicFuncOP ang_func  (new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_2_over_3(), 0.02 ) );
	CircularHarmonicFuncOP ang_func2 (new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_over_3(), 0.02 ) );
	CircularHarmonicFuncOP dih_func  (new CircularHarmonicFunc( numeric::NumericTraits<float>::pi(), 0.02 ) );
	CircularHarmonicFuncOP dih_func2  (new CircularHarmonicFunc( (0-numeric::NumericTraits<float>::pi_over_3()), 0.02 ) );

	AtomID aidCYH( pose.residue( hbs_pre_position ).atom_index("CYH"), hbs_pre_position );
	AtomID aidCZH( pose.residue( hbs_pre_position+2 ).atom_index("CZH"), hbs_pre_position+2 );
	AtomID aidN  ( pose.residue( hbs_pre_position+2 ).atom_index("N"), hbs_pre_position+2 );
	AtomID aidCA  ( pose.residue( hbs_pre_position+2 ).atom_index("CA"), hbs_pre_position+2 );
	AtomID aidVZH( pose.residue( hbs_pre_position ).atom_index("VZH"), hbs_pre_position );
	AtomID aidVYH( pose.residue( hbs_pre_position+2 ).atom_index("VYH"), hbs_pre_position+2 );
	AtomID aidCY3( pose.residue( hbs_pre_position ).atom_index("CY3"), hbs_pre_position );
	AtomID aidCY2( pose.residue( hbs_pre_position ).atom_index("CY2"), hbs_pre_position );

	ConstraintCOP atompair ( new AtomPairConstraint( aidCYH, aidCZH, harm_func ) );
	ConstraintCOP atompair2 ( new AtomPairConstraint( aidCYH, aidVYH, harm_func_0 ) );
	ConstraintCOP atompair3 ( new AtomPairConstraint( aidCZH, aidVZH, harm_func_0 ) );
	ConstraintCOP angle ( new AngleConstraint( aidCZH, aidCYH, aidCY3, ang_func2 ) );
	ConstraintCOP angle2 ( new AngleConstraint( aidN, aidCZH, aidCYH, ang_func2 ) );
	ConstraintCOP dihedral ( new DihedralConstraint( aidCZH, aidCYH, aidCY3, aidCY2, dih_func ) );
	ConstraintCOP dihedral2 ( new DihedralConstraint( aidCA, aidN, aidCZH, aidCYH, dih_func2 ) );


	pose.add_constraint( atompair );
	pose.add_constraint( atompair2 );
	pose.add_constraint( atompair3 );
	pose.add_constraint( angle );
	pose.add_constraint( angle2 );
	pose.add_constraint( dihedral );
	pose.add_constraint( dihedral2 );

	TR << "added atom pair constraint to a3b hbs with distance: " << distance << " and std: "<< std << std::endl;
	TR << "and atom pair constraints with the virtual atoms" << std::endl;

}

void
constrain_ring_atoms( core::pose::Pose & pose, utility::vector1< core::id::AtomID > ids ) {

	using namespace core;
	using namespace scoring;
	using namespace func;
	using namespace constraints;

	using namespace core::id;

	CircularHarmonicFuncOP chf( new CircularHarmonicFunc( 0, 0.01 ) );
	CircularHarmonicFuncOP ahf( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi() * ( ids.size() - 2 ) / ids.size(), 0.01 ) );
	for ( Size i = 1; i <= ids.size(); ++i ) {
		pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
			ids[ ( i ) % ids.size() + 1 ],
			ids[ (i+1) % ids.size() + 1 ],
			ids[ (i+2) % ids.size() + 1 ],
			ids[ (i+3) % ids.size() + 1 ], chf ) ) );

		pose.add_constraint( AngleConstraintOP( new AngleConstraint(
			ids[ ( i ) % ids.size() + 1 ],
			ids[ (i+1) % ids.size() + 1 ],
			ids[ (i+2) % ids.size() + 1 ], ahf ) ) );
	}
}

void add_triazole_constraint( core::pose::Pose & pose, core::Size i /* triazole_position */)
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;
	runtime_assert_msg( pose.residue( i ).has_variant_type(core::chemical::TRIAZOLAMERC), "residue must have TRIAZOLAMERC variant type" );
	runtime_assert_msg( pose.residue( i+1 ).has_variant_type(core::chemical::TRIAZOLAMERN), "next residue must have TRIAZOLAMERN variant type" );

	HarmonicFuncOP hf( new HarmonicFunc( 1.347, 0.05 ) );
	HarmonicFuncOP zf( new HarmonicFunc( 0.000, 0.01 ) );
	CircularHarmonicFuncOP chf( new CircularHarmonicFunc( 0, 0.01 ) );
	CircularHarmonicFuncOP chf180( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi(), 0.01 ) );
	CircularHarmonicFuncOP ahf1( new CircularHarmonicFunc( 130.1*numeric::NumericTraits<float>::pi()/180.0, 0.01 ) );
	CircularHarmonicFuncOP ahf2( new CircularHarmonicFunc( 122.6*numeric::NumericTraits<float>::pi()/180.0, 0.01 ) );

	std::string cai   = pose.residue(  i  ).is_protein() ? "CA" : "CH3";
	std::string caip1 = pose.residue( i+1 ).is_protein() ? "CA" : "CH3";

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
		AtomID( pose.residue( i   ).atom_index( "CT1" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT3" ), i+1 ), hf ) ) );

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
		AtomID( pose.residue( i   ).atom_index( "CT2" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ), hf ) ) );

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
		AtomID( pose.residue( i   ).atom_index(  "CT1" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "VCT1" ), i+1 ), zf ) ) );

	pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint(
		AtomID( pose.residue( i   ).atom_index( "VNT3" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index(  "NT3" ), i+1 ), zf ) ) );

	// exo bonds to triazolamers
	pose.add_constraint( AngleConstraintOP( new AngleConstraint(
		AtomID( pose.residue( i   ).atom_index(  cai  ), i   ),
		AtomID( pose.residue( i   ).atom_index( "CT1" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT3" ), i+1 ), ahf2 ) ) );

	pose.add_constraint( AngleConstraintOP( new AngleConstraint(
		AtomID( pose.residue( i   ).atom_index( "CT2" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ),
		AtomID( pose.residue( i+1 ).atom_index( caip1 ), i+1 ), ahf1 ) ) );

	pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
		AtomID( pose.residue( i+1 ).atom_index( caip1 ), i+1 ),
		AtomID( pose.residue( i   ).atom_index( "CT2" ), i  ),
		AtomID( pose.residue( i+1 ).atom_index( "NT2" ), i+1 ),
		AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ), chf ) ) );

	pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
		AtomID( pose.residue( i   ).atom_index(  cai  ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT3" ), i+1 ),
		AtomID( pose.residue( i   ).atom_index( "CT2" ), i   ),
		AtomID( pose.residue( i   ).atom_index( "CT1" ), i   ), chf ) ) );

	pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
		AtomID( pose.residue( i+1 ).atom_index( caip1 ), i+1 ),
		AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ),
		AtomID( pose.residue( i   ).atom_index( "CT2" ), i   ),
		AtomID( pose.residue( i   ).atom_index( "CT1" ), i   ), chf180 ) ) );

	pose.add_constraint( DihedralConstraintOP( new DihedralConstraint(
		AtomID( pose.residue( i   ).atom_index(  cai  ), i   ),
		AtomID( pose.residue( i   ).atom_index( "CT1" ), i   ),
		AtomID( pose.residue( i+1 ).atom_index( "NT3" ), i+1 ),
		AtomID( pose.residue( i+1 ).atom_index( "NT2" ), i+1 ), chf180 ) ) );

	utility::vector1< AtomID > ids;
	ids.push_back( AtomID( pose.residue( i+1 ).atom_index( "NT2" ), i+1 ) );
	ids.push_back( AtomID( pose.residue( i+1 ).atom_index( "NT3" ), i+1 ) );
	ids.push_back( AtomID( pose.residue( i   ).atom_index( "CT1" ), i   ) );
	ids.push_back( AtomID( pose.residue( i   ).atom_index( "CT2" ), i   ) );
	ids.push_back( AtomID( pose.residue( i+1 ).atom_index( "NT1" ), i+1 ) );
	constrain_ring_atoms( pose, ids );

	TR << "added atom pair constraint to triazole with distance: 1.347 and std 0.05" << std::endl;
	TR << "and atom pair constraints to the virtual atoms" << std::endl;

}

void add_hbs_constraint( core::pose::Pose & pose, core::Size hbs_pre_position )
{
	add_hbs_constraint( pose, hbs_pre_position, 1.479871, 0.05 );
}

void add_hbs_constraint( core::pose::Pose & pose, core::Size hbs_pre_position, core::Real distance, core::Real std )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;

	//kdrew: add constraint
	HarmonicFuncOP harm_func( new HarmonicFunc( distance, std ) );
	HarmonicFuncOP harm_func_0( new HarmonicFunc( 0, std ) );
	CircularHarmonicFuncOP ang_func( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_2_over_3(), 0.02 ) );
	CircularHarmonicFuncOP ang_func2( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_over_3(), 0.02 ) );
	CircularHarmonicFuncOP ang_func3( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi()/180*109.5, 0.02 ) );
	CircularHarmonicFuncOP dih_func( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi(), 0.02 ) );
	CircularHarmonicFuncOP dih_func2( new CircularHarmonicFunc( 0, 0.02 ) );

	AtomID aidCYH( pose.residue( hbs_pre_position ).atom_index("CYH"), hbs_pre_position );
	AtomID aidCZH( pose.residue( hbs_pre_position+2 ).atom_index("CZH"), hbs_pre_position+2 );
	AtomID aidN  ( pose.residue( hbs_pre_position+2 ).atom_index("N"), hbs_pre_position+2 );
	AtomID aidVZH( pose.residue( hbs_pre_position ).atom_index("VZH"), hbs_pre_position );
	AtomID aidVYH( pose.residue( hbs_pre_position+2 ).atom_index("VYH"), hbs_pre_position+2 );
	AtomID aidCY2( pose.residue( hbs_pre_position ).atom_index("CY2"), hbs_pre_position );
	AtomID aidCY1( pose.residue( hbs_pre_position ).atom_index("CY1"), hbs_pre_position );

	ConstraintCOP atompair( ConstraintOP( new AtomPairConstraint( aidCYH, aidCZH, harm_func ) ) );
	ConstraintCOP atompair2( ConstraintOP( new AtomPairConstraint( aidCYH, aidVYH, harm_func_0 ) ) );
	ConstraintCOP atompair3( ConstraintOP( new AtomPairConstraint( aidCZH, aidVZH, harm_func_0 ) ) );
	ConstraintCOP angle( ConstraintOP( new AngleConstraint( aidCZH, aidCYH, aidCY2, ang_func2 ) ) );
	ConstraintCOP angle2( ConstraintOP( new AngleConstraint( aidN, aidCZH, aidCYH, ang_func3 ) ) );
	ConstraintCOP dihedral( ConstraintOP( new DihedralConstraint( aidCZH, aidCYH, aidCY2, aidCY1, dih_func ) ) );
	//ConstraintCOP dihedral2( ConstraintOP( new DihedralConstraint( aidCZH, aidCYH, aidCY2, aidCY1, dih_func2 ) ) );

	pose.add_constraint( atompair );
	pose.add_constraint( atompair2 );
	pose.add_constraint( atompair3 );
	pose.add_constraint( angle );
	pose.add_constraint( angle2 );
	pose.add_constraint( dihedral );
	//pose.add_constraint( dihedral2 );

	TR << "added atom pair constraint to hbs with distance: " << distance << " and std: "<< std << std::endl;
	TR << "and atom pair constraints with the virtual atoms" << std::endl;

}

/// @details loops through residues searching for oop variants and sets up constraints
utility::vector1< core::Size >
initialize_oops(
	Pose & pose
) {
	utility::vector1< core::Size > oop_seq_positions;
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).has_variant_type(chemical::OOP_PRE) == 1 ) {
			oop_seq_positions.push_back( i );
			core::pose::ncbb::add_oop_constraint(pose, i);
		}
	}
	return oop_seq_positions;
}

// Add constraints to keep oligooxopiperazine (oop) ring closed, default values (distance = 1.5, std = 0.05)
// ??? oop should be 1.518685 according to params and QM optimized oop_dimer_llll!
/// @details Overloaded function which defines default values, calls more general add_oop_constraint function
void add_oop_constraint( core::pose::Pose & pose, core::Size oop_seq_position )
{
	add_oop_constraint( pose, oop_seq_position, 1.518685, 0.05 );
}

// Add constraints to keep oligooxopiperazine (oop) ring closed
/// @details General function to add atom pair constraints to CYP and CZP for given distance \n
/// Other constraints are setup to keep virtual atoms VYP fixed near CYP and VZP fixed near CZP
void add_oop_constraint( core::pose::Pose & pose, core::Size oop_seq_position, core::Real distance, core::Real std )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	runtime_assert_msg( pose.residue(oop_seq_position).has_variant_type(chemical::OOP_PRE) == 1, "residue must have OOP_PRE variant type" );

	//kdrew: add constraint
	core::scoring::func::HarmonicFuncOP harm_func( new core::scoring::func::HarmonicFunc( distance, std ) );

	//kdrew: constrain: CYP VYP  and CZP VZP, hard coded to have zero distance
	core::scoring::func::HarmonicFuncOP virtual_atom_overlap_harm_func( new core::scoring::func::HarmonicFunc( 0.0, 0.1 ) );

	AtomID aidCYP( pose.residue( oop_seq_position ).atom_index("CYP"), oop_seq_position );
	AtomID aidCZP( pose.residue( oop_seq_position+1 ).atom_index("CZP"), oop_seq_position+1 );
	AtomID aidVYP( pose.residue( oop_seq_position+1 ).atom_index("VYP"), oop_seq_position+1 );
	AtomID aidVZP( pose.residue( oop_seq_position ).atom_index("VZP"), oop_seq_position );

	AtomPairConstraintCOP CYP_CZP_atompair( AtomPairConstraintOP( new AtomPairConstraint( aidCYP, aidCZP, harm_func ) ) );
	//kdrew: setup virtual atoms constraints
	AtomPairConstraintCOP CYP_VYP_atompair( AtomPairConstraintOP( new AtomPairConstraint( aidCYP, aidVYP, virtual_atom_overlap_harm_func ) ) );
	AtomPairConstraintCOP CZP_VZP_atompair( AtomPairConstraintOP( new AtomPairConstraint( aidCZP, aidVZP, virtual_atom_overlap_harm_func ) ) );

	//kdrew: remove old constraints that are identical
	core::scoring::constraints::ConstraintCOPs cs = pose.constraint_set()->get_all_constraints();
	for ( Size i = 1; i <= cs.size(); i++ ) {
		Constraint const & other_cst = *cs[i];
		if ( !dynamic_cast< AtomPairConstraint const * > ( &other_cst ) ) {
			continue;
		}

		AtomPairConstraint const & constraint_i( static_cast< AtomPairConstraint const & > (other_cst) );

		if ( (constraint_i.atom(1) == CYP_CZP_atompair->atom(1) && constraint_i.atom(2) == CYP_CZP_atompair->atom(2))
				|| (constraint_i.atom(1) == CYP_VYP_atompair->atom(1) && constraint_i.atom(2) == CYP_VYP_atompair->atom(2))
				|| (constraint_i.atom(1) == CZP_VZP_atompair->atom(1) && constraint_i.atom(2) == CZP_VZP_atompair->atom(2))  ) {
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
