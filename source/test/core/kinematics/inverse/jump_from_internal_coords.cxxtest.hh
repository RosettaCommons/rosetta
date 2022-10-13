// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/inverse/jump_from_internal_coords.cxxtest.hh
/// @brief  Utility functions for calculating jumps by knowing desired jump_from_internal_coords measurements
/// @author Jack Maguire

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/kinematics/inverse/jump.hh>
#include <core/kinematics/inverse/jump_from_internal_coords.hh>
#include <core/kinematics/inverse/util.hh>
#include <core/kinematics/inverse/AlignmentAtom.hh>

#include <core/kinematics/Jump.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <test/util/pose_funcs.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>

// Utility headers

// C/C++
#include <stdexcept>

//Auto Headers

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "test.core.kinematics.inverse.jump_from_internal_coords" );

using namespace core::kinematics;
using namespace core::kinematics::inverse;
using numeric::constants::d::pi;

using XYZ = numeric::xyzVector< core::Real >;

//helper function to reduce code duplication
void
test_calc_new_atom_location(
	XYZ const & a, XYZ const & b, XYZ const & c,
	core::Real const dist, core::Real const ba_rad, core::Real const tor_rad
){
	XYZ const d = calc_new_atom_location( a, b, c, dist, ba_rad, tor_rad );

	TR << "d x " << d.x() << std::endl;
	TR << "d y " << d.y() << std::endl;
	TR << "d z " << d.z() << std::endl;

	core::Real const final_dist = c.distance( d );
	TR << "dist: " << final_dist << std::endl;
	TS_ASSERT_DELTA( final_dist, dist, 0.01 );

	core::Real const final_ba_rad = numeric::angle_of( c-d, b-c );
	TR << "bond angle: " << final_ba_rad << std::endl;
	//gotta subtract from pi for some reason
	TS_ASSERT_DELTA( final_ba_rad, ba_rad, 0.01 );

	core::Real const final_tor_rad = numeric::dihedral_radians( a, b, c, d );
	TR << "torsion angle: " << final_tor_rad << std::endl;
	TS_ASSERT_DELTA( final_tor_rad, tor_rad, 0.01 );
}

class InverseJumpFromInternalCoordsTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	void test_angles_with_zero_torsion(){
		TR << "STARTING test_angles_with_zero_torsion" << std::endl;

		XYZ const a( 0, 0, 0 );
		XYZ const b( 0, 0, 1 );
		XYZ const c( 0, 1, 1 );

		core::Real const dist = 1;
		core::Real const ba_rad = 3.14/2;
		core::Real const tor_rad = 0.0;

		test_calc_new_atom_location( a, b, c, dist, ba_rad, tor_rad );
	}

	void test_angles_with_zero_torsion2(){
		TR << "STARTING test_angles_with_zero_torsion2" << std::endl;

		XYZ const a( 0, 0, 0 );
		XYZ const b( 0, 0, 1 );
		XYZ const c( 0, 1, 1 );

		core::Real const dist = 5;
		core::Real const ba_rad = 0.2; //found (2.9415 != 0.2000)
		core::Real const tor_rad = 0.0;

		test_calc_new_atom_location( a, b, c, dist, ba_rad, tor_rad );
	}

	void test_angles_with_zero_torsion3(){
		TR << "STARTING test_angles_with_zero_torsion3" << std::endl;

		XYZ const a( 0, 0, 0 );
		XYZ const b( 0, 1, 1 );
		XYZ const c( 0, 2, 1 );

		core::Real const dist = 5;
		core::Real const ba_rad = 3.14/2;
		core::Real const tor_rad = 0.0;

		test_calc_new_atom_location( a, b, c, dist, ba_rad, tor_rad );
	}

	void test_angles_w_torsion(){
		TR << "STARTING test_angles_w_torsion" << std::endl;

		XYZ const a( 0, 0, 0 );
		XYZ const b( 0, 1, 1 );
		XYZ const c( 0, 2, 1 );

		core::Real const dist = 5;
		core::Real const ba_rad = 3.10; //found (0.0415 != 3.1000)
		core::Real const tor_rad = 3.14;

		test_calc_new_atom_location( a, b, c, dist, ba_rad, tor_rad );
	}

	void test_angles_w_torsion2(){
		TR << "STARTING test_angles_w_torsion2" << std::endl;

		XYZ const a( 0, 0, 0 );
		XYZ const b( 0, 1, 1 );
		XYZ const c( 0, 2, 1 );

		core::Real const dist = 3;
		core::Real const ba_rad = 2.10; //found (1.0415 != 2.1000)
		core::Real const tor_rad = -1.14;

		test_calc_new_atom_location( a, b, c, dist, ba_rad, tor_rad );
	}

	void test_bind_to_CO(){
		TR << "STARTING test_bind_to_CO" << std::endl;

		//MAKE POSE
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/multistage_rosetta_scripts/3U3B_A.pdb" );
		core::pose::Pose peptide;
		core::pose::make_pose_from_sequence( peptide, "GGG", "fa_standard" );
		core::Size const original_pose_size = pose.size();
		core::Size const original_num_jump = pose.num_jump();
		pose.append_pose_by_jump( peptide, 1 );

		//DEFINE DESIRES
		core::Real const CA_C_O_H_torsion_angle_rad = 0.0;
		core::Real const C_O_H_N_torsion_angle_rad  = 3.0;
		core::Real const O_H_N_CA_torsion_angle_rad = 3.0;
		core::Real const C_O_H_bond_angle_rad = 3.10;
		core::Real const O_H_N_bond_angle_rad = 2.80;
		core::Real const O_H_dist_A = 2.0;

		//DEFINE ACCEPTOR ATOMS
		core::Size const acc_resid = 30;
		core::conformation::Residue const & acc_res = pose.residue( acc_resid );
		InternalCoordAtoms acc_atoms;
		acc_atoms.child.set( acc_res.atom_index( "O" ), acc_resid );
		acc_atoms.parent.set( acc_res.atom_index( "C" ), acc_resid );
		acc_atoms.grandparent.set( acc_res.atom_index( "CA" ), acc_resid );

		//DEFINE DONOR ATOMS
		core::Size const don_resid = original_pose_size+2;//middle peptide res
		core::conformation::Residue const & don_res = pose.residue( don_resid );
		InternalCoordAtoms don_atoms;
		don_atoms.child.set( don_res.atom_index( "H" ), don_resid );
		don_atoms.parent.set( don_res.atom_index( "N" ), don_resid );
		don_atoms.grandparent.set( don_res.atom_index( "CA" ), don_resid );

		//COLLECT GEOMETRY
		InternalCoordGeometry geom;
		geom.A_B_C_D_torsion_angle_rad = CA_C_O_H_torsion_angle_rad;
		geom.B_C_D_E_torsion_angle_rad = C_O_H_N_torsion_angle_rad;
		geom.C_D_E_F_torsion_angle_rad = O_H_N_CA_torsion_angle_rad;
		geom.B_C_D_bond_angle_rad = C_O_H_bond_angle_rad;
		geom.C_D_E_bond_angle_rad = O_H_N_bond_angle_rad;
		geom.C_D_dist_Ang = O_H_dist_A;

		//Calc new jump
		core::Size const jump_id = original_num_jump+1;
		Jump const j = jump_from_internal_coords( pose.conformation(), acc_atoms, don_atoms, geom, jump_id );
		pose.set_jump( jump_id, j );

		XYZ const new_acc_Oxyz = pose.xyz( acc_atoms.child );
		XYZ const new_acc_Cxyz = pose.xyz( acc_atoms.parent );
		XYZ const new_acc_CAxyz = pose.xyz( acc_atoms.grandparent );

		XYZ const new_don_Hxyz = pose.xyz( don_atoms.child );
		XYZ const new_don_Nxyz = pose.xyz( don_atoms.parent );
		XYZ const new_don_CAxyz = pose.xyz( don_atoms.grandparent );

		TS_ASSERT_DELTA( new_don_Hxyz.distance(new_acc_Oxyz), O_H_dist_A, 0.01 );

		TS_ASSERT_DELTA( O_H_N_bond_angle_rad, numeric::angle_of(
			new_don_Hxyz - new_acc_Oxyz,
			new_don_Hxyz - new_don_Nxyz
			), 0.01 );

		TS_ASSERT_DELTA( C_O_H_bond_angle_rad, numeric::angle_of(
			new_acc_Oxyz - new_acc_Cxyz,
			new_acc_Oxyz - new_don_Hxyz
			), 0.01 );

		TS_ASSERT_DELTA( CA_C_O_H_torsion_angle_rad,
			numeric::dihedral_radians(
			new_acc_CAxyz,
			new_acc_Cxyz,
			new_acc_Oxyz,
			new_don_Hxyz
			), 0.01 );

		TS_ASSERT_DELTA( C_O_H_N_torsion_angle_rad,
			numeric::dihedral_radians(
			new_acc_Cxyz,
			new_acc_Oxyz,
			new_don_Hxyz,
			new_don_Nxyz
			), 0.01 );

		TS_ASSERT_DELTA( O_H_N_CA_torsion_angle_rad,
			numeric::dihedral_radians(
			new_acc_Oxyz,
			new_don_Hxyz,
			new_don_Nxyz,
			new_don_CAxyz
			), 0.01 );

		TS_ASSERT_DELTA( new_don_Hxyz.distance(new_acc_Oxyz), O_H_dist_A, 0.01 );

	}

	void test_alignment_to_known_position(){
		TR << "STARTING test_alignment_to_known_position" << std::endl;

		//MAKE POSE
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/kinematics/inverse/3u3b.pdb.gz" );

		//Pardon the indenting - Rosetta's beutifier doesn't handle lamdas well
		auto const && get_resid_for_str = [&]( std::string const & s ) -> core::Size {
			TR << "!!!" << s << std::endl;
			core::select::residue_selector::ResidueIndexSelector selector( s );
			auto const sele = selector.apply( pose );
			for( core::Size resid = 1; resid <= sele.size(); ++resid ){
			if( sele[resid] ) return resid;
			}
			TS_ASSERT( false );
			return 0;
			};

		InternalCoordAtoms chain1_atoms;
		chain1_atoms.grandparent.set( 1, get_resid_for_str( "37A" ) ); //A
		chain1_atoms.parent.set( 1, get_resid_for_str( "24A" ) );      //B
		chain1_atoms.child.set( 1, get_resid_for_str( "45A" ) );       //C

		InternalCoordAtoms chain2_atoms;
		chain2_atoms.child.set( 1, get_resid_for_str( "41B" ) );       //D
		chain2_atoms.parent.set( 1, get_resid_for_str( "17B" ) );      //E
		chain2_atoms.grandparent.set( 1, get_resid_for_str( "52B" ) ); //F

		InternalCoordGeometry geom;
		geom.A_B_C_D_torsion_angle_rad = numeric::conversions::radians(134.8);//134.8 deg
		geom.B_C_D_E_torsion_angle_rad = numeric::conversions::radians(-98.2);//-98.2 deg
		geom.C_D_E_F_torsion_angle_rad = numeric::conversions::radians(-34.9);//-34.9 deg
		geom.B_C_D_bond_angle_rad = numeric::conversions::radians(158.8); //158.8 deg
		geom.C_D_E_bond_angle_rad = numeric::conversions::radians(93.3); //93.3 deg
		geom.C_D_dist_Ang = 11.7; //angstroms
		//

		auto const original_pose = pose.clone();
		auto && assert_trivial_deviation = [&](){
			// We should have reconstructed the original setup, so check for deviation
			for ( core::Size resid = 1; resid <= pose.size(); ++resid ) {
				core::Real const dist = pose.residue(resid).xyz(1).distance( original_pose->residue(resid).xyz(1) );
				TS_ASSERT_DELTA( dist, 0, 0.1 );
			}
		};

		TR << "NO CHANGE YET" << std::endl;
		assert_trivial_deviation();

		TS_ASSERT_EQUALS( pose.num_jump(), 1 );
		Jump const j = jump_from_internal_coords( pose.conformation(), chain1_atoms, chain2_atoms, geom, 1 );
		pose.set_jump( 1, j );

		TR << "NEW JUMP" << std::endl;
		assert_trivial_deviation();


		//Test init_from_current
		InternalCoordGeometry geom2;
		geom2.init_from_current( pose.conformation(), chain1_atoms, chain2_atoms );
		TS_ASSERT_DELTA( geom.A_B_C_D_torsion_angle_rad, geom2.A_B_C_D_torsion_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.B_C_D_bond_angle_rad, geom2.B_C_D_bond_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.C_D_dist_Ang, geom2.C_D_dist_Ang, 0.01 );
		TS_ASSERT_DELTA( geom.B_C_D_E_torsion_angle_rad, geom2.B_C_D_E_torsion_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.C_D_E_bond_angle_rad, geom2.C_D_E_bond_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.C_D_E_F_torsion_angle_rad, geom2.C_D_E_F_torsion_angle_rad, 0.01 );
	}

	void
	helper_assert_small_deviation(
		core::pose::Pose const & poseA,
		core::pose::Pose const & poseB,
		core::Real const expected_max_dist
	){
		core::Real max_dist = 0;
		for ( core::Size resid = 1; resid <= poseA.size(); ++resid ) {
			core::Real const dist = poseA.residue(resid).xyz(1).distance( poseB.residue(resid).xyz(1) );
			if ( dist > max_dist ) max_dist = dist;
		}
		//TR << "Max Dist: " << max_dist << std::endl;
		TS_ASSERT_DELTA( max_dist, expected_max_dist, 0.001 );
	};

	void test_smooth_angle_sampling(){
		//Sometimes negative angles can result in drastically different outcomes
		//I pushed a torsion-angle consideration for this and this unit test ensures it works
		//The strategy is simply to iteratively sample the bond angles and ensure that the pose doesn't make any big changes
		TR << "STARTING test_smooth_angle_sampling" << std::endl;

		//MAKE POSE
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/kinematics/inverse/3u3b.pdb.gz" );

		//Pardon the indenting - Rosetta's beautifier doesn't handle lambdas well
		auto const && get_resid_for_str = [&]( std::string const & s ) -> core::Size {
			core::select::residue_selector::ResidueIndexSelector selector( s );
			auto const sele = selector.apply( pose );
			for( core::Size resid = 1; resid <= sele.size(); ++resid ){
			if( sele[resid] ) return resid;
			}
			TS_ASSERT( false );
			return 0;
			};

		InternalCoordAtoms chain1_atoms;
		chain1_atoms.grandparent.set( 1, get_resid_for_str( "37A" ) ); //A
		chain1_atoms.parent.set( 1, get_resid_for_str( "24A" ) );      //B
		chain1_atoms.child.set( 1, get_resid_for_str( "45A" ) );       //C

		InternalCoordAtoms chain2_atoms;
		chain2_atoms.child.set( 1, get_resid_for_str( "41B" ) );       //D
		chain2_atoms.parent.set( 1, get_resid_for_str( "17B" ) );      //E
		chain2_atoms.grandparent.set( 1, get_resid_for_str( "52B" ) ); //F

		InternalCoordGeometry geom;
		geom.init_from_current( pose.conformation(), chain1_atoms, chain2_atoms );

		TR << "Test B_C_D_bond_angle_rad" << std::endl;
		{
			core::pose::PoseOP prev_pose = nullptr; //pose.clone();
			for ( core::Real angle_rad = -3.95; angle_rad <= 4.0; angle_rad += 0.1 ) {
				InternalCoordGeometry temp_geom = geom;
				temp_geom.B_C_D_bond_angle_rad = angle_rad;

				constexpr core::Size jump_id = 1;
				Jump const j = jump_from_internal_coords( pose.conformation(), chain1_atoms, chain2_atoms, temp_geom, jump_id );
				pose.set_jump( jump_id, j );

				if ( prev_pose != nullptr ) {
					helper_assert_small_deviation( pose, *prev_pose, 2.91921 );
				}
				prev_pose = pose.clone();
			}
		}

		TR << "Test C_D_E_bond_angle_rad" << std::endl;
		{
			core::pose::PoseOP prev_pose = nullptr;
			for ( core::Real angle_rad = -3.95; angle_rad <= 4.0; angle_rad += 0.1 ) {
				InternalCoordGeometry temp_geom = geom;
				temp_geom.C_D_E_bond_angle_rad = angle_rad;

				constexpr core::Size jump_id = 1;
				Jump const j = jump_from_internal_coords( pose.conformation(), chain1_atoms, chain2_atoms, temp_geom, jump_id );
				pose.set_jump( jump_id, j );

				if ( prev_pose != nullptr ) {
					helper_assert_small_deviation( pose, *prev_pose, 3.31759 );
				}
				prev_pose = pose.clone();
			}
		}


		TR << "Test A_B_C_D_torsion_angle_rad" << std::endl;
		{
			core::pose::PoseOP prev_pose = nullptr;
			for ( core::Real angle_rad = -3.95; angle_rad <= 4.0; angle_rad += 0.1 ) {
				InternalCoordGeometry temp_geom = geom;
				temp_geom.A_B_C_D_torsion_angle_rad = angle_rad;

				constexpr core::Size jump_id = 1;
				Jump const j = jump_from_internal_coords( pose.conformation(), chain1_atoms, chain2_atoms, temp_geom, jump_id );
				pose.set_jump( jump_id, j );

				if ( prev_pose != nullptr ) {
					helper_assert_small_deviation( pose, *prev_pose, 3.18903 );
				}
				prev_pose = pose.clone();
			}
		}

		TR << "Test B_C_D_E_torsion_angle_rad" << std::endl;
		{
			core::pose::PoseOP prev_pose = nullptr;
			for ( core::Real angle_rad = -3.95; angle_rad <= 4.0; angle_rad += 0.1 ) {
				InternalCoordGeometry temp_geom = geom;
				temp_geom.B_C_D_E_torsion_angle_rad = angle_rad;

				constexpr core::Size jump_id = 1;
				Jump const j = jump_from_internal_coords( pose.conformation(), chain1_atoms, chain2_atoms, temp_geom, jump_id );
				pose.set_jump( jump_id, j );

				if ( prev_pose != nullptr ) {
					helper_assert_small_deviation( pose, *prev_pose, 3.22914 );
				}
				prev_pose = pose.clone();
			}
		}

		TR << "Test C_D_E_F_torsion_angle_rad" << std::endl;
		{
			core::pose::PoseOP prev_pose = nullptr;
			for ( core::Real angle_rad = -3.95; angle_rad <= 4.0; angle_rad += 0.1 ) {
				InternalCoordGeometry temp_geom = geom;
				temp_geom.C_D_E_F_torsion_angle_rad = angle_rad;

				constexpr core::Size jump_id = 1;
				Jump const j = jump_from_internal_coords( pose.conformation(), chain1_atoms, chain2_atoms, temp_geom, jump_id );
				pose.set_jump( jump_id, j );

				if ( prev_pose != nullptr ) {
					helper_assert_small_deviation( pose, *prev_pose, 1.91758 );
				}
				prev_pose = pose.clone();
			}
		}

	}


	void test_agreement_between_forwards_and_backwards_calculations(){
		TR << "STARTING test_agreement_between_forwards_and_backwards_calculations" << std::endl;

		//MAKE POSE
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/kinematics/inverse/3u3b.pdb.gz" );

		//Pardon the indenting - Rosetta's beautifier doesn't handle lambdas well
		auto const && get_resid_for_str = [&]( std::string const & s ) -> core::Size {
			TR << "!!!" << s << std::endl;
			core::select::residue_selector::ResidueIndexSelector selector( s );
			auto const sele = selector.apply( pose );
			for( core::Size resid = 1; resid <= sele.size(); ++resid ){
			if( sele[resid] ) return resid;
			}
			TS_ASSERT( false );
			return 0;
			};

		InternalCoordAtoms chain1_atoms;
		chain1_atoms.grandparent.set( 1, get_resid_for_str( "32A" ) ); //A
		chain1_atoms.parent.set( 1, get_resid_for_str( "33A" ) );      //B
		chain1_atoms.child.set( 1, get_resid_for_str( "34A" ) );       //C

		InternalCoordAtoms chain2_atoms;
		chain2_atoms.child.set( 1, get_resid_for_str( "16B" ) );       //D
		chain2_atoms.parent.set( 1, get_resid_for_str( "17B" ) );      //E
		chain2_atoms.grandparent.set( 1, get_resid_for_str( "18B" ) ); //F

		InternalCoordGeometry geom;
		geom.init_from_current( pose.conformation(), chain1_atoms, chain2_atoms );

		Jump const j = jump_from_internal_coords( pose.conformation(), chain1_atoms, chain2_atoms, geom, 1 );
		pose.set_jump( 1, j );

		InternalCoordGeometry geom2;
		geom2.init_from_current( pose.conformation(), chain1_atoms, chain2_atoms );

		//Test that the cycle is mostly lossless
		TS_ASSERT_DELTA( geom.A_B_C_D_torsion_angle_rad, geom2.A_B_C_D_torsion_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.B_C_D_E_torsion_angle_rad, geom2.B_C_D_E_torsion_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.C_D_E_F_torsion_angle_rad, geom2.C_D_E_F_torsion_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.B_C_D_bond_angle_rad, geom2.B_C_D_bond_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.C_D_E_bond_angle_rad, geom2.C_D_E_bond_angle_rad, 0.01 );
		TS_ASSERT_DELTA( geom.C_D_dist_Ang, geom2.C_D_dist_Ang, 0.01 );
		//If any of these fail, then the method for calculating these values differs from the method for applying them
	}

};
