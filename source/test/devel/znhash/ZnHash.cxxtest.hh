// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/devel/zn_hash/ZnHash.cxxtest.hh
/// @brief  test suite for devel::zn_hash::ZnHash
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <devel/znhash/ZnHash.hh>

// Utility headers
//#include <utility/fixedsizearray1.hh>
//#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/residue_io.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/Constraint.hh>

//#include <protocols/match/SixDHasher.hh>
//#include <protocols/match/MatchSet.hh>

// C++ headers
//#include <string>
//#include <iostream>

//#include <numeric/HomogeneousTransform.hh>

// Boost headers
//#include <boost/unordered_map.hpp>
//#include <boost/algorithm/string/predicate.hpp>

//Auto Headers
//#include <protocols/match/Hit.hh>
//#include <utility/vector1.hh>
//#include <boost/unordered/unordered_map_fwd.hpp>

#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/enzdes/EnzdesMovers.hh>


// --------------- Test Class --------------- //

class ZnHashTests : public CxxTest::TestSuite {

public:

	typedef core::Vector Vector;
	typedef core::Size   Size;
	typedef core::Real   Real;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
	typedef numeric::HomogeneousTransform< Real >    HTReal;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-run::preserve_header -extra_res_fa devel/znhash/ZNX.params" );
		//core_init_with_additional_options( "-run::preserve_header -symmetry:symmetry_definition devel/znhash/symm_def_file_cn4_generic -symmetry:initialize_rigid_body_dofs" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_ZnCoordinateConstraint_transform_coordinates() { // empty test for the moment


		core::pose::Pose match_pose, match1;
		core::import_pose::pose_from_file( match_pose, "devel/znhash/1EER_A.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( match1, "devel/znhash/UM_1_H5E11_95_1EER_A_ZNX_1.pdb" , core::import_pose::PDB_file);

		core::id::AtomID zn_atom_id( match1.residue(3).atom_index( "ZN" ), 3);


		protocols::enzdes::AddOrRemoveMatchCsts adder_1H_2DE;
		adder_1H_2DE.cstfile( "devel/znhash/ZNX.cst" );
		adder_1H_2DE.set_cst_action( protocols::enzdes::ADD_NEW );

		adder_1H_2DE.apply( match1 );

		protocols::enzdes::EnzdesConstraintReporter reporter;
		reporter.ligand_resno( 3 );
		reporter.find_constraints_to_ligand( match1 );

		utility::vector1< std::pair< core::id::AtomID, core::Real > > closest_atoms(2, std::make_pair( core::id::AtomID(), -1 ));


		for ( core::Size ii = 1; ii <= reporter.constrained_nonligand_atoms().size(); ++ii ) {
			core::id::AtomID iiatid = reporter.constrained_nonligand_atoms()[ ii ];
			TS_ASSERT( iiatid.rsd() == 1 || iiatid.rsd() == 2 );
			//std::cout << "Constraint to ligand #" << ii << " Res " <<
			// iiatid.rsd() << " atom " <<
			// match1.residue( iiatid.rsd() ).atom_name( iiatid.atomno() ) <<
			// std::endl;
			core::Real const d2 = match1.xyz(zn_atom_id).distance_squared( match1.xyz( iiatid ) );
			if ( closest_atoms[ iiatid.rsd() ].second < 0 || d2 < closest_atoms[ iiatid.rsd() ].second ) {
				closest_atoms[ iiatid.rsd() ].first = iiatid;
				closest_atoms[ iiatid.rsd() ].second = d2;
			}
		}
		for ( core::Size ii = 1; ii <= 2; ++ii ) {
			//std::cout << "Closest atom rsd: " << closest_atoms[ii].first.rsd() << " " <<
			// match1.residue( closest_atoms[ii].first.rsd() ).atom_name( closest_atoms[ii].first.atomno() ) <<
			// " with distance " << std::sqrt( closest_atoms[ii].second ) << std::endl;
		}


		core::chemical::ResidueType const & znx_restype =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->
			name_map( "ZNX" );

		core::Vector znpos = znx_restype.atom( znx_restype.atom_index( "ZN" )).ideal_xyz();
		core::Vector v1pos = znx_restype.atom( znx_restype.atom_index( "V1" )).ideal_xyz();
		core::Vector v2pos = znx_restype.atom( znx_restype.atom_index( "V2" )).ideal_xyz();
		core::Vector v12_midpoint = v1pos + v2pos / 2;

		numeric::HomogeneousTransform< core::Real > znframe( v1pos, v12_midpoint, znpos );

		/// Now measure the coordinates of the four virtual atoms in this idealized frame
		numeric::HomogeneousTransform< core::Real > invznframe = znframe.inverse();
		utility::vector1< core::Vector > znx_ideal_coords( znx_restype.natoms() );
		for ( core::Size ii = 1; ii <= znx_restype.natoms(); ++ii ) {
			znx_ideal_coords[ ii ] = invznframe * znx_restype.atom( ii ).ideal_xyz();
			if ( ii != 1 ) {
				znx_ideal_coords[ ii ].normalize();
			}
		}

		core::Vector r1coord_atom = match1.xyz( closest_atoms[1].first );
		core::Vector r2coord_atom = match1.xyz( closest_atoms[2].first );
		//core::Vector r12_halfway = r1coord_atom + r2coord_atom / 2;

		core::Vector zn_atom = match1.xyz( core::id::AtomID( 1, 3 ) );
		core::Vector zn_vect = ((r1coord_atom-zn_atom).normalize() + (r2coord_atom-zn_atom).normalize()) / 2;

		core::Vector walk_from_zn = zn_atom + zn_vect;

		numeric::HomogeneousTransform< core::Real > actual_frame( r1coord_atom, walk_from_zn, zn_atom );

		utility::vector1< core::Vector > znx_match_coords( znx_restype.natoms() );
		for ( Size ii = 1; ii <= znx_restype.natoms(); ++ii ) {
			core::Vector iipos = actual_frame * znx_ideal_coords[ ii ];
			znx_match_coords[ ii ] = iipos;
			//std::cout << " ZNX coordinates " << znx_restype.atom_name( ii ) << " "
			// << iipos.x() << " "
			// << iipos.y() << " "
			// << iipos.z() << std::endl;
		}

		/// Now lets consider the transformation to place these coordinates into the frame of the original
		/// 1EER_A.pdb given that they are being built on chain B of a

		core::pose::Pose c4_centroid_pose;
		core::import_pose::pose_from_file( c4_centroid_pose,
			* core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ),
			"devel/znhash/1EER_cn4_0043.pdb", core::import_pose::PDB_file );
		//core::pose::symmetry::make_symmetric_pose( c4_centroid_pose );

		core::id::AtomID aunit_n(1,  1), aunit_ca(2,  1), aunit_c(3,  1);
		core::id::AtomID bunit_n(1,167), bunit_ca(2,167), bunit_c(3,167);

		numeric::HomogeneousTransform< core::Real > aframe(
			c4_centroid_pose.xyz( aunit_n ),
			c4_centroid_pose.xyz( aunit_ca),
			c4_centroid_pose.xyz( aunit_c ));

		numeric::HomogeneousTransform< core::Real > bframe(
			c4_centroid_pose.xyz( bunit_n ),
			c4_centroid_pose.xyz( bunit_ca),
			c4_centroid_pose.xyz( bunit_c ));

		numeric::HomogeneousTransform< core::Real > matchframe(
			match_pose.xyz( aunit_n ),
			match_pose.xyz( aunit_ca),
			match_pose.xyz( aunit_c ));

		numeric::HomogeneousTransform< core::Real > inv_matchframe = matchframe.inverse();
		utility::vector1< core::Vector > zn_match_coords_in_reference_frame( znx_match_coords.size() );
		for ( Size ii = 1; ii <= znx_restype.natoms(); ++ii ) {
			zn_match_coords_in_reference_frame[ ii ] = inv_matchframe * znx_match_coords[ ii ];
		}

		core::conformation::ResidueOP znA = core::conformation::ResidueFactory::create_residue( znx_restype );
		for ( Size ii = 1; ii <= znA->natoms(); ++ii ) {
			znA->set_xyz( ii, aframe * zn_match_coords_in_reference_frame[ii] );
		}

		core::conformation::ResidueOP znB = core::conformation::ResidueFactory::create_residue( znx_restype );
		for ( Size ii = 1; ii <= znB->natoms(); ++ii ) {
			znB->set_xyz( ii, bframe * zn_match_coords_in_reference_frame[ii] );
		}

		//c4_centroid_pose.append_residue_by_jump( *znA, 1 );
		//c4_centroid_pose.append_residue_by_jump( *znB, 1 );
		//c4_centroid_pose.dump_pdb( "1EER_c4_w_2ZNX.pdb" );

		devel::znhash::ZnCoordinationScorer znscore;
		znscore.set_symm_resid( 167 );
		znscore.set_reference_pdb( "devel/znhash/1EER_A.pdb" );
		znscore.set_matcher_constraint_file_name( "devel/znhash/ZNX.cst" );
		znscore.add_match_from_file( "devel/znhash/UM_1_H5E11_95_1EER_A_ZNX_1.pdb" );

		devel::znhash::ZnCoord first_zn = znscore.original_frame_coordinate_for_match( 1, c4_centroid_pose );
		numeric::HomogeneousTransform< core::Real > inv_aframe = aframe.inverse();
		for ( Size ii = 1; ii <= 5; ++ii ) {
			TS_ASSERT( first_zn[ii].distance( inv_aframe*znB->xyz(ii) ) < 1e-6 );
			//std::cout << ii << " " << first_zn[ii].x() << " " << first_zn[ii].y() << " " << first_zn[ii].z() << " vs" << std::endl;
			//std::cout << ii << " " << (inv_aframe*znB->xyz(ii)).x() << " " << (inv_aframe*znB->xyz(ii)).y() << " " << (inv_aframe*znB->xyz(ii)).z() << std::endl;
		}

		znscore.finalize_after_all_matches_added();
		TS_ASSERT( znscore.score( c4_centroid_pose ) == 0.0 );

	}

	void test_ZnCoordConstraint_score_two_matches() {

		core::pose::Pose c4_centroid_pose;
		core::import_pose::pose_from_file( c4_centroid_pose,
			* core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ),
			"devel/znhash/1EER_cn4_0043.pdb", core::import_pose::PDB_file );

		devel::znhash::ZnCoordinationScorer znscore;
		znscore.set_symm_resid( 167 );
		znscore.set_reference_pdb( "devel/znhash/1EER_A.pdb" );

		std::ifstream match_cst_pairs( "devel/znhash/match_cst_pairs.txt" );
		if ( ! match_cst_pairs.good() ) {
			std::cerr << "Error opening devel/znhash/match_cst_pairs.txt" << std::endl;
			TS_ASSERT( match_cst_pairs.good() );
			return;
		}
		//while ( match_cst_pairs ) {
		// std::string fname1, fname2;
		// match_cst_pairs >> fname1 >> fname2;
		// if ( fname1 == "" || fname2 == "" ) break;
		// //std::cout << "reading " << fname1 << " and " << fname2 << std::endl;
		// znscore.add_match_from_file( fname1, fname2 );
		//}
		znscore.set_matcher_constraint_file_name( "devel/znhash/ZNX.cst" );

		znscore.add_match_from_file( "devel/znhash/UM_39_H104D100_28_1EER_A_ZNX_1.pdb" );
		znscore.add_match_from_file( "devel/znhash/UM_50_H14E10_32_1EER_A_ZNX_1.pdb" );

		znscore.add_match_from_file( "devel/znhash/UM_26_H100D103_1_1EER_A_ZNX_1.pdb" );
		znscore.add_match_from_file( "devel/znhash/UM_41_H13E17_3_1EER_A_ZNX_1.pdb" );

		znscore.finalize_after_all_matches_added();
		//clock_t starttime = clock();
		//for ( Size ii = 1; ii <= 10000; ++ii ) {
		// znscore.score( c4_centroid_pose );
		//}
		//clock_t stoptime = clock();
		//std::cout << "timing: " << ((double) stoptime-starttime ) / ( CLOCKS_PER_SEC * 10000 ) << std::endl;

		//Size temp_precision( std::cout.precision());
		//std::cout.precision(12);
		//std::cout << "SCORE: " << znscore.score( c4_centroid_pose ) << std::endl;
		//std::cout.precision( temp_precision );

		znscore.set_clash_weight( 0.0 );
		devel::znhash::ZnCoordinationScorer::CoordinationData result =
			znscore.score_and_index_for_best_match( c4_centroid_pose.residue(1), c4_centroid_pose.residue(167) );
		TS_ASSERT( result.first.first  == 2 );
		TS_ASSERT( result.first.second == 1 );
		TS_ASSERT_DELTA( result.second.first, -2.87761791297, 1e-9 );

		znscore.set_clash_weight( 0.5 );
		result = znscore.score_and_index_for_best_match( c4_centroid_pose.residue(1), c4_centroid_pose.residue(167) );
		TS_ASSERT( result.first.first  == 2 );
		TS_ASSERT( result.first.second == 1 );
		//int start_precision( std::cout.precision() ); std::cout.precision( 12 ); std::cout << "w/ clash score: " << result.second.first << std::endl; std::cout.precision(start_precision);
		TS_ASSERT_DELTA( result.second.first, -0.894529044997, 1e-9 );


	}
	void test_score_zn_constraint_in_pose() {

		core::pose::Pose c4_centroid_pose;
		core::import_pose::pose_from_file( c4_centroid_pose,
			* core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ),
			"devel/znhash/1EER_cn4_0043.pdb", core::import_pose::PDB_file );

		devel::znhash::ZnCoordinationScorerOP znscore( new devel::znhash::ZnCoordinationScorer );
		znscore->set_symm_resid( 167 );
		znscore->set_reference_pdb( "devel/znhash/1EER_A.pdb" );
		znscore->set_matcher_constraint_file_name( "devel/znhash/ZNX.cst" );
		znscore->add_match_from_file( "devel/znhash/UM_39_H104D100_28_1EER_A_ZNX_1.pdb" );
		znscore->add_match_from_file( "devel/znhash/UM_50_H14E10_32_1EER_A_ZNX_1.pdb" );
		znscore->add_match_from_file( "devel/znhash/UM_26_H100D103_1_1EER_A_ZNX_1.pdb" );
		znscore->add_match_from_file( "devel/znhash/UM_41_H13E17_3_1EER_A_ZNX_1.pdb" );

		znscore->finalize_after_all_matches_added();
		znscore->set_clash_weight( 0.0 );
		devel::znhash::ZnCoordinationConstraintOP zncst( new devel::znhash::ZnCoordinationConstraint( znscore ) );
		c4_centroid_pose.add_constraint( zncst );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::metalhash_constraint, 1.0 );

		core::Real score = sfxn( c4_centroid_pose );
		TS_ASSERT_DELTA( score, -2.87761791297, 1e-9 );

	}
};
