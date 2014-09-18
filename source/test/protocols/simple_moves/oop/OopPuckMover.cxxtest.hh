// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/oop/OopPuckMover.cxxtest.hh
/// @brief  test suite for protocols::simple_moves::oop::OopPuckMover.cc
/// @author Kevin Drew


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Unit headers
#include <protocols/simple_moves/oop/OopPuckMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>

// Package headers
#include <core/io/pdb/pose_io.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/database/open.hh>

#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <utility/tools/make_vector1.hh>

static basic::Tracer TR("protocols.simple_moves.oop.OopPuckMover.cxxtest");

// using declarations
using namespace core;
using namespace core::pose;
using namespace protocols;
using namespace protocols::moves;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// @name OopPuckMoverTest
/// @brief: test conformation switches from puck plus and puck minus
/// @detailed
///
/// @author Kevin Drew October 17, 2011
///////////////////////////////////////////////////////////////////////////
class OopPuckMoverTest : public CxxTest::TestSuite {

public:
	//chemical::ResidueTypeSetCAP residue_set;

	pose::Pose pose;
	OopPuckMoverTest() {}

	void setUp() {
		core_init();
		basic::options::option[ basic::options::OptionKeys::chemical::include_patches](utility::tools::make_vector1( std::string("patches/oop_pre.txt"), std::string("patches/oop_post.txt") ) );

		// create score function
		scoring::ScoreFunctionOP score_fxn( get_score_function() );
		//scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

		//kdrew: create new residue type set with oop patches included, cannot use chemical manager residue type set singleton because already initialized without oop patches
		std::string const directory( basic::database::full_name( "chemical/residue_type_sets/fa_standard/" ) );
		std::string const tag( "fa_standard" );
		chemical::ResidueTypeSetCAP residue_set = new chemical::ResidueTypeSet( tag, directory );

		// read pdb file
		core::import_pose::pose_from_pdb( pose, *residue_set, "protocols/simple_moves/oop/oop_test.pdb" );
		//scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

		core::init::init_random_generators(1000, "mt19937");
	}

	void tearDown() {

	}

	void test_OopPuckMoverPlus() {
		core::Real correct_phi ( -131.27 );
		core::Real correct_psi ( -10.38 );
		core::Real correct_CZP_CYP_N_CA ( -57.4 );
		core::Real correct_CYP_CZP_N_CA ( 160.7 );
		core::Real correct_HYP1_CYP_N_CA ( 63.5 );
		core::Real correct_HYP2_CYP_N_HYP1( 116.2 );
		core::Real correct_HZP1_CZP_NZ_CAZ ( 38.5 );
		core::Real correct_HZP2_CZP_NZ_HZP1( -116.0 );

		//kdrew: positions of oop in mdm2/p53 oop mimic (based on pdbid: 1ycr.pdb)
		//utility::vector1< core::Size > oop_seq_pos(1);
		core::Size oop_seq_pos = 88;

		//simple_moves::oop::OopPuckMoverOP opm( new simple_moves::oop::OopPuckMover( oop_seq_pos , true, false, false, false ) );
		simple_moves::oop::OopPuckPlusMoverOP opm( new simple_moves::oop::OopPuckPlusMover( oop_seq_pos ) );
		//moves::OopPuckMoverOP opm( new moves::OopPuckMover( oop_seq_pos ) );
		opm->apply(pose);

		core::Real phi = pose.phi( oop_seq_pos );
		core::Real psi = pose.psi( oop_seq_pos );
		TR << "phi and psi: " << phi << ' ' << psi << std::endl;

		//kdrew: get xyz coords of other angles to test
		Vector const& CA_xyz ( pose.residue(oop_seq_pos).xyz("CA") );
		Vector const& N_xyz ( pose.residue(oop_seq_pos).xyz("N") );
		Vector const& CYP_xyz ( pose.residue(oop_seq_pos).xyz("CYP") );
		Vector const& HYP1_xyz ( pose.residue(oop_seq_pos).xyz("1HYP") );
		Vector const& HYP2_xyz ( pose.residue(oop_seq_pos).xyz("2HYP") );

		Vector const& CAZ_xyz ( pose.residue(oop_seq_pos+1).xyz("CA") );
		Vector const& NZ_xyz ( pose.residue(oop_seq_pos+1).xyz("N") );
		Vector const& CZP_xyz ( pose.residue(oop_seq_pos+1).xyz("CZP") );
		Vector const& HZP1_xyz ( pose.residue(oop_seq_pos+1).xyz("1HZP") );
		Vector const& HZP2_xyz ( pose.residue(oop_seq_pos+1).xyz("2HZP") );

		core::Real CZP_CYP_N_CA = numeric::dihedral_degrees( CZP_xyz, CYP_xyz, N_xyz, CA_xyz ) ;
		core::Real CYP_CZP_N_CA = numeric::dihedral_degrees( CYP_xyz, CZP_xyz, NZ_xyz, CAZ_xyz) ;
		TR<< "CZP_CYP_N_CA and CYP_CZP_N_CA: " << CZP_CYP_N_CA << ' ' << CYP_CZP_N_CA << std::endl;

		core::Real HYP1_CYP_N_CA = numeric::dihedral_degrees( HYP1_xyz, CYP_xyz, N_xyz, CA_xyz ) ;
		core::Real HYP2_CYP_N_HYP1 = numeric::dihedral_degrees( HYP2_xyz, CYP_xyz, N_xyz, HYP1_xyz ) ;
		TR<< "HYP1_CYP_N_CA and HYP1_CYP_N_HYP1: " << HYP1_CYP_N_CA << ' ' << HYP2_CYP_N_HYP1 << std::endl;

		core::Real HZP1_CZP_NZ_CAZ = numeric::dihedral_degrees( HZP1_xyz, CZP_xyz, NZ_xyz, CAZ_xyz ) ;
		core::Real HZP2_CZP_NZ_HZP1 = numeric::dihedral_degrees( HZP2_xyz, CZP_xyz, NZ_xyz, HZP1_xyz ) ;
		TR<< "HZP1_CZP_NZ_CAZ and HZP1_CZP_NZ_HZP1: " << HZP1_CZP_NZ_CAZ << ' ' << HZP2_CZP_NZ_HZP1 << std::endl;

		// compare to the correct answer
		float TOLERATED_ERROR = 0.01;

		TS_ASSERT_DELTA( phi, correct_phi, TOLERATED_ERROR );
		TS_ASSERT_DELTA( psi, correct_psi, TOLERATED_ERROR );

		TOLERATED_ERROR = 1.00;

		TS_ASSERT_DELTA( CZP_CYP_N_CA, correct_CZP_CYP_N_CA, TOLERATED_ERROR );
		TS_ASSERT_DELTA( CYP_CZP_N_CA, correct_CYP_CZP_N_CA, TOLERATED_ERROR );

		//kdrew: test hydrogen placement
		TS_ASSERT_DELTA( HYP1_CYP_N_CA, correct_HYP1_CYP_N_CA, TOLERATED_ERROR );
		TS_ASSERT_DELTA( HYP2_CYP_N_HYP1, correct_HYP2_CYP_N_HYP1, TOLERATED_ERROR );
		TS_ASSERT_DELTA( HZP1_CZP_NZ_CAZ, correct_HZP1_CZP_NZ_CAZ, TOLERATED_ERROR );
		TS_ASSERT_DELTA( HZP2_CZP_NZ_HZP1, correct_HZP2_CZP_NZ_HZP1, TOLERATED_ERROR );

		TR << "test_OopPuckMoverPlus completed!! " << std::endl;
	}

	void test_OopPuckMoverMinus() {
		core::Real correct_phi ( -147.05 );
		core::Real correct_psi ( -36.90 );
		core::Real correct_CZP_CYP_N_CA ( 28.85 );
		core::Real correct_CYP_CZP_N_CA ( -136.31 );
		core::Real correct_HYP1_CYP_N_CA ( 150.5 );
		core::Real correct_HYP2_CYP_N_HYP1( 116.2 );
		core::Real correct_HZP1_CZP_NZ_CAZ ( 102.5 );
		core::Real correct_HZP2_CZP_NZ_HZP1( -116.0 );

		//kdrew: positions of oop in mdm2/p53 oop mimic (based on pdbid: 1ycr.pdb)
		core::Size oop_seq_pos = 88;

		//simple_moves::oop::OopPuckMoverOP opm( new simple_moves::oop::OopPuckMover( oop_seq_pos , false, true, false, false ) );
		simple_moves::oop::OopPuckMinusMoverOP opm( new simple_moves::oop::OopPuckMinusMover( oop_seq_pos ) );
		//simple_moves::oop::OopPuckMoverOP opm( new simple_moves::oop::OopPuckMover( oop_seq_pos , false, true, false, false ) );
		//moves::OopPuckMoverOP opm( new moves::OopPuckMover( oop_seq_pos ) );
		opm->apply(pose);

		core::Real phi = pose.phi( oop_seq_pos );
		core::Real psi = pose.psi( oop_seq_pos );
		TR << "phi and psi: " << phi << ' ' << psi << std::endl;

		//kdrew: get xyz coords of other angles to test
		Vector const& CA_xyz ( pose.residue(oop_seq_pos).xyz("CA") );
		Vector const& N_xyz ( pose.residue(oop_seq_pos).xyz("N") );
		Vector const& CYP_xyz ( pose.residue(oop_seq_pos).xyz("CYP") );
		Vector const& HYP1_xyz ( pose.residue(oop_seq_pos).xyz("1HYP") );
		Vector const& HYP2_xyz ( pose.residue(oop_seq_pos).xyz("2HYP") );

		Vector const& CAZ_xyz ( pose.residue(oop_seq_pos+1).xyz("CA") );
		Vector const& NZ_xyz ( pose.residue(oop_seq_pos+1).xyz("N") );
		Vector const& CZP_xyz ( pose.residue(oop_seq_pos+1).xyz("CZP") );
		Vector const& HZP1_xyz ( pose.residue(oop_seq_pos+1).xyz("1HZP") );
		Vector const& HZP2_xyz ( pose.residue(oop_seq_pos+1).xyz("2HZP") );

		core::Real CZP_CYP_N_CA = numeric::dihedral_degrees( CZP_xyz, CYP_xyz, N_xyz, CA_xyz ) ;
		core::Real CYP_CZP_N_CA = numeric::dihedral_degrees( CYP_xyz, CZP_xyz, NZ_xyz, CAZ_xyz) ;
		TR<< "CZP_CYP_N_CA and CYP_CZP_N_CA: " << CZP_CYP_N_CA << ' ' << CYP_CZP_N_CA << std::endl;

		core::Real HYP1_CYP_N_CA = numeric::dihedral_degrees( HYP1_xyz, CYP_xyz, N_xyz, CA_xyz ) ;
		core::Real HYP2_CYP_N_HYP1 = numeric::dihedral_degrees( HYP2_xyz, CYP_xyz, N_xyz, HYP1_xyz ) ;
		TR<< "HYP1_CYP_N_CA and HYP1_CYP_N_HYP1: " << HYP1_CYP_N_CA << ' ' << HYP2_CYP_N_HYP1 << std::endl;

		core::Real HZP1_CZP_NZ_CAZ = numeric::dihedral_degrees( HZP1_xyz, CZP_xyz, NZ_xyz, CAZ_xyz ) ;
		core::Real HZP2_CZP_NZ_HZP1 = numeric::dihedral_degrees( HZP2_xyz, CZP_xyz, NZ_xyz, HZP1_xyz ) ;
		TR<< "HZP1_CZP_NZ_CAZ and HZP1_CZP_NZ_HZP1: " << HZP1_CZP_NZ_CAZ << ' ' << HZP2_CZP_NZ_HZP1 << std::endl;

		// compare to the correct answer
		float TOLERATED_ERROR = 0.01;

		TS_ASSERT_DELTA( phi, correct_phi, TOLERATED_ERROR );
		TS_ASSERT_DELTA( psi, correct_psi, TOLERATED_ERROR );

		TOLERATED_ERROR = 1.00;

		TS_ASSERT_DELTA( CZP_CYP_N_CA, correct_CZP_CYP_N_CA, TOLERATED_ERROR );
		TS_ASSERT_DELTA( CYP_CZP_N_CA, correct_CYP_CZP_N_CA, TOLERATED_ERROR );

		//kdrew: test hydrogen placement
		TS_ASSERT_DELTA( HYP1_CYP_N_CA, correct_HYP1_CYP_N_CA, TOLERATED_ERROR );
		TS_ASSERT_DELTA( HYP2_CYP_N_HYP1, correct_HYP2_CYP_N_HYP1, TOLERATED_ERROR );
		TS_ASSERT_DELTA( HZP1_CZP_NZ_CAZ, correct_HZP1_CZP_NZ_CAZ, TOLERATED_ERROR );
		TS_ASSERT_DELTA( HZP2_CZP_NZ_HZP1, correct_HZP2_CZP_NZ_HZP1, TOLERATED_ERROR );


		TR << "test_OopPuckMoverMinus completed!! " << std::endl;
	}

	void test_OopPuckMoverSmall() {
		//kdrew: positions of oop in mdm2/p53 oop mimic (based on pdbid: 1ycr.pdb)
		utility::vector1< core::Size > oop_seq_pos(1);
		oop_seq_pos[1] = 88;

		core::Real old_phi = pose.phi( oop_seq_pos[1] );
		core::Real old_psi = pose.psi( oop_seq_pos[1] );
		TR << "old phi and old psi: " << old_phi << ' ' << old_psi << std::endl;

		core::Real small_angle_max = 2.0;
		//simple_moves::oop::OopPuckMoverOP opm_small( new simple_moves::oop::OopPuckMover( oop_seq_pos, false, false, true, false, true, small_angle_max ) );
		simple_moves::oop::OopRandomSmallMoverOP opm_small( new simple_moves::oop::OopRandomSmallMover( oop_seq_pos, small_angle_max ) );
		opm_small->apply(pose);

		core::Real phi = pose.phi( oop_seq_pos[1] );
		core::Real psi = pose.psi( oop_seq_pos[1] );
		TR << "phi and psi: " << phi << ' ' << psi << std::endl;

		//kdrew: get xyz coords of other angles to test
		Vector const& CA_xyz ( pose.residue(oop_seq_pos[1]).xyz("CA") );
		Vector const& N_xyz ( pose.residue(oop_seq_pos[1]).xyz("N") );
		Vector const& CYP_xyz ( pose.residue(oop_seq_pos[1]).xyz("CYP") );
		Vector const& CAZ_xyz ( pose.residue(oop_seq_pos[1]+1).xyz("CA") );
		Vector const& NZ_xyz ( pose.residue(oop_seq_pos[1]+1).xyz("N") );
		Vector const& CZP_xyz ( pose.residue(oop_seq_pos[1]+1).xyz("CZP") );

		core::Real CZP_CYP_N_CA = numeric::dihedral_degrees( CZP_xyz, CYP_xyz, N_xyz, CA_xyz ) ;
		core::Real CYP_CZP_N_CA = numeric::dihedral_degrees( CYP_xyz, CZP_xyz, NZ_xyz, CAZ_xyz) ;
		TR<< "CZP_CYP_N_CA and CYP_CZP_N_CA: " << CZP_CYP_N_CA << ' ' << CYP_CZP_N_CA << std::endl;

		// compare to the correct answer
		//float TOLERATED_ERROR = 0.01;

		TS_ASSERT_LESS_THAN_EQUALS( phi, old_phi + ( small_angle_max / 2.0 ) );
		TS_ASSERT_LESS_THAN_EQUALS( old_phi - ( small_angle_max / 2.0 ), phi );

		TS_ASSERT_LESS_THAN_EQUALS( psi, old_psi + ( small_angle_max / 2.0 ) );
		TS_ASSERT_LESS_THAN_EQUALS( old_psi - ( small_angle_max / 2.0 ), psi );

		//kdrew: test hydrogen placement

		TR << "test_OopPuckMoverSmall completed!! " << std::endl;
	}
};
