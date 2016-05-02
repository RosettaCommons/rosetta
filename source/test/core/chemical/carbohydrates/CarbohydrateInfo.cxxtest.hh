// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/chemical/carbohydrates/CarbohydrateInfo.cxxtest.hh
/// @brief   Test suite for CarbohydrateInfo
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("core.chemical.carbohydrates.CarbohydrateInfo.cxxtest");

class CarbohydrateInfoTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace core::pose;
		using namespace core::import_pose;
		using namespace basic::options;

		core_init_with_additional_options( "-include_sugars -out:level 500" );

		// Test that oligosaccharides are loaded correctly.
		pose_from_file( maltotriose_, "core/chemical/carbohydrates/maltotriose.pdb" , core::import_pose::PDB_file);
		pose_from_file( isomaltose_, "core/chemical/carbohydrates/isomaltose.pdb" , core::import_pose::PDB_file);

		// Test branched oligosaccharide.
		pose_from_file( branched_fragment_, "core/chemical/carbohydrates/amylopectin_fragment.pdb" , core::import_pose::PDB_file);

		// Test N-linked glycosylation.
		pose_from_file( N_linked_, "core/chemical/carbohydrates/glycosylated_peptide.pdb" , core::import_pose::PDB_file);

		// Test modified sugar patch system.
		pose_from_file( glucosamine_, "core/chemical/carbohydrates/GlcN.pdb" , core::import_pose::PDB_file);

		// Test a combination of the above.
		pose_from_file( Lex_, "core/chemical/carbohydrates/Lex.pdb" , core::import_pose::PDB_file);

		// Test that oligosaccharides can be created from a given sequence.
		//make_pose_from_saccharide_sequence( lactose_, "beta-D-Galp-(1->4)-Glcp" );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Confirm that CarbohydrateInfo.short_name_ is assigned correctly.
	void test_Pose_chain_sequence_w_polysaccharide()
	{
		TR << "Testing chain_sequence() method of Pose with polysaccharide chains."  << std::endl;
		TS_ASSERT_EQUALS( maltotriose_.chain_sequence( 1 ), "alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp" );
		TS_ASSERT_EQUALS( isomaltose_.chain_sequence( 1 ), "alpha-D-Glcp-(1->6)-D-Glcp" );
		//TS_ASSERT_EQUALS( lactose_.chain_sequence( 1 ), "beta-D-Galp-(1->4)-D-Glcp" );
		TS_ASSERT_EQUALS( glucosamine_.chain_sequence( 1 ), "D-GlcpN" );
	}

	// Confirm that backbone torsion angles are assigned correctly.
	void test_Pose_phi_psi_omega_w_polysaccharide()
	{
		TR << "Testing phi(), psi(), and omega() methods of Pose with polysaccharide chains." << std::endl;
		TS_ASSERT_DELTA( maltotriose_.phi( 1 ), 0.000, 0.02 );

		TS_ASSERT_DELTA( isomaltose_.phi( 2 ), 44.3268, 0.02 );
		TS_ASSERT_DELTA( isomaltose_.psi( 2 ), -170.869, 0.02 );
		TS_ASSERT_DELTA( isomaltose_.omega( 2 ), 49.383, 0.02 );

		TS_ASSERT_DELTA( branched_fragment_.phi( 5 ), 111.187, 0.02 );

		TS_ASSERT_DELTA( N_linked_.phi( 6 ), -103.691, 0.02 );

		TS_ASSERT_DELTA( Lex_.phi( 1 ), 0.0, 0.1 );  // undefined torsion
		TS_ASSERT_DELTA( Lex_.psi( 1 ), 0.0, 0.1 );  // undefined torsion
		TS_ASSERT_DELTA( Lex_.omega( 1 ), 0.0, 0.1 );  // undefined torsion
		TS_ASSERT_DELTA( Lex_.phi( 2 ), -85.8, 0.1 );
		TS_ASSERT_DELTA( Lex_.psi( 2 ), 135.6, 0.1 );
		TS_ASSERT_DELTA( Lex_.omega( 2 ), 0.0, 0.1 );  // undefined torsion
		TS_ASSERT_DELTA( Lex_.phi( 3 ), -76.9, 0.1 );
		TS_ASSERT_DELTA( Lex_.psi( 3 ), -97.0, 0.1 );
		TS_ASSERT_DELTA( Lex_.omega( 3 ), 0.0, 0.1 );  // undefined torsion

		Lex_.set_phi( 1, 10.0 );  // undefined torsion; should be ignored
		Lex_.set_psi( 1, 20.0 );  // undefined torsion; should be ignored
		Lex_.set_omega( 1, 30.0 );  // undefined torsion; should be ignored
		Lex_.set_phi( 2, 40.0 );
		Lex_.set_psi( 2, 50.0 );
		Lex_.set_omega( 2, 60.0 );  // undefined torsion; should be ignored
		Lex_.set_phi( 3, 70.0 );
		Lex_.set_psi( 3, 80.0 );
		Lex_.set_omega( 3, 90.0 );  // undefined torsion; should be ignored

		TS_ASSERT_DELTA( Lex_.phi( 1 ), 0.0, 0.1 );  // undefined torsion
		TS_ASSERT_DELTA( Lex_.psi( 1 ), 0.0, 0.1 );  // undefined torsion
		TS_ASSERT_DELTA( Lex_.omega( 1 ), 0.0, 0.1 );  // undefined torsion
		TS_ASSERT_DELTA( Lex_.phi( 2 ), 40.0, 0.1 );
		TS_ASSERT_DELTA( Lex_.psi( 2 ), 50.0, 0.1 );
		TS_ASSERT_DELTA( Lex_.omega( 2 ), 0.0, 0.1 );  // undefined torsion
		TS_ASSERT_DELTA( Lex_.phi( 3 ), 70.0, 0.1 );
		TS_ASSERT_DELTA( Lex_.psi( 3 ), 80.0, 0.1 );
		TS_ASSERT_DELTA( Lex_.omega( 3 ), 0.0, 0.1 );  // undefined torsion
	}

	// Confirm that side-chain torsion angles are assigned correctly.
	void test_Pose_chi_w_polysaccharide()
	{
		TR << "Testing chi() method of Pose with polysaccharide chains."  << std::endl;
		//TS_ASSERT_DELTA( maltotriose_.chi( 1, 2 ), 0.000, 0.02 );
		TS_ASSERT_DELTA( maltotriose_.chi( 2, 2 ), -179.959, 0.02 );
		TS_ASSERT_DELTA( maltotriose_.chi( 3, 2 ), 175.924, 0.02 );
		TS_ASSERT_DELTA( maltotriose_.chi( 4, 2 ), maltotriose_.psi( 3 ), 0.02 );
		TS_ASSERT_DELTA( maltotriose_.chi( 5, 2 ), -161.7, 0.02 );
		TS_ASSERT_DELTA( maltotriose_.chi( 6, 2 ), -178.781, 0.02 );
	}

	// Confirm that branches are handled properly.
	void test_CarbohydrateInfo_branch_point()
	{
		using namespace core;
		using namespace conformation;

		TR << "Testing branch_point() method of CarbohydrateInfo."  << std::endl;
		Residue const & res2( branched_fragment_.residue( 2 ) );
		TS_ASSERT_EQUALS( res2.carbohydrate_info()->branch_point( 1 ), 6 );
	}


private:  // Private data /////////////////////////////////////////////////////
	core::pose::Pose maltotriose_;  // a (1alpha->4) trisaccharide of D-glucose
	core::pose::Pose isomaltose_;  // a (1alpha->6) disaccharide of D-glucose
	core::pose::Pose branched_fragment_;  // a branched fragment of amylopectin
	core::pose::Pose N_linked_;  // a 5-mer peptide with an N-linked glycan
	core::pose::Pose lactose_;  // a (1beta->4) disaccharide of D-glucose and D-galactose
	core::pose::Pose glucosamine_;  // 2-amino-2-deoxy-D-glucopyranose
	core::pose::Pose Lex_;  // Lewisx: beta-D-Galp-(1->4)-[alpha-D-Fucp-(1->3)]-D-GlcpNAc

};  // class CarbohydrateInfoTests
