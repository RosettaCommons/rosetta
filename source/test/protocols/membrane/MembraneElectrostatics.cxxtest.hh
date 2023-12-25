// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/MembraneElectrostatics.cxxtest.hh
/// @brief  Tests the electrostatic features of membrane
/// @author Rituparnasamanta (rituparna@utexas.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/AddMembraneMoverCreator.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/energy_methods/ImplicitMembraneCoulomb.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>

static basic::Tracer TR("MembraneElectrostatics");


class MembraneElectrostatics : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();

	}

	void tearDown() {

	}


	/// @brief Calculate the energy of a protein at the center and normal to the membrane.
	void test_energy_components() {

		using namespace core;
		using namespace core::scoring;
		using namespace protocols::membrane;
		using namespace core::import_pose;
		using namespace core::pose;
		using namespace core::conformation;
		using namespace core::conformation::membrane;
		using namespace protocols::membrane_benchmark;

		// 1. TM domain of the M2 proton channel (single helix): 1mp6
		PoseOP m2_pose ( utility::pointer::make_shared< Pose >() );
		pose_from_file( *m2_pose, "protocols/membrane/1mp6_transformed.pdb" , core::import_pose::PDB_file);

		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1mp6.span" );
		add_memb1->apply( *m2_pose );

		MembraneEnergyLandscapeSamplerOP test_memb = utility::pointer::make_shared< MembraneEnergyLandscapeSampler >();
		TR << "This is a unit test for: "<< test_memb->get_name()<< std::endl;
		//2. Set up the score function
		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "franklin2023" );

		//3. align the protein to the center of the membrane and perpendicular
		core::Real limit1( 0 );
		core::Size membrane_jump( m2_pose->conformation().membrane_info()->membrane_jump() );

		// Apply an initial translation of the membrane center
		// Perform an initial transformation of the pose into the membrane
		TransformIntoMembraneMoverOP transform_into_memb( utility::pointer::make_shared< TransformIntoMembraneMover >() );
		transform_into_memb->apply( *m2_pose );

		Vector peptide_center( 0,0,0 );
		core::Real flag_for_center( 0.0 );
		peptide_center = test_memb->getcenter( *m2_pose, flag_for_center );

		// Apply an initial translation of the membrane center
		Vector initial_move( 0, 0, limit1 );
		Vector difference_in_center = initial_move - peptide_center;
		TranslationMoverOP adjusting_center( utility::pointer::make_shared< TranslationMover >( difference_in_center, membrane_jump ) );
		adjusting_center->apply( *m2_pose );
		peptide_center = test_memb->getcenter( *m2_pose, flag_for_center );

		//Align the protein with the z-axis as the initial position
		Vector peptide_normal( 0,0,0 );
		peptide_normal = test_memb->getaxis( *m2_pose, flag_for_center );
		Vector axis( 0,0,1 );//rotation about YZ
		RotationMoverOP align_axis( utility::pointer::make_shared< RotationMover >( peptide_normal, axis, peptide_center, membrane_jump ) );
		align_axis->apply( *m2_pose );
		//Vector rot_center = pose_tm_com( m2_pose );

		// core::Real total_score =
		scorefxn->score(*m2_pose);
		TR << "name of residue 3: " << m2_pose->residue(3).name() << std::endl;
		TR << "position: " << m2_pose->residue(3).xyz(3).z() << std::endl;
		TR << " fa_water_bilayer: " << m2_pose->energies().residue_total_energies(3)[fa_water_to_bilayer] << std::endl;
		TR << " fa_imm_elec: " << m2_pose->energies().residue_total_energies(3)[fa_imm_elec] << std::endl;
		TR << " fa_elec_lipidlayer " << m2_pose->energies().residue_total_energies(3)[f_elec_lipidlayer] << std::endl;

		TS_ASSERT_DELTA( m2_pose->energies().residue_total_energies(3)[fa_water_to_bilayer], 0.95 , 0.5 );
		TS_ASSERT_DELTA( m2_pose->energies().residue_total_energies(3)[fa_imm_elec], -0.762, 0.50);
		TS_ASSERT_DELTA( m2_pose->energies().residue_total_energies(3)[f_elec_lipidlayer], -0.815, 0.50);

		core::Vector translation_delta( 0, 0, -10 ); // Angstroms
		TranslationMoverOP translate_memb( utility::pointer::make_shared< TranslationMover >( translation_delta, membrane_jump ) );
		translate_memb->apply( *m2_pose );
		peptide_normal = test_memb->getaxis( *m2_pose, flag_for_center);
		peptide_center = test_memb->getcenter( *m2_pose, flag_for_center );
		TR << " peptide normal after alignment is ::" << peptide_normal.x() << "i+ " << peptide_normal.y() << "j +" << peptide_normal.z() << "k" << std::endl;
		TR << " peptide center after alignment is ::" << peptide_center.x() << "i+ " << peptide_center.y() << "j +" << peptide_center.z() << "k" << std::endl;
		TS_ASSERT_DELTA( std::abs(peptide_center.z()), 10.0 , 0.1 );
		TS_ASSERT_DELTA( translation_delta, peptide_center, 0.1);

		TS_ASSERT(true);

	}

	/// @brief Calculate the electrostatic term due to low dielectric constant in the membrane.
	void test_dielectric_constant_calc() {

		using namespace core;
		using namespace protocols::membrane;
		using namespace core::import_pose;
		using namespace core::pose;

		using namespace protocols::membrane;
		using namespace core::energy_methods;

		ImplicitMembraneCoulombOP imc( new ImplicitMembraneCoulomb() );
		//verifying dielectric constant
		TS_ASSERT_DELTA( imc->compute_depth_and_bilayer_dep_dielectric( 1.0, 1.0, 10.0 ), 62.38, 0.5 );
		TS_ASSERT_DELTA( imc->compute_depth_and_bilayer_dep_dielectric( 0.9, 0.9, 6.4 ), 37.43, 0.5 );
		TS_ASSERT_DELTA( imc->compute_depth_and_bilayer_dep_dielectric( 0.0, 0.0, 0.0 ), 3.0, 0.5 );
		TS_ASSERT_DELTA( imc->compute_depth_and_bilayer_dep_dielectric( 0.5, 0.5, 6.4 ), 23.59, 0.5 );

		// Test electrostatics calculation
		core::Vector i_xyz( 0, 0, 0 );
		core::Vector j_xyz_close( 0, 0, 2.0 );
		core::Vector j_xyz_far( 0, 0, 5.0 );

		core::Real pos_charge( 1.0 );
		core::Real neg_charge( -1.0 );

		core::Real in_memb( 0.5 );
		core::Real in_aqueous( 1.0 );

		// Big Charge
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_aqueous, j_xyz_close, pos_charge, in_aqueous ) , 0.0 , 0.01 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_aqueous, j_xyz_close, neg_charge, in_aqueous ) , 0.0, 0.01 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_aqueous, j_xyz_far, pos_charge, in_aqueous ) , 0.0, 0.01 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_aqueous, j_xyz_far, neg_charge, in_aqueous ) , 0.0, 0.01 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_memb, j_xyz_close, pos_charge, in_memb ) , 6.91, 0.1 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_memb, j_xyz_close, neg_charge, in_memb ) , -6.91, 0.1 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_memb, j_xyz_far, pos_charge, in_memb ) , 0.215, 0.1 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_memb, j_xyz_far, neg_charge, in_memb ) , -0.215, 0.1 );

		// Small Chrge
		pos_charge = 0.090;
		neg_charge = -0.604;

		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_aqueous, j_xyz_close, neg_charge, in_aqueous ), 0.0 , 0.01 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_aqueous, j_xyz_far, neg_charge, in_aqueous ), 0.0 , 0.01 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_memb, j_xyz_close, neg_charge, in_memb ), -0.387 , 0.1 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_memb, j_xyz_close, neg_charge, in_aqueous ), -0.191 , 0.1 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_memb, j_xyz_far, neg_charge, in_memb ), -0.02607, 0.1 );
		TS_ASSERT_DELTA( imc->eval_atom_atom_fa_elecE( i_xyz, pos_charge, in_memb, j_xyz_far, neg_charge, in_aqueous ), -0.01203, 0.1 );

		TS_ASSERT(true)
	}


};
