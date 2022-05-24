// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/crosslinker/Thioether_Helper_Tests.cxxtest.hh
/// @brief  Unit tests for thioether functionality in Rosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/cyclic_peptide/CrosslinkerMover.hh>
#include <protocols/backbone_moves/RandomizeBBByRamaPrePro.hh>
#include <protocols/simple_moves/DeclareBond.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("Thioether_Helper_Tests");


class Thioether_Helper_Tests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init_with_additional_options( "-symmetric_gly_tables true" );
	}

	void tearDown() {

	}

	protocols::relax::FastRelaxOP get_frlx() {
		return utility::pointer::make_shared< protocols::relax::FastRelax >(
			core::scoring::ScoreFunctionFactory::get_instance()->create_score_function( "ref2015_cst.wts" ),
			1
		);
	}

	void configure_genkic_for_alpha_aa( protocols::generalized_kinematic_closure::GeneralizedKIC & genkic ) {
		genkic.set_selector_scorefunction( core::scoring::get_score_function() );
		genkic.set_selector_type( protocols::generalized_kinematic_closure::selector::selector_type::lowest_energy_selector );
		genkic.set_closure_attempts( 1000 );
		genkic.set_min_solution_count( 1 );
		genkic.set_correct_polymer_dependent_atoms( true );
		genkic.add_loop_residue( 5 );
		genkic.add_loop_residue( 6 );
		genkic.add_loop_residue( 7 );
		genkic.add_loop_residue( 8 );
		genkic.add_loop_residue( 1 );
		genkic.add_loop_residue( 2 );
		genkic.add_loop_residue( 3 );
		genkic.add_tail_residue( 9 );
		genkic.add_tail_residue( 10 );
		genkic.set_pivot_atoms( 5, "CA", 8, "CB", 3, "CA" );
		genkic.close_bond( 8, "SG", 1, "CP2", 0, "", 0, "", 1.826, 102.028, 112.284, 180.0, true, false );
		genkic.add_perturber( protocols::generalized_kinematic_closure::perturber::perturber_effect::randomize_backbone_by_rama_prepro );
		genkic.add_residue_to_perturber_residue_list( 5 );
		genkic.add_residue_to_perturber_residue_list( 6 );
		genkic.add_residue_to_perturber_residue_list( 7 );
		genkic.add_residue_to_perturber_residue_list( 2 );
		genkic.add_residue_to_perturber_residue_list( 3 );
		genkic.add_perturber( protocols::generalized_kinematic_closure::perturber::perturber_effect::randomize_dihedral );
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "N", 8 ),
			core::id::NamedAtomID( "CA", 8 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CA", 8 ),
			core::id::NamedAtomID( "CB", 8 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CB", 8 ),
			core::id::NamedAtomID( "SG", 8 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "SG", 8 ),
			core::id::NamedAtomID( "CP2", 1 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CP2", 1 ),
			core::id::NamedAtomID( "CO", 1 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "N", 1 ),
			core::id::NamedAtomID( "CA", 1 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CA", 1 ),
			core::id::NamedAtomID( "C", 1 )
			}
		);
		genkic.add_filter( protocols::generalized_kinematic_closure::filter::filter_type::loop_bump_check );
	}

	void configure_genkic_for_beta_aa( protocols::generalized_kinematic_closure::GeneralizedKIC & genkic ) {
		genkic.set_selector_scorefunction( core::scoring::get_score_function() );
		genkic.set_selector_type( protocols::generalized_kinematic_closure::selector::selector_type::lowest_energy_selector );
		genkic.set_closure_attempts( 1000 );
		genkic.set_min_solution_count( 1 );
		genkic.set_correct_polymer_dependent_atoms( true );
		genkic.add_loop_residue( 5 );
		genkic.add_loop_residue( 6 );
		genkic.add_loop_residue( 7 );
		genkic.add_loop_residue( 8 );
		genkic.add_loop_residue( 1 );
		genkic.add_loop_residue( 2 );
		genkic.add_loop_residue( 3 );
		genkic.add_tail_residue( 9 );
		genkic.add_tail_residue( 10 );
		genkic.set_pivot_atoms( 5, "CA", 8, "CB", 3, "CA" );
		genkic.close_bond( 8, "SG", 1, "CP2", 0, "", 0, "", 1.826, 102.028, 112.284, 180.0, true, false  );
		genkic.add_perturber( protocols::generalized_kinematic_closure::perturber::perturber_effect::randomize_backbone_by_rama_prepro );
		genkic.add_residue_to_perturber_residue_list( 5 );
		genkic.add_residue_to_perturber_residue_list( 6 );
		genkic.add_residue_to_perturber_residue_list( 7 );
		genkic.add_residue_to_perturber_residue_list( 2 );
		genkic.add_residue_to_perturber_residue_list( 3 );
		genkic.add_perturber( protocols::generalized_kinematic_closure::perturber::perturber_effect::randomize_dihedral );
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "N", 8 ),
			core::id::NamedAtomID( "CA", 8 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CA", 8 ),
			core::id::NamedAtomID( "CB", 8 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CB", 8 ),
			core::id::NamedAtomID( "SG", 8 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "SG", 8 ),
			core::id::NamedAtomID( "CP2", 1 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CP2", 1 ),
			core::id::NamedAtomID( "CO", 1 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "N", 1 ),
			core::id::NamedAtomID( "CA", 1 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CA", 1 ),
			core::id::NamedAtomID( "CM", 1 )
			}
		);
		genkic.add_atomset_to_perturber_atomset_list(
			utility::vector1< core::id::NamedAtomID >{
			core::id::NamedAtomID( "CM", 1 ),
			core::id::NamedAtomID( "C", 1 )
			}
		);
		genkic.add_filter( protocols::generalized_kinematic_closure::filter::filter_type::loop_bump_check );
	}

	void check_geometry( core::pose::Pose const & pose, bool const check_bondlength ) {
		core::Real const bondangle_N(
			numeric::angle_degrees(
			( check_bondlength ? pose.xyz( core::id::NamedAtomID( "SG", 8 ) ) : pose.xyz( core::id::NamedAtomID( "VTH", 1 ) ) ),
			pose.xyz( core::id::NamedAtomID( "CP2", 1 ) ),
			pose.xyz( core::id::NamedAtomID( "CO", 1 ) )
			)
		);
		core::Real const bondangle_C(
			numeric::angle_degrees(
			pose.xyz( core::id::NamedAtomID( "CB", 8 ) ),
			pose.xyz( core::id::NamedAtomID( "SG", 8 ) ),
			( check_bondlength ? pose.xyz( core::id::NamedAtomID( "CP2", 1 ) ) : pose.xyz( core::id::NamedAtomID( "V1", 8 ) ) )
			)
		);
		TS_ASSERT_DELTA( bondangle_C, 102.028, 3.0 );
		TS_ASSERT_DELTA( bondangle_N, 112.284, 3.0 );

		if ( !check_bondlength ) return;

		core::Real const bondlength(
			pose.xyz( core::id::NamedAtomID( "SG", 8 ) ).distance(
			pose.xyz( core::id::NamedAtomID( "CP2", 1 ) )
			)
		);
		core::Real const v1dist(
			pose.xyz( core::id::NamedAtomID( "V1", 8 ) ).distance(
			pose.xyz( core::id::NamedAtomID( "CP2", 1 ) )
			)
		);
		core::Real const vthdist(
			pose.xyz( core::id::NamedAtomID( "VTH", 1 ) ).distance(
			pose.xyz( core::id::NamedAtomID( "SG", 8 ) )
			)
		);

		TS_ASSERT_DELTA( bondlength, 1.827, 0.1 );
		TS_ASSERT_DELTA( v1dist, 0.0, 0.1 );
		TS_ASSERT_DELTA( vthdist, 0.0, 0.1 );
	}

	void set_up_geometry( core::pose::PoseOP pose, std::string const & resname, std::string const & cysresname ) {

		{
			//Build peptide
			protocols::cyclic_peptide::PeptideStubMover stubmover;
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 1, true, "", 1, 1, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 2, false, "", 1, 1, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 3, false, "", 1, 2, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 4, false, "", 1, 3, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 5, false, "", 1, 4, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 6, false, "", 1, 5, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 7, false, "", 1, 6, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, cysresname, 8, false, "", 1, 7, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 9, false, "", 1, 8, nullptr, "" );
			stubmover.add_residue( protocols::cyclic_peptide::PSM_StubMode::PSM_append, resname, 10, false, "", 1, 9, nullptr, "" );
			stubmover.apply( *pose );
		}

		{
			//Add termini
			protocols::simple_moves::DeclareBond declbond;
			declbond.set( 1, "C", 2, "N", true );
			declbond.apply( *pose );
		}

		{
			//Set omegas
			for ( core::Size i(1); i<pose->total_residue(); ++i ) {
				pose->set_omega(i, 180.0);
			}
			pose->update_residue_neighbors();
		}

		{
			//Randomize rama
			protocols::backbone_moves::RandomizeBBByRamaPrePro randrama;
			randrama.apply( *pose );
		}

		{
			//Crosslink
			core::select::residue_selector::ResidueIndexSelectorOP ressel(
				utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( "1,8" )
			);
			protocols::cyclic_peptide::CrosslinkerMover xlink;
			xlink.set_scorefxn( core::scoring::ScoreFunctionFactory::get_instance()->create_score_function( "ref2015_cst" ) );
			xlink.set_behaviour( true, true, false, false );
			xlink.set_filter_behaviour( false, false, false, 1.0, 1.0, 1.0 );
			xlink.set_linker_name( "thioether" );
			xlink.set_residue_selector( ressel );
			xlink.apply( *pose );
		}
	}

	void test_build_gly_thioether() {
		for ( core::Size i(1); i<=3; ++i ) {
			core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );

			set_up_geometry( pose, "GLY", "CYS" );

			check_geometry( *pose, false );

			{
				//Close with GenKIC
				protocols::generalized_kinematic_closure::GeneralizedKIC genkic;
				configure_genkic_for_alpha_aa(genkic);
				genkic.apply( *pose );
			}

			check_geometry( *pose, true );

			{
				//Relax
				protocols::relax::FastRelaxOP frlx( get_frlx() );
				frlx->apply( *pose );
			}

			check_geometry( *pose, true );
		}
	}

	void test_build_ala_thioether() {
		for ( core::Size i(1); i<=3; ++i ) {
			core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );

			set_up_geometry( pose, "ALA", "DCYS" );

			check_geometry( *pose, false );

			{
				//Close with GenKIC
				protocols::generalized_kinematic_closure::GeneralizedKIC genkic;
				configure_genkic_for_alpha_aa(genkic);
				genkic.apply( *pose );
			}

			check_geometry( *pose, true );

			{
				//Relax
				protocols::relax::FastRelaxOP frlx( get_frlx() );
				frlx->apply( *pose );
			}

			check_geometry( *pose, true );
		}
	}

	void test_build_b3a_thioether() {
		for ( core::Size i(1); i<=3; ++i ) {
			core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );

			set_up_geometry( pose, "B3A", "DCYS" );

			check_geometry( *pose, false );

			{
				//Close with GenKIC
				protocols::generalized_kinematic_closure::GeneralizedKIC genkic;
				configure_genkic_for_beta_aa(genkic);
				genkic.apply( *pose );
			}

			check_geometry( *pose, true );

			{
				//Relax
				protocols::relax::FastRelaxOP frlx( get_frlx() );
				frlx->apply( *pose );
			}

			check_geometry( *pose, true );
		}
	}


};
