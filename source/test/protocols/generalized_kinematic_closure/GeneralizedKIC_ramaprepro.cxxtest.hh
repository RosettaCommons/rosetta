// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/generalized_kinematic_closure/GeneralizedKIC.cxxtest.hh
/// @brief  Unit tests for the GeneralizedKIC mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers:
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pdb1ubq.hh>

// GeneralizedKIC headers:
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>
#include <protocols/generalized_kinematic_closure/util.hh>

// Other Rosetta libraries:
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>

// Basic Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("GeneralizedKIC_Tests_ramaprepro");


// --------------- Test Class --------------- //

class GeneralizedKIC_Tests_ramaprepro : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
	core::scoring::ScoreFunctionOP scorefxn_;

public:

	void setUp() {
		// amw: no longer have to load Ds
		core_init(); // d-caa/DALA.params" );

		scorefxn_ = core::scoring::get_score_function();
		testpose_ = pdb1ubq5to13_poseop();

	}

	void tearDown() {
	}

	/// @brief Test protocol used by the next eight or so sub-tests.
	void do_ramaprepro_test(
		std::string const &seq,
		core::Size const sampled_pos,
		bool const n_methylate,
		bool const flip_chirality
	) {
		using namespace protocols::generalized_kinematic_closure;

		core::Size const ntrials(50);

		//Make the pose:
		core::pose::PoseOP pose( new core::pose::Pose );
		core::pose::make_pose_from_sequence(*pose, seq, "fa_standard", false);

		if ( n_methylate ) {
			core::chemical::ResidueTypeSetCOP rsd_set( pose->residue_type_set_for_pose( pose->residue_type(sampled_pos).mode() ) );
			core::chemical::ResidueTypeCOP rsd_type( pose->residue_type_ptr(sampled_pos) );
			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_residue_type_with_variant_added( *rsd_type,
				core::chemical::ResidueProperties::get_variant_from_string( "N_METHYLATION" ) ).get_self_ptr() );
			core::pose::replace_pose_residue_copying_existing_coordinates( *pose, sampled_pos, *new_rsd_type );
			//pose->dump_pdb("vtemp_genkic_ramaprepro.pdb"); //DELETE ME
		}

		if ( flip_chirality ) {
			protocols::cyclic_peptide::FlipChiralityMover flipper;
			flipper.apply( *pose );
		}

		TR << "\nPHI\tPSI\n";
		core::Size leftcount(0);

		for ( core::Size repeats=1; repeats<=ntrials; ++repeats ) {

			//Fully randomize mainchain torsions:
			for ( core::Size ir=1, irmax=pose->total_residue(); ir<=irmax; ++ir ) {
				pose->set_omega(ir, 180);
				pose->set_phi(ir, numeric::random::rg().uniform()*360.0-180.0);
				pose->set_psi(ir, numeric::random::rg().uniform()*360.0-180.0);
			}

			GeneralizedKICOP genkic( new GeneralizedKIC ); //Create the mover.

			//Add the loop residues:
			for ( core::Size i=2; i<pose->total_residue(); ++i ) genkic->add_loop_residue(i);

			//Set the pivots:
			genkic->set_pivot_atoms( 2, "CA", 9, "CA", pose->total_residue()-1, "CA" );

			//Add a rama_prepro filter on the pivots:
			utility::vector1 <core::Size> pivresidues(3);
			pivresidues[1] = 2; pivresidues[2] = 9; pivresidues[3] = pose->total_residue() - 1;
			for ( core::Size i=1; i<=3; ++i ) {
				genkic->add_filter( "rama_prepro_check" );
				genkic->set_filter_resnum(pivresidues[i]);
				genkic->set_filter_rama_cutoff_energy( 100.0 ); //I don't actually want to reject anything; I just want to make sure that the filter doesn't crash.
			}

			//Add a randomizing-by-rama-prepro perturber:
			genkic->add_perturber( "randomize_backbone_by_rama_prepro" );
			genkic->add_residue_to_perturber_residue_list( sampled_pos );

			//Add full randomization for the other residues:
			genkic->add_perturber( "randomize_dihedral" );
			for ( core::Size ir=2, irmax=pose->total_residue() -1; ir<=irmax; ++ir ) {
				if ( ir == sampled_pos ) continue;
				utility::vector1< core::id::NamedAtomID > tors, tors2;
				tors.push_back( core::id::NamedAtomID("N", ir) );
				tors.push_back( core::id::NamedAtomID("CA", ir) );
				tors2.push_back( core::id::NamedAtomID("CA", ir) );
				tors2.push_back( core::id::NamedAtomID("C", ir) );
				genkic->add_atomset_to_perturber_atomset_list( tors );
				genkic->add_atomset_to_perturber_atomset_list( tors2 );
			}

			//Set options:
			genkic->set_closure_attempts(100); //Try a maximum of 10 times.
			genkic->set_min_solution_count(1); //Stop when a solution is found.

			//Add a lowest_energy selector:
			genkic->set_selector_type("random_selector");
			genkic->apply(*pose);

			if ( !genkic->last_run_successful() ) {
				--repeats;
				continue;
			}

			TR << pose->phi(sampled_pos) << "\t" << pose->psi(sampled_pos) << "\n";
			if ( pose->phi(sampled_pos) <= 0 ) ++leftcount;
		}
		TR << "LEFT: " << leftcount << "\tRIGHT: " << ntrials-leftcount << std::endl;
		if ( flip_chirality ) {
			TS_ASSERT( leftcount < static_cast<core::Size>(static_cast<core::Real>(ntrials)/2.0) );
		} else {
			TS_ASSERT( leftcount > static_cast<core::Size>(static_cast<core::Real>(ntrials)/2.0) );
		}
	}

	/// @brief Test GeneralizedKIC randomizing using RamaPrePro tables, for a canonical L-amino
	/// acid that is not adjacent to a proline.
	void test_GeneralizedKIC_ramaprepro_perturbation_canonical() {
		do_ramaprepro_test( "AAAAAAAAAAAAAAA", 5, false, false );
	}

	/// @brief Test GeneralizedKIC randomizing using RamaPrePro tables, for a canonical D-amino
	/// acid that is not adjacent to a proline.
	void test_GeneralizedKIC_ramaprepro_perturbation_canonical_D() {
		do_ramaprepro_test( "AAAAAAAAAAAAAAA", 5, false, true );
	}

	/// @brief Test GeneralizedKIC randomizing using RamaPrePro tables, for a canonical L-amino
	/// acid that is adjacent to a proline.
	void test_GeneralizedKIC_ramaprepro_perturbation_canonical_beforepro() {
		do_ramaprepro_test( "AAAAAPAAAAAAAAA", 5, false, false );
	}

	/// @brief Test GeneralizedKIC randomizing using RamaPrePro tables, for a canonical D-amino
	/// acid that is adjacent to a proline.
	void test_GeneralizedKIC_ramaprepro_perturbation_canonical_D_beforepro() {
		do_ramaprepro_test( "AAAAAPAAAAAAAAA", 5, false, true );
	}

}; //class GeneralizedKIC_Tests
