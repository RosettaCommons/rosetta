// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/RamaPrePro.cxxtest.hh
/// @brief  Unit tests for the RamaPrePro energy.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <test/util/pdb1ubq.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/scoring/RamaPrePro.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

// Protocol Headers
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <sstream>

static THREAD_LOCAL basic::Tracer TR("RamaPreProTests");


class RamaPreProTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options("-output_virtual true -write_all_connect_info true");
	}

	void tearDown(){
	}

	void do_test(
		std::string const &seq,
		bool const add_nmethyl,
		bool const flip_chirality
	) {
		core::pose::PoseOP pose( new core::pose::Pose );
		core::pose::make_pose_from_sequence(*pose, seq, "fa_standard", false);

		if ( add_nmethyl ) {
			core::chemical::ResidueTypeSetCOP rsd_set(pose->residue_type_set_for_pose( pose->residue_type(2).mode() ));
			core::chemical::ResidueTypeCOP rsd_type( pose->residue_type_ptr(2) );
			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_residue_type_with_variant_added( *rsd_type,
				core::chemical::ResidueProperties::get_variant_from_string( "N_METHYLATION" ) ).get_self_ptr() );
			core::pose::replace_pose_residue_copying_existing_coordinates( *pose, 2, *new_rsd_type );
			//pose->dump_pdb("vtemp_ramaprepro.pdb"); //DELETE ME
		}

		if ( flip_chirality ) {
			protocols::cyclic_peptide::FlipChiralityMover flipper;
			flipper.apply( *pose );
		}

		core::scoring::RamaPrePro const & rama( core::scoring::ScoringManager::get_instance()->get_RamaPrePro() );
		TR << "\nPHI\tPSI\n";
		core::Size leftcount(0);
		for ( core::Size i=1; i<=1000; ++i ) {
			utility::vector1 < core::Real > phipsi;
			rama.random_mainchain_torsions( pose->conformation(), pose->residue_type_ptr(2), pose->residue_type_ptr(3), phipsi);
			TR << phipsi[1] << "\t" << phipsi[2] << "\n";
			if ( phipsi[1] <= 0 ) ++leftcount;
		}
		TR << "LEFT: " << leftcount << "\tRIGHT: " << 1000-leftcount << std::endl;
		if ( flip_chirality ) {
			TS_ASSERT( leftcount < 500 );
		} else {
			TS_ASSERT( leftcount > 500 );
		}
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution.
	void test_random_phipsi_canonical() {
		do_test("AAAA", false, false);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a pre-proline amino acid.
	void test_random_phipsi_canonical_prepro() {
		do_test("AAPA", false, false);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a D-amino acid.
	void test_random_phipsi_canonical_d_aa() {
		do_test("AAAA", false, true);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a D-pre-proline amino acid.
	void test_random_phipsi_canonical_d_aa_prepro() {
		do_test("AAPA", false, true);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a noncanonical (N-methyl-trp).
	void test_random_phipsi_noncanonical() {
		do_test("AWAA", true, false);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a noncanonical (N-methyl-trp)
	/// using the pre-proline map.
	void test_random_phipsi_noncanonical_prepro() {
		do_test("AWPA", true, false);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a noncanonical D-amino acid (D-N-methyl-trp).
	void test_random_phipsi_noncanonical_d_aa() {
		do_test("AWAA", true, true);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a noncanonical D-amino acid (D-N-methyl-trp)
	/// using the pre-proline map.
	void test_random_phipsi_noncanonical_d_aa_prepro() {
		do_test("AWPA", true, true);
	}

	/// @brief Checks that the rama_prepro score term is storing the value that it has calculated in the
	/// energies object of the pose proprely.
	void test_correct_score_stored() {
		core::pose::Pose pose( pdb1ubq5to13_pose() );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::rama_prepro, 1.0 );
		(sfxn)(pose);

		core::scoring::ScoringManager const &score_man( *(core::scoring::ScoringManager::get_instance()) );
		core::scoring::RamaPrePro const &rama( score_man.get_RamaPrePro() );
		utility::vector1< core::Real > gradient; //Unused, but needed below.

		core::Real old_score2(0), old_score3(0);

		for ( core::Size ir=2, irmax=pose.total_residue(); ir<irmax; ++ir ) {
			core::Real direct_calc(0.0), direct_calc_2(0.0);

			utility::vector1< core::Real > mainchain_tors(2);
			mainchain_tors[1] = pose.phi(ir);
			mainchain_tors[2] = pose.psi(ir);
			rama.eval_rpp_rama_score( pose.conformation(), pose.residue_type_ptr(ir), pose.residue_type_ptr(ir+1), mainchain_tors, direct_calc, gradient, false);
			if ( ir > 2 ) {
				utility::vector1< core::Real > mainchain_tors2(2);
				mainchain_tors2[1] = pose.phi(ir-1);
				mainchain_tors2[2] = pose.psi(ir-1);
				rama.eval_rpp_rama_score( pose.conformation(), pose.residue_type_ptr(ir-1), pose.residue_type_ptr(ir), mainchain_tors2, direct_calc_2, gradient, false);
			}

			// Note that there is a subtlety, here: because the rama_prepro energy is a two-body energy computed between residue i and i+1, the score is the sum of the rama_prepro score for
			// residue i and residue i-1.
			TS_ASSERT_DELTA( pose.energies().residue_total_energies(ir)[ core::scoring::rama_prepro ], (direct_calc+direct_calc_2) / 2.0, 0.0001 );

			//Store scores for residue 2 and 3
			if ( ir==2 ) {
				old_score2 = pose.energies().residue_total_energies(ir)[core::scoring::rama_prepro];
			} else if ( ir==3 ) {
				old_score3 = pose.energies().residue_total_energies(ir)[core::scoring::rama_prepro];
			}
		}

		//Now check that residues 2 and 3 still have the same score if we add a cutpoint variant.
		//pose.dump_scored_pdb( "rama_pdb_before.pdb", sfxn ); //DELETE ME
		core::Real const old_psi2(pose.psi(2)), old_omega2(pose.omega(2)), old_phi3(pose.phi(3));
		core::pose::correctly_add_cutpoint_variants(pose, 2, false, 3);
		core::kinematics::FoldTree new_foldtree;
		std::istringstream foldtree_setup;
		foldtree_setup.clear();
		foldtree_setup.str("FOLD_TREE EDGE 2 1 -1 EDGE 2 3 1 EDGE 3 9 -1");
		foldtree_setup >> new_foldtree;
		pose.conformation().fold_tree(new_foldtree);
		pose.set_psi(2, old_psi2);
		pose.set_omega(2, old_omega2);
		pose.set_phi(3, old_phi3);
		(sfxn)(pose);
		TS_ASSERT_DELTA( old_score2, pose.energies().residue_total_energies(2)[ core::scoring::rama_prepro ], 0.001);
		TS_ASSERT_DELTA( old_score3, pose.energies().residue_total_energies(3)[ core::scoring::rama_prepro ], 0.001);
		//pose.dump_scored_pdb( "rama_pdb_after.pdb", sfxn ); //DELETE ME
	}


};



