// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/cyclic_peptide/AchiralResidueTypeTestHeaders.hh
/// @brief  Code for achiral residue type tests.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_test_protocols_cyclic_peptide_AchiralResidueTypeTestHeaders_HH
#define INCLUDED_test_protocols_cyclic_peptide_AchiralResidueTypeTestHeaders_HH

// Test headers
#include <test/util/pdb1ubq.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>

namespace test {
namespace protocols {
namespace cyclic_peptide {

//static THREAD_LOCAL basic::Tracer TR("protocols.cyclic_peptide.AchiralResidueTypeTestHeaders");

class AchiralResidueTypeTestHelper {

public:

	void test_import(
		basic::Tracer &TR,
		std::string const &filename,
		core::Size const res_index,
		std::string const &expected_3letter_code,
		std::string const &expected_full_name
	) {
		core::pose::PoseOP pose( core::import_pose::pose_from_file( filename, false, core::import_pose::PDB_file) );
		TR << "\nRES" << res_index << "-name3\tRES" << res_index << "-name\tEXP-name3\tEXP-name\n";
		TR << pose->residue(res_index).name3() << "\t" << pose->residue(res_index).name() << "\t" << expected_3letter_code << "\t" << expected_full_name << std::endl;
		TS_ASSERT_EQUALS( pose->residue(res_index).name3(), expected_3letter_code );
		TS_ASSERT_EQUALS( pose->residue(res_index).name(), expected_full_name );
	}

	void test_mirror_symmetry(
		basic::Tracer &TR,
		std::string const &residue_type_name,
		bool const prepro,
		bool const n_methylate
	) {
		core::pose::PoseOP pose( core::import_pose::pose_from_file( "protocols/cyclic_peptide/c4m_test_pose.pdb", false, core::import_pose::PDB_file) );

		::protocols::simple_moves::MutateResidue mut3( 2, residue_type_name );
		mut3.apply(*pose);
		if ( n_methylate ) {
			::protocols::simple_moves::ModifyVariantTypeMover modvartype;
			::core::select::residue_selector::ResidueIndexSelectorOP index( new ::core::select::residue_selector::ResidueIndexSelector );
			index->set_index("2");
			modvartype.set_residue_selector(index);
			modvartype.set_additional_type_to_add( "N_METHYLATION" );
			modvartype.apply(*pose);
		}
		if ( prepro ) {
			::protocols::simple_moves::MutateResidue mut4( 3, "PRO" );
			mut4.apply(*pose);
		}

		core::pose::PoseOP pose2( pose->clone() );

		::protocols::cyclic_peptide::FlipChiralityMover flip;
		flip.set_center( 0, 0, 0 );
		flip.set_normal( 0, 0, 1 );
		flip.apply( *pose2 );

		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() ); //Get whatever the current scorefunction is.
		core::scoring::symmetry::SymmetricScoreFunctionOP sfxn_symm( core::scoring::symmetry::symmetrize_scorefunction( *sfxn ) ); //Set up a symmetric version
		core::Real const asymm_energy_1( (*sfxn)(*pose) );
		core::Real const asymm_energy_2( (*sfxn)(*pose) );
		TR << "E1=" << asymm_energy_1 << std::endl;
		TR << "E2=" << asymm_energy_2 << std::endl;
		TS_ASSERT_DELTA( asymm_energy_1, asymm_energy_2, 0.00001 ); //Make sure that the asymmetric pose and its mirror image score identically

		::protocols::simple_moves::symmetry::SetupForSymmetryMover symm_setup("protocols/cyclic_peptide/c4m.symm");
		symm_setup.apply(*pose);
		symm_setup.apply(*pose2);

		::protocols::cyclic_peptide::DeclareBond bond;
		bond.set( 4, "C", 5, "N", false );
		bond.apply(*pose);
		bond.apply(*pose2);
		core::Real const symm_energy_1( (*sfxn_symm)(*pose) );
		core::Real const symm_energy_2( (*sfxn_symm)(*pose2) );
		TR << "E1_symm=" << symm_energy_1 << std::endl;
		TR << "E2_symm=" << symm_energy_2 << std::endl;
		TS_ASSERT_DELTA( symm_energy_1, symm_energy_2, 0.00001 ); //Make sure that the symmetric pose and its mirror image score identically

		//DELETE ME -- just for testing:
		//pose->dump_pdb( "symmtemp1.pdb");
		//pose2->dump_pdb( "symmtemp2.pdb");
	}

	/// @brief Given a residue type, cycle through its rama map and ensure that it is symmetric.
	///
	void test_symmetric_rama_prepro_scoring(
		basic::Tracer &TR,
		std::string const &residue_type_name,
		bool const prepro,
		bool const n_methylate,
		core::Real const &threshold
	) {
		core::pose::PoseOP pose( pdb1ubq5to13_poseop() );
		::protocols::simple_moves::MutateResidue mut2( 2, residue_type_name ); //Position 2 is the residue, not before a pro
		::protocols::simple_moves::MutateResidue mut4( 4, residue_type_name ); //Position 4 is the residue, before a pro
		::protocols::simple_moves::MutateResidue mut5( 5, "PRO" );
		mut2.apply(*pose);
		mut4.apply(*pose);
		mut5.apply(*pose);

		if ( n_methylate ) {
			::protocols::simple_moves::ModifyVariantTypeMover modvartype;
			::core::select::residue_selector::ResidueIndexSelectorOP index( new ::core::select::residue_selector::ResidueIndexSelector );
			index->set_index("2,4");
			modvartype.set_residue_selector(index);
			modvartype.set_additional_type_to_add( "N_METHYLATION" );
			modvartype.apply(*pose);
		}

		core::scoring::ScoringManager const &score_man( *(core::scoring::ScoringManager::get_instance()) );
		core::scoring::RamaPrePro const &rama( score_man.get_RamaPrePro() );


		TR << "Testing rama_prepro scoring of residue type " << residue_type_name << " when it is" << (prepro == false ? " not " : " ") <<  "before a proline.\n";
		TR << "PHI\tPSI\tE\tE_mirror\n";
		utility::vector1< core::Real > gradient; //Unused but needed for eval_rpp_rama_score().

		core::Size const curpos( prepro == false ? 2 : 4 );

		for ( core::Real phi(-179.0), phimax(180.0); phi<=phimax; phi+=3.0 ) {
			for ( core::Real psi(-179.0), psimax(180.0); psi<=psimax; psi+=3.0 ) {
				core::Real energy1(0.0), energy2(0.0);
				utility::vector1< core::Real > mainchain_tors(2, 0.0);
				mainchain_tors[1] = phi; mainchain_tors[2] = psi;
				rama.eval_rpp_rama_score( pose->conformation(), pose->residue_type(curpos).get_self_ptr(), pose->residue_type(curpos+1).get_self_ptr(), mainchain_tors, energy1, gradient, false); //Get energy
				mainchain_tors[1] = -1.0*phi; mainchain_tors[2] = -1.0*psi;
				rama.eval_rpp_rama_score( pose->conformation(), pose->residue_type(curpos).get_self_ptr(), pose->residue_type(curpos+1).get_self_ptr(), mainchain_tors, energy2, gradient, false); //Get energy of mirror-image conformation.
				TR << phi << "\t" << psi << "\t" << energy1 << "\t" << energy2 << "\n";
				TS_ASSERT_DELTA( energy1, energy2, threshold);
			}
		}
		TR.flush();


	}

};

} //namespace cyclic_peptide
} //namespace protocols
} //namespace test

#endif //INCLUDED_test_protocols_cyclic_peptide_NMethylation_functions_HH
