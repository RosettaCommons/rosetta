// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/chemical/CtermAmidationTests.cxxtest.hh
/// @brief  Test whether the Cterm_amidation patch works as advertised.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/chemical/ResidueType.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("CtermAmidationTests");


class CtermAmidationTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();
	}

	void tearDown(){
	}

	void test_cterm_amidation_on_each_aa(){
		TR << "Starting CtermAmidationTests::test_cterm_amidation_on_each_aa()." << std::endl;

		for ( core::Size i( static_cast<core::Size>(core::chemical::first_l_aa) );
				i<=static_cast<core::Size>( core::chemical::num_canonical_aas );
				++i
				) {
			std::string const aaname( core::chemical::name_from_aa( static_cast< core::chemical::AA >(i) ) );
			TR << "\tTrying CTERM_AMIDATION on restype " << aaname << ":CtermProteinFull..." << std::endl;

			core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
			{ //Generate the pose.
				protocols::cyclic_peptide::PeptideStubMover stubmover;
				stubmover.set_reset_mode(true);
				stubmover.add_residue( protocols::cyclic_peptide::PSM_append, "GLY:NtermProteinFull", 1, true, "", 0, 0, nullptr, "" );
				stubmover.add_residue( protocols::cyclic_peptide::PSM_append, "GLY", 2, false, "", 0, 0, nullptr, "" );
				stubmover.add_residue( protocols::cyclic_peptide::PSM_append, aaname+":CtermProteinFull", 3, false, "", 0, 0, nullptr, "" );
				stubmover.apply(*pose);
			}

			{ //Add the variant type.
				protocols::simple_moves::ModifyVariantTypeMover modvartype;
				core::select::residue_selector::ResidueIndexSelectorOP select_three(
					utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector>()
				);
				select_three->append_index(3);
				modvartype.set_additional_type_to_add( "CTERM_AMIDATION" );
				modvartype.set_residue_selector(select_three);
				TS_ASSERT_THROWS_NOTHING( modvartype.apply(*pose) );
				TS_ASSERT( pose->residue_type(3).has("NT") );
			}

			for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
				pose->set_phi( ir, -61 );
				pose->set_psi( ir, -41 );
				pose->set_omega( ir, 179 );
			}

			//pose->dump_pdb( "VTEMP_" + aaname + ".pdb"  ); //DELETE ME
		}
		TR << "Finished CtermAmidationTests::test_cterm_amidation_on_each_aa()." << std::endl;
	}

	void test_cterm_amidation_on_each_aa_notermini(){
		TR << "Starting CtermAmidationTests::test_cterm_amidation_on_each_aa_notermini()." << std::endl;

		for ( core::Size i( static_cast<core::Size>(core::chemical::first_l_aa) );
				i<=static_cast<core::Size>( core::chemical::num_canonical_aas );
				++i
				) {
			std::string const aaname( core::chemical::name_from_aa( static_cast< core::chemical::AA >(i) ) );
			TR << "\tTrying CTERM_AMIDATION on restype " << aaname << "..." << std::endl;

			core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
			{ //Generate the pose.
				protocols::cyclic_peptide::PeptideStubMover stubmover;
				stubmover.set_reset_mode(true);
				stubmover.add_residue( protocols::cyclic_peptide::PSM_append, "GLY", 1, true, "", 0, 0, nullptr, "" );
				stubmover.add_residue( protocols::cyclic_peptide::PSM_append, "GLY", 2, false, "", 0, 0, nullptr, "" );
				stubmover.add_residue( protocols::cyclic_peptide::PSM_append, aaname, 3, false, "", 0, 0, nullptr, "" );
				stubmover.apply(*pose);
			}

			{ //Add the variant type.
				protocols::simple_moves::ModifyVariantTypeMover modvartype;
				core::select::residue_selector::ResidueIndexSelectorOP select_three(
					utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector>()
				);
				select_three->append_index(3);
				modvartype.set_additional_type_to_add( "CTERM_AMIDATION" );
				modvartype.set_residue_selector(select_three);
				TS_ASSERT_THROWS_NOTHING( modvartype.apply(*pose) );
				TS_ASSERT( pose->residue_type(3).has("NT") );
			}

			for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
				pose->set_phi( ir, -61 );
				pose->set_psi( ir, -41 );
				pose->set_omega( ir, 179 );
			}

			//pose->dump_pdb( "VTEMP_" + aaname + ".pdb"  ); //DELETE ME
		}
		TR << "Finished CtermAmidationTests::test_cterm_amidation_on_each_aa_notermini()." << std::endl;
	}



};
