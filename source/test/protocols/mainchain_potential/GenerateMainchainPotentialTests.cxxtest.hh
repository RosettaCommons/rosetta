// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/mainchain_potential/GenerateMainchainPotentialTests.cxxtest.hh
/// @brief  Unit tests for the GenerateMainchainPotential class, which underlies the make_mainchain_potential application.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/mainchain_potential/GenerateMainchainPotential.hh>
#include <protocols/mainchain_potential/GenerateMainchainPotentialOptions.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/mainchain_potential/MainchainScoreTable.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("GenerateMainchainPotentialTests");


class GenerateMainchainPotentialTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-ex1 -extrachi_cutoff 0");
	}

	void tearDown(){

	}

	/// @brief Test the initial creation of the one-residue pose.
	void test_pose_generation() {
		using namespace protocols::mainchain_potential;

		GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) );
		options->set_residue_name( "ALA" );
		GenerateMainchainPotential generator(options);
		core::pose::PoseOP pose( generator.generate_pose() );

		TS_ASSERT_EQUALS( pose->total_residue(), 1 );
		core::chemical::ResidueType const & restype( pose->residue_type(1) );
		TR << "Residue name: " << restype.name() << std::endl;
		TS_ASSERT_EQUALS( restype.base_name(), "ALA" );

		TS_ASSERT( restype.has_variant_type( core::chemical::ACETYLATED_NTERMINUS_VARIANT ) );
		TS_ASSERT( restype.has_variant_type( core::chemical::METHYLATED_CTERMINUS_VARIANT ) );
		TS_ASSERT( !restype.has_variant_type( core::chemical::DIMETHYLATED_CTERMINUS_VARIANT ) );

		//pose->dump_pdb( "V_GENMP_TEMP1.pdb" ); //DELETE ME
	}

	/// @brief Test the initial creation of the one-residue pose at a pre-proline position.
	void test_pose_generation_prepro() {
		using namespace protocols::mainchain_potential;

		GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) );
		options->set_residue_name( "ALA" );
		options->set_make_pre_proline_potential( true );
		GenerateMainchainPotential generator(options);
		core::pose::PoseOP pose( generator.generate_pose() );

		TS_ASSERT_EQUALS( pose->total_residue(), 1 );
		core::chemical::ResidueType const & restype( pose->residue_type(1) );
		TR << "Residue name: " << restype.name() << std::endl;
		TS_ASSERT_EQUALS( restype.base_name(), "ALA" );

		TS_ASSERT( restype.has_variant_type( core::chemical::ACETYLATED_NTERMINUS_VARIANT ) );
		TS_ASSERT( !restype.has_variant_type( core::chemical::METHYLATED_CTERMINUS_VARIANT ) );
		TS_ASSERT( restype.has_variant_type( core::chemical::DIMETHYLATED_CTERMINUS_VARIANT ) );

		//pose->dump_pdb( "V_GENMP_TEMP2.pdb" ); //DELETE ME
	}

	/// @brief Actually run the protocol and make sure that it does what we think it should.
	/// @details This version tests glycine, and makes an asymmetric Rama map.
	void test_full_run_GLY_asymm() {
		using namespace protocols::mainchain_potential;

		GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) );
		options->set_residue_name( "GLY" );
		options->set_make_pre_proline_potential(false);
		options->set_dimensions( utility::vector1<int>{ 18, 9 } ); //Low-res potential.
		options->set_output_filename("generated_gly_mainchain_potential_asymmetric.txt");

		GenerateMainchainPotential generator(options);
		TS_ASSERT_THROWS_NOTHING( generator.run() );
		TS_ASSERT_THROWS_NOTHING( generator.write_last_generated_to_disk() );
		//TODO: check size of mainchain potential, etc.
	}

	/// @brief Actually run the protocol and make sure that it does what we think it should.
	/// @details This version tests glycine, and makes a symmetric Rama map.
	void test_full_run_GLY_symm() {
		using namespace protocols::mainchain_potential;

		GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) );
		options->set_residue_name( "GLY" );
		options->set_make_pre_proline_potential(false);
		options->set_dimensions( utility::vector1<int>{ 18, 9 } ); //Low-res potential.
		options->set_symmetrize_output(true);
		options->set_output_filename("generated_gly_mainchain_potential_symmetric.txt");

		GenerateMainchainPotential generator(options);
		TS_ASSERT_THROWS_NOTHING( generator.run() );
		TS_ASSERT_THROWS_NOTHING( generator.write_last_generated_to_disk() );
		//TODO: check size of mainchain potential, etc.

		//Check that the generated potential is symmetric:
		for ( core::Size i(0); i<18; ++i ) {
			for ( core::Size j(0); j<9; ++j ) {
				utility::vector1 < core::Size > const coords1({ i, j });
				utility::vector1 < core::Size > const coords2({ 17-i, 8-j });
				TS_ASSERT_DELTA( generator.last_generated_scoretable_->energy_tensor(coords1), generator.last_generated_scoretable_->energy_tensor(coords2), 1e-7 );
			}
		}
	}

	/// @brief Actually run the protocol and make sure that it does what we think it should.
	/// @details This version tests alanine.
	void test_full_run_ALA() {
		using namespace protocols::mainchain_potential;

		GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) );
		options->set_residue_name( "ALA" );
		options->set_make_pre_proline_potential(false);
		options->set_dimensions( utility::vector1<int>{ 18, 9 } ); //Low-res potential.
		options->set_output_filename("generated_ala_mainchain_potential.txt");

		GenerateMainchainPotential generator(options);
		TS_ASSERT_THROWS_NOTHING( generator.run() );
		TS_ASSERT_THROWS_NOTHING( generator.write_last_generated_to_disk() );
		//TODO: check size of mainchain potential, etc.
	}

	/// @brief Actually run the protocol and make sure that it does what we think it should.
	/// @details This version tests valine.
	void test_full_run_VAL() {
		using namespace protocols::mainchain_potential;

		GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) );
		options->set_residue_name( "VAL" );
		options->set_make_pre_proline_potential(false);
		options->set_dimensions( utility::vector1<int>{ 18, 9 } ); //Low-res potential.
		options->set_output_filename("generated_val_mainchain_potential.txt");

		GenerateMainchainPotential generator(options);
		TS_ASSERT_THROWS_NOTHING( generator.run() );
		TS_ASSERT_THROWS_NOTHING( generator.write_last_generated_to_disk() );
		//TODO: check size of mainchain potential, etc.
	}

	/// @brief Actually run the protocol and make sure that it does what we think it should.
	/// @details This version tests valine with a pre-proline potential.
	void test_full_run_VAL_prepro() {
		using namespace protocols::mainchain_potential;

		GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) );
		options->set_residue_name( "VAL" );
		options->set_make_pre_proline_potential(true);
		options->set_dimensions( utility::vector1<int>{ 18, 9 } ); //Low-res potential.
		options->set_output_filename("generated_val_mainchain_potential_prepro.txt");

		GenerateMainchainPotential generator(options);
		TS_ASSERT_THROWS_NOTHING( generator.run() );
		TS_ASSERT_THROWS_NOTHING( generator.write_last_generated_to_disk() );
		//TODO: check size of mainchain potential, etc.
	}

	/// @brief Run the protocol on a beta-amino acid to confirm that it works for higher-dimensional potentials.
	/// @details Note that the beta-3-Asn rotamer library is already packaged with Rosetta for unit testing.
	void test_full_run_B3N() {
		using namespace protocols::mainchain_potential;

		GenerateMainchainPotentialOptionsOP options( utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) );
		options->set_residue_name( "B3N" );
		options->set_make_pre_proline_potential(false);
		options->set_dimensions( utility::vector1<int>{ 4,2,4 } ); //Very low-res potential.
		options->set_output_filename("generated_B3A_mainchain_potential.txt");

		GenerateMainchainPotential generator(options);
		TS_ASSERT_THROWS_NOTHING( generator.run() );
		TS_ASSERT_THROWS_NOTHING( generator.write_last_generated_to_disk() );
		//TODO: check size of mainchain potential, etc.
	}




};
