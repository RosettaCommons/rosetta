// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/relax/FastRelax.cxxtest.hh
/// @brief  test for FastRelax
/// @author Jack Maguire, jackmaguire1444@gmail.com

// I created this file but am testing only my contribution to the mover. There is still a lot to be tested. -Jack

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit header
#include <protocols/relax/FastRelax.hh>
#include <protocols/denovo_design/movers/FastDesign.hh>

// project headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/variant_util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/parser/ScoreFunctionLoader.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

static basic::Tracer TR( "FastRelaxTests" );

// --------------- Test Class --------------- //

class FastRelaxTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	// ------------- Helper Functions ------------- //

	void check_fa_rep_weights_in_script(
		std::vector< protocols::relax::RelaxScriptCommand > const & script,
		std::vector< core::Real > const & expected_weights
	) {
		core::Size counter = 0;
		for ( protocols::relax::RelaxScriptCommand const & cmd : script ) {
			if ( cmd.command == "ramp_repack_min" ) {
				core::Real const fa_rep_weight = cmd.param1;
				TS_ASSERT( counter <= expected_weights.size() );
				TS_ASSERT_LESS_THAN( std::abs( fa_rep_weight - expected_weights[ counter ] ), 0.01 );
				++counter;
			}
		}
		TS_ASSERT_EQUALS( counter, expected_weights.size() );
	}


	// --------------- Test Cases --------------- //

	void test_legacy_script() {
		TR << "test_legacy_script" << std::endl;

		using namespace basic::options;
		option[ OptionKeys::corrections::beta_nov16 ]( false );

		impl_test_legacy_script< protocols::relax::FastRelax >();
		impl_test_legacy_script< protocols::denovo_design::movers::FastDesign >();
	}

	template < typename T >
	void impl_test_legacy_script(){
		std::stringstream ss;
		ss << "<" << T::mover_name() << " name=\"name\" scorefxn=\"dummy\" relaxscript=\"legacy\"/>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::Pose pose;

		std::stringstream ss2;
		ss2 << "<SCOREFXNS>\n<ScoreFunction name=\"dummy\"/>\n</SCOREFXNS>" << std::endl; //use legacy weights
		utility::tag::TagOP sfxn_tag( new utility::tag::Tag() );
		sfxn_tag->read( ss2 );
		protocols::parser::ScoreFunctionLoader loader;
		TS_ASSERT_THROWS_NOTHING( loader.load_data( pose, sfxn_tag, dm ) );

		//protocols::relax::FastRelax mover;
		T mover;
		mover.parse_my_tag( tag, dm, fm, mm, pose );

		std::vector< core::Real > expected_weights = { 0.02, 0.25, 0.55, 1.0 };
		check_fa_rep_weights_in_script( mover.script_, expected_weights );
	}

	void test_legacy_dualspace_script() {
		TR << "test_dualspace_script" << std::endl;

		using namespace basic::options;
		option[ OptionKeys::corrections::beta_nov16 ]( false );

		impl_test_legacy_dualspace_script< protocols::relax::FastRelax >();
		impl_test_legacy_dualspace_script< protocols::denovo_design::movers::FastDesign >();
	}

	template < typename T >
	void impl_test_legacy_dualspace_script(){
		std::stringstream ss;
		ss << "<" << T::mover_name() << " name=\"name\" scorefxn=\"dummy\" dualspace=\"true\" relaxscript=\"legacy\"/>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::Pose pose;

		std::stringstream ss2;
		ss2 << "<SCOREFXNS>\n<ScoreFunction name=\"dummy\"/>\n</SCOREFXNS>" << std::endl; //use default weights
		utility::tag::TagOP sfxn_tag( new utility::tag::Tag() );
		sfxn_tag->read( ss2 );
		protocols::parser::ScoreFunctionLoader loader;
		TS_ASSERT_THROWS_NOTHING( loader.load_data( pose, sfxn_tag, dm ) );

		//protocols::relax::FastRelax mover;
		T mover;
		mover.parse_my_tag( tag, dm, fm, mm, pose );

		std::vector< core::Real > expected_weights = { 0.02, 0.25, 0.55, 1.0, 0.02, 0.25, 0.55, 1.0 };
		check_fa_rep_weights_in_script( mover.script_, expected_weights );
	}


	void test_beta_custom_script() {
		TR << "test_beta_custom_script" << std::endl;

		protocols_init_with_additional_options( "-beta_nov16" );

		impl_test_beta_custom_script< protocols::relax::FastRelax >();
		impl_test_beta_custom_script< protocols::denovo_design::movers::FastDesign >();
	}

	template < typename T >
	void impl_test_beta_custom_script(){
		std::stringstream ss;
		ss << "<" << T::mover_name() << " name=\"name\" scorefxn=\"dummy\" relaxscript=\"rosettacon2018\"/>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::Pose pose;

		std::stringstream ss2;
		ss2 << "<SCOREFXNS>\n<ScoreFunction name=\"dummy\"/>\n</SCOREFXNS>" << std::endl; //use default weights
		utility::tag::TagOP sfxn_tag( new utility::tag::Tag() );
		sfxn_tag->read( ss2 );
		protocols::parser::ScoreFunctionLoader loader;
		TS_ASSERT_THROWS_NOTHING( loader.load_data( pose, sfxn_tag, dm ) );

		//protocols::relax::FastRelax mover;
		T mover;
		mover.parse_my_tag( tag, dm, fm, mm, pose );

		std::vector< core::Real > expected_weights = { 0.079, 0.295, 0.577, 1 };
		check_fa_rep_weights_in_script( mover.script_, expected_weights );
	}





	void test_local_script() {
		TR << "test_local_script" << std::endl;
		impl_test_local_script< protocols::relax::FastRelax >();
		impl_test_local_script< protocols::denovo_design::movers::FastDesign >();
	}

	template < typename T >
	void impl_test_local_script(){
		std::stringstream ss;
		ss << "<" << T::mover_name() << " name=\"name\" scorefxn=\"dummy\" relaxscript=\"protocols/relax/fast_relax_test_script.dat\"/>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::Pose pose;

		std::stringstream ss2;
		ss2 << "<SCOREFXNS>\n<ScoreFunction name=\"dummy\"/>\n</SCOREFXNS>" << std::endl; //use default weights
		utility::tag::TagOP sfxn_tag( new utility::tag::Tag() );
		sfxn_tag->read( ss2 );
		protocols::parser::ScoreFunctionLoader loader;
		TS_ASSERT_THROWS_NOTHING( loader.load_data( pose, sfxn_tag, dm ) );

		//protocols::relax::FastRelax mover;
		T mover;
		mover.parse_my_tag( tag, dm, fm, mm, pose );

		std::vector< core::Real > expected_weights = { 0.1, 0.2, 0.5, 1 };
		check_fa_rep_weights_in_script( mover.script_, expected_weights );
	}



	void test_nonexistent_local_script() {
		TR << "test_nonexistent_local_script" << std::endl;
		impl_test_nonexistent_local_script< protocols::relax::FastRelax >();
		impl_test_nonexistent_local_script< protocols::denovo_design::movers::FastDesign >();
	}

	template < typename T >
	void impl_test_nonexistent_local_script(){
		std::stringstream ss;
		ss << "<" << T::mover_name() << " name=\"name\" scorefxn=\"dummy\" relaxscript=\"DOES.NOT.EXIST.nope\"/>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::Pose pose;

		std::stringstream ss2;
		ss2 << "<SCOREFXNS>\n<ScoreFunction name=\"dummy\"/>\n</SCOREFXNS>" << std::endl; //use default weights
		utility::tag::TagOP sfxn_tag( new utility::tag::Tag() );
		sfxn_tag->read( ss2 );
		protocols::parser::ScoreFunctionLoader loader;
		TS_ASSERT_THROWS_NOTHING( loader.load_data( pose, sfxn_tag, dm ) );

		//protocols::relax::FastRelax mover;
		T mover;
		set_throw_on_next_assertion_failure(); // Suppress backtrace printing.
		TS_ASSERT_THROWS_ANYTHING( mover.parse_my_tag( tag, dm, fm, mm, pose ) );
	}

	/// @brief Check that this works for a noncanonical
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	void test_fastrelax_B3N() {
		TR << "Testing with beta-3-asparagine..." << std::endl;
		impl_test_fastrelax_B3N< protocols::relax::FastRelax >();
		//Note: this isn't expected to work with FastDesign at the moment.
		//impl_test_fastrelax_B3N< protocols::denovo_design::movers::FastDesign >();
	}

	/// @brief Check that this works for a noncanonical
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	template < typename T >
	void impl_test_fastrelax_B3N() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "GX[B3N]G", "fa_standard");

		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
		T mover( sfxn, 3 );
		TS_ASSERT_THROWS_NOTHING( mover.apply(pose) ); //Just run and make sure we don't crash.
	}

	/// @brief Check that this works for a noncanonical with modified termini.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	void test_fastrelax_B3N_mod_termini() {
		TR << "Testing with beta-3-asparagine with acetylated N-terminus and N-methylated C-terminus..." << std::endl;
		impl_test_fastrelax_B3N_mod_termini< protocols::relax::FastRelax >();
		//Note: this isn't expected to work with FastDesign at the moment.
		//impl_test_fastrelax_B3N_mod_termini< protocols::denovo_design::movers::FastDesign >();
	}

	/// @brief Check that this works for a noncanonical
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	template < typename T >
	void impl_test_fastrelax_B3N_mod_termini() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "X[B3N]", "fa_standard");

		// Strip off terminal types:
		core::pose::remove_lower_terminus_type_from_pose_residue( pose, 1 );
		core::pose::remove_upper_terminus_type_from_pose_residue( pose, 1 );

		// Add new terminal types:
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::METHYLATED_CTERMINUS_VARIANT, 1 );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::ACETYLATED_NTERMINUS_VARIANT, 1 );

		// Check that they're on there:
		TS_ASSERT( pose.residue_type(1).has_variant_type( core::chemical::METHYLATED_CTERMINUS_VARIANT ) );
		TS_ASSERT( pose.residue_type(1).has_variant_type( core::chemical::ACETYLATED_NTERMINUS_VARIANT ) );

		//pose.dump_pdb("FRLX_B3N_MOD_TERMINI.pdb"); //DELETE ME

		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
		T mover( sfxn, 3 );
		TS_ASSERT_THROWS_NOTHING( mover.apply(pose) ); //Just run and make sure we don't crash.
	}

};//end class
