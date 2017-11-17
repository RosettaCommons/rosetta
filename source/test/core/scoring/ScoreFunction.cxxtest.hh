// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction.cxxtest.hh
/// @brief  unit testing for core/scoring/ScoreFunction.*
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/UTracer.hh>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// Package headers
#include <utility/vector1.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKey.hh>
#include <basic/options/util.hh>

// option key includes
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.ScoreFunction.cxxtest");

// using declarations
using namespace core;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// ScoreFunctionUtilityTest
/// moved to ScoreFunctionFactory.cxxtest.hh
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
/// @name ScoreFunctionTest
/// @brief: unified tests for ScoreFunction object
///////////////////////////////////////////////////////////////////////////
class ScoreFunctionTest : public CxxTest::TestSuite {

public:
	typedef utility::keys::VariantKey< utility::options::OptionKey > VariantOptionKey;

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_add_weights_from_file() {

		// Test local loading
		ScoreFunction sfxn;
		sfxn.add_weights_from_file("core/scoring/test");
		TS_ASSERT_DELTA( sfxn.get_weight(fa_rep) , 3.1416, 0.0001 );

		// Test database loading
		ScoreFunction sfxn2;
		sfxn2.add_weights_from_file("pre_talaris_2013_standard");
		TS_ASSERT_DELTA( sfxn2.get_weight(fa_atr) , 0.8, 0.0001 );
	}

	void test_apply_patch_from_file() {

		// Test local loading
		ScoreFunction sfxn;
		sfxn.add_weights_from_file("core/scoring/test");
		sfxn.apply_patch_from_file("core/scoring/test");
		TS_ASSERT_DELTA( sfxn.get_weight(fa_rep) , 3.1416, 0.0001 );
		TS_ASSERT_DELTA( sfxn.get_weight(fa_atr) , 2.7183, 0.0001 );

		// Test database loading
		ScoreFunction sfxn2;

		sfxn2.add_weights_from_file("pre_talaris_2013_standard");
		sfxn2.apply_patch_from_file("score12");
		TS_ASSERT_DELTA( sfxn2.get_weight(fa_atr) , 0.8, 0.0001 );
		TS_ASSERT_DELTA( sfxn2.get_weight(p_aa_pp) , 0.32, 0.0001 );
	}


	void test_get_sub_score_exclude_res() {

		ScoreFunction scfxn;
		core::pose::Pose pose = create_trpcage_ideal_pose();
		Real sc(scfxn(pose));

		utility::vector1< Size > exclude_list;
		Real sc_exc1(scfxn.get_sub_score_exclude_res(pose, exclude_list));
		TS_ASSERT_DELTA(sc_exc1, sc, .0000001);

		utility::vector1< bool > residue_mask(pose.size(), true);
		Real sc_exc2(scfxn.get_sub_score(pose, residue_mask));
		TS_ASSERT_DELTA(sc_exc2, sc, .0000001);
	}

	/// @detail the long range terms maintain their own neighbor maps, so
	///require iterating over the energy methods then over the residue
	///pairs rather then the other way around.
	void test_get_sub_score_exclude_res_long_range_terms() {

		ScoreFunction scfxn;

		//CartesianBondedEnergy is a ContextIndependentLRTwoBodyEnergy
		scfxn.set_weight(core::scoring::cart_bonded, 1);

		core::pose::Pose pose = create_trpcage_ideal_pose();
		Real sc(scfxn(pose));

		utility::vector1< Size > exclude_list;
		Real sc_exc1(scfxn.get_sub_score_exclude_res(pose, exclude_list));
		TS_ASSERT_DELTA(sc_exc1, sc, .0000001);

		utility::vector1< bool > residue_mask(pose.size(), true);
		Real sc_exc2(scfxn.get_sub_score(pose, residue_mask));
		TS_ASSERT_DELTA(sc_exc2, sc, .0000001);
	}

	void test_ScoreFunction_list_options_read_in_sync() {
		using namespace utility::keys;
		TS_ASSERT( true );

		utility::options::OptionKeyList sfxn_ctor_options;
		ScoreFunction::list_options_read( sfxn_ctor_options );

		TS_ASSERT_EQUALS( std::count( sfxn_ctor_options.begin(), sfxn_ctor_options.end(), VariantOptionKey( basic::options::OptionKeys::score::elec_die )), 1 );
		TS_ASSERT_EQUALS( std::count( sfxn_ctor_options.begin(), sfxn_ctor_options.end(), VariantOptionKey( basic::options::OptionKeys::score::grp_cpfxn )), 1 );
		TS_ASSERT_EQUALS( std::count( sfxn_ctor_options.begin(), sfxn_ctor_options.end(), VariantOptionKey( basic::options::OptionKeys::corrections::score::lj_hbond_hdis )), 1 );

		// cursory check that not all options are in here, because that would be weird
		TS_ASSERT_EQUALS( std::count( sfxn_ctor_options.begin(), sfxn_ctor_options.end(), VariantOptionKey( basic::options::OptionKeys::in::file::s )), 0 );

		utility::options::OptionCollectionCOP sfxn_option_collection =
			basic::options::read_subset_of_global_option_collection( sfxn_ctor_options );

		// Now, try to call get_scoref_fxn using the local option collection
		try {
			set_throw_on_next_assertion_failure(); // just in case
			ScoreFunction sfxn( *sfxn_option_collection );
		} catch (utility::excn::Exception const & e ) {
			std::cerr << e.msg() << std::endl;
			TS_ASSERT( false ); // we screwed the pooch
		}

	}

	void test_ScoreFunction_ctor_actually_reads_option_collection()
	{

		using namespace utility::keys;
		TS_ASSERT( true );

		utility::options::OptionKeyList sfxn_ctor_options;
		ScoreFunction::list_options_read( sfxn_ctor_options );

		// Now drop one of the options from the list, one which always gets read, and make sure that
		// when we call the initialize_from_options function, that an assertion failure occurs
		sfxn_ctor_options.remove( VariantOptionKey( basic::options::OptionKeys::corrections::score::lj_hbond_hdis ));

		utility::options::OptionCollectionCOP sfxn_option_collection =
			basic::options::read_subset_of_global_option_collection( sfxn_ctor_options );

		// Now, try to create a ScoreFunction from the local option collection
		try {
			set_throw_on_next_assertion_failure();
			ScoreFunction sfxn( *sfxn_option_collection );
			TS_ASSERT( false ); // we screwed the pooch
		} catch ( ... ) {
			// good -- if we don't list an option that we're going to read, then
			// an exception will be thrown / an assertion failure will get triggered
			TS_ASSERT( true );
		}

	}

};
