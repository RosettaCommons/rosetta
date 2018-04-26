// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/features/InterfaceDdGMoverTests.cxxtest.hh
/// @brief  Test the InterfaceDdGMover class
/// @author Kyle Barlow (kb@kylebarlow.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb3fod.hh>
#include <test/util/pdb3fod_mutated.hh>
#include <test/protocols/init_util.hh>

#include <protocols/features/InterfaceDdGMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("InterfaceDdGMoverTests");

class InterfaceDdGMoverTests : public CxxTest::TestSuite {

private:
	basic::datacache::DataMap data_;
	Filters_map filters_;
	Movers_map movers_;
	protocols::features::InterfaceDdGMoverOP test_instantiation_;
	core::pose::PoseOP test_8mer_pose_;
	core::pose::PoseOP test_dimer_pose_;

public:

	void setUp(){
		using namespace protocols::features;
		core_init();
		prime_Movers( movers_ ); // Adds "null" to movers map
		prime_Filters( filters_ ); // Adds true_filter, false_filter

		// Create dummy "commandline" score function
		core::scoring::ScoreFunctionOP dummy_commandline_sfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
		data_.add( "scorefxns", "commandline", dummy_commandline_sfxn );

		test_instantiation_ = InterfaceDdGMoverOP( new InterfaceDdGMover() );
		test_8mer_pose_ = pdb3fod_poseop();
		test_dimer_pose_ = create_2res_1ten_2res_trp_cage_poseop();
	}

	void tearDown(){
		// Reset RosettaScripts global data
		data_ = basic::datacache::DataMap();
		filters_ = Filters_map();
		movers_ = Movers_map();
	}

	bool is_dimer_pose_bound( core::pose::PoseCOP pose ) const {
		return pose->residue( 1 ).xyz( 1 ).distance( pose->residue( 3 ).xyz( 1 ) ) < 55;
	}

	bool is_dimer_pose_unbound( core::pose::PoseCOP pose ) const {
		return pose->residue( 1 ).xyz( 1 ).distance( pose->residue( 3 ).xyz( 1 ) ) > 900;
	}

	void assert_member_dimer_is_bound() const {
		TS_ASSERT( is_dimer_pose_bound( test_dimer_pose_ ) );
	}

	void assert_member_dimer_is_unbound() const {
		TS_ASSERT( is_dimer_pose_unbound( test_dimer_pose_ ) );
	}

	bool is_8mer_pose_bound( core::pose::PoseCOP pose, core::Size chain1, core::Size chain2 ) const {
		core::Size res1 = pose->conformation().chain_begin(chain1);
		core::Size res2 = pose->conformation().chain_begin(chain2);
		return pose->residue( res1 ).xyz( 1 ).distance( pose->residue( res2 ).xyz( 1 ) ) < 55;
	}

	bool is_8mer_pose_unbound( core::pose::PoseCOP pose, core::Size chain1, core::Size chain2 ) const {
		core::Size res1 = pose->conformation().chain_begin(chain1);
		core::Size res2 = pose->conformation().chain_begin(chain2);
		return pose->residue( res1 ).xyz( 1 ).distance( pose->residue( res2 ).xyz( 1 ) ) > 900;
	}

	void assert_member_8mer_is_bound( core::Size chain1, core::Size chain2 ) const {
		TS_ASSERT( is_8mer_pose_bound( test_8mer_pose_, chain1, chain2 ) );
	}

	void assert_member_8mer_is_unbound( core::Size chain1, core::Size chain2 ) const {
		TS_ASSERT( is_8mer_pose_unbound( test_8mer_pose_, chain1, chain2 ) );
	}

	void test_uninitialized() {
		// Tests that exception is thrown if which chain to move is not defined; via c++ interface
		TS_ASSERT_THROWS( test_instantiation_->unbind( *test_dimer_pose_ ), utility::excn::BadInput & );
	}

	void test_dimer_unbind_jump1() {
		test_dimer_pose_->dump_pdb("dimer_bound.pdb");
		assert_member_dimer_is_bound();
		test_instantiation_->add_jump_id( 1, *test_dimer_pose_ );
		test_instantiation_->unbind( *test_dimer_pose_ );
		assert_member_dimer_is_unbound();
		test_dimer_pose_->dump_pdb("dimer_unbound.pdb");
	}

	void test_dimer_unbind_jump2() {
		// Tests that failure is caught with out of range jump_id
		TS_ASSERT_THROWS( test_instantiation_->add_jump_id( 2, *test_dimer_pose_ ), utility::excn::BadInput & );
	}

	void test_dimer_unbind_chain1() {
		assert_member_dimer_is_bound();
		test_instantiation_->add_chain_id( 1, *test_dimer_pose_ );
		test_instantiation_->unbind( *test_dimer_pose_ );
		assert_member_dimer_is_unbound();
	}

	void test_dimer_unbind_chain2() {
		assert_member_dimer_is_bound();
		test_instantiation_->add_chain_id( 2, *test_dimer_pose_ );
		test_instantiation_->unbind( *test_dimer_pose_ );
		assert_member_dimer_is_unbound();
	}

	void test_dimer_unbind_chain3() {
		// Tests that failure is caught with out of range chain_id
		TS_ASSERT_THROWS( test_instantiation_->add_chain_id( 3, *test_dimer_pose_ ), utility::excn::BadInput & );
	}

	void test_dimer_unbind_chainA() {
		assert_member_dimer_is_bound();
		test_instantiation_->add_chain_name( "A", *test_dimer_pose_ );
		test_instantiation_->unbind( *test_dimer_pose_ );
		assert_member_dimer_is_unbound();
	}

	void test_dimer_unbind_chainB() {
		assert_member_dimer_is_bound();
		test_instantiation_->add_chain_name( "B", *test_dimer_pose_ );
		test_instantiation_->unbind( *test_dimer_pose_ );
		assert_member_dimer_is_unbound();
	}

	void test_dimer_unbind_chainC() {
		// Tests that failure is caught with out of range chain name
		TS_ASSERT_THROWS( test_instantiation_->add_chain_name( "C", *test_dimer_pose_ ), utility::excn::BadInput & );
	}

	void test_8mer_unbind_by_single_chain_id() {
		std::vector<core::Size> all_chains = {1, 2, 3, 4, 5, 6, 7, 8};
		std::vector<core::Size> chains_to_move = {1, 2, 5, 8}; // test a subset of chains for faster run

		for ( auto & outer_chain : chains_to_move ) {
			test_8mer_pose_->dump_pdb("8mer_bound_" + utility::to_string( outer_chain ) + ".pdb");

			// Test that all chains start off bound to all other chains
			for ( auto & all_chain : all_chains ) {
				for ( auto & inner_chain : all_chains ) {
					if ( all_chain != inner_chain ) {
						assert_member_8mer_is_bound( inner_chain, all_chain );
					}
				}
			}

			// Unbind chain
			test_instantiation_->add_chain_id( outer_chain, *test_8mer_pose_ );
			test_instantiation_->unbind( *test_8mer_pose_ );

			// Test that the chain we meant to unbind is unbound compared to all other chains
			for ( auto & inner_chain : all_chains ) {
				if ( inner_chain != outer_chain ) {
					assert_member_8mer_is_unbound(inner_chain, outer_chain);
				}
			}

			// Test that the chains we didn't want to move are all still next to each other
			for ( auto & all_chain : all_chains ) {
				for ( auto & inner_chain : all_chains ) {
					if ( all_chain != inner_chain && all_chain != outer_chain && inner_chain != outer_chain ) {
						assert_member_8mer_is_bound( inner_chain, all_chain );
					}
				}
			}

			test_8mer_pose_->dump_pdb("8mer_unbound_" + utility::to_string( outer_chain ) + ".pdb");

			// Run only necessary bits of setUp each time through the for loop
			test_instantiation_ = protocols::features::InterfaceDdGMoverOP( new protocols::features::InterfaceDdGMover() );
			test_8mer_pose_ = pdb3fod_poseop();
		}
	}

	void test_parse_my_tag() {
		// Test that emptiest tag possible works and defaults to move jump 1
		TagCOP empty_tag = tagptr_from_string( "<InterfaceDdGMover name=test_interface_ddg_filter/>\n" );
		test_instantiation_->parse_my_tag( empty_tag, data_, filters_, movers_, *test_8mer_pose_ );
		TS_ASSERT( test_instantiation_->get_chain_ids().size() == 1 );
		TS_ASSERT( test_instantiation_->get_chain_ids()[1] == 1 );

		// Can specify either chain_name, jump, or chain_num, but not multiples of those
		std::vector<std::string> good_tags = {
			"<InterfaceDdGMover name=test_interface_ddg_filter chain_name=\"A\" />\n",
			"<InterfaceDdGMover name=test_interface_ddg_filter chain_name=\"B\" />\n",
			"<InterfaceDdGMover name=test_interface_ddg_filter chain_num=1 />\n",
			"<InterfaceDdGMover name=test_interface_ddg_filter chain_num=2 />\n",
			"<InterfaceDdGMover name=test_interface_ddg_filter jump=1 />\n"
			};
		for ( auto tag_str : good_tags ) {
			TagCOP tag = tagptr_from_string(tag_str);
			test_instantiation_->parse_my_tag(tag, data_, filters_, movers_, *test_dimer_pose_);
		}

		// Test failure for conflicting tag options
		std::vector<std::string> bad_tags = {
			"<InterfaceDdGMover name=test_interface_ddg_filter chain_name=\"A\" chain_num=1 jump=1/>\n",
			"<InterfaceDdGMover name=test_interface_ddg_filter chain_name=\"A\" chain_num=1/>\n",
			"<InterfaceDdGMover name=test_interface_ddg_filter chain_name=\"A\" jump=1/>\n",
			"<InterfaceDdGMover name=test_interface_ddg_filter wt_ref_savepose_mover=\"whatever\" mut_ref_savepose_mover=\"nihilism\"/>\n",
			"<InterfaceDdGMover name=test_interface_ddg_filter chain_num=1 jump=1/>\n"
			};
		for ( auto tag_str : bad_tags ) {
			TagCOP tag = tagptr_from_string(tag_str);
			TS_ASSERT_THROWS( test_instantiation_->parse_my_tag(tag, data_, filters_, movers_, *test_dimer_pose_), utility::excn::RosettaScriptsOptionError & );
		}
	}

	void test_wt_pose_setter_and_getter() {
		// Tests that the setter sets and unbinds the given pose
		test_instantiation_->add_jump_id( 1, *test_dimer_pose_ );
		assert_member_dimer_is_bound();
		test_instantiation_->set_wt_pose( *test_dimer_pose_ );
		assert_member_dimer_is_bound(); // member dimer pose should be cloned and remain unchanged
		TS_ASSERT( is_dimer_pose_bound( test_instantiation_->get_bound_wt_pose() ) );
		TS_ASSERT( is_dimer_pose_unbound( test_instantiation_->get_unbound_wt_pose() ) );
	}

	void test_mut_pose_setter_and_getter() {
		// Tests that the setter sets and unbinds the given pose
		test_instantiation_->add_jump_id( 1, *test_dimer_pose_ );
		assert_member_dimer_is_bound();
		test_instantiation_->set_mut_pose( *test_dimer_pose_ );
		assert_member_dimer_is_bound(); // member dimer pose should be cloned and remain unchanged
		TS_ASSERT( is_dimer_pose_bound( test_instantiation_->get_bound_mut_pose() ) );
		TS_ASSERT( is_dimer_pose_unbound( test_instantiation_->get_unbound_mut_pose() ) );
	}

	void test_mutation_list() {
		core::pose::PoseOP test_8mer_mut_pose = pdb3fod_mutated_poseop();
		test_instantiation_->add_chain_id( 2, *test_8mer_pose_ );
		test_instantiation_->set_wt_pose( *test_8mer_pose_ );
		test_instantiation_->set_mut_pose( *test_8mer_mut_pose );
		utility::vector1< std::tuple<std::string, core::Size, std::string> > mutation_list = test_instantiation_->mutation_list();

		TS_ASSERT( std::get<0>(mutation_list[1]) == "SER" );
		TS_ASSERT( std::get<1>(mutation_list[1]) == 16 );
		TS_ASSERT( std::get<2>(mutation_list[1]) == "ARG" );

		TS_ASSERT( std::get<0>(mutation_list[2]) == "ILE" );
		TS_ASSERT( std::get<1>(mutation_list[2]) == 26 );
		TS_ASSERT( std::get<2>(mutation_list[2]) == "ALA" );

	}

};
