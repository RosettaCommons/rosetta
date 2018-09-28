// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief test suite for protocols::simple_moves::DisulfideInsertionMover
/// @author Yuval Sadan (yuval.sedan@mail.huji.ac.il)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/DisulfideInsertionMover.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// Basic headers
#include <basic/options/option.hh>

class DisulfideInsertionMoverTests : public CxxTest::TestSuite {
public:

	core::pose::PoseOP test_edge_pose_;
	core::pose::PoseOP test_shorter_pose_;
	protocols::simple_moves::DisulfideInsertionMoverOP inserter_edges_;

	DisulfideInsertionMoverTests() {}


	void setUp() {
		protocols_init_with_additional_options("-max_dslf_energy 0");
		test_edge_pose_ = core::pose::PoseOP( new core::pose::Pose() );
		core::import_pose::pose_from_file( *test_edge_pose_, "protocols/simple_moves/1A2A_A.pdb" , core::import_pose::PDB_file);


		core::Size const PEP_CHAIN = 1;
		core::Size const N_TER = test_edge_pose_->conformation().chain_begin(PEP_CHAIN);
		core::Size const C_TER = test_edge_pose_->conformation().chain_end(PEP_CHAIN);
		core::Size const N_TER_SHORT_PEP = N_TER+1;

		test_shorter_pose_ = core::pose::PoseOP( new core::pose::Pose(*test_edge_pose_, N_TER_SHORT_PEP, C_TER) );

		inserter_edges_ = protocols::simple_moves::DisulfideInsertionMoverOP( new protocols::simple_moves::DisulfideInsertionMover() );
		inserter_edges_->set_peptide_chain(PEP_CHAIN);

	}

	void test_cyclization_viability_function() {
		try {
			core::Size const PEP_CHAIN = 1;
			core::Size const N_TER = test_edge_pose_->conformation().chain_begin(PEP_CHAIN);
			core::Size const C_TER = test_edge_pose_->conformation().chain_end(PEP_CHAIN);

			// check that the function reports that the two residues could not be mutated to cys to form disulfide
			// since the C-ter is already part of another one
			protocols::simple_moves::DisulfideCyclizationViability term_cyclization_viable (inserter_edges_->determine_cyclization_viability(*test_edge_pose_, N_TER, C_TER));
			TS_ASSERT_EQUALS(term_cyclization_viable, protocols::simple_moves::DCV_NOT_CYCLIZABLE);

			// check that the existing disulfide is identified as one
			protocols::simple_moves::DisulfideCyclizationViability existing_disulfide_maintained (inserter_edges_->determine_cyclization_viability(*test_shorter_pose_, N_TER, C_TER-1));
			TS_ASSERT_EQUALS(existing_disulfide_maintained, protocols::simple_moves::DCV_ALREADY_CYCLIZED);
		}
catch (utility::excn::Exception & e ) {
	std::cerr << "Raised exception: " << e.msg() << std::endl;
	TS_ASSERT( false );
}
	}

	void test_apply_function () {
		try {
			// check that apply function returns and the status is fail retry for DCV_NOT_CYCLIZABLE
			inserter_edges_->apply(*test_edge_pose_);
			TS_ASSERT_EQUALS(inserter_edges_->get_last_move_status(), protocols::moves::FAIL_RETRY);

			// check the cyclic existing disulfide is successfully modeled within the mover
			inserter_edges_->apply(*test_shorter_pose_);
			TS_ASSERT_EQUALS(inserter_edges_->get_last_move_status(), protocols::moves::MS_SUCCESS);
		}
catch (utility::excn::Exception & e ) {
	std::cerr << "Raised exception: " << e.msg() << std::endl;
	TS_ASSERT( false );
}
	}


	void tearDown() {
		test_edge_pose_ = NULL;
		test_shorter_pose_ = NULL;
	}

	void test_rosettascripts_options() {

		std::ifstream xml_stream( "protocols/simple_moves/DisulfideInsertionMover.xml" );

		utility::tag::TagCOP tag =  utility::tag::Tag::create( xml_stream );
		try {
			protocols::rosetta_scripts::RosettaScriptsParser parser;
			// NOTE : the following assumes that RosettaScriptParser::parse_protocol_tag() returns a pointer ParsedProtocol instance, even though the interface is to a MoverOP.
			protocols::rosetta_scripts::ParsedProtocolOP protocol(utility::pointer::dynamic_pointer_cast<protocols::rosetta_scripts::ParsedProtocol> (parser.parse_protocol_tag(tag,basic::options::option)) );
			TS_ASSERT_EQUALS(protocol->size(), 1);
			protocols::rosetta_scripts::ParsedProtocol::MoverFilterPair pair = protocol->get_mover_filter_pair(1);
			TS_ASSERT_EQUALS(protocol->get_mover(1)->get_name(), "DisulfideInsertion");
			protocols::simple_moves::DisulfideInsertionMoverCOP const mover( utility::pointer::dynamic_pointer_cast<protocols::simple_moves::DisulfideInsertionMover const> (protocol->get_mover(1)) );

			TS_ASSERT_EQUALS(mover->get_peptide_chain(), 2);
			TS_ASSERT_EQUALS(mover->get_n_cyd_seqpos(), 2);
			TS_ASSERT_EQUALS(mover->get_c_cyd_seqpos(), 6);
			TS_ASSERT_EQUALS(mover->get_scorefxn()->get_name(), "ref2015_soft");
			TS_ASSERT_DELTA(mover->get_constraint_weight(), 0.1, 1e-6);
		}

catch (utility::excn::Exception & e ) {
	std::cerr << "Raised exception: " << e.msg() << std::endl;
	TS_ASSERT( false );
}

	}

}; // DisulfideInsertionMoverTests
