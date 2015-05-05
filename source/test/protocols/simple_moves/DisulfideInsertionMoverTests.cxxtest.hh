// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/simple_moves/DisulfideInsertionMover.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

class DisulfideInsertionMoverTests : public CxxTest::TestSuite {
public:
	// TODO : setup a pose
	DisulfideInsertionMoverTests() {}

	void setUp() {
		protocols_init();
	}

	void tearDown() {
	}

	void test_rosettascripts_options() {

		std::ifstream xml_stream( "protocols/simple_moves/DisulfideInsertionMover.xml" );

		utility::tag::TagCOP tag =  utility::tag::Tag::create( xml_stream );
		try {
			protocols::rosetta_scripts::RosettaScriptsParser parser;
			// NOTE : the following assumes that RosettaScriptParser::parse_protocol_tag() returns a pointer ParsedProtocol instance, even though the interface is to a MoverOP.
			protocols::rosetta_scripts::ParsedProtocolOP protocol(utility::pointer::dynamic_pointer_cast<protocols::rosetta_scripts::ParsedProtocol> (parser.parse_protocol_tag(tag)) );
			TS_ASSERT_EQUALS(protocol->size(), 1);
			protocols::rosetta_scripts::ParsedProtocol::MoverFilterPair pair = protocol->get_mover_filter_pair(1);
			TS_ASSERT_EQUALS(protocol->get_mover(1)->get_name(), "DisulfideInsertionMover");
			protocols::simple_moves::DisulfideInsertionMoverCOP const mover( utility::pointer::dynamic_pointer_cast<protocols::simple_moves::DisulfideInsertionMover const> (protocol->get_mover(1)) );

			TS_ASSERT_EQUALS(mover->get_peptide_chain(), 2);
			TS_ASSERT_EQUALS(mover->get_n_cyd_seqpos(), 2);
			TS_ASSERT_EQUALS(mover->get_c_cyd_seqpos(), 6);
			TS_ASSERT_EQUALS(mover->get_is_cyd_res_at_termini(), false);
			TS_ASSERT_EQUALS(mover->get_scorefxn()->get_name(), "soft_rep");
			TS_ASSERT_DELTA(mover->get_constraint_weight(), 0.1, 1e-6);

		}
		catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

}; // DisulfideInsertionMoverTests
