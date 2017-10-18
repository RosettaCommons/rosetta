// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/residue_selectors/LigandMetalContactSelector.cxxtest.hh
/// @brief test suite for protocols::residue_selectors::LigandMetalContactSelector
/// @author


// Test headers
#include <cxxtest/TestSuite.h>
#include <protocols/residue_selectors/LigandMetalContactSelector.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <test/protocols/init_util.hh>
#include <core/import_pose/import_pose.hh>




static THREAD_LOCAL basic::Tracer TR("protocols.residue_selectors.LigandMetalContactSelectorTests");

class LigandMetalContactSelectorTests : public CxxTest::TestSuite {

private:
	core::pose::Pose pose_;

public:

	void setUp(){
		protocols_init();
		const std::string original_file_name("core/util/2c9v_stripped.pdb");
		core::import_pose::pose_from_file(pose_, original_file_name, core::import_pose::PDB_file);
	}

	void tearDown(){
	}

	void test_parse_tag_invalid(){
		basic::datacache::DataMap datamap;
		protocols::residue_selectors::LigandMetalContactSelector selector;

		//Trying to give resnums and an embedded selector
		std::stringstream ss_string_and_selector;
		ss_string_and_selector << "<LigandMetalContactSelector name=\"metal\" resnums=\"154,155\"><ResidueName name=\"select_metal\" residue_name3=\" ZN\"/></LigandMetalContactSelector>"<< std::endl;
		utility::tag::TagOP tag(new utility::tag::Tag());
		tag->read(ss_string_and_selector);
		TR << "Tag with embedded selector and resnums" << std::endl;
		TS_ASSERT_THROWS_ANYTHING(selector.parse_my_tag(tag, datamap));

		//Trying to give two embedded selectors
		std::stringstream ss_two_embedded;
		ss_two_embedded << "<LigandMetalContactSelector name=\"metal\"><Index resnums=\"154,155\" /><Chain chains=\"A\" /></LigandMetalContactSelector>" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_two_embedded );
		TR << "Tag with two embedded selectors" << std::endl;
		TS_ASSERT_THROWS_ANYTHING( selector.parse_my_tag(tag, datamap));

		//Provide a selector that isn't in the data map
		std::stringstream ss_bad_selector;
		ss_bad_selector << "<LigandMetalContactSelector name=\"metal\" residue_selector=\"dummy\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_bad_selector );
		TR << "Tag with undefined residue selector" << std::endl;
		TS_ASSERT_THROWS_ANYTHING( selector.parse_my_tag( tag, datamap) );

		//Add selector to the data map
		core::select::residue_selector::ResidueIndexSelectorOP dummy( new core::select::residue_selector::ResidueIndexSelector );
		dummy->set_index( "154,155" );
		datamap.add( "ResidueSelector", "dummy", dummy );
		//Confirm that the above tag would work once the selector was added
		std::stringstream ss_good_selector;
		ss_good_selector << "<LigandMetalContactSelector name=\"metal\" residue_selector=\"dummy\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_good_selector );
		TR << "Tag with defined residue selector" << std::endl;
		TS_ASSERT_THROWS_NOTHING( selector.parse_my_tag( tag, datamap) );
		//Try to provide a selector string and an embedded selector
		std::stringstream ss_selector_two_ways;
		ss_selector_two_ways << "<LigandMetalContactSelector name=\"metal\" residue_selector=\"dummy\" ><Index resnums=\"154,155\" /></LigandMetalContactSelector>" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_selector_two_ways );
		TR << "Tag with selector two ways" << std::endl;
		TS_ASSERT_THROWS_ANYTHING( selector.parse_my_tag( tag, datamap ) );
		//Try to provide a selector string and a resnum string
		std::stringstream ss_select_string_and_resnum_string;
		ss_select_string_and_resnum_string << "<LigandMetalContactSelector name=\"metal\" resnums=\"154,155\" residue_selector=\"dummy\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_select_string_and_resnum_string );
		TR << "Tag with selector string and resnum string" << std::endl;
		TS_ASSERT_THROWS_ANYTHING( selector.parse_my_tag( tag, datamap) );
	}

	void test_parse_tag_embedded_selector(){
		//This should only set up metal ions specified by the selector
		//Select all metals
		protocols::residue_selectors::LigandMetalContactSelector selector;
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_all;
		ss_all << "<LigandMetalContactSelector name=\"metal\"><Index resnums=\"154,155,309,310\" /></LigandMetalContactSelector>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_all );
		TR << "Tag with embedded selector to select all metals" << std::endl;
		basic::datacache::DataMap datamap;
		TS_ASSERT_THROWS_NOTHING( selector.parse_my_tag( tag, datamap ) );
		//Apply the mover
		core::select::residue_selector::ResidueSubset res_sub = selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( res_sub[46] );
		TS_ASSERT( res_sub[48] );
		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[120] );

		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[71] );
		TS_ASSERT( res_sub[80] );
		TS_ASSERT( res_sub[83] );

		TS_ASSERT( res_sub[201] );
		TS_ASSERT( res_sub[203] );
		TS_ASSERT( res_sub[218] );
		TS_ASSERT( res_sub[275] );

		TS_ASSERT( res_sub[218] );
		TS_ASSERT( res_sub[226] );
		TS_ASSERT( res_sub[235] );
		TS_ASSERT( res_sub[238] );

		//Select metals only on chain A
		selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_a;
		ss_a << "<LigandMetalContactSelector name=\"metal\"  ><Chain chains=\"A\" /></LigandMetalContactSelector>" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_a );
		TR << "Tag with embedded selector to select chain A" << std::endl;
		TS_ASSERT_THROWS_NOTHING( selector.parse_my_tag( tag, datamap ) );

		//Apply the mover
		res_sub = selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( res_sub[46] );
		TS_ASSERT( res_sub[48] );
		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[120] );

		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[71] );
		TS_ASSERT( res_sub[80] );
		TS_ASSERT( res_sub[83] );

		TS_ASSERT( !res_sub[201] );
		TS_ASSERT( !res_sub[203] );
		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[275] );

		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[226] );
		TS_ASSERT( !res_sub[235] );
		TS_ASSERT( !res_sub[238] );

		//Select only zinc ions
		selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_zn;
		ss_zn << "<LigandMetalContactSelector name=\"metal\" ><ResidueName residue_name3=\" ZN\" /></LigandMetalContactSelector>" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_zn );
		TR << "Tag with embedded selector to select zinc ions" << std::endl;
		TS_ASSERT_THROWS_NOTHING( selector.parse_my_tag( tag, datamap) );

		//Apply the mover
		res_sub = selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT(!res_sub[46]);
		TS_ASSERT(!res_sub[48]);
		// TS_ASSERT(!res_sub[63]);
		TS_ASSERT(!res_sub[120]);

		TS_ASSERT(res_sub[63]);
		TS_ASSERT(res_sub[71]);
		TS_ASSERT(res_sub[80]);
		TS_ASSERT(res_sub[83]);

		TS_ASSERT(!res_sub[201]);
		TS_ASSERT(!res_sub[203]);
		// TS_ASSERT(!res_sub[218]);
		TS_ASSERT(!res_sub[275]);

		TS_ASSERT(res_sub[218]);
		TS_ASSERT(res_sub[226]);
		TS_ASSERT(res_sub[235]);
		TS_ASSERT(res_sub[238]);

		//Select none of the metals (none of the metals are in the selection)
		selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_none;
		ss_none << "<LigandMetalContactSelector name=\"metal\" ><Index resnums=\"1,2,3\" /></LigandMetalContactSelector>" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_none );
		TR << "Tag with embedded selector that should not include any of the metal ions" << std::endl;
		TS_ASSERT_THROWS_NOTHING( selector.parse_my_tag( tag, datamap) );

		//Apply the mover
		res_sub = selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT(!res_sub[46]);
		TS_ASSERT(!res_sub[48]);
		TS_ASSERT(!res_sub[63]);
		TS_ASSERT(!res_sub[120]);

		TS_ASSERT(!res_sub[63]);
		TS_ASSERT(!res_sub[71]);
		TS_ASSERT(!res_sub[80]);
		TS_ASSERT(!res_sub[83]);

		TS_ASSERT(!res_sub[201]);
		TS_ASSERT(!res_sub[203]);
		TS_ASSERT(!res_sub[218]);
		TS_ASSERT(!res_sub[275]);

		TS_ASSERT(!res_sub[218]);
		TS_ASSERT(!res_sub[226]);
		TS_ASSERT(!res_sub[235]);
		TS_ASSERT(!res_sub[238]);

	}

	void test_parse_tag_named_selector(){
		//This should only set up metal ions specified by the selector
		basic::datacache::DataMap datamap;
		//Add the selectors to the data map
		//Names will be all, chain, zinc, and none
		core::select::residue_selector::ResidueIndexSelectorOP all( new core::select::residue_selector::ResidueIndexSelector );
		all->set_index( "154,155,309,310" );
		datamap.add( "ResidueSelector", "all", all );
		core::select::residue_selector::ResidueIndexSelectorOP none( new core::select::residue_selector::ResidueIndexSelector );
		none->set_index( "1,2,3" );
		datamap.add( "ResidueSelector", "none", none );
		core::select::residue_selector::ResidueNameSelectorOP zinc( new core::select::residue_selector::ResidueNameSelector );
		zinc->set_residue_name3( " ZN" );
		datamap.add( "ResidueSelector", "zinc", zinc );
		core::select::residue_selector::ChainSelectorOP chain( new core::select::residue_selector::ChainSelector );
		utility::vector1< std::string > chain_ids;
		chain_ids.push_back( "A" );
		chain->set_chain_strings( chain_ids );
		datamap.add( "ResidueSelector", "chain", chain );

		//Select all metals
		protocols::residue_selectors::LigandMetalContactSelector new_selector;
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_all;
		ss_all << "<LigandMetalContactSelector name=\"metal\" residue_selector=\"all\" />" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_all );
		TR << "Tag with embedded selector to select all metals" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap ) );

		//Apply the mover
		core::select::residue_selector::ResidueSubset res_sub = new_selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( res_sub[46] );
		TS_ASSERT( res_sub[48] );
		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[120] );

		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[71] );
		TS_ASSERT( res_sub[80] );
		TS_ASSERT( res_sub[83] );

		TS_ASSERT( res_sub[201] );
		TS_ASSERT( res_sub[203] );
		TS_ASSERT( res_sub[218] );
		TS_ASSERT( res_sub[275] );

		TS_ASSERT( res_sub[218] );
		TS_ASSERT( res_sub[226] );
		TS_ASSERT( res_sub[235] );
		TS_ASSERT( res_sub[238] );

		//Select metals only on chain A
		new_selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_a;
		ss_a << "<LigandMetalContactSelector name=\"metal\" residue_selector=\"chain\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_a );
		TR << "Tag with embedded selector to select chain A" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap ) );

		//Apply the mover
		res_sub = new_selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( res_sub[46] );
		TS_ASSERT( res_sub[48] );
		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[120] );

		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[71] );
		TS_ASSERT( res_sub[80] );
		TS_ASSERT( res_sub[83] );

		TS_ASSERT( !res_sub[201] );
		TS_ASSERT( !res_sub[203] );
		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[275] );

		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[226] );
		TS_ASSERT( !res_sub[235] );
		TS_ASSERT( !res_sub[238] );

		//Select only zinc ions
		new_selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_zn;
		ss_zn << "<LigandMetalContactSelector name=\"metal\" residue_selector=\"zinc\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_zn );
		TR << "Tag with embedded selector to select zinc ions" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap ) );

		//Apply the mover
		res_sub = new_selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( !res_sub[46] );
		TS_ASSERT( !res_sub[48] );
		//TS_ASSERT( !res_sub[63] );
		TS_ASSERT( !res_sub[120] );

		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[71] );
		TS_ASSERT( res_sub[80] );
		TS_ASSERT( res_sub[83] );

		TS_ASSERT( !res_sub[201] );
		TS_ASSERT( !res_sub[203] );
		//TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[275] );

		TS_ASSERT( res_sub[218] );
		TS_ASSERT( res_sub[226] );
		TS_ASSERT( res_sub[235] );
		TS_ASSERT( res_sub[238] );

		//Select none of the metals (none of the metals are in the selection)
		new_selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_none;
		ss_none << "<LigandMetalContactSelector name=\"metal\" residue_selector=\"none\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_none );
		TR << "Tag with embedded selector that should not include any of the metal ions" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap) );

		//Apply the mover
		res_sub = new_selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( !res_sub[46] );
		TS_ASSERT( !res_sub[48] );
		TS_ASSERT( !res_sub[63] );
		TS_ASSERT( !res_sub[120] );

		TS_ASSERT( !res_sub[63] );
		TS_ASSERT( !res_sub[71] );
		TS_ASSERT( !res_sub[80] );
		TS_ASSERT( !res_sub[83] );

		TS_ASSERT( !res_sub[201] );
		TS_ASSERT( !res_sub[203] );
		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[275] );

		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[226] );
		TS_ASSERT( !res_sub[235] );
		TS_ASSERT( !res_sub[238] );
	}

	void test_parse_tag_resnum_list(){
		basic::datacache::DataMap datamap;
		//Select all metals
		protocols::residue_selectors::LigandMetalContactSelector new_selector;
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_all;
		ss_all << "<LigandMetalContactSelector name=\"metal\" resnums=\"154,155,309,310\" />" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss_all );
		TR << "Tag with embedded selector to select all metals" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap) );

		//Apply the mover
		core::select::residue_selector::ResidueSubset res_sub = new_selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( res_sub[46] );
		TS_ASSERT( res_sub[48] );
		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[120] );

		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[71] );
		TS_ASSERT( res_sub[80] );
		TS_ASSERT( res_sub[83] );

		TS_ASSERT( res_sub[201] );
		TS_ASSERT( res_sub[203] );
		TS_ASSERT( res_sub[218] );
		TS_ASSERT( res_sub[275] );

		TS_ASSERT( res_sub[218] );
		TS_ASSERT( res_sub[226] );
		TS_ASSERT( res_sub[235] );
		TS_ASSERT( res_sub[238] );

		//Select metals only on chain A
		new_selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_a;
		ss_a << "<LigandMetalContactSelector name=\"metal\" resnums=\"1-155\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_a );
		TR << "Tag with embedded selector to select chain A" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap) );

		//Apply the mover
		res_sub = new_selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( res_sub[46] );
		TS_ASSERT( res_sub[48] );
		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[120] );

		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[71] );
		TS_ASSERT( res_sub[80] );
		TS_ASSERT( res_sub[83] );

		TS_ASSERT( !res_sub[201] );
		TS_ASSERT( !res_sub[203] );
		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[275] );

		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[226] );
		TS_ASSERT( !res_sub[235] );
		TS_ASSERT( !res_sub[238] );

		//Select only zinc ions
		new_selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_zn;
		ss_zn << "<LigandMetalContactSelector name=\"metal\" resnums=\"155,310\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_zn );
		TR << "Tag with embedded selector to select zinc ions" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap ) );

		//Apply the mover
		res_sub = new_selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( !res_sub[46] );
		TS_ASSERT( !res_sub[48] );
		//TS_ASSERT( !res_sub[63] );
		TS_ASSERT( !res_sub[120] );

		TS_ASSERT( res_sub[63] );
		TS_ASSERT( res_sub[71] );
		TS_ASSERT( res_sub[80] );
		TS_ASSERT( res_sub[83] );

		TS_ASSERT( !res_sub[201] );
		TS_ASSERT( !res_sub[203] );
		//TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[275] );

		TS_ASSERT( res_sub[218] );
		TS_ASSERT( res_sub[226] );
		TS_ASSERT( res_sub[235] );
		TS_ASSERT( res_sub[238] );

		//Select none of the metals (none of the metals are in the selection)
		new_selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_none;
		ss_none << "<LigandMetalContactSelector name=\"metal\"  resnums=\"1,2,3\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_none );
		TR << "Tag with embedded selector that should not include any of the metal ions" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap) );

		//Apply the mover
		res_sub = new_selector.apply( pose_ );
		//Check that all bonds were created
		TS_ASSERT( !res_sub[46] );
		TS_ASSERT( !res_sub[48] );
		TS_ASSERT( !res_sub[63] );
		TS_ASSERT( !res_sub[120] );

		TS_ASSERT( !res_sub[63] );
		TS_ASSERT( !res_sub[71] );
		TS_ASSERT( !res_sub[80] );
		TS_ASSERT( !res_sub[83] );

		TS_ASSERT( !res_sub[201] );
		TS_ASSERT( !res_sub[203] );
		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[275] );

		TS_ASSERT( !res_sub[218] );
		TS_ASSERT( !res_sub[226] );
		TS_ASSERT( !res_sub[235] );
		TS_ASSERT( !res_sub[238] );

	}

	void test_dist_cutoff_multiplier(){

		basic::datacache::DataMap datamap;
		//Testing the distance cutoff multiplier providing a number
		protocols::residue_selectors::LigandMetalContactSelector new_selector;
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_none;
		ss_none << "<LigandMetalContactSelector name=\"metal\" resnums=\"1,2,3\" dist_cutoff_multiplier=\"3\" />" << std::endl;
		utility::tag::TagOP tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_none );
		TR << "Tag with embedded selector that should not include any of the metal ions" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap) );

		//Test the values are equal
		TS_ASSERT_EQUALS(new_selector.get_dist_cutoff_multiplier(), 3);


		//Testing the distance cutoff multiplier default
		new_selector = protocols::residue_selectors::LigandMetalContactSelector();
		//This should behave like above except it should override default values, including command line options
		std::stringstream ss_none2;
		ss_none2 << "<LigandMetalContactSelector name=\"metal\" resnums=\"1,2,3\" />" << std::endl;
		tag = utility::tag::TagOP( new utility::tag::Tag() );
		tag->read( ss_none2 );
		TR << "Tag with embedded selector that should not include any of the metal ions" << std::endl;
		TS_ASSERT_THROWS_NOTHING( new_selector.parse_my_tag( tag, datamap) );

		//Test the values are equal
		TS_ASSERT_EQUALS(new_selector.get_dist_cutoff_multiplier(), 1);


	}



};
