// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/FilterFactory.cxxtest.hh
/// @brief  test suite for protocols::filters::FilterFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/types.hh>

// Unit headers
#include <protocols/filters/FilterFactory.hh>

// Package headers
#include <protocols/filters/Filter.hh>
#include <protocols/filters/filter_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Test utility headers
#include <test/util/schema_utilities.hh>

using namespace protocols::filters;
using namespace utility::tag;

class FilterFactoryTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_FilterFactory_all_Filter_complexTypes_have_descriptions() {
		ensure_all_cts_for_creators_have_documentation_strings(
			FilterFactory::get_instance()->filter_creator_map(),
			"Filter",
			& complex_type_name_for_filter );
	}

	void test_FilterFactory_all_Filter_all_attributes_have_descriptions() {
		FilterFactory * mf = FilterFactory::get_instance();
		FilterFactory::FilterMap const & creator_map = mf->filter_creator_map();

		for ( auto iter : creator_map ) {
			XMLSchemaDefinition xsd;
			iter.second->provide_xml_schema( xsd );
			std::string full_def = xsd.full_definition();
			//std::cout << "full def: " << iter.first << "\n" << full_def << std::endl;
			TagCOP tag( Tag::create( full_def ) );
			recurse_through_subtags_for_attribute_descriptions( tag, iter.first );
		}
	}

	void test_all_filters_define_valid_xsds() {
		XMLSchemaDefinition xsd;
		FilterFactory::get_instance()->define_filter_xml_schema( xsd );
		ensure_rosetta_scripts_like_XSD_validates_w_group( xsd, FilterFactory::filter_xml_schema_group_name() );
	}


};

