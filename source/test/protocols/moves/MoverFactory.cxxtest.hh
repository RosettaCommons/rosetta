// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MoverFactory.cxxtest.hh
/// @brief  test suite for protocols::moves::MoverFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/types.hh>

// Unit headers
#include <protocols/moves/MoverFactory.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/mover_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Test utility headers
#include <test/util/schema_utilities.hh>

using namespace protocols::moves;
using namespace utility::tag;

class MoverFactoryTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_MoverFactory_all_Mover_complexTypes_have_descriptions() {
		ensure_all_cts_for_creators_have_documentation_strings(
			MoverFactory::get_instance()->mover_creator_map(),
			"Mover",
			& complex_type_name_for_mover );
	}

	void test_MoverFactory_all_Mover_all_attributes_have_descriptions() {
		MoverFactory * mf = MoverFactory::get_instance();
		MoverFactory::MoverMap const & creator_map = mf->mover_creator_map();

		for ( auto iter : creator_map ) {
			XMLSchemaDefinition xsd;
			iter.second->provide_xml_schema( xsd );
			std::string full_def = xsd.full_definition();
			//std::cout << "full def: " << iter.first << "\n" << full_def << std::endl;
			TagCOP tag( Tag::create( full_def ) );
			recurse_through_subtags_for_attribute_descriptions( tag, iter.first );
		}
	}

	void test_all_movers_define_valid_xsds() {
		XMLSchemaDefinition xsd;
		MoverFactory::get_instance()->define_mover_xml_schema( xsd );
		ensure_rosetta_scripts_like_XSD_validates_w_group( xsd, MoverFactory::mover_xml_schema_group_name() );
	}

};

