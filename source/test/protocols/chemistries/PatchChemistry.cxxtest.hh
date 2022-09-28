// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/PatchChemistry.cxxtest.hh
/// @brief  test for PatchChemistry
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/schema_utilities.hh>

// Project Headers
#include <protocols/chemistries/PatchChemistry.hh>

#include <protocols/chemistries/ChemistryFactory.hh>
#include <protocols/chemistries/util.hh>

#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <string>

static basic::Tracer TR("protocols.chemistries.PatchChemistry.cxxtest.hh");

// --------------- Test Class --------------- //

class PatchChemistryTests : public CxxTest::TestSuite {

private:
public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_apply_patch() {
		protocols::chemistries::PatchChemistry patchchem;
		patchchem.patch_file("core/chemical/1pqc_test.patch");
		patchchem.add_patch_operation_line("SET_INTERCHANGEABILITY_GROUP FOO");
		patchchem.add_patch_operation_line("NBR_RADIUS 10.0");

		core::chemical::ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		core::chemical::MutableResidueTypeOP restype = core::chemical::read_topology_file("core/chemical/1pqc.params", rts);

		TS_ASSERT( restype->has( "H12" ) );
		TS_ASSERT( !restype->has( "OP1" ) );
		TS_ASSERT_EQUALS( restype->interchangeability_group(), "QC1" );
		TS_ASSERT_EQUALS( restype->nbr_radius(), 8.924279 );

		TR << "Restype loaded" << std::endl;

		patchchem.apply( *restype );

		TR << "Restype patched" << std::endl;

		TS_ASSERT( !restype->has( "H12" ) );
		TS_ASSERT( restype->has( "OP1" ) );
		TS_ASSERT_EQUALS( restype->interchangeability_group(), "FOO" );
		TS_ASSERT_EQUALS( restype->nbr_radius(), 10.0 );

	}

	void test_xml() {
		basic::datacache::DataMap data;
		std::istringstream tag_stream(R"raw_string(
			<PatchChemistry name="patch" patch_file="core/chemical/1pqc_test.patch" >
				<Op line="SET_INTERCHANGEABILITY_GROUP FOO" />
				<Op line="NBR_RADIUS 10.0" />
			</PatchChemistry>
		)raw_string");

		utility::tag::TagCOP tag = utility::tag::Tag::create(tag_stream);

		protocols::chemistries::ChemistryOP patchchem = protocols::chemistries::ChemistryFactory::get_instance()->new_chemistry( tag, data );

		core::chemical::ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		core::chemical::MutableResidueTypeOP restype = core::chemical::read_topology_file("core/chemical/1pqc.params", rts);

		TS_ASSERT( restype->has( "H12" ) );
		TS_ASSERT( !restype->has( "OP1" ) );
		TS_ASSERT_EQUALS( restype->interchangeability_group(), "QC1" );
		TS_ASSERT_EQUALS( restype->nbr_radius(), 8.924279 );

		TR << "Restype loaded" << std::endl;

		patchchem->apply( *restype );

		TR << "Restype patched" << std::endl;

		TS_ASSERT( !restype->has( "H12" ) );
		TS_ASSERT( restype->has( "OP1" ) );
		TS_ASSERT_EQUALS( restype->interchangeability_group(), "FOO" );
		TS_ASSERT_EQUALS( restype->nbr_radius(), 10.0 );
	}


	void test_xsd() {
		using namespace protocols::chemistries;

		std::string const tag = R"raw_string(
			<PatchChemistry name="patch" patch_file="core/chemical/1pqc_test.patch" >
				<Op line="SET_INTERCHANGEABILITY_GROUP FOO" />
				<Op line="NBR_RADIUS 10.0" />
			</PatchChemistry>
		)raw_string";

		utility::tag::XMLSchemaDefinition xsd;
		check_if_tag_validates<PatchChemistry>( tag, xsd, PatchChemistry::class_name(), complex_type_name_for_chemistry( PatchChemistry::class_name() ) );
	}
};
