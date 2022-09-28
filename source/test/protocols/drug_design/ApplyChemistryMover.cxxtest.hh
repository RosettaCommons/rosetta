// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/drug_design/ApplyChemistryMover.cxxtest.hh
/// @brief  test for ApplyChemistryMover
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/schema_utilities.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/drug_design/ApplyChemistryMover.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/residue_io.hh>
#include <protocols/chemistries/WrappedChemistries.hh>
#include <protocols/chemistries/ChemistryFactory.hh>
#include <protocols/chemistries/PatchChemistry.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>
#include <core/chemical/sdf/MolFileIOData.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <string>

static basic::Tracer TR("protocols.drug_design.ApplyChemistryMover.cxxtest.hh");

// --------------- Test Class --------------- //

class ApplyChemistryMoverTests : public CxxTest::TestSuite {

private:
public:

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	void test_apply_chemistry() {
		protocols::drug_design::ApplyChemistryMover appchemmover;
		protocols::chemistries::PatchChemistryOP chem = utility::pointer::make_shared< protocols::chemistries::PatchChemistry >();

		chem->add_patch_operation_line("SET_INTERCHANGEABILITY_GROUP XYZ" );
		chem->add_patch_operation_line("NBR_RADIUS 31.415" );

		appchemmover.add_chemistry( chem );
		appchemmover.residue_id( 1 );

		core::chemical::sdf::MolFileIOReader molfile_reader;

		utility::vector1< core::chemical::sdf::MolFileIOMoleculeOP > molfile_data = molfile_reader.parse_file( "protocols/drug_design/carboxypyrimidine.sdf" );
		utility::vector1< core::chemical::MutableResidueTypeOP > alltypes = convert_to_ResidueTypes( molfile_data, false );
		TS_ASSERT_EQUALS( alltypes.size(), 1 );
		core::chemical::ResidueTypeCOP restype( core::chemical::ResidueType::make(*alltypes[1]) );

		core::conformation::ResidueOP res( core::conformation::ResidueFactory::create_residue( *restype ) );
		core::pose::Pose pose;
		pose.append_residue_by_jump( *res, 0 );

		TS_ASSERT_DELTA( pose.residue(1).nbr_radius(), 3.8398, 0.001 );
		TS_ASSERT_EQUALS( pose.residue_type(1).interchangeability_group(), "CPD"  );

		appchemmover.apply( pose );

		TS_ASSERT_EQUALS( pose.residue(1).nbr_radius(), 31.415 );
		TS_ASSERT_EQUALS( pose.residue_type(1).interchangeability_group(), "XYZ"  );
	}

	void test_xml() {
		std::istringstream chem_stream(R"raw_string("
				<PatchChemistry >
					<Op line="SET_INTERCHANGEABILITY_GROUP ABC" />
					<Op line="NBR_RADIUS 4.13" />
				</PatchChemistry>
		)raw_string");
		utility::tag::TagCOP chem_tag = utility::tag::Tag::create(chem_stream);

		basic::datacache::DataMap data;

		protocols::chemistries::ChemistryOP reprot = protocols::chemistries::ChemistryFactory::get_instance()->new_chemistry( chem_tag, data );
		data.add( "chemistry", "reprot", reprot );

		std::istringstream tag_stream(R"raw_string(<ApplyChemistryMover name="chem" chemistry="reprot" residue="1" />)raw_string");
		utility::tag::TagCOP tag = utility::tag::Tag::create(tag_stream);

		protocols::moves::MoverOP mover = protocols::moves::MoverFactory::get_instance()->newMover( tag, data );

		protocols::drug_design::ApplyChemistryMoverOP acm = utility::pointer::dynamic_pointer_cast<protocols::drug_design::ApplyChemistryMover>(mover);

		TS_ASSERT_DIFFERS( acm, nullptr );

		utility::vector1< protocols::chemistries::ChemistryCOP > const & chemistries = acm->chemistries();
		TS_ASSERT_EQUALS( chemistries.size(), 1 );

		auto subchem = utility::pointer::dynamic_pointer_cast< protocols::chemistries::PatchChemistry const >( chemistries[1] );

		TS_ASSERT_DIFFERS( subchem, nullptr );
	}

	void test_xml_nested() {
		std::string const tag = R"raw_string(
			<ApplyChemistryMover name="chem" residue="1A" >
				<PatchChemistry >
					<Op line="SET_INTERCHANGEABILITY_GROUP FOO" />
					<Op line="NBR_RADIUS 10.0" />
				</PatchChemistry>
			</ApplyChemistryMover>
		)raw_string";

		// //Nice idea, but the limitations of check_if_mover_tag_validates() means we can't use it with the ApplyChemistryMover with its subtags.
		// // -- you get an invalid schema error.
		// check_if_mover_tag_validates< protocols::drug_design::ApplyChemistryMover >(tag);
	}
};
