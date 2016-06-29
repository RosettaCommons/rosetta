// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/denovo_design/components/StructureDataFactoryTests.cxxtest.hh
/// @brief  tests for the StructureDataFactory, which handles creation/storage of StructureData objects in poses
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Protocol Headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>

// Core Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

// Basic/Utility Headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>

// Boost Headers
#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>

static THREAD_LOCAL basic::Tracer TR("StructureDataFactoryTests");

using namespace protocols::denovo_design::components;

class StructureDataFactoryTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}

	void test_structuredata_io()
	{
		StructureDataFactory const & factory = *( StructureDataFactory::get_instance() );

		core::pose::Pose pdbpose;
		core::io::pdb::build_pose_from_pdb_as_is( pdbpose, "protocols/denovo_design/test_pdbcomp.pdb" );

		core::pose::Pose const input_pose( pdbpose );

		StructureDataCOP const_read_perm = factory.create_from_pose( input_pose, "UnitTest" );
		StructureData read_perm = factory.get_from_pose( pdbpose, "UnitTest" );

		TS_ASSERT_THROWS_NOTHING( read_perm.check_pose_consistency( pdbpose ) );
		TS_ASSERT( input_pose.pdb_info() );
		TS_ASSERT( pdbpose.pdb_info() );
		core::io::Remarks const & rem1 = input_pose.pdb_info()->remarks();
		core::io::Remarks const & rem2 = pdbpose.pdb_info()->remarks();

		// StructureDataFactory::REMARK_NUM should not be present in rem2
		// But otherwise, order should be preserved
		core::io::Remarks::const_iterator r1 = rem1.begin();
		core::io::Remarks::const_iterator r2 = rem2.begin();

		while ( ( r1 != rem1.end() ) && ( r2 != rem2.end() ) ) {
			TS_ASSERT_DIFFERS( r2->num, StructureDataFactory::REMARK_NUM );
			if ( r1->num == StructureDataFactory::REMARK_NUM ) {
				++r1;
				continue;
			}
			TS_ASSERT_EQUALS( r1->num, r2->num );
			TS_ASSERT_EQUALS( r1->value, r2->value );
			++r1;
			++r2;
		}

		// save into pose
		std::string const new_id = "NewID";
		read_perm.set_id( new_id );
		TS_ASSERT_EQUALS( read_perm.id(), new_id );

		factory.save_into_pose( pdbpose, read_perm );
		StructureData const & saved = factory.get_from_pose( pdbpose, "UnitTest" );
		TS_ASSERT_EQUALS( saved.id(), read_perm.id() );

		// test stored version vs original
		std::stringstream orig;
		orig << read_perm;
		std::stringstream stored;
		stored << saved;
		TS_ASSERT_EQUALS( orig.str(), stored.str() );
	}

	void test_silent_io()
	{
		StructureDataFactory const & factory = *( StructureDataFactory::get_instance() );

		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/test_pdbcomp.pdb" );

		std::string const new_id = "StructureData_SilentFile_IO";
		StructureData read_perm = factory.get_from_pose( input_pose, new_id );
		read_perm.set_id( new_id );

		// now write data to the pose and test
		factory.save_into_pose( input_pose, read_perm );

		std::stringstream origss;
		origss << read_perm;

		core::io::silent::SilentFileData sfd;
		for ( core::Size z=1; z<=1; ++z ) {
			core::io::silent::SilentStructOP silent =
				core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( input_pose );
			silent->fill_struct( input_pose, '_' + boost::lexical_cast< std::string >( z ) );
			sfd.add_structure( *silent );
		}
		sfd.write_all( "test.silent" );

		core::io::silent::SilentFileData sftest;
		sftest.read_file( "test.silent" );
		typedef utility::vector1< core::io::silent::SilentStructOP > SilentStructOPs;
		SilentStructOPs read_structures = sftest.structure_list();
		TR << "Read " << read_structures.size() << " structures" << std::endl;
		for ( SilentStructOPs::const_iterator s=read_structures.begin(); s!=read_structures.end(); ++s ) {
			core::pose::Pose pose_from_silent;
			(*s)->set_residue_numbers( utility::vector1< core::Size >() );
			(*s)->fill_pose( pose_from_silent );
			StructureData const & newperm = factory.get_from_pose( pose_from_silent, "UnitTest" );
			TS_ASSERT( pose_from_silent.pdb_info() );
			TS_ASSERT( factory.observer_attached( pose_from_silent ) );
			TS_ASSERT_EQUALS( newperm.id(), new_id );

			// test vs original
			std::stringstream newss;
			newss << newperm;
			TS_ASSERT_EQUALS( newss.str(), origss.str() );
		}

		TS_ASSERT( !remove( "test.silent" ) );
	}

	void test_xml_io()
	{
		StructureDataFactory const & factory = *( StructureDataFactory::get_instance() );

		utility::io::izstream tag_xml( "protocols/denovo_design/test_sd.xml" );
		utility::tag::TagCOP tag1 = utility::tag::Tag::create( tag_xml );

		utility::io::izstream input_xml( "protocols/denovo_design/test_sd.xml" );
		StructureDataOP sd = factory.create_from_xml( input_xml );
		TS_ASSERT( sd );
		TS_ASSERT_EQUALS( tag1->getName(), "StructureData" );
		TS_ASSERT_EQUALS( tag1->getOption< core::Size >( "length" ), sd->length() );
		TS_ASSERT_EQUALS( tag1->getOption< core::Size >( "pose_length" ), sd->pose_length() );

		// replace 'UnitTest' name with cat_sheet (the original name) so that the strings can be directly compared
		test::UTracer UT( "protocols/denovo_design/test_sd.xml" );
		UT << *sd;
	}

};

