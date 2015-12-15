// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/denovo_design/StructureData.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::components::StructureData
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>

// Protocol headers
#include <protocols/denovo_design/util.hh>

// Core headers
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/sequence/ABEGOManager.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Boost
#include <boost/algorithm/string/replace.hpp>
#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>

// C++ headers

static THREAD_LOCAL basic::Tracer TR("protocols.denovo_design.StructureData.cxxtest");

// --------------- Test Class --------------- //
class StructureDataTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init_with_additional_options( "-extra_res_fa protocols/denovo_design/QTS.params -extra_patch_fa protocols/denovo_design/QTS_connectC6.txt -extra_patch_fa protocols/denovo_design/CYS_connectSG.txt" );

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_structuredata_io() {
		using namespace protocols::denovo_design::components;
		core::pose::Pose input_pose;
		core::pose::PoseOP pdbpose = core::pose::PoseOP( new core::pose::Pose(input_pose) );
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/test_pdbcomp.pdb" );
		core::io::pdb::build_pose_from_pdb_as_is( *pdbpose, "protocols/denovo_design/test_pdbcomp.pdb" );
		StructureDataOP read_perm = StructureData::create_from_pose( *pdbpose, "UnitTest" );
		TS_ASSERT( read_perm );
		TS_ASSERT( input_pose.pdb_info() );
		TS_ASSERT( pdbpose->pdb_info() );
		core::pose::Remarks const & rem1 = input_pose.pdb_info()->remarks();
		core::pose::Remarks const & rem2 = pdbpose->pdb_info()->remarks();
		TS_ASSERT_EQUALS( rem1.size(), rem2.size() );
		// start at 0 because Remarks derives from std::vector for some reason
		for ( core::Size i=0; i<rem1.size(); ++i ) {
			if ( rem1[i].num == StructureData::REMARK_NUM ) {
				TS_ASSERT_EQUALS( rem1[i].num, rem2[i].num );
				TS_ASSERT_EQUALS( rem1[i].value, rem2[i].value );
			}
		}

		// now write data to the pose and test
		// data includes silent files
		read_perm->save_into_pose( *pdbpose );
		core::pose::Pose testpose = *(read_perm->pose());
		TS_ASSERT( StructureData::has_cached_string( *(read_perm->pose()) ) );

		std::stringstream silentfile;
		core::io::silent::SilentFileData sfd;
		for ( core::Size z=1; z<=1; ++z ) {
			core::io::silent::SilentStructOP silent =
				core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( *(read_perm->pose()) );
			silent->fill_struct( *(read_perm->pose()), '_' + boost::lexical_cast< std::string >( z ) );
			sfd.add_structure( *silent );
		}
		sfd.write_all( "test.silent" );
		std::stringstream origstr;
		origstr << *read_perm;
		TS_ASSERT_EQUALS( StructureData::cached_string( *(read_perm->pose()) ), origstr.str() );

		core::io::silent::SilentFileData sftest;
		sftest.read_file( "test.silent" );
		utility::vector1< core::io::silent::SilentStructOP > read_structures = sftest.structure_list();
		TR << "Read " << read_structures.size() << " structures" << std::endl;
		for ( core::Size i=1, endi=read_structures.size(); i<=endi; ++i ) {
			read_structures[i]->set_residue_numbers( utility::vector1< core::Size >() );
			core::pose::Pose newpose;
			read_structures[i]->fill_pose( newpose );
			TR << "Cached string : " << StructureData::cached_string( newpose ) << std::endl;
			StructureDataOP newperm = StructureData::create_from_pose( newpose, "UnitTest" );
			std::stringstream newss;
			newss << *newperm;
			TS_ASSERT_EQUALS( newss.str(), origstr.str() );
		}

		TS_ASSERT( !remove( "test.silent" ) );

		core::pose::Remarks const & rem3 = pdbpose->pdb_info()->remarks();

		// extract tags for both
		utility::vector1< std::string > lines;
		for ( core::pose::Remarks::const_iterator it_rem=rem1.begin(), it_end=rem1.end(); it_rem != it_end; ++it_rem ) {
			if ( it_rem->num != StructureData::REMARK_NUM ) {
				continue;
			}
			lines.push_back( protocols::denovo_design::get_remark_line( it_rem, it_end ) );
		}
		// create full xml tag
		std::stringstream xmltag1;
		xmltag1 << utility::join( lines, "\n" );
		utility::tag::TagOP tag1 = utility::tag::Tag::create( xmltag1 );

		lines.clear();
		for ( core::pose::Remarks::const_iterator it_rem=rem3.begin(), it_end=rem3.end(); it_rem != it_end; ++it_rem ) {
			if ( it_rem->num != StructureData::REMARK_NUM ) {
				continue;
			}
			lines.push_back( protocols::denovo_design::get_remark_line( it_rem, it_end ) );
		}
		// create full xml tag
		std::stringstream xmltag3;
		xmltag3 << utility::join( lines, "\n" );
		utility::tag::TagOP tag3 = utility::tag::Tag::create( xmltag3 );

		/// test permutation tag options
		TS_ASSERT_EQUALS( tag1->getName(), "StructureData" );
		TS_ASSERT_EQUALS( tag1->getName(), tag3->getName() );
		TS_ASSERT_EQUALS( tag3->getOption< std::string >( "name" ), "UnitTest" );
		TS_ASSERT_EQUALS( tag1->getOption< bool >( "multi" ), tag3->getOption< bool >( "multi" ) );
		TS_ASSERT_EQUALS( tag1->getOption< core::Size >( "length" ), tag3->getOption< core::Size >( "length" ) );
		TS_ASSERT_EQUALS( tag1->getOption< core::Size >( "pose_length" ), tag3->getOption< core::Size >( "pose_length" ) );

		utility::vector0< utility::tag::TagCOP >::const_iterator t1 = tag1->getTags().begin(), end = tag1->getTags().end(),
			t3 = tag3->getTags().begin(), end3 = tag3->getTags().end();
		// same # of tags
		TS_ASSERT_EQUALS( std::distance( t3, end3 ), std::distance( t1, end ) );
		std::stringstream stream1, stream3;
		stream1 << *tag1;
		stream3 << *tag3;
		// replace 'UnitTest' name with cat_sheet (the original name) so that the strings can be directly compared
		std::string string3 = stream3.str();
		boost::replace_all( string3, "UnitTest", "cat_sheet" );
		TS_ASSERT_EQUALS( stream1.str(), string3 );
	}

	void test_delete_segment()
	{
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		std::ifstream in_xml( "protocols/denovo_design/test_sd.xml" );

		StructureDataOP perm = StructureData::create_from_xml( in_xml, "UnitTest" );
		StructureDataOP orig = perm->clone();

		for ( StringList::const_iterator c = perm->segments_begin(); c != perm->segments_end(); ++c ) {
			TS_ASSERT( perm->segment( *c ).stop() >= perm->segment( *c ).start() );
			TS_ASSERT( perm->segment( *c ).safe() >= perm->segment( *c ).start() );
		}

		// delete bb.h1_sheet -- it shouldn't mess up any residues
		perm->delete_segment( "bb.h1_sheet" );
		for ( StringList::const_iterator c = perm->segments_begin(); c != perm->segments_end(); ++c ) {
			TS_ASSERT( orig->has_segment( *c ) );
			TS_ASSERT( perm->segment( *c ).stop() >= perm->segment( *c ).start() );
			TS_ASSERT( perm->segment( *c ).safe() >= perm->segment( *c ).start() );
			TS_ASSERT_EQUALS( perm->segment( *c ).stop() - perm->segment( *c ).start(),
				orig->segment( *c ).stop() - orig->segment( *c ).start() );
			TS_ASSERT_EQUALS( perm->segment( *c ).safe() - perm->segment( *c ).start(),
				orig->segment( *c ).safe() - orig->segment( *c ).start() );
		}
	}

	void test_non_peptidic_bonds()
	{
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/test_structuredata_nonpeptidic.pdb" );
		StructureDataOP sd = StructureData::create_from_pose( input_pose, "" );

		core::conformation::Residue const prev_cys = input_pose.residue( sd->alias_resnum( "bb.cys" ) );
		core::conformation::Residue const prev_qts = input_pose.residue( sd->alias_resnum( "bb.ligand" ) );
		sd->consolidate_movable_groups( boost::assign::list_of("1") );
		TS_ASSERT( sd->pose()->fold_tree().check_fold_tree() );
		input_pose.fold_tree( sd->pose()->fold_tree() );

		// rotate a chi angle in QTS
		core::Size const cysres = sd->alias_resnum( "bb.cys" );
		core::Size const ligres = sd->alias_resnum( "bb.ligand" );
		for ( core::Size chi=1; chi<=6; ++chi ) {
			input_pose.set_chi( chi, ligres, 50.0 );
		}

		TS_ASSERT_DELTA( prev_cys.xyz( "CB" ).distance( input_pose.residue( cysres ).xyz( "CB" ) ), 0.0, 1e-4 );
		TS_ASSERT_DELTA( prev_cys.xyz( "SG" ).distance( input_pose.residue( cysres ).xyz( "SG" ) ), 0.0, 1e-4 );
		TS_ASSERT_DELTA( prev_qts.xyz( "C6" ).distance( input_pose.residue( ligres ).xyz( "C6" ) ), 0.0, 1e-4 );
		TS_ASSERT_DELTA( prev_qts.xyz( "O2" ).distance( input_pose.residue( ligres ).xyz( "O2" ) ), 0.0, 1e-4 );
		input_pose.dump_pdb( "test.pdb" );

		// rotate chi in cys
		input_pose.set_chi( 1, cysres, 50.0 );

		TS_ASSERT_DELTA( prev_cys.xyz( "CB" ).distance( input_pose.residue( cysres ).xyz( "CB" ) ), 0.0, 1e-4 );
		TS_ASSERT( prev_cys.xyz( "SG" ).distance( input_pose.residue( cysres ).xyz( "SG" ) ) > 1e-4 );
		TS_ASSERT( prev_qts.xyz( "C6" ).distance( input_pose.residue( ligres ).xyz( "C6" ) ) > 1e-4 );
		TS_ASSERT( prev_qts.xyz( "O2" ).distance( input_pose.residue( ligres ).xyz( "O2" ) ) > 1e-4 );
		TR << input_pose.fold_tree() << std::endl;
		input_pose.dump_pdb( "test2.pdb" );
	}
};
