// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/denovo_design/components/ExtendedPoseBuilderTests.cxxtest.hh
/// @brief  Tests for Pose Builders
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
// this should be the only protocols file that includes this cc file
#include <test/protocols/denovo_design/test_utils.cc>

// Unit Headers
#include <protocols/denovo_design/components/ExtendedPoseBuilder.hh>

// Project Headers
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>

// Core Headers
#include <core/conformation/Conformation.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

// Boost Headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR("ExtendedPoseBuilderTests");

using namespace protocols::denovo_design;
using namespace protocols::denovo_design::architects;
using namespace protocols::denovo_design::components;

// this is a copy of the function in core::pose::util.hh, except it uses distances between equivalent atoms
// instead of std::floor() -- this is better for floating point errors, where std::floor( 191.999999 ) is not
// the same as std::floor( 1.920000 ), even though those numbers are within 6 decimal places.
// I made this because the comments indicate there is a reason for using std::floor() over
// std::abs( x1 - x2 ) < epsilon in the original function in core::pose::util.cc, so I didn't want to change
// it.
bool
compare_atom_coordinates(
	core::pose::Pose const & lhs,
	core::pose::Pose const & rhs,
	core::Size const n_dec_places );

void
check_residue_dihedrals( core::pose::Pose const & pose, core::pose::Pose const & orig, core::Real const tol );

class ExtendedPoseBuilderTests : public CxxTest::TestSuite {
public:
	void setUp()
	{
		protocols_init();
	}

	void tearDown(){

	}

	void test_add_chain_from_segment() {
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/twohelix_structuredata.pdb" );
		input_pose = core::pose::Pose( input_pose, 1, 35 );

		StructureData sd = *StructureDataFactory::get_instance()->create_from_pose( input_pose );
		sd.set_template_pose( "H01", input_pose, 2, 34 );

		core::pose::Pose pose;
		append_new_chain_from_template_segment( pose, sd.segment( "H01" ) );

		TS_ASSERT_EQUALS( pose.sequence(), input_pose.sequence() );
		check_residue_dihedrals( pose, input_pose, 1e-5 );
		TS_ASSERT( compare_atom_coordinates( pose, input_pose, 3 ) );
	}

	void test_build_pose_simple() {
		StructureData sd( "test_build_pose" );
		sd.add_segment( "h1", Segment( "LHHHHHHHHHH", "XAAAAAAAAAA", false, true ) );
		Segment l1_seg( "LLLL", "BAGB", true, true );
		l1_seg.set_movable_group( 2 );
		sd.add_segment( "l1", l1_seg );
		sd.add_segment( "h2", Segment( "HHHHHHHHHL", "AAAAAAAAAO", true, false ) );
		sd.connect_segments( "h1", "l1" );
		sd.connect_segments( "l1", "h2" );

		ExtendedPoseBuilder builder;
		core::pose::PoseOP pose = builder.apply( sd );
		TS_ASSERT( pose );
		sd.check_pose_consistency( *pose );
		TS_ASSERT( StructureDataFactory::get_instance()->observer_attached( *pose ) );

		// one chain
		TS_ASSERT_EQUALS( pose->conformation().num_chains(), 1 );

		// proper termini variants
		TS_ASSERT( pose->residue( 1 ).is_lower_terminus() );
		TS_ASSERT( pose->residue( pose->total_residue() ).is_upper_terminus() );
	}

	void test_build_pose_with_lower_template() {
		core::pose::Pose helix15;
		core::io::pdb::build_pose_from_pdb_as_is( helix15, "protocols/denovo_design/components/helix15.pdb" );

		// let's extend c-terminally
		StructureData sd( "test_build_pose" );
		sd.add_segment( "h1", Segment( "LHHHHHHHHHHHHHHH", "XAAAAAAAAAAAAAAA", false, true ) );
		Segment l1_seg( "LLLL", "BAGB", true, true );
		l1_seg.set_movable_group( 2 );
		sd.add_segment( "l1", l1_seg );
		sd.add_segment( "h2", Segment( "HHHHHHHHHL", "AAAAAAAAAO", true, false ) );
		sd.connect_segments( "h1", "l1" );
		sd.connect_segments( "l1", "h2" );
		sd.set_template_pose( "h1", helix15, sd.segment( "h1" ).start(), sd.segment( "h1" ).stop() );

		ExtendedPoseBuilder builder;
		core::pose::PoseOP poseptr = builder.apply( sd );
		TS_ASSERT( poseptr );

		core::pose::Pose const & pose = *poseptr;
		sd.check_pose_consistency( pose );
		TS_ASSERT( StructureDataFactory::get_instance()->observer_attached( pose ) );

		// one chain
		TS_ASSERT_EQUALS( pose.conformation().num_chains(), 1 );

		// proper termini variants
		TS_ASSERT( pose.residue( 1 ).is_lower_terminus() );
		TS_ASSERT( pose.residue( pose.total_residue() ).is_upper_terminus() );

		// phi/psi should match for residues in helix15
		for ( core::Size resid=2; resid<=sd.segment( "h1" ).stop(); ++resid ) {
			TS_ASSERT_DELTA( pose.phi( resid ), helix15.phi( resid ), 1e-5 );
			TS_ASSERT_DELTA( pose.psi( resid ), helix15.psi( resid ), 1e-5 );
			TS_ASSERT_DELTA( pose.omega( resid ), helix15.omega( resid ), 1e-5 );
		}
	}

	void test_build_pose_with_upper_template() {
		core::pose::Pose helix15;
		core::io::pdb::build_pose_from_pdb_as_is( helix15, "protocols/denovo_design/components/helix15.pdb" );

		// let's extend c-terminally
		StructureData sd( "test_build_pose" );
		sd.add_segment( "h2", Segment( "LHHHHHHHHH", "XAAAAAAAAA", false, true ) );
		Segment l1_seg( "LLLL", "BAGB", true, true );
		l1_seg.set_movable_group( 2 );
		sd.add_segment( "l1", l1_seg );
		sd.add_segment( "h1", Segment( "HHHHHHHHHHHHHHHL", "AAAAAAAAAAAAAAAX", true, false ) );
		sd.connect_segments( "h2", "l1" );
		sd.connect_segments( "l1", "h1" );
		sd.set_template_pose( "h1", helix15, 2, 16 );

		ExtendedPoseBuilder builder;
		core::pose::PoseOP poseptr = builder.apply( sd );
		TS_ASSERT( poseptr );

		core::pose::Pose const & pose = *poseptr;
		sd.check_pose_consistency( pose );
		TS_ASSERT( StructureDataFactory::get_instance()->observer_attached( pose ) );

		// one chain
		TS_ASSERT_EQUALS( pose.conformation().num_chains(), 1 );

		// proper termini variants
		TS_ASSERT( pose.residue( 1 ).is_lower_terminus() );
		TS_ASSERT( pose.residue( pose.total_residue() ).is_upper_terminus() );

		// phi/psi should match for residues in helix15
		for ( core::Size resid=2, pose_resid=sd.segment( "h1" ).lower(); resid<=helix15.total_residue()-1; ++resid, ++pose_resid ) {
			TR.Debug << "Testing residue " << pose_resid << " vs " << resid << std::endl;
			TS_ASSERT_DELTA( pose.psi( pose_resid ), helix15.psi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.phi( pose_resid ), helix15.phi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.omega( pose_resid ), helix15.omega( resid ), 1e-6 );
		}
	}

	void test_build_pose_with_middle_template() {
		core::pose::Pose helix15;
		core::io::pdb::build_pose_from_pdb_as_is( helix15, "protocols/denovo_design/components/helix15.pdb" );
		core::pose::Pose helix14;
		core::io::pdb::build_pose_from_pdb_as_is( helix14, "protocols/denovo_design/components/helix14.pdb" );

		// let's extend c-terminally
		StructureData sd( "test_build_pose" );
		sd.add_segment( "h1", Segment( "HHHHHHHHHHHHHHHL", "AAAAAAAAAAAAAAAX", false, true ) );
		Segment l1_seg( "LLLL", "BAGB", true, true );
		l1_seg.set_movable_group( 2 );
		sd.add_segment( "l1", l1_seg );
		sd.add_segment( "h2", Segment( "HHHHHHHHHHHHHHL",  "AAAAAAAAAAAAAAX", true, false ) );
		sd.connect_segments( "h1", "l1" );
		sd.connect_segments( "l1", "h2" );
		sd.set_template_pose( "h1", helix15, 2, 16 );
		sd.set_template_pose( "h2", helix14, 2, 15 );

		ExtendedPoseBuilder builder;
		core::pose::PoseOP poseptr = builder.apply( sd );
		TS_ASSERT( poseptr );

		core::pose::Pose const & pose = *poseptr;
		sd.check_pose_consistency( pose );
		TS_ASSERT( StructureDataFactory::get_instance()->observer_attached( pose ) );

		// this one is two chains
		TS_ASSERT_EQUALS( pose.conformation().num_chains(), 1 );

		// proper termini variants
		TS_ASSERT( pose.residue( 1 ).is_lower_terminus() );
		TS_ASSERT( pose.residue( pose.total_residue() ).is_upper_terminus() );

		// phi/psi should match for residues in helix15
		for ( core::Size resid=2; resid<helix15.total_residue()-1; ++resid ) {
			TR.Debug << "Testing residue " << resid << std::endl;
			TS_ASSERT_DELTA( pose.psi( resid ), helix15.psi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.phi( resid ), helix15.phi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.omega( resid ), helix15.omega( resid ), 1e-6 );
		}

		// phi/psi should also match for helix14
		for ( core::Size resid=2, pose_resid=sd.segment( "h2" ).lower(); resid<=helix14.total_residue()-1; ++resid, ++pose_resid ) {
			TR.Debug << "Testing residue " << pose_resid << " vs " << resid << std::endl;
			TS_ASSERT_DELTA( pose.psi( pose_resid ), helix14.psi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.phi( pose_resid ), helix14.phi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.omega( pose_resid ), helix14.omega( resid ), 1e-6 );
		}
	}

	void test_build_pose_with_middle_template_and_cut() {
		core::pose::Pose helix15;
		core::io::pdb::build_pose_from_pdb_as_is( helix15, "protocols/denovo_design/components/helix15.pdb" );
		core::pose::Pose helix14;
		core::io::pdb::build_pose_from_pdb_as_is( helix14, "protocols/denovo_design/components/helix14.pdb" );

		// let's extend c-terminally
		StructureData sd( "test_build_pose" );
		sd.add_segment( "h1", Segment( "HHHHHHHHHHHHHHHL", "AAAAAAAAAAAAAAAX", false, true ) );
		Segment l1_seg( "LLLL", "BAGB", true, true );
		l1_seg.set_movable_group( 2 );
		sd.add_segment( "l1", l1_seg );
		sd.add_segment( "h2", Segment( "HHHHHHHHHHHHHHL",  "AAAAAAAAAAAAAAX", true, false ) );
		sd.connect_segments( "h1", "l1" );
		sd.connect_segments( "l1", "h2" );
		sd.set_template_pose( "h1", helix15, 2, 16 );
		sd.set_template_pose( "h2", helix14, 2, 15 );
		sd.set_cutpoint( "l1", 2 );

		ExtendedPoseBuilder builder;
		core::pose::PoseOP poseptr = builder.apply( sd );
		TS_ASSERT( poseptr );

		core::pose::Pose const & pose = *poseptr;
		sd.check_pose_consistency( pose );
		TS_ASSERT( StructureDataFactory::get_instance()->observer_attached( pose ) );

		// this one is two chains
		TS_ASSERT_EQUALS( pose.conformation().num_chains(), 1 );

		// proper termini variants
		TS_ASSERT( pose.residue( 1 ).is_lower_terminus() );
		TS_ASSERT( pose.residue( pose.total_residue() ).is_upper_terminus() );

		// phi/psi should match for residues in helix15
		for ( core::Size resid=2; resid<helix15.total_residue()-1; ++resid ) {
			TR.Debug << "Testing residue " << resid << std::endl;
			TS_ASSERT_DELTA( pose.psi( resid ), helix15.psi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.phi( resid ), helix15.phi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.omega( resid ), helix15.omega( resid ), 1e-6 );
		}

		// phi/psi should also match for helix14
		for ( core::Size resid=2, pose_resid=sd.segment( "h2" ).lower(); resid<=helix14.total_residue()-1; ++resid, ++pose_resid ) {
			TR.Debug << "Testing residue " << pose_resid << " vs " << resid << std::endl;
			TS_ASSERT_DELTA( pose.psi( pose_resid ), helix14.psi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.phi( pose_resid ), helix14.phi( resid ), 1e-6 );
			TS_ASSERT_DELTA( pose.omega( pose_resid ), helix14.omega( resid ), 1e-6 );
		}
	}

	void test_build_copy_simple() {
		core::pose::Pose const input_pose = create_trpcage_ideal_pose();

		PoseArchitect architect( "pose" );
		ExtendedPoseBuilder builder;

		StructureDataOP sd = architect.apply( input_pose );
		TS_ASSERT( sd );
		core::pose::PoseOP newpose = builder.apply( *sd );
		TS_ASSERT( newpose );
		TS_ASSERT( StructureDataFactory::get_instance()->observer_attached( *newpose ) );

		check_residue_dihedrals( *newpose, input_pose, 1e-4 );
		// poses should be equal
		TS_ASSERT( compare_atom_coordinates( *newpose, input_pose, 3 ) );
	}

	void test_build_copy_multi_chain()
	{
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/twohelix_structuredata.pdb" );

		PoseArchitect architect( "pose" );
		ExtendedPoseBuilder builder;

		StructureDataOP sd = architect.apply( input_pose );
		TS_ASSERT( sd );
		core::pose::PoseOP newpose = builder.apply( *sd );
		TS_ASSERT( newpose );

		TS_ASSERT_EQUALS( newpose->total_residue(), input_pose.total_residue() );
		for ( core::Size resid=1; resid<=newpose->total_residue(); ++resid ) {
			TR.Debug << "Checking dihedrals for residue " << resid << " vs "
				<< input_pose.phi( resid ) << " "
				<< input_pose.psi( resid ) << " "
				<< input_pose.omega( resid ) << std::endl;
			TS_ASSERT_DELTA( newpose->phi( resid ), input_pose.phi( resid ), 1e-5 );
			TS_ASSERT_DELTA( newpose->psi( resid ), input_pose.psi( resid ), 1e-5 );
			TS_ASSERT_DELTA( newpose->omega( resid ), input_pose.omega( resid ), 1e-5 );
		}
		// poses should be equal
		TS_ASSERT( compare_atom_coordinates( *newpose, input_pose, 3 ) );
	}

	void test_build_complicated() {
		using protocols::denovo_design::SegmentNameList;

		StructureDataFactory const & factory = *StructureDataFactory::get_instance();

		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/test_pdbcomp.pdb" );

		std::stringstream xml;
		xml << "<StructureData name=\"test_pdbcomp\" length=37 pose_length=47 >";
		xml << "<ResidueRange name=E1 start=2 safe=0 nterm=0 cterm=1 group=1 cutpoint=0 ss=LLEE abego=EBBB upper_segment=E1_2 />";
		xml << "<ResidueRange name=E1_2 start=5 safe=0 nterm=1 cterm=1 group=1 cutpoint=0 ss=EEE abego=BBB lower_segment=E1 upper_segment=cat1 />";
		xml << "<ResidueRange name=cat1 start=8 safe=0 nterm=1 cterm=0 group=2 cutpoint=0 ss=LLLLLLL abego=BBBAAAX lower_segment=E1_2 />";
		xml << "<ResidueRange name=E2 start=16 safe=0 nterm=0 cterm=0 group=1 cutpoint=0 ss=LEEEEEELL abego=EBBBBBBBX />";
		xml << "<ResidueRange name=E3 start=25 safe=0 nterm=0 cterm=0 group=1 cutpoint=0 ss=LEEEEEEL abego=EBBBBBBX />";
		xml << "<ResidueRange name=E4 start=33 safe=0 nterm=0 cterm=0 group=1 cutpoint=0 ss=LLEEEEL abego=EBBBBBX />";
		xml << "<ResidueRange name=cat2 start=40 safe=0 nterm=0 cterm=0 group=2 cutpoint=0 ss=LLLLLLLLL abego=EBAAGABBX />";
		xml << "</StructureData>";
		StructureDataOP sd = factory.create_from_xml( xml );
		TS_ASSERT( sd );
		sd->check_pose_consistency( input_pose );

		for ( SegmentNameList::const_iterator s=sd->segments_begin(); s!=sd->segments_end(); ++s ) {
			core::Size const start_resid = sd->segment( *s ).start();
			core::Size const stop_resid = sd->segment( *s ).stop();
			sd->set_template_pose( *s, input_pose, start_resid, stop_resid );
			for ( core::Size r=1; r<=sd->segment( *s ).template_pose()->total_residue(); ++r ) {
				TR.Debug << *s << " " << r << " " << sd->segment( *s ).template_pose()->residue( r ).name() << " " << input_pose.residue( start_resid + r - 1 ).name() << std::endl;
			}
		}

		ExtendedPoseBuilder builder;
		core::pose::PoseOP newpose = builder.apply( *sd );
		TS_ASSERT( newpose );

		TR.Debug << "Resid Orig New" << std::endl;
		for ( core::Size r=1; r<=newpose->total_residue(); ++r ) {
			TR.Debug << r << " " << input_pose.residue( r ).name() << " " << newpose->residue( r ).name() << std::endl;
		}

		TS_ASSERT_EQUALS( newpose->total_residue(), input_pose.total_residue() );
		TS_ASSERT_EQUALS( newpose->conformation().num_chains(), input_pose.conformation().num_chains() );
		TS_ASSERT_EQUALS( newpose->sequence(), input_pose.sequence() );
		check_residue_dihedrals( *newpose, input_pose, 1e-5 );
		TS_ASSERT( compare_atom_coordinates( *newpose, input_pose, 3 ) );
	}

	void test_bridgechains_add_template_segment_functions() {
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/test_pdbcomp_BridgeChains.pdb" );

		PoseArchitect architect( "" );
		StructureData sd = *architect.apply( input_pose );
		TS_ASSERT_THROWS_NOTHING( sd.check_pose_consistency( input_pose ) );

		TR << sd << std::endl;
		StructureData const orig = sd;
		sd.move_segment( "L03", "L06" );
		core::pose::Pose pose;

		// append segments one by one
		append_new_chain_from_template_segment( pose, sd.segment( "L01" ) );
		append_residues_from_template_segment( pose, sd.segment( "L01" ), sd.segment( "E01" ) );
		append_residues_from_template_segment( pose, sd.segment( "E01" ), sd.segment( "L02" ) );

		append_new_chain_from_template_segment( pose, sd.segment( "E02" ) );
		append_residues_from_template_segment( pose, sd.segment( "E02" ), sd.segment( "L03" ) );

		append_new_chain_from_template_segment( pose, sd.segment( "L06" ) );

		append_new_chain_from_template_segment( pose, sd.segment( "L04" ) );
		append_residues_from_template_segment( pose, sd.segment( "L04" ), sd.segment( "E03" ) );

		append_new_chain_from_template_segment( pose, sd.segment( "L05" ) );
		append_residues_from_template_segment( pose, sd.segment( "L05" ), sd.segment( "E04" ) );

		check_unwanted_movement( orig, input_pose, sd, pose );
	}

	void test_bridgechains_case() {
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/test_pdbcomp_BridgeChains.pdb" );

		PoseArchitect architect( "pose" );
		StructureData sd = *architect.apply( input_pose );
		TS_ASSERT_THROWS_NOTHING( sd.check_pose_consistency( input_pose ) );

		StructureData const orig = sd;
		TR << "Orig=" << sd << std::endl;

		Segment bridge_seg( "LLLLLL", "XXXXXX", false, false );
		bridge_seg.set_movable_group( 2 );
		bridge_seg.set_cutpoint( 2 );

		std::string const conn_id = "bridge";
		std::string const s1 = "pose.L03";
		std::string const s2 = "pose.L06";

		sd.add_segment( conn_id, bridge_seg );
		sd.move_segment( s1, conn_id );
		sd.move_segment( conn_id, s2 );
		sd.delete_trailing_residues( s1 );
		sd.delete_leading_residues( conn_id );
		sd.delete_trailing_residues( conn_id );
		sd.delete_leading_residues( s2 );
		sd.connect_segments( s1, conn_id );
		sd.connect_segments( conn_id, s2 );

		TR << "Conn=" << sd << std::endl;

		ExtendedPoseBuilder builder;
		core::pose::PoseOP newpose = builder.apply( sd );
		TS_ASSERT( newpose );
		sd.check_pose_consistency( *newpose );
		std::string new_seq = "";
		for ( SegmentNameList::const_iterator s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
			if ( ! orig.has_segment( *s ) ) new_seq += "VVVV";
			else if ( *s == s1 ) new_seq += input_pose.sequence().substr( orig.segment( s1 ).lower() - 1, orig.segment( s1 ).length() - 1 );
			else if ( *s == s2 ) new_seq += input_pose.sequence().substr( orig.segment( s2 ).start() - 1, orig.segment( s2 ).length() - 1 );
			else new_seq += input_pose.sequence().substr( orig.segment( *s ).lower() - 1, orig.segment( *s ).length() );
		}

		TS_ASSERT_EQUALS( newpose->total_residue(), 2 + input_pose.total_residue() );
		TS_ASSERT_EQUALS( newpose->conformation().num_chains() + 1, input_pose.conformation().num_chains() );
		TS_ASSERT_EQUALS( newpose->sequence(), new_seq );
		check_unwanted_movement( orig, input_pose, sd, *newpose );
	}

	void test_gerard_case() {
		// three sidechains need to be fixed relative to one another
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/components/test_fixed_trimer.pdb" );

		std::stringstream orig_xml;
		orig_xml << "<StructureData name=\"compound\" length=\"3\" pose_length=\"9\" >" << std::endl
			<< "<ResidueRange name=\"L01\" start=\"2\" stop=\"2\" safe=\"0\" nterm=\"0\" cterm=\"0\" group=\"1\" cutpoint=\"0\" ss=\"LLL\" abego=\"GAX\" />" << std::endl
			<< "<ResidueRange name=\"L02\" start=\"5\" stop=\"5\" safe=\"0\" nterm=\"0\" cterm=\"0\" group=\"1\" cutpoint=\"0\" ss=\"LLL\" abego=\"GAX\" />" << std::endl
			<< "<ResidueRange name=\"L03\" start=\"8\" stop=\"8\" safe=\"0\" nterm=\"0\" cterm=\"0\" group=\"1\" cutpoint=\"0\" ss=\"LLL\" abego=\"GAX\" />" << std::endl
			<< "</StructureData>" << std::endl;
		StructureData const orig = *StructureDataFactory::get_instance()->create_from_xml( orig_xml );
		orig.check_pose_consistency( input_pose );

		std::stringstream xml;
		xml << "<StructureData name=\"compound\" length=\"25\" pose_length=\"31\" >" << std::endl;
		xml << "<ResidueRange name=\"L02\" start=\"2\" stop=\"2\" safe=\"0\" nterm=\"0\" cterm=\"0\" group=\"1\" cutpoint=\"0\" ss=\"LLL\" abego=\"GAX\" />" << std::endl;
		xml << "<ResidueRange name=\"L03\" start=\"5\" stop=\"5\" safe=\"0\" nterm=\"0\" cterm=\"0\" group=\"1\" cutpoint=\"0\" ss=\"LLL\" abego=\"GAX\" />" << std::endl;
		xml << "<ResidueRange name=\"compound.h1\" start=\"8\" stop=\"17\" safe=\"0\" nterm=\"0\" cterm=\"1\" group=\"2\" cutpoint=\"0\" ss=\"LHHHHHHHHHH\" abego=\"XAAAAAAAAAA\" upper_segment=\"compound.staple1\" />" << std::endl;
		xml << "<ResidueRange name=\"compound.staple1\" start=\"18\" stop=\"18\" safe=\"0\" nterm=\"1\" cterm=\"1\" group=\"4\" cutpoint=\"0\" ss=\"L\" abego=\"X\" lower_segment=\"compound.h1\" upper_segment=\"L01\" />" << std::endl;
		xml << "<ResidueRange name=\"L01\" start=\"19\" stop=\"19\" safe=\"0\" nterm=\"1\" cterm=\"1\" group=\"1\" cutpoint=\"0\" ss=\"L\" abego=\"A\" lower_segment=\"compound.staple1\" upper_segment=\"compound.staple2\" />" << std::endl;
		xml << "<ResidueRange name=\"compound.staple2\" start=\"20\" stop=\"20\" safe=\"0\" nterm=\"1\" cterm=\"1\" group=\"5\" cutpoint=\"0\" ss=\"L\" abego=\"X\" lower_segment=\"L01\" upper_segment=\"compound.h2\" />" << std::endl;
		xml << "<ResidueRange name=\"compound.h2\" start=\"21\" stop=\"30\" safe=\"0\" nterm=\"1\" cterm=\"0\" group=\"3\" cutpoint=\"0\" ss=\"HHHHHHHHHHL\" abego=\"AAAAAAAAAAX\" lower_segment=\"compound.staple2\" />" << std::endl;
		xml << "<Str name=\"compound.staple1#segment1\" value=\"compound.h1\" />" << std::endl;
		xml << "<Str name=\"compound.staple1#segment2\" value=\"L01\" />" << std::endl;
		xml << "<Str name=\"compound.staple2#segment1\" value=\"L01\" />" << std::endl;
		xml << "<Str name=\"compound.staple2#segment2\" value=\"compound.h2\" />" << std::endl;
		xml << "</StructureData>";

		// Create SD and add template segments
		StructureData sd = *StructureDataFactory::get_instance()->create_from_xml( xml );
		StructureData pose_sd = *StructureDataFactory::get_instance()->create_from_pose( input_pose );
		for ( SegmentNameList::const_iterator s=pose_sd.segments_begin(); s!=pose_sd.segments_end(); ++s ) {
			sd.set_template_pose( *s, input_pose, pose_sd.segment( *s ).start(), pose_sd.segment( *s ).stop() );
		}

		ExtendedPoseBuilder builder;
		core::pose::PoseOP pose_ptr = builder.apply( sd );
		TS_ASSERT( pose_ptr );
		core::pose::Pose const & pose = *pose_ptr;
		check_unwanted_movement( orig, input_pose, sd, pose );
	}
};

void
check_residue_dihedrals( core::pose::Pose const & pose, core::pose::Pose const & orig, core::Real const tol )
{
	TS_ASSERT_EQUALS( pose.total_residue(), orig.total_residue() );
	for ( core::Size resid=1; resid<=pose.total_residue(); ++resid ) {
		TR.Debug << resid << " phi1 " << orig.phi( resid ) << " phi2 " << pose.phi( resid ) << std::endl;
		TR.Debug << resid << " psi1 " << orig.psi( resid ) << " psi2 " << pose.psi( resid ) << std::endl;
		TR.Debug << resid << " omega1 " << orig.omega( resid ) << " omega2 " << pose.omega( resid ) << std::endl;
		TS_ASSERT_EQUALS( orig.residue( resid ).name(), pose.residue( resid ).name() );
		TS_ASSERT_DELTA( orig.phi( resid ), pose.phi( resid ), tol );
		TS_ASSERT_DELTA( orig.psi( resid ), pose.psi( resid ), tol );
		TS_ASSERT_DELTA( orig.omega( resid ), pose.omega( resid ), tol );
	}
}

// this is a copy of the function in core::pose::util.hh, except it uses distances between equivalent atoms
// instead of std::floor() -- this is better for floating point errors, where std::floor( 191.999999 ) is not
// the same as std::floor( 1.920000 ), even though those numbers are essentially equivalent
bool
compare_atom_coordinates(
	core::pose::Pose const & lhs,
	core::pose::Pose const & rhs,
	core::Size const n_dec_places )
{
	//number of decimal places of precision - 3 (1000) is equivalent to PDB precision.
	core::Real const n_dec( pow( static_cast< core::Real >(10), static_cast< int >( n_dec_places ) ) );

	//first compare pose sizes; prerequisite to iterating through length
	core::Size const lhssize(lhs.total_residue()), rhssize(rhs.total_residue());

	if ( lhssize != rhssize ) {
		TR.Warning << "poses of different length in compare_atom_coordinates; doomed to fail!" << std::endl;
		return false;
	}

	//now iterate through residues and make comparisons
	for ( core::Size resid(1); resid<=lhssize; ++resid ) {
		//check equality of residue types
		if ( lhs.residue( resid ).is_terminus() ) {
			if ( !rhs.residue( resid ).is_terminus() ) return false;
			continue;
		}
		core::chemical::ResidueType const & lhstype(lhs.residue_type(resid)), & rhstype(rhs.residue_type(resid));
		if ( lhstype.name() != rhstype.name() ) { //string matching is sufficient because ResidueType objects have unique names
			TR.Warning << "nonmatching ResidueTypes at " << resid << " in compare_atom_coordinates -- skipping" << std::endl;
			TR.Warning << lhstype.name() << " vs " << rhstype.name() << std::endl;
			return false;
		}

		//get atoms vectors to compare
		core::conformation::Atoms const & lhsatoms(lhs.residue(resid).atoms()), rhsatoms(rhs.residue(resid).atoms());
		core::Size const lhsatmsize(lhsatoms.size()), rhsatmsize(rhsatoms.size());
		if ( lhsatmsize != rhsatmsize ) { //check vector length equality
			TR.Warning << "nonmatching numbers of atoms at residue " << resid << " in compare_atom_coordinates" << std::endl;
			TR.Warning << "How did we even get here?  ResidueType comparison should have failed!" << std::endl;
			return false;
		}

		//iterate through atoms vector
		for ( core::Size atm(1); atm <= lhsatmsize; ++atm ) {
			core::Real const xdist = n_dec * std::abs( lhsatoms[atm].xyz().x() - rhsatoms[atm].xyz().x() );
			core::Real const ydist = n_dec * std::abs( lhsatoms[atm].xyz().y() - rhsatoms[atm].xyz().y() );
			core::Real const zdist = n_dec * std::abs( lhsatoms[atm].xyz().z() - rhsatoms[atm].xyz().z() );
			if ( ( xdist >= 1.0 ) || ( ydist >= 1.0 ) || ( zdist >= 1.0 ) ) {
				TR.Debug << "Nonmatching coordinates at residue " << resid << " atom " << atm << " in compare_atom_coordinates" << std::endl;
				TR.Debug << "Raw coords: " << lhsatoms[atm].xyz().x() << " " << lhsatoms[atm].xyz().y() << " " << lhsatoms[atm].xyz().z()
					<< " vs "
					<< rhsatoms[atm].xyz().x() << " " << rhsatoms[atm].xyz().y() << " " << rhsatoms[atm].xyz().z() << std::endl;
				return false; //no warning messages, this is the "expected" failure
			}
		}//iterate through atoms vector
	}//over all residues

	return true; //whoo! we made it!
}//compare_atom_coordinates

