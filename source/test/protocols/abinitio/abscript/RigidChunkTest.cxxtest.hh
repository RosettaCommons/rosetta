// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <core/id/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/pack/task/residue_selector/ResidueIndexSelector.hh>

#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvExcn.hh>

#include <protocols/filters/Filter.hh>

#include <protocols/abinitio/abscript/RigidChunkCM.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <test/core/init_util.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <utility/tag/Tag.hh>
// #include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <iostream>

using namespace core;
using namespace protocols::environment;

std::string const TEST_CHUNK_PDB = "protocols/abinitio/abscript/ubq_frag.pdb";
core::Real const TOLERANCE = 1e-6;

class RigidChunkTest : public CxxTest::TestSuite {
public:

	// Shared data elements go here.
	core::pose::Pose pose;

	core::pose::Pose make_seq( std::string const& str ) const {
		core::pose::Pose pose;

		core::pose::make_pose_from_sequence(pose, str, "fa_standard");

		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			pose.set_phi( i, -65 );
			pose.set_psi( i, -41 );
			pose.set_omega( i, 180 );
		}

		return pose;
	}

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();

		pose = make_seq( "FRMQIFVYTLTGNDSS" );

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void ANGLE_DELTA( Real x, Real y, Real d ) const {
		TS_ASSERT_DELTA( std::cos( x ) , std::cos( y ), d );
	}

	void test_chunk() {
		TS_TRACE( "Starting test_chunk" );

		using namespace protocols::environment;
		using namespace core::pack::task::residue_selector;
		using namespace protocols::abinitio::abscript;

		basic::datacache::DataMap datamap;

		pose = make_seq( "MQIFV" );

		std::string const SELECTOR_NAME = "sele";
		core::Size const BEGIN_CHUNK = 1;
		core::Size const END_CHUNK = 5;

		std::ostringstream os;
		os << BEGIN_CHUNK << "-" << END_CHUNK;
		datamap.add( "ResidueSelector", SELECTOR_NAME, ResidueSelectorOP( new ResidueIndexSelector( os.str() ) ) );

		std::stringstream ss;
		ss << "<RigidChunkCM name=chunk template=\"" << TEST_CHUNK_PDB << "\" selector=\"" << SELECTOR_NAME << "\" />";
		utility::tag::TagPtr tag( new utility::tag::Tag );
		tag->read( ss );

		RigidChunkCMOP rigid_chunk( new RigidChunkCM() );
		TS_ASSERT_THROWS_NOTHING( rigid_chunk->parse_my_tag( tag , datamap, protocols::filters::Filters_map(), protocols::moves::Movers_map(), core::pose::Pose() ) );
		TS_ASSERT_THROWS_NOTHING( tag->die_for_unaccessed_options() );

		EnvironmentOP env( new Environment( "env" ) );
		env->register_mover( rigid_chunk );

		core::pose::Pose ppose;
		TS_ASSERT_THROWS_NOTHING( ppose = env->start( pose ) );

		for ( core::Size i = 1; i <= ppose.total_residue(); ++i ) {
			if ( i >= BEGIN_CHUNK && i <= END_CHUNK ) {
				core::Size const seqpos = i - BEGIN_CHUNK + 1;
				if ( seqpos != 1 ) {
					ANGLE_DELTA( ppose.phi( i ), rigid_chunk->templ().phi( seqpos ), TOLERANCE );
				}
				if ( seqpos != rigid_chunk->templ().total_residue() &&
						seqpos != ppose.total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), rigid_chunk->templ().psi( seqpos ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ) , rigid_chunk->templ().omega( seqpos ), TOLERANCE );
				}
			} else if ( ( i == BEGIN_CHUNK - 1 ) ||
					( i == END_CHUNK + 1 ) ) {
				// inconsistent contributions.
			} else {
				if ( i != 1 ) ANGLE_DELTA( ppose.phi( i ), pose.phi( i ), TOLERANCE );
				if ( i != pose.total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), pose.psi( i ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ), pose.omega( i ), TOLERANCE );
				}
			}
		}

		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env->end( ppose ) );
	}

	void test_two_chunk() {
		TS_TRACE( "Starting test_two_chunk" );

		using namespace protocols::environment;
		using namespace core::pack::task::residue_selector;
		using namespace protocols::abinitio::abscript;

		// We need an extra long sequence in this one.
		pose = make_seq( "FRMQIFVYTLTGMQIFVSS" );

		basic::datacache::DataMap datamap;

		std::string const SELECTOR_NAME_1 = "sele1";
		core::Size const BEGIN_CHUNK_1 = 3;
		core::Size const END_CHUNK_1 = 7;

		std::ostringstream os_1;
		os_1 << BEGIN_CHUNK_1 << "-" << END_CHUNK_1;
		datamap.add( "ResidueSelector", SELECTOR_NAME_1, ResidueSelectorOP( new ResidueIndexSelector( os_1.str() ) ) );
		std::stringstream ss_1;
		ss_1 << "<RigidChunkCM name=chunk template=\"" << TEST_CHUNK_PDB << "\" selector=\"" << SELECTOR_NAME_1 << "\" />";
		utility::tag::TagPtr tag( new utility::tag::Tag );
		tag->read( ss_1 );

		RigidChunkCMOP rigid_chunk_1( new RigidChunkCM() );
		TS_ASSERT_THROWS_NOTHING( rigid_chunk_1->parse_my_tag( tag , datamap, protocols::filters::Filters_map(), protocols::moves::Movers_map(), core::pose::Pose() ) );
		TS_ASSERT_THROWS_NOTHING( tag->die_for_unaccessed_options() );

		std::string const SELECTOR_NAME_2 = "sele2";
		core::Size const BEGIN_CHUNK_2 = BEGIN_CHUNK_1 + END_CHUNK_1 + 3;
		core::Size const END_CHUNK_2 = BEGIN_CHUNK_2 + ( END_CHUNK_1 - BEGIN_CHUNK_1 );

		std::ostringstream os_2;
		os_2 << BEGIN_CHUNK_2 << "-" << END_CHUNK_2;
		datamap.add( "ResidueSelector", SELECTOR_NAME_2, ResidueSelectorOP( new ResidueIndexSelector( os_2.str() ) ) );
		std::stringstream ss_2;
		ss_2 << "<RigidChunkCM name=chunk template=\"" << TEST_CHUNK_PDB << "\" selector=\"" << SELECTOR_NAME_2 << "\" />";
		tag = utility::tag::TagPtr( new utility::tag::Tag() );
		tag->read( ss_2 );

		RigidChunkCMOP rigid_chunk_2( new RigidChunkCM() );
		TS_ASSERT_THROWS_NOTHING( rigid_chunk_2->parse_my_tag( tag , datamap, protocols::filters::Filters_map(), protocols::moves::Movers_map(), core::pose::Pose() ) );
		TS_ASSERT_THROWS_NOTHING( tag->die_for_unaccessed_options() );


		EnvironmentOP env( new Environment( "env" ) );
		env->register_mover( rigid_chunk_1 );
		env->register_mover( rigid_chunk_2 );

		core::pose::Pose ppose;
		TS_ASSERT_THROWS_NOTHING( ppose = env->start( pose ) );

		for ( core::Size i = 1; i <= ppose.total_residue(); ++i ) {
			if ( i >= BEGIN_CHUNK_1 && i <= END_CHUNK_1 ) {
				core::Size const seqpos = i - BEGIN_CHUNK_1 + 1;
				if ( seqpos != 1 ) {
					ANGLE_DELTA( ppose.phi( i ), rigid_chunk_1->templ().phi( seqpos ), TOLERANCE );
				}
				if ( seqpos != rigid_chunk_1->templ().total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), rigid_chunk_1->templ().psi( seqpos ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ) , rigid_chunk_1->templ().omega( seqpos ), TOLERANCE );
				}
			} else if ( ( i == BEGIN_CHUNK_1 - 1 ) ||
					( i == END_CHUNK_1 + 1 ) ||
					( i == BEGIN_CHUNK_2 - 1 ) ||
					( i == END_CHUNK_2 + 1 ) ) {
				// inconsistent contributions.
			} else if ( i >= BEGIN_CHUNK_2 && i <= END_CHUNK_2 ) {
				core::Size const seqpos = i - BEGIN_CHUNK_2 + 1;
				if ( seqpos != 1 ) {
					ANGLE_DELTA( ppose.phi( i ), rigid_chunk_1->templ().phi( seqpos ), TOLERANCE );
				}
				if ( seqpos != rigid_chunk_1->templ().total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), rigid_chunk_1->templ().psi( seqpos ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ) , rigid_chunk_1->templ().omega( seqpos ), TOLERANCE );
				}
			} else {
				if ( i != 1 ) ANGLE_DELTA( ppose.phi( i ), pose.phi( i ), TOLERANCE );
				if ( i != pose.total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), pose.psi( i ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ), pose.omega( i ), TOLERANCE );
				}
			}
		}

		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env->end( ppose ) );
	}

	void test_apply_to_template() {
		TS_TRACE( "Starting test_apply_to_template" );

		using namespace protocols;
		using namespace protocols::environment;
		using namespace core::pack::task::residue_selector;
		using namespace protocols::abinitio::abscript;

		basic::datacache::DataMap datamap;
		moves::Movers_map movers_map;

		std::string const SELECTOR_NAME = "sele";
		std::string const APPLY_MOVER = "centroid";

		movers_map[ APPLY_MOVER ] = moves::MoverOP( new simple_moves::SwitchResidueTypeSetMover( "centroid" ) );

		core::Size const BEGIN_CHUNK = 3;
		core::Size const END_CHUNK = 7;
		std::ostringstream os;
		os << BEGIN_CHUNK << "-" << END_CHUNK;
		datamap.add( "ResidueSelector", SELECTOR_NAME, ResidueSelectorOP( new ResidueIndexSelector( os.str() ) ) );

		std::stringstream ss;
		ss << "<RigidChunkCM name=\"name"
			<< "\" template=\"" << TEST_CHUNK_PDB
			<< "\" selector=\"" << SELECTOR_NAME
			<< "\" apply_to_template=\"" << APPLY_MOVER
			<< "\" />";
		utility::tag::TagPtr tag( new utility::tag::Tag );
		tag->read( ss );

		RigidChunkCMOP rigid_chunk( new RigidChunkCM() );
		TS_ASSERT_THROWS_NOTHING( rigid_chunk->parse_my_tag( tag , datamap, filters::Filters_map(), movers_map, core::pose::Pose() ) );
		TS_ASSERT_THROWS_NOTHING( tag->die_for_unaccessed_options() );

		EnvironmentOP env( new Environment( "env" ) );
		env->register_mover( rigid_chunk );

		core::pose::Pose ppose;
		TS_ASSERT_THROWS( ppose = env->start( pose ), utility::excn::EXCN_BadInput );

		simple_moves::SwitchResidueTypeSetMover( "centroid" ).apply( pose );
		TS_ASSERT_THROWS_NOTHING( ppose = env->start( pose ) );

		for ( core::Size i = 1; i <= ppose.total_residue(); ++i ) {
			if ( i >= BEGIN_CHUNK && i <= END_CHUNK ) {
				core::Size const seqpos = i - BEGIN_CHUNK + 1;
				if ( seqpos != 1 ) {
					ANGLE_DELTA( ppose.phi( i ), rigid_chunk->templ().phi( seqpos ), TOLERANCE );
				}
				if ( seqpos != rigid_chunk->templ().total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), rigid_chunk->templ().psi( seqpos ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ) , rigid_chunk->templ().omega( seqpos ), TOLERANCE );
				}
			} else if ( ( i == BEGIN_CHUNK - 1 ) ||
					( i == END_CHUNK + 1 ) ) {
				// inconsistent contributions.
			} else {
				if ( i != 1 ) ANGLE_DELTA( ppose.phi( i ), pose.phi( i ), TOLERANCE );
				if ( i != pose.total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), pose.psi( i ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ), pose.omega( i ), TOLERANCE );
				}
			}
		}

		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env->end( ppose ) );
	}

	void test_two_part_chunk() {
		TS_TRACE( "Starting test_two_part_chunk" );

		using namespace protocols::environment;
		using namespace core::pack::task::residue_selector;
		using namespace protocols::abinitio::abscript;

		pose = make_seq( "FRMQIFVYKTLTGKSSSS" );

		basic::datacache::DataMap datamap;

		std::string const SELECTOR_NAME = "sele";

		core::Size const BEGIN_CHUNK_1 = 3;
		core::Size const END_CHUNK_1 = 7;
		core::Size const BEGIN_CHUNK_2 = 9;
		core::Size const END_CHUNK_2 = 12;
		std::ostringstream os;
		os << BEGIN_CHUNK_1 << "-" << END_CHUNK_1 << "," << BEGIN_CHUNK_2 << "-" << END_CHUNK_2;
		datamap.add( "ResidueSelector", SELECTOR_NAME, ResidueSelectorOP( new ResidueIndexSelector( os.str() ) ) );

		std::stringstream ss;
		ss << "<RigidChunkCM name=chunk template=\"" << TEST_CHUNK_PDB << "\" selector=\"" << SELECTOR_NAME << "\" />";
		utility::tag::TagPtr tag( new utility::tag::Tag );
		tag->read( ss );

		RigidChunkCMOP rigid_chunk( new RigidChunkCM() );
		TS_ASSERT_THROWS_NOTHING( rigid_chunk->parse_my_tag( tag , datamap, protocols::filters::Filters_map(), protocols::moves::Movers_map(), core::pose::Pose() ) );
		TS_ASSERT_THROWS_NOTHING( tag->die_for_unaccessed_options() );

		EnvironmentOP env( new Environment( "env" ) );
		env->auto_cut( true );
		env->register_mover( rigid_chunk );

		core::pose::Pose ppose;
		TS_ASSERT_THROWS_NOTHING( ppose = env->start( pose ) );

		assert( std::abs( ppose.psi( 9 ) - rigid_chunk->templ().psi( 6 ) ) < 1e-6 );
		ANGLE_DELTA( ppose.psi( 9 ), rigid_chunk->templ().psi( 6 ), TOLERANCE );

		TS_ASSERT_EQUALS( ppose.fold_tree().num_cutpoint(), 1 );

		for ( core::Size i = 1; i <= ppose.total_residue(); ++i ) {
			if ( ( i >= BEGIN_CHUNK_1 && i <= END_CHUNK_1 ) ||
					( i >= BEGIN_CHUNK_2 && i <= END_CHUNK_2 ) ) {
				core::Size const seqpos = rigid_chunk->sim_origin().at( i );
				if ( seqpos != 1 && !ppose.fold_tree().is_cutpoint( i - 1 ) ) {
					ANGLE_DELTA( ppose.phi( i ), rigid_chunk->templ().phi( seqpos ), TOLERANCE );
				}
				if ( seqpos != rigid_chunk->templ().total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), rigid_chunk->templ().psi( seqpos ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ) , rigid_chunk->templ().omega( seqpos ), TOLERANCE );
				}
			} else if ( i == ( BEGIN_CHUNK_1 - 1 ) ||
					i == ( END_CHUNK_1 + 1 ) ||
					i == ( BEGIN_CHUNK_2 - 1 ) ||
					i == ( END_CHUNK_2 + 1 ) ) {
				// inconsistent contributions.
			} else {
				if ( i != 1 ) ANGLE_DELTA( ppose.phi( i ), pose.phi( i ), TOLERANCE );
				if ( i != pose.total_residue() ) {
					ANGLE_DELTA( ppose.psi( i ), pose.psi( i ), TOLERANCE );
					ANGLE_DELTA( ppose.omega( i ), pose.omega( i ), TOLERANCE );
				}
			}
		}

		// Compute the distance between two points in the rigid chunk and verify they're the same in
		// simulation and template poses
		core::Size const templ_pos_1 = rigid_chunk->sim_origin().at( BEGIN_CHUNK_1 );
		core::Size const templ_pos_2 = rigid_chunk->sim_origin().at( BEGIN_CHUNK_2 );

		numeric::xyzVector< core::Length > sim_vect = ppose.residue( BEGIN_CHUNK_1 ).xyz( 1 ) - ppose.residue( BEGIN_CHUNK_2 ).xyz( 1 );
		numeric::xyzVector< core::Length > templ_vect = rigid_chunk->templ().residue( templ_pos_1 ).xyz( 1 ) - rigid_chunk->templ().residue( templ_pos_2 ).xyz( 1 );

		std::cout << sim_vect.length() << "," << templ_vect.length() << std::endl;
		TS_ASSERT_DELTA( sim_vect.length(), templ_vect.length(), TOLERANCE );

		// Compute the above for a different pair of points, this guarantees similarity within a coordinate transform.
		sim_vect = ppose.residue( BEGIN_CHUNK_1 ).xyz( 2 ) - ppose.residue( BEGIN_CHUNK_2 ).xyz( 2 );
		templ_vect = rigid_chunk->templ().residue( templ_pos_1 ).xyz( 2 ) - rigid_chunk->templ().residue( templ_pos_2 ).xyz( 2 );
		TS_ASSERT_DELTA( sim_vect.length(), templ_vect.length(), TOLERANCE );

		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env->end( ppose ) );
	}
};
