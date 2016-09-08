// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_daemon/PackDaemon.cxxtest.hh
/// @brief  test suite for PackDaemon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/pack_daemon/EntityCorrespondence.hh>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Core headers

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <sstream>

//Auto Headers
#include <utility/vector1.hh>


//static basic::Tracer TR("EntityCorrespondenceTest.cxxtest");

using namespace core;
using namespace protocols::pack_daemon;

class EntityCorrespondenceTests : public CxxTest::TestSuite
{
public:
	typedef core::pose::PoseOP PoseOP;

private:
	bool initialized_;

public:
	EntityCorrespondenceTests() : initialized_( false ) {}

	void setUp() {
		if ( ! initialized_ ) {
			core_init();
			initialized_ = true;
		}
	}

	EntityCorrespondence create_default_correspondence()
	{
		PoseOP trpcage = create_trpcage_ideal_poseop();

		EntityCorrespondence corr;
		corr.set_pose( trpcage );
		corr.set_num_entities( 2 );
		corr.add_resid_to_entity_list( 1, 10 );
		corr.add_resid_to_entity_list( 1, 11 );
		corr.add_resid_to_entity_list( 2, 12 );
		return corr;
	}

	void test_entity_correspondence_ctor() {
		EntityCorrespondence corr( create_default_correspondence() );

		TS_ASSERT( corr.num_entities() == 2 );
		TS_ASSERT( corr.num_residues() == 20 );
		TS_ASSERT( corr.n_residues_for_entity( 1 ) == 2 );
		TS_ASSERT( corr.n_residues_for_entity( 2 ) == 1 );

		EntityCorrespondence::ResIDListConstIter iter, iter_end;

		iter     = corr.residues_for_entity_begin( 1 );
		iter_end = corr.residues_for_entity_end(   1 );

		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 10 );
		++iter;
		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 11 );
		++iter;
		TS_ASSERT( iter == iter_end );

		iter     = corr.residues_for_entity_begin( 2 );
		iter_end = corr.residues_for_entity_end(   2 );

		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 12 );

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii >= 10 && ii <= 12 ) continue;
			TS_ASSERT( corr.entity_for_residue( ii ) == 0 );
		}
		TS_ASSERT( corr.entity_for_residue( 10 ) == 1 );
		TS_ASSERT( corr.entity_for_residue( 11 ) == 1 );
		TS_ASSERT( corr.entity_for_residue( 12 ) == 2 );
	}


	void test_correspondence_entity_index_out_of_bounds() {

		EntityCorrespondence corr( create_default_correspondence() );
		try {
			corr.n_residues_for_entity( 3 );
			TS_ASSERT( false ); // should never reach here.
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			//std::cout << excn.msg() << std::endl;
			TS_ASSERT( excn.msg() == "EntityCorrespondence::n_residues_for_entity attempted to access information for entity 3 (only 2 entities exist)" );
		}
	}

	void test_correspondence_residue_index_out_of_bounds() {

		EntityCorrespondence corr( create_default_correspondence() );
		try {
			corr.entity_for_residue( 21 );
			TS_ASSERT( false ); // should never reach here.
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			//std::cout << excn.msg() << std::endl;
			TS_ASSERT( excn.msg() == "EntityCorrespondence::entity_for_residue attempted to access information for residue 21 (only 20 residues exist)" );
		}
	}

	void test_correspondence_from_instream() {

		EntityCorrespondence corr;
		PoseOP trpcage = create_trpcage_ideal_poseop();
		corr.set_pose( trpcage );
		corr.set_num_entities( 2 );

		std::string corr_file( "1 10 A\n1 11 A\n2 12 A\n" );
		std::istringstream corr_stream( corr_file );


		try {
			corr.initialize_from_correspondence_file( corr_stream );
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			std::cout << excn.msg() << std::endl;
			TS_ASSERT( false );
		}


		EntityCorrespondence::ResIDListConstIter iter, iter_end;
		iter     = corr.residues_for_entity_begin( 1 );
		iter_end = corr.residues_for_entity_end(   1 );

		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 10 );
		++iter;
		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 11 );
		++iter;
		TS_ASSERT( iter == iter_end );

		iter     = corr.residues_for_entity_begin( 2 );
		iter_end = corr.residues_for_entity_end(   2 );

		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 12 );

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii >= 10 && ii <= 12 ) continue;
			TS_ASSERT( corr.entity_for_residue( ii ) == 0 );
		}
		TS_ASSERT( corr.entity_for_residue( 10 ) == 1 );
		TS_ASSERT( corr.entity_for_residue( 11 ) == 1 );
		TS_ASSERT( corr.entity_for_residue( 12 ) == 2 );
	}

	void test_correspondence_from_instream_ignoring_comments() {

		EntityCorrespondence corr;
		PoseOP trpcage = create_trpcage_ideal_poseop();
		corr.set_pose( trpcage );
		corr.set_num_entities( 2 );

		std::string corr_file( "1 10 A anything after the first two numbers and chain character on a line\n1 11 A is ignored\n2 12 A\n" );
		std::istringstream corr_stream( corr_file );


		try {
			corr.initialize_from_correspondence_file( corr_stream );
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			std::cout << excn.msg() << std::endl;
			TS_ASSERT( false );
		}


		EntityCorrespondence::ResIDListConstIter iter, iter_end;
		iter     = corr.residues_for_entity_begin( 1 );
		iter_end = corr.residues_for_entity_end(   1 );

		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 10 );
		++iter;
		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 11 );
		++iter;
		TS_ASSERT( iter == iter_end );

		iter     = corr.residues_for_entity_begin( 2 );
		iter_end = corr.residues_for_entity_end(   2 );

		TS_ASSERT( iter != iter_end );
		if ( iter == iter_end ) return;
		TS_ASSERT( *iter == 12 );

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii >= 10 && ii <= 12 ) continue;
			TS_ASSERT( corr.entity_for_residue( ii ) == 0 );
		}
		TS_ASSERT( corr.entity_for_residue( 10 ) == 1 );
		TS_ASSERT( corr.entity_for_residue( 11 ) == 1 );
		TS_ASSERT( corr.entity_for_residue( 12 ) == 2 );
	}

	void test_correspondence_bad_instream_entity_type() {

		EntityCorrespondence corr;
		PoseOP trpcage = create_trpcage_ideal_poseop();
		corr.set_pose( trpcage );
		corr.set_num_entities( 2 );

		std::string corr_file( "one 10\n1 11\n2 12\n" );
		std::istringstream corr_stream( corr_file );


		try {
			corr.initialize_from_correspondence_file( corr_stream );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			//std::cout << excn.msg() << std::endl;
			TS_ASSERT( excn.msg() == "Failed to read entity id on line 1 of the EntityCorrespondence file:\none 10" );
		}
	}

	void test_correspondence_bad_instream_resid_type() {

		EntityCorrespondence corr;
		PoseOP trpcage = create_trpcage_ideal_poseop();
		corr.set_pose( trpcage );
		corr.set_num_entities( 2 );

		std::string corr_file( "1 ten A\n1 11 A\n2 12 A\n" );
		std::istringstream corr_stream( corr_file );


		try {
			corr.initialize_from_correspondence_file( corr_stream );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			//std::cout << excn.msg() << std::endl;
			//TS_ASSERT( excn.msg() == "Failed to read residue id on line 1 of the EntityCorrespondence file:\n1 ten" );
			TS_ASSERT( excn.msg() == "Failed to read PDB residue id on line 1 of the EntityCorrespondence file:\n1 ten A\nCharacter1 of 'ten' is not a digit");
		}
	}

	void test_correspondence_bad_instream_entity_range() {

		EntityCorrespondence corr;
		PoseOP trpcage = create_trpcage_ideal_poseop();
		corr.set_pose( trpcage );
		corr.set_num_entities( 2 );

		std::string corr_file( "3 10 A\n1 11 A\n2 12 A\n" );
		std::istringstream corr_stream( corr_file );


		try {
			corr.initialize_from_correspondence_file( corr_stream );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			//std::cout << excn.msg() << std::endl;
			TS_ASSERT( excn.msg() == "Entity ID read on line 1 of the EntityCorrespondence file exceeds the number of entities: 3 vs 2" );
		}
	}

	void test_correspondence_bad_instream_residue_range() {

		EntityCorrespondence corr;
		PoseOP trpcage = create_trpcage_ideal_poseop();
		corr.set_pose( trpcage );
		corr.set_num_entities( 2 );

		std::string corr_file( "1 21 A\n1 11 A\n2 12 A\n" );
		std::istringstream corr_stream( corr_file );


		try {
			corr.initialize_from_correspondence_file( corr_stream );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & excn ) {
			//std::cout << excn.msg() << std::endl;
			//TS_ASSERT( excn.msg() == "Residue ID read on line 1 of the EntityCorrespondence file exceeds the number of residues: 21 vs 20" );
			TS_ASSERT( excn.msg() == "Residue ID read on line 1 of the EntityCorrespondence file is not present in the pose: 21 ch: A vs pose.size()= 20");
		}
	}


};


