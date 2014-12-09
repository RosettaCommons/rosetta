// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/MultistateFitnessFunction.cxxtest.hh
/// @brief  test suite for MultistateFitnessFunction
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/pack_daemon/MultistateFitnessFunction.hh>

#include <test/core/init_util.hh>

// Core headers
#include <core/types.hh>

// Protocols headers
#include <protocols/genetic_algorithm/Entity.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <list>

//static basic::Tracer TR("MultistateFitnessFunctionTest.cxxtest");

using namespace protocols::genetic_algorithm;
using namespace protocols::pack_daemon;

class MultistateFitnessFunctionTest : public CxxTest::TestSuite
{
public:
	void setUp() {
		core_init();
	}

	void test_top_entity_set_keep_four()
	{

		Entity ent1( "traits AA:1:A fitness 4.0" );
		Entity ent2( "traits AA:1:C fitness 3.0" );
		Entity ent3( "traits AA:1:D fitness 2.0" );
		Entity ent4( "traits AA:1:E fitness 1.0" );
		Entity ent5( "traits AA:1:F fitness 0.0" );

		TopEntitySet::StateEnergiesAndNPDs seanpds;
		bool was_added( false );
		std::list< EntityOP > drop_list;

		{
			TopEntitySet tes;
			tes.desired_entity_history_size( 4 );
			drop_list = tes.update_entity_history( ent1, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent2, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent3, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent4, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );

			drop_list = tes.update_entity_history( ent5, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( ! drop_list.empty() );
			TS_ASSERT( drop_list.size() == 1 );
			EntityOP dropped = *drop_list.begin();
			TS_ASSERT( dropped->fitness() == 4.0 );
		}

		{
			TopEntitySet tes;
			tes.desired_entity_history_size( 4 );
			drop_list = tes.update_entity_history( ent5, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent3, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent1, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent4, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );

			drop_list = tes.update_entity_history( ent2, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( ! drop_list.empty() );
			TS_ASSERT( drop_list.size() == 1 );
			EntityOP dropped = *drop_list.begin();
			TS_ASSERT( dropped->fitness() == 4.0 );
		}

		{
			TopEntitySet tes;
			tes.desired_entity_history_size( 4 );
			drop_list = tes.update_entity_history( ent4, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent2, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent3, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent5, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );

			drop_list = tes.update_entity_history( ent1, seanpds, was_added );
			TS_ASSERT( ! was_added );
			TS_ASSERT( drop_list.empty() );
		}


	}

	void test_top_entity_set_handle_ties()
	{

		Entity ent1( "traits AA:1:A fitness 4.0" );
		Entity ent2( "traits AA:1:C fitness 4.0" );
		Entity ent3( "traits AA:1:D fitness 4.0" );
		Entity ent4( "traits AA:1:E fitness 1.0" );
		Entity ent5( "traits AA:1:F fitness 0.0" );
		Entity ent6( "traits AA:1:G fitness -1.0" );
		Entity ent7( "traits AA:1:H fitness 5.0" );

		TopEntitySet::StateEnergiesAndNPDs seanpds;
		bool was_added( false );
		std::list< EntityOP > drop_list;

		{
			TopEntitySet tes;
			tes.desired_entity_history_size( 2 );
			drop_list = tes.update_entity_history( ent1, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent2, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent3, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent4, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );

			drop_list = tes.update_entity_history( ent5, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( ! drop_list.empty() );
			TS_ASSERT( drop_list.size() == 3 );
			for ( std::list< EntityOP >::const_iterator iter = drop_list.begin(); iter != drop_list.end(); ++iter ) {
				EntityOP dropped = *iter;
				TS_ASSERT( dropped->fitness() == 4.0 );
			}
		}

		{
			TopEntitySet tes;
			tes.desired_entity_history_size( 2 );
			drop_list = tes.update_entity_history( ent7, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent1, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent2, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( ! drop_list.empty() );
			TS_ASSERT( (*drop_list.begin())->fitness() == 5.0 );

			drop_list = tes.update_entity_history( ent3, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );

			drop_list = tes.update_entity_history( ent4, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );

			drop_list = tes.update_entity_history( ent5, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( ! drop_list.empty() );
			TS_ASSERT( drop_list.size() == 3 );
			for ( std::list< EntityOP >::const_iterator iter = drop_list.begin(); iter != drop_list.end(); ++iter ) {
				EntityOP dropped = *iter;
				TS_ASSERT( dropped->fitness() == 4.0 );
			}
		}


		{
			TopEntitySet tes;
			tes.desired_entity_history_size( 2 );
			drop_list = tes.update_entity_history( ent1, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent7, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );
			drop_list = tes.update_entity_history( ent2, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( ! drop_list.empty() );
			TS_ASSERT( (*drop_list.begin())->fitness() == 5.0 );

			drop_list = tes.update_entity_history( ent4, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );

			drop_list = tes.update_entity_history( ent3, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( drop_list.empty() );

			drop_list = tes.update_entity_history( ent5, seanpds, was_added );
			TS_ASSERT( was_added );
			TS_ASSERT( ! drop_list.empty() );
			TS_ASSERT( drop_list.size() == 3 );
			for ( std::list< EntityOP >::const_iterator iter = drop_list.begin(); iter != drop_list.end(); ++iter ) {
				EntityOP dropped = *iter;
				TS_ASSERT( dropped->fitness() == 4.0 );
			}
		}


	}


};


