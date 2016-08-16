// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/EnergyGraph.cxxtest.hh
/// @brief  test suite for core::scoring::EnergyGraph.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/types.hh>
#include <core/graph/Graph.hh>
#include <core/scoring/EnergyGraph.hh>

// C++ headers, for debugging your tests
#include <sstream>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

// --------------- Test Class --------------- //

class EnergyGraphTests : public CxxTest::TestSuite {
public:

	void test_serialize_energy_graph() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		using namespace core;
		using namespace core::scoring;
		EnergyGraph egraph( 5 );
		ScoreTypes sts;
		sts.push_back( fa_atr );
		sts.push_back( fa_rep );
		egraph.active_score_types( sts );

		std::set< std::pair< Size, Size > > edges;
		edges.insert( std::make_pair( 1, 2 ));
		edges.insert( std::make_pair( 2, 3 ));
		edges.insert( std::make_pair( 3, 4 ));
		edges.insert( std::make_pair( 4, 5 ));
		EnergyMap emap;
		for ( std::set< std::pair< Size, Size > >::const_iterator iter = edges.begin(), iter_end = edges.end(); iter != iter_end; ++iter ) {
			egraph.add_edge( iter->first, iter->second );
			EnergyEdge * eeptr = egraph.find_energy_edge( iter->first, iter->second );
			emap[ fa_atr ] = -1.125*iter->first; emap[ fa_rep ] = 1.125*iter->second;
			eeptr->store_active_energies( emap );
			eeptr->square_distance( iter->first * iter->second );
			if ( iter->first % 2 ) {
				eeptr->mark_energies_computed();
			} else {
				eeptr->mark_energies_uncomputed();
			}
		}

		for ( Size ii = 1; ii <= 5; ++ii ) {
			egraph.get_energy_node( ii )->moved( ii % 2 );
		}

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arch( oss );
			arch( egraph );
		}

		EnergyGraph egraph2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arch( iss );
			arch( egraph2 );
		}

		TS_ASSERT_EQUALS( egraph2.num_nodes(), 5 );
		for ( Size ii = 1; ii <= 5; ++ii ) {
			for ( Size jj = ii+1; jj <= 5; ++jj ) {
				if ( std::find( edges.begin(), edges.end(), std::make_pair( ii, jj ) ) == edges.end() ) {
					TS_ASSERT( ! egraph2.get_edge_exists( ii, jj ) );
				} else {
					if ( ! egraph2.get_edge_exists( ii, jj ) ) {
						TS_ASSERT( egraph2.get_edge_exists( ii, jj ) );
						continue;
					}
					EnergyEdge const * ee = egraph2.find_energy_edge( ii, jj );
					TS_ASSERT_EQUALS( (*ee)[ fa_atr ], -1.125 * ii );
					TS_ASSERT_EQUALS( (*ee)[ fa_rep ],  1.125 * jj );
					TS_ASSERT_EQUALS( ee->square_distance(), ii * jj );
					TS_ASSERT_EQUALS( ee->energies_not_yet_computed(), ii % 2 == 0 )
				}
			}
		}

		for ( Size ii = 1; ii <= 5; ++ii ) {
			TS_ASSERT_EQUALS( egraph2.get_energy_node( ii )->moved(), ii % 2 == 1 );
		}

#endif
	}

	void test_serialize_energy_graph_from_graph_op() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		using namespace core;
		using namespace core::graph;
		using namespace core::scoring;
		EnergyGraphOP egraph( new EnergyGraph( 5 ) );
		ScoreTypes sts;
		sts.push_back( fa_atr );
		sts.push_back( fa_rep );
		egraph->active_score_types( sts );

		std::set< std::pair< Size, Size > > edges;
		edges.insert( std::make_pair( 1, 2 ));
		edges.insert( std::make_pair( 2, 3 ));
		edges.insert( std::make_pair( 3, 4 ));
		edges.insert( std::make_pair( 4, 5 ));
		EnergyMap emap;
		for ( std::set< std::pair< Size, Size > >::const_iterator iter = edges.begin(), iter_end = edges.end(); iter != iter_end; ++iter ) {
			egraph->add_edge( iter->first, iter->second );
			EnergyEdge * eeptr = egraph->find_energy_edge( iter->first, iter->second );
			emap[ fa_atr ] = -1.125*iter->first; emap[ fa_rep ] = 1.125*iter->second;
			eeptr->store_active_energies( emap );
			eeptr->square_distance( iter->first * iter->second );
			if ( iter->first % 2 ) {
				eeptr->mark_energies_computed();
			} else {
				eeptr->mark_energies_uncomputed();
			}
		}

		for ( Size ii = 1; ii <= 5; ++ii ) {
			egraph->get_energy_node( ii )->moved( ii % 2 );
		}

		GraphOP g( egraph );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arch( oss );
			arch( g );
		}

		GraphOP g2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arch( iss );
			arch( g2 );
		}

		EnergyGraphOP egraph2 = utility::pointer::dynamic_pointer_cast< EnergyGraph > ( g2 );
		TS_ASSERT( egraph2 );
		if ( ! egraph2 ) {
			std::cerr << "Failed to extract an EnergyGraph from a pointer to the base class, a GraphOP." << std::endl;
			return;
		}

		TS_ASSERT_EQUALS( egraph2->num_nodes(), 5 );
		for ( Size ii = 1; ii <= 5; ++ii ) {
			for ( Size jj = ii+1; jj <= 5; ++jj ) {
				if ( std::find( edges.begin(), edges.end(), std::make_pair( ii, jj ) ) == edges.end() ) {
					TS_ASSERT( ! egraph2->get_edge_exists( ii, jj ) );
				} else {
					if ( ! egraph2->get_edge_exists( ii, jj ) ) {
						TS_ASSERT( egraph2->get_edge_exists( ii, jj ) );
						continue;
					}
					EnergyEdge const * ee = egraph2->find_energy_edge( ii, jj );
					TS_ASSERT_EQUALS( (*ee)[ fa_atr ], -1.125 * ii );
					TS_ASSERT_EQUALS( (*ee)[ fa_rep ],  1.125 * jj );
					TS_ASSERT_EQUALS( ee->square_distance(), ii * jj );
					TS_ASSERT_EQUALS( ee->energies_not_yet_computed(), ii % 2 == 0 )
				}
			}
		}

		for ( Size ii = 1; ii <= 5; ++ii ) {
			TS_ASSERT_EQUALS( egraph2->get_energy_node( ii )->moved(), ii % 2 == 1 );
		}

#endif
	}

};
