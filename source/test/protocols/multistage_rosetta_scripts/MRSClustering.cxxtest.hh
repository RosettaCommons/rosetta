// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/cluster/MRSClustering.cxxtest.hh
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Package headers
#include <protocols/multistage_rosetta_scripts/cluster/KMedoids.hh>

#include <protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetric.hh>
#include <protocols/multistage_rosetta_scripts/cluster/metrics/SequenceMetric.hh>
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetricFactory.hh>

#include <basic/Tracer.hh>
#include <cxxtest/TestSuite.h>

#include <test/protocols/init_util.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>

#ifdef SERIALIZATION
// Utility headers
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

#endif

static basic::Tracer TR("test.protocols.multistage_rosetta_scripts.MRSClustering");

using namespace protocols::multistage_rosetta_scripts;
using namespace protocols::multistage_rosetta_scripts::cluster;
using namespace protocols::multistage_rosetta_scripts::cluster::metrics;

class MRSClusteringTests : public CxxTest::TestSuite {

public:
	void setUp() override {
		protocols_init();
	}

	void test_factory(){
		ClusterMetricFactory::get_instance()->has_type( "JumpMetric" );
		ClusterMetricFactory::get_instance()->has_type( "SequenceMetric" );
	}

#ifdef SERIALIZATION
	template < class T >
	std::string
	serialized_T( T const & data ) {
		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( data );
		}
		return oss.str();
	}


	template < class T >
	T
	deserialized_T( std::string const & str ) {
		std::istringstream iss( str );
		T return_object;
		cereal::BinaryInputArchive arc( iss );
		arc( return_object );
		return return_object;
	}
#endif

	void test_sequence_serialization(){
#ifdef SERIALIZATION
		TR << "Running test_sequence_serialization" << std::endl;
		SequenceMetric simple( "Jack Maguire Is Very Handsome");
		std::string data = serialized_T< SequenceMetric >( simple );
		SequenceMetric simple2 = deserialized_T< SequenceMetric >( data );
		TS_ASSERT( ! simple.sequence().compare( simple2.sequence() ) );
		TS_ASSERT_EQUALS( simple.distance( simple2 ), 0 );
#endif
		TS_ASSERT( true );
	}

	void test_jump_serialization(){
#ifdef SERIALIZATION
		TR << "Running test_jump_serialization" << std::endl;
		utility::vector1< core::Real > dofs( 6, 8.5 );
		JumpMetric simple( dofs );
		std::string data = serialized_T< JumpMetric >( simple );
		JumpMetric simple2 = deserialized_T< JumpMetric >( data );

		for ( int i=1; i<=6; ++i ) {
			TS_ASSERT_EQUALS( simple2.dofs()[ i ], simple.dofs()[ i ] );
			TS_ASSERT_EQUALS( simple2.dofs()[ i ], 8.5 );
		}
		TS_ASSERT_EQUALS( simple.distance( simple2 ), 0 );
#endif
		TS_ASSERT( true );
	}

};


