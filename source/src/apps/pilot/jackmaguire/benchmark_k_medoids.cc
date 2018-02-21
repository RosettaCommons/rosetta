// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jackmaguire/benchmark_k_medoids.cc
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/multistage_rosetta_scripts/cluster/KMedoids.hh>
#include <protocols/multistage_rosetta_scripts/cluster/KMedoidsOnTheFly.hh>

#include <protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetric.hh>
#include <protocols/multistage_rosetta_scripts/cluster/metrics/SequenceMetric.hh>

#include <numeric/random/random.hh>

#include <string>
#include <chrono>

static basic::Tracer TR( "apps.pilot.jackmaguire.BenchmarkKMedoids" );

using namespace protocols::multistage_rosetta_scripts::cluster;
using namespace protocols::multistage_rosetta_scripts::cluster::metrics;

void test(){
	utility::vector1< SequenceMetricOP > points;

	points.emplace_back( utility::pointer::make_shared< SequenceMetric >( "AAAAAAAAAAAAAAAAAAAAA" ) );
	points.emplace_back( utility::pointer::make_shared< SequenceMetric >( "AAAAAASAAASAAAAAAAAAA" ) );
	points.emplace_back( utility::pointer::make_shared< SequenceMetric >( "AAAAAASAAAAAAAAAAAAAA" ) );
	points.emplace_back( utility::pointer::make_shared< SequenceMetric >( "AAAAAAAAAAAAAAAAAAAAA" ) );
	points.emplace_back( utility::pointer::make_shared< SequenceMetric >( "AAAAAAAAASSSAAAAAAAAA" ) );
	points.emplace_back( utility::pointer::make_shared< SequenceMetric >( "AAAASAAAASSAAAAAAAAAA" ) );
	points.emplace_back( utility::pointer::make_shared< SequenceMetric >( "RRRRRRRRRSSSRRRRRRRRR" ) );
	points.emplace_back( utility::pointer::make_shared< SequenceMetric >( "RRRRRRRRRRRRRRRRRSRRR" ) );

	for ( auto i = 1; i <= 8; ++i ) {
		for ( auto j = 1; j <= 8; ++j ) {
			TR << "\t" << points[ i ]->distance( *points[ j ] );
		}
		TR << std::endl;
	}

	utility::vector1< bool > results =
		k_medoids_with_edge_precalculation( points, 2 );

	for ( bool b : results ) {
		TR << b << std::endl;
	}

}

SequenceMetricOP
generate_sequence_binary_alphabet( core::Size length ){
	std::string sequence( length, 'A' );

	for ( core::Size ii = 0; ii < length; ++ii ) {
		if ( numeric::random::rg().random_range( 0, 2 ) == 1 ) {
			sequence[ ii ] = 'R';
		}
	}

	return SequenceMetricOP( new SequenceMetric( sequence ) );
}


int main2( int argc, char* argv[] ){

	devel::init( argc, argv );
	test();
	return 0;

}

int main( int argc, char* argv[] ){

	try {

		devel::init( argc, argv );

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		//repurposing existing option so that I do not have to create a new one
		bool on_the_fly = option[ jd3::do_not_archive_on_node0 ].user();

		core::Size nloop = 2500;
		if ( option[ jd3::n_archive_nodes ].user() ) {
			//repurposing existing option so that I do not have to create a new one
			nloop = option[ jd3::n_archive_nodes ]();
		}

		core::Size num_clusters = 0;

		for ( core::Size size = 100; size <= nloop; size += 100 ) {
			utility::vector1< SequenceMetricOP > sequences( size );
			for ( core::Size index = 1; index <= size; ++index ) {
				sequences[ index ] = generate_sequence_binary_alphabet( 250 );
			}

			auto start = std::chrono::steady_clock::now();
			//utility::vector1< bool > results = k_medoids_with_edge_precalculation( sequences, ++num_clusters );
			num_clusters += 10;
			if ( on_the_fly ) {
				k_medoids_on_the_fly( sequences, num_clusters );
			} else {
				k_medoids_with_edge_precalculation( sequences, num_clusters );
			}
			auto end = std::chrono::steady_clock::now();
			auto time_elapsed = end - start;
			TR << size << "\t" << std::chrono::duration < double, std::nano > ( time_elapsed ).count() << std::endl;
		}

		return 0;

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
