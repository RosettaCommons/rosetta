// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file cluster_silent.cc
/// @brief
/// @author James Thompson

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector1.hh>

#include <string>
#include <iostream>

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <numeric/random/random.fwd.hh>

#include <utility/excn/Exceptions.hh>

utility::vector1< unsigned int >
get_dm(
	core::pose::Pose & pose,
	double const dist_threshold = 9,
	std::string const & atom_name = "CA"
) {
	using utility::vector1;

	double const dist_threshold_sq( dist_threshold * dist_threshold );
	unsigned int const N( pose.size() );

	vector1< unsigned int > distances( (unsigned int)(N*N-N)/2, 0 );
	unsigned int feat_idx(1);
	for ( unsigned int ii = 1; ii <= pose.size(); ++ii ) {
		for ( unsigned int jj = ii+1; jj <= pose.size(); ++jj ) {
			double const dist_sq(
				pose.residue(ii).xyz(atom_name).distance_squared(pose.residue(jj).xyz(atom_name))
			);

			if ( dist_sq < dist_threshold_sq ) {
				distances[feat_idx] = 1;
			}
			feat_idx++;
		} // jj
	} // ii

	return distances;
}

void print_debug_information(
	//utility::vector1< utility::vector1< double > > const & clusters,
	utility::vector1< utility::vector1< double > > const & resp
) {
	for ( unsigned int ii = 1; ii <= resp.size(); ++ii ) {
		for ( unsigned int jj = 1; jj <= resp.front().size(); ++jj ) {
			std::cout << " " << resp[ii][jj] << std::endl;
		}
		std::cout << std::endl;
	}
}

void print_dm(
	utility::vector1< unsigned int > dm
) {
	for ( unsigned int ii = 1; ii <= dm.size(); ++ii ) {
		std::cout << dm[ii];
	}
	std::cout << std::endl;
}

int
main( int argc, char* argv [] ) {
	try {

		devel::init( argc, argv );

		using core::pose::Pose;
		using utility::vector1;
		using namespace basic;
		using namespace core::import_pose::pose_stream;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		core::chemical::ResidueTypeSetCOP rsd_set(
			core::chemical::ChemicalManager::get_instance()->residue_type_set(
			option[ in::file::residue_type_set ]()
			)
		);

		MetaPoseInputStream input = streams_from_cmd_line();
		vector1< core::pose::PoseOP > poses( input.get_all_poses( *rsd_set ) );

		vector1< vector1< unsigned int > > data_points;
		for ( unsigned int ii = 1; ii <= poses.size(); ++ii ) {
			data_points.push_back( get_dm(*poses[ii]) );
			//std::cout << "point " << ii << ":";
			//print_dm( data_points[ii] );
		}
		//std::cout << "have " << data_points.size() << " points." << std::endl;

		// make starting clusters
		unsigned int const K(10);
		unsigned int const n_features( data_points.front().size() );
		//double const distance_cutoff(9);  // unused ~Labonte
		vector1< vector1< double > > clusters( K, vector1< double >( n_features, 0.0 ) );
		runtime_assert( data_points.size() > K );
		for ( unsigned int ii = 1; ii <= clusters.size(); ++ii ) {
			for ( unsigned int kk = 1; kk <= n_features; ++kk ) {
				clusters[ii][kk] = numeric::random::uniform();
			}
		}
		std::cout << "have " << clusters.size() << " clusters." << std::endl;

		// compute responsibility of each cluster for each data point
		vector1< vector1< double > > responsibilities(
			data_points.size(), vector1< double >( clusters.size(), 0.0 )
		);
		for ( unsigned int ii = 1; ii <= data_points.size(); ++ii ) {
			double total_ii(0.0);
			for ( unsigned int jj = 1; jj <= clusters.size(); ++jj ) {
				for ( unsigned int kk = 1; kk <= data_points[ii].size(); ++kk ) {
					if ( data_points[ii][kk] > 0 ) {
						total_ii += clusters[ii][kk];
						responsibilities[ii][jj] += clusters[ii][kk];
					}
				} // kk
			} // jj

			// normalize according to total_ii
			for ( unsigned int jj = 1; jj <= clusters.size(); ++jj ) {
				responsibilities[ii][jj] /= total_ii;
				std::cout << "resp(" << ii << "," << jj << ") = " << responsibilities[ii][jj] << std::endl;
			}
		} // ii

		// update cluster probabilities according to the responsibilities
		for ( unsigned int ii = 1; ii <= data_points.size(); ++ii ) {
			for ( unsigned int jj = 1; jj <= clusters.size(); ++jj ) {
				for ( unsigned int kk = 1; kk <= n_features; ++kk ) {
					if ( data_points[ii][kk] > 0 ) {
						double const & resp( responsibilities[ii][jj] );
						clusters[ii][kk] += resp;
					}
				}
			}
		}

		//print_debug_information( responsibilities );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
