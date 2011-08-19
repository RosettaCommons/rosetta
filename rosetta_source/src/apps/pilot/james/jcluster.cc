// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file cluster_silent.cc
/// @brief
/// @author James Thompson

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/rms_util.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector1.hh>

#include <string>
#include <iostream>

#include <numeric/random/random.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>

template< typename T >
void vec_show( utility::vector1< T > const & values, std::ostream & out ) {
	typedef typename utility::vector1< T >::const_iterator iter;
	for ( iter it = values.begin(), end = values.end();
			it != end; ++it
	) {
		it->show( out );
		out << std::endl;
	}
}

template< typename T >
void output_vec( utility::vector1< T > const & values, std::ostream & out ) {
	typedef typename utility::vector1< T >::const_iterator iter;
	for ( iter it = values.begin(), end = values.end();
			it != end; ++it
	) {
		out << *it << " ";
		out << std::endl;
	}
}

/// @brief class for representing a square distance matrix
class DistanceMatrix : public utility::pointer::ReferenceCount {
public:
	DistanceMatrix( core::Size const max_dim )
		: distances_( max_dim, utility::vector1< core::Real >( max_dim, 0.0 ) )
	{}

	~DistanceMatrix() {}

	core::Size size() const {
		return distances_.size();
	}

	core::Real distance( core::Size const dim1, core::Size const dim2 ) const {
		if ( dim1 == dim2 ) {
			return 0;
		}

		core::Size x_dim( std::min( dim1, dim2 ) );
		core::Size y_dim( std::max( dim2, dim2 ) );
		check_bounds_( x_dim, y_dim );
		return distances_[y_dim][x_dim];
	}

	void distance( core::Size const dim1, core::Size const dim2, core::Real const val ) {
		core::Size x_dim( std::min( dim1, dim2 ) );
		core::Size y_dim( std::max( dim2, dim2 ) );
		check_bounds_( x_dim, y_dim );
		distances_[y_dim][x_dim] = val;
	}

	void show( std::ostream & out ) {
		for ( Size ii = 1; ii <= size(); ++ii ) {
			out << ii;
			for ( Size jj = 1; jj <= size(); ++jj ) {
				out << " " << distance( ii, jj );
			}
			out << std::endl;
		}
	}

private:
	void check_bounds_( core::Size const x_dim, core::Size const y_dim ) const {
		if ( x_dim > size() || y_dim > size() ) {
			utility_exit_with_message( "Error: dimension out of bounds!" );
		}
	}

	utility::vector1< utility::vector1< core::Real > > distances_;
}; // DistanceMatrix

/// assumes that first value is centroid!
class Cluster {
public:
	Cluster() {}
	Cluster( core::Size idx ) : vals_( 1, idx ) {}
	Cluster( utility::vector1< core::Size > vals ) : vals_( vals ) {}

	core::Size operator[] ( core::Size idx ) {
		assert( idx <= vals_.size() );
		return vals_[idx];
	}

	void add( core::Size const idx ) {
		vals_.push_back( idx );
	}

	core::Size size() const {
		return vals_.size();
	}

	core::Size centroid() const {
		runtime_assert( vals_.size() > 0 );
		return centroid_;
		//return vals_[1];
	}

	void centroid( core::Size new_cen ) {
		centroid_ = new_cen;
	}

	core::Size operator[]( core::Size const pos ) const {
		runtime_assert( pos <= vals_.size() );
		return vals_[pos];
	}

	void show( std::ostream & out ) const {
		out << "CLUSTER " << size() << " points, centered at " << centroid()
			<< ". values: ";

		typedef utility::vector1< core::Size >::const_iterator iter;
		for ( iter it = vals_.begin(), end = vals_.end(); it != end; ++it ) {
			out << " " << *it;
		}
	}

private:
	core::Size centroid_;
	utility::vector1< core::Size > vals_;
};

vector1< Cluster > associate_points_with_centroids(
	DistanceMatrix const & matrix,
	vector1< core::Size > const & centroids
) {
	using core::Size;
	using core::Real;
	using utility::vector1;
	vector1< Cluster > clusters;

	vector1< Size > unassigned( matrix.size(), 0 );
	for ( Size ii = 1; ii <= unassigned.size(); ++ii ) unassigned[ii] = ii;

	for ( vector1< Size >::const_iterator cen_it = centroids.begin(),
				cen_end = centroids.end(); cen_it != cen_end; ++cen_it
	) {
		Cluster this_cluster( *cen_it );
		this_cluster.centroid( *cen_it );
		clusters.push_back( this_cluster );

		//std::cout << "looking for " << *cen_it << std::endl;
		vector1< Size >::iterator it(
			find( unassigned.begin(), unassigned.end(), *cen_it )
		);

		if ( it != unassigned.end() ) {
			//std::cout << "erasing (" << *cen_it << "," << *it << ")" << std::endl;
			//std::cout << "before erase unassigned is:";
			//for ( core::Size ii = 1; ii <= unassigned.size(); ++ii ) {
			//	std::cout << " " << unassigned[ii];
			//}
			//std::cout << std::endl;

			unassigned.erase( it );
			//std::cout << "after erase unassigned is:";
			//for ( core::Size ii = 1; ii <= unassigned.size(); ++ii ) {
			//	std::cout << " " << unassigned[ii];
			//}
			//std::cout << std::endl;
		} else {
			//std::cout << "can't find value " << *cen_it << std::endl;
			//std::cout << "unassigned is:";
			//for ( core::Size ii = 1; ii <= unassigned.size(); ++ii ) {
			//	std::cout << " " << unassigned[ii];
			//}
			//std::cout << std::endl;
		}
	}

	//std::cout << "before loop!" << std::endl;
	//vec_show< Cluster >( clusters, std::cout );

	// associate each pose with a cluster
	typedef vector1< Size >::iterator iter;
	typedef vector1< Cluster >::iterator clust_iter;
	for ( iter it = unassigned.begin(), end = unassigned.end();
				it != end; ++it
	) {
		clust_iter closest_cluster = clusters.begin();
		Real min_dist( matrix.distance( closest_cluster->centroid(), *it ) );

		for ( clust_iter c_it = clusters.begin(), c_end = clusters.end();
					c_it != c_end; ++c_it
		) {
			Real const this_dist( matrix.distance( *it, c_it->centroid() ) );
			if ( this_dist <= min_dist ) {
				min_dist        = this_dist;
				closest_cluster = c_it;
			}
		}

		closest_cluster->add( *it );
	}
	//std::cout << "after loop!" << std::endl;
	//vec_show< Cluster >( clusters, std::cout );

	return clusters;
} // associate_points_with_centroids

utility::vector1< core::Size > update_centroids(
	DistanceMatrix const & matrix,
	utility::vector1< Cluster > clusters
) {
	// for each cluster, find the point that is on average nearest
	// to all cluster members.
	using core::Size;
	using core::Real;
	using utility::vector1;

	vector1< Size > new_centroids;
	typedef utility::vector1< Cluster >::const_iterator iter;
	for ( iter outer = clusters.begin(), outer_end = clusters.end();
				outer != outer_end; ++outer
	) {
		//std::cout << "reassigning center " << outer->centroid() << std::endl;
		Real best_total_dist( 1e20 );
		Size most_average_member( clusters.size() + 1 );
		for ( Size ii = 1; ii <= outer->size(); ++ii ) {
			Real total_dist( 0.0 );
			for ( Size jj = 1; jj <= outer->size(); ++jj ) {
				total_dist += matrix.distance( (*outer)[ii], (*outer)[jj] );
			}
			if ( total_dist <= best_total_dist ) {
				best_total_dist = total_dist;
				most_average_member = (*outer)[ii];
			}
		}
		//std::cout << "reassigned to " << most_average_member << std::endl;
		new_centroids.push_back( most_average_member );
	} // outer
	//std::cout << "DONE!" << std::endl;
	return new_centroids;
} // update_centroids

template< typename T >
bool
vectors_are_equal(
	utility::vector1< T > vec1,
	utility::vector1< T > vec2
) {
	if ( vec1.size() != vec2.size() ) return false;
	typedef typename utility::vector1< T >::const_iterator iter;
	for ( iter it1 = vec1.begin(), end1 = vec1.end(),
				it2 = vec2.begin(), end2 = vec2.end();
				it1 != end1 && it2 != end2;
				++it1, ++it2
	) {
		if ( *it1 != *it2 ) return false;
	}

	return true;
}

int
main( int argc, char* argv [] ) {
	// options, random initialization
	devel::init( argc, argv );

	using std::find;
	using core::Size;
	using core::Real;
	using core::pose::Pose;
	using utility::vector1;
	using numeric::random::random_range;
	using namespace basic;
	using namespace core::import_pose::pose_stream;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	MetaPoseInputStream input = streams_from_cmd_line();

	// make starting clusters
	vector1< Pose > poses( input.get_all_poses( *rsd_set ) );

	// calculate distances between all Poses
	DistanceMatrix matrix( poses.size() );
	std::cout << "calculating distances for " << matrix.size() << " objects."
		<< std::endl;
	PROF_START( KDTREE_CONSTRUCT );
	for ( Size ii = 1; ii <= poses.size(); ++ii ) {
		for ( Size jj = ii + 1; jj <= poses.size(); ++jj ) {
			core::Real rmsd = core::scoring::CA_rmsd( poses[ii], poses[jj] );
			matrix.distance( ii, jj, rmsd );
		}
	}
	PROF_STOP( KDTREE_CONSTRUCT );
	prof_show();

	Size const n_compare( ( poses.size() * poses.size() - poses.size() ) / 2 );
	std::cout << "finished with " << n_compare << " comparisons." << std::endl;

	Size const K( 5 ); // number of cluster centers
	runtime_assert( K < matrix.size() );

	vector1< Size > unassigned( matrix.size(), 0 );
	for ( Size ii = 1; ii <= unassigned.size(); ++ii ) unassigned[ii] = ii;

	// randomly guess K cluster centroids
	vector1< core::Size > centroids;
	for ( Size ii = 1; ii <= K; ++ii ) {
		Size random_pt = numeric::random::random_range( 1, unassigned.size() );
		while( find( centroids.begin(), centroids.end(), random_pt ) != centroids.end() ) {
			random_pt = numeric::random::random_range( 1, unassigned.size() );
		} // pick a point without replacement

		vector1< Size >::iterator it(
			unassigned.begin() + random_pt - 1
		);
		centroids.push_back( *it );
	}

	// iteratively update clusters and assign points to nearest clusters
	vector1< Cluster > clusters;
	vector1< Size > last_centroids;
	Size iter_count( 1 );
	Size const max_iter( 200 );
	while ( !vectors_are_equal< Size >( centroids, last_centroids ) ) {
	 	last_centroids = centroids;

		clusters  = associate_points_with_centroids( matrix, centroids );
		//std::cout << "centroids are: ";
		//for ( Size ii = 1; ii <= centroids.size(); ++ii ) {
		//	std::cout << " " << centroids[ii] << std::endl;
		//}
		centroids = update_centroids( matrix, clusters );
		//std::cout << "iteration " << iter_count << std::endl;
		//vec_show< Cluster >( clusters, std::cout );
		//std::cout << std::endl;
		iter_count++;

		std::sort( centroids.begin(), centroids.end() );

		if ( iter_count >= max_iter ) break;
	} // update clustering
	std::cout << "done in " << iter_count << " iterations!" << std::endl;
	vec_show< Cluster >( clusters, std::cout );

	return 0;
}
