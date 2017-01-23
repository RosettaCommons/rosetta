// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// the incorrect documentation below is an embarassment.
/// @file   protocols/jd2/Job.hh
/// @brief  header file for ThreadingJob classes, part of August 2008 job distributor as planned at RosettaCon08.  This file is responsible for three ideas: "inner" jobs, "outer" jobs (with which the job distributor works) and job container (currently just typdefed in the .fwd.hh)
/// @author Steven Lewis smlewi@gmail.com


#include <protocols/toolbox/Cluster.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {

using namespace ObjexxFCL;
using namespace format;

static THREAD_LOCAL basic::Tracer tr( "protocols.cluster", basic::t_info );


using namespace core;

ClusterOptions::ClusterOptions( bool assign_new_cluster_tag_in )
: assign_new_cluster_tag( assign_new_cluster_tag_in ),
	limit_cluster_size( basic::options::option[ basic::options::OptionKeys::cluster::limit_cluster_size ] ),
	cluster_radius( basic::options::option[ basic::options::OptionKeys::cluster::radius ] ),
	keep_center( true )
{ }

void ClusterBase::print_cluster_assignment( std::ostream & out ) const {
	for ( Size i=1; i <= clusterlist_.size(); i++ ) {
		for ( Size j=1; j<= clusterlist_[i].size(); j++ ) {
			out << RJ( 3,clusterlist_[i] [ j-1 ]) << " " << RJ( 5,i) << RJ(5, j)
				<< std::endl;
		}
	}
}


std::ostream & operator<< ( std::ostream & out, ClusterBase const & cl) {
	cl.show( out );
	return out;
}

std::ostream & operator<< (
	std::ostream & out,
	ClusterBase::ClusterList const & cl
) {
	Size ncl( 1 );
	for ( auto const & it : cl ) {
		out << "CLUSTER " << RJ(3, ncl++ ) << " " << RJ(4, it.size()) << " ";
		for ( core::Size cit : it ) {
			out << RJ( 4, cit ) << " ";
		}
	}
	return out;
}

std::istream & operator>>( std::istream & in, ClusterBase & cl ) {
	cl.read( in );
	return in;
}


std::istream & operator>>( std::istream & in, ClusterBase::ClusterList & cl ) {
	cl.clear();
	while ( in ) {
		std::string tag;
		Size size,nr;
		in >> tag >> nr >> size;
		if ( !in.good() ) break;
		if ( tag != "CLUSTER" ) {
			in.setstate( std::ios_base::failbit );
			break;
		}
		ClusterBase::Cluster new_clust;
		for ( Size i = 1; i<=size && in.good(); i++ ) {
			Size new_elem;
			in >> new_elem;
			if ( in.good() ) new_clust.push_back( new_elem );
		}
		if ( new_clust.size() == size ) {
			cl.push_back( new_clust );
		}
	}
	return in;
}

void ClusterBase::show( std::ostream & out ) const {
	out << clusterlist_;
}

void ClusterBase::read( std::istream & in ) {
	in >> clusterlist_;
}

void ClusterBase::print_summary( utility::vector1< std::string > tags, utility::vector1< core::Real > all_energies ) {
	tr.Info << "---------- Summary ---------------------------------" << std::endl;

	int count = 0;

	for ( Size i = 1; i <= clusterlist_.size(); i++ ) {
		tr.Info << "Cluster:  " << i << "  N: " << clusterlist_[i].size() << "   c." << i << ".*.pdb " ;
		tr.Info << std::endl;
		count += clusterlist_[i].size();
		for ( Size j=1; j <= clusterlist_[i].size(); j++ ) {
			tr.Info << tags[ clusterlist_[i][j-1] ] << "    " << all_energies[ clusterlist_[i][j-1] ]
				<< " " << dist( clusterlist_[i][j-1], clusterlist_[i][0] )
				<< std::endl;
		}
		tr.Info << std::endl;
	}
	tr.Info << "----------------------------------------------------" << std::endl;
	tr.Info << "  Clusters: " << clusterlist_.size() << std::endl;
	tr.Info << "  Structures: " << count << std::endl;
	tr.Info << "----------------------------------------------------" << std::endl;
}

void ClusterPhilStyle::compute() {
	utility::vector1< Size > neighbors ( dim(), 0 );
	utility::vector1< Size > clusternr ( dim(), 0 );
	utility::vector1< Size > clustercenter;
	Size    mostneighbors;
	Size    nclusters = 0;

	Size listsize = dim();

	utility::vector1<Size> clustercentre;
	// now assign groupings
	while ( true ) {
		// count each's neighbors
		for ( Size i=1; i <= listsize; i++ ) {
			neighbors[i] = 0;
			if ( clusternr[i] > 0 ) continue; // ignore ones already taken
			for ( Size j=1; j <= listsize; j++ ) {
				if ( clusternr[j] > 0 ) continue; // ignore ones already taken
				if ( dist( i, j ) < cluster_radius_ ) neighbors[i]++;
			}
		}

		mostneighbors = 1;
		for ( Size i=1; i <= listsize; i++ ) {
			if ( neighbors[i] > neighbors[mostneighbors] ) mostneighbors=i;
		}
		if ( neighbors[ mostneighbors ] == 0 ) break;  // finished!
		for ( Size i=1; i <= listsize; i++ ) {
			if ( clusternr[i] > 0 ) continue; // ignore ones already taken
			if ( dist( i, mostneighbors ) < cluster_radius_ ) {
				clusternr[i] = mostneighbors;
			}
		}

		clustercentre.push_back(mostneighbors);
		nclusters++;

		if ( nclusters > n_max_cluster_ ) break;  // currently fixed but ought to be a paraemter
	}


	for ( Size i=1; i <= clustercentre.size(); i++ ) {
		std::deque< Size > newlist;
		for ( Size j=1; j <=listsize; j++ ) {
			// if that struture belongs to a given cluster
			if ( clusternr[j] == clustercentre[i] ) {
				//add it
				//newlist.push_back(j);
				if ( j!= clustercentre[i] ) newlist.push_back(j);         // add structure
				else                     newlist.push_front(j);        // add cluster centre at the front
			}
		}

		// add the new list to the clusterlist_
		clusterlist_.push_back(newlist);
	}

	// redistribute groups - i.e. take each structure, calculate the rms to each cluster centre.
	for ( Size i = 1; i <= clusterlist_.size(); i++ ) {
		for ( Size j = 1; j <= clusterlist_[i].size(); j++ ) {
			Size lowrmsi=i;
			Real lowrms=10000.0;
			for ( Size m=1; m <= clusterlist_.size(); m++ ) {
				Real rms;
				rms = dist( clusterlist_[i][j - 1],  // current structure
					clusterlist_[m][0]); // clustercentre m
				//std::cout << "C: " << i << " S: " << j << " rms " << m << " --> " << rms << " " << std::endl;

				if ( rms < lowrms ) {
					lowrms  = rms;
					lowrmsi = m;
				}
			}
			if ( lowrmsi != i ) { // is a different cluster centre closer to us than our current centre ?
				// switch over;
				clusterlist_[ lowrmsi ].push_back( clusterlist_[i][j-1] );
				auto it = clusterlist_[i].begin();
				it += j-1;
				clusterlist_[ i ].erase( it );
			}
		} // j
	} // i
} // compute()

bool compareIndexEnergyPair(
	const std::pair< int, float > & p1, const std::pair< int, float > & p2 )
{
	return p1.second  < p2.second;
}


void ClusterBase::sort_each_group_by_energy( utility::vector1< core::Real > all_energies, bool keep_center ) {
	Size const offset( keep_center ? 1 : 0 );
	for ( Size i=1; i<= clusterlist_.size(); i++ ) {
		utility::vector1< std::pair< int, float > > cluster_energies;
		for ( Size j=1; j<= clusterlist_[i].size(); j++ ) {
			Real score = all_energies[ clusterlist_[i][j-1] ];
			cluster_energies.push_back( std::pair< int, Real > ( clusterlist_[i][j-1], score ) );
		}
		// dont sort in the first member (the cluster center!!)
		std::sort( cluster_energies.begin() + offset, cluster_energies.end(), compareIndexEnergyPair );
		clusterlist_[i].clear();
		for ( Size j=1; j <= cluster_energies.size(); j++ ) clusterlist_[i].push_back( cluster_energies[j].first );
	}
}

// take first 'limit' from each cluster
void ClusterBase::limit_groupsize( Size limit ) {
	if ( limit == 0 ) clusterlist_.clear();
	for ( Size i=1; i <= clusterlist_.size(); i++ ) {
		Cluster temp = clusterlist_[i];
		clusterlist_[i].clear();
		for ( Size j=1; ( j<= temp.size() ) && ( j<=limit ); j++ ) {
			clusterlist_[i].push_back( temp[j-1] );
		}
	}
}


} // toolbox
} // protocols
