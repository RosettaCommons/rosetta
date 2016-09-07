// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ShortestPathInFoldTree.cc
/// @brief helper class to FoldTree: allows to measure distance in fold-trees
/// @details This class provides a fast mechanism to determine the distance between residues
/// according to a given fold-tree
/// instead of storing a full NxN matrix with N number of residues
/// we store only MxM matrix for distances between jump_residue (at most M=2*J J nr of jumps)
/// a table with N entries gives for each peptide the distance to the next jump-node.
/// thus memory requirement is low and still,
/// a single dist-evaluation requires to check only 4 possible pathways
///
/// @author Oliver Lange


// Unit Headers
#include <core/kinematics/ShortestPathInFoldTree.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>

// Project Headers
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers

//// C++ headers
#include <string>

#include <basic/Tracer.hh>

//Auto using namespaces
#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/serialization/ObjexxFCL/FArray2D.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace ObjexxFCL {
namespace format {
}
}
using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace core {
namespace kinematics {

/// @details Auto-generated virtual destructor
ShortestPathInFoldTree::~ShortestPathInFoldTree() = default;

static THREAD_LOCAL basic::Tracer tr( "core.kinematics.ShortestPathInFoldTree", basic::t_info );


/// @detail cs-tor
ShortestPathInFoldTree::ShortestPathInFoldTree(
	core::kinematics::FoldTree const& f
) : nres_( f.nres() ), max_dist_( 0 )
{
	simple_fold_tree_ = ( f.num_jump() == 0 );
	if ( simple_fold_tree_ ) max_dist_=nres_;
	if ( !simple_fold_tree_ ) {
		build_jumpres_distmap( f );
		build_peptide_table( f );
	}
}

/// @detail cs-tor
ShortestPathInFoldTree::ShortestPathInFoldTree(
	ShortestPathInFoldTree const & src
) :
	jump_res_( src.jump_res_ ),
	node_dist_( src.node_dist_ ),
	res2jumps_( src.res2jumps_ ),
	nres_( src.nres_ ),
	simple_fold_tree_( src.simple_fold_tree_ ),
	max_dist_( src.max_dist_ )
{}


/// @detail the core of the distance cache is build here: node_dist_
/// a 2D array that knows distances for each pair of jump-residues
/// i.e., a 10 100 1 /  100 120 -1 / 120 150 2 -- fold-tree would have
/// jump_residues 10, 100, 120 and 150.
/// distances would be 10, 100   1
///                    120 150   1
//                     10  150   22
///                    100 150   21
///                    ....
/// for the jump_residues we use internal numbering, each  seqpos that is a jump-res in one or more jumps will get a
/// individual number, starting at 1 and counting continously.
/// thus the distance of jump_res-pair i,j is found as node_dist(i,j)
///
/// first we go through the fold-tree to find all jump_residues, and build the map: jump_res_
/// to keep track between mapping    "seqpos<-->internal_numbering"
///
/// then assign distance 1 to each pair of jump_residues connected by a jump
/// and using Warshall algorithm to build up full distance matrix
/// diagonal is distance 0
void
ShortestPathInFoldTree::build_jumpres_distmap( core::kinematics::FoldTree const& f ) {
	using namespace core::kinematics;
	// go through jumps and find the jump-residues start distance table
	// maintain an EdgeList to memorize the pairs until full number of unique jump-residues is known

	//  typedef utility::vector1< Edge > EdgeList;
	EdgeList edges;


	/// cycle through list of jumps --> store jump_residues and assign running number (ct) to them
	int ct = 1; /// for giving individual numbers to jump_residues
	for ( Size jump = 1; jump <= f.num_jump(); jump++ ) {
		Size const start ( f.jump_edge( jump ).start() );
		Size const stop ( f.jump_edge( jump ).stop() );
		tr.Trace << "add jump " << start << "-" << stop << std::endl;


		// search in our list of jump-residues, assign to my_start/my_stop if we got it already
		//
		// setup
		int my_start = -1; //-1 denotes: "not found"
		int my_stop = -1;
		std::map< Size, Size>::const_iterator fit;

		// look for start residue
		fit =  jump_res_.find( start );
		if ( fit != jump_res_.end() ) {
			my_start = fit->second;
		} else {
			my_start = ct;
			jump_res_[ start ] = ct++;
		}

		// look for stop residue
		fit = jump_res_.find( stop );
		if ( fit != jump_res_.end() ) {
			my_stop = fit->second;
		} else {
			my_stop = ct;
			jump_res_[ stop ] = ct++;
		}
		edges.push_back( Edge(my_start, my_stop, jump ) );
	};

	//if root is not residue 1
	//add the root of the tree to the list of jumps as a jump-onto-itself
	if ( f.root() > 1 ) {
		int nr = get_jump( f.root() );
		if ( nr < 0 ) { //haven't got root already in my list of jumps
			jump_res_[ f.root() ] = ct++;
		}
	}

	// create some debug output
	if ( tr.Trace.visible() ) {
		tr.Trace << " jump_res_nr -  seqpos  \n";
		for ( std::map< Size, Size>::const_iterator it=jump_res_.begin(), eit=jump_res_.end();
				it!=eit; ++it ) {
			tr.Trace << it->second << " - " << it->first << std::endl;
		}

		tr.Trace << " jump-edges -- internal enumeration \n" ;
		for ( EdgeList::const_iterator it=edges.begin(), eit=edges.end(); it!=eit; ++it ) {
			tr.Trace << it->start() <<  " -- " << it->stop() << std::endl;
		}
	}

	init_dist_map( edges );
	compute_dist_map( f );
}


/// @detail initialize dist map with dist 1 for each pair of residues connected by jumps
/// as stored in the EdgeList
void
ShortestPathInFoldTree::init_dist_map( EdgeList const& edges ) {

	//  Warshall algorithm
	Size const inf( 12345678 ); //assumption: fewer than 12 million nodes in the graph.
	debug_assert( jump_res_.size() < inf );
	node_dist_.dimension( jump_res_.size(), jump_res_.size(), inf );

	// initialize distance array with jump-edges
	for (const auto & edge : edges) {
		node_dist_( edge.start(), edge.stop() ) = 1;
		node_dist_( edge.stop(), edge.start() ) = 1;
		node_dist_( edge.start(), edge.start() ) = 0;
		node_dist_( edge.stop(), edge.stop() ) = 0;
	}
}


/// @detail
// to compute the full dist map we go through 2 steps:
//
// (1) get distanes on "simple paths" i.e., length of peptide edges that connect two jumps
// (2) use warshall algo to get all distances by combining distances via jumps and simple peptide edges
void
ShortestPathInFoldTree::compute_dist_map( FoldTree const& f ) {


	// look for peptid edges that connect two jumps
	for (const auto & it : f) {
		if ( it.is_jump() ) continue; // only look at peptide edges
		std::map< Size, Size>::const_iterator fit;

		// do we have start and stop reside listed as jump_residues?
		int my_start = -1; //-1 denotes "not found"
		int my_stop = -1;

		// look for start residue
		fit =  jump_res_.find( it.start() );
		if ( fit != jump_res_.end() ) {
			my_start = fit->second;
		}

		// look for stop residue
		fit =  jump_res_.find( it.stop() );
		if ( fit != jump_res_.end() ) {
			my_stop = fit->second;
		}

		// if start and stop are jump-residues this is an internal peptide edge!
		if ( my_start > 0 && my_stop > 0 ) {
			Size dd = node_dist_( my_start, my_stop ) = std::abs( it.start() - it.stop() );
			node_dist_( my_stop, my_start ) = dd;
		};
	} // for fold-tree edges


	// Warshall algorithm
	// symmetry makes this marginally inefficient, but easy to read
	// if this shows up in a hotspot, it can be made more efficient
	for ( Size ii = 1; ii <= jump_res_.size(); ++ii ) {
		for ( Size jj = 1; jj <= jump_res_.size(); ++jj ) {
			for ( Size kk = 1; kk <= jump_res_.size(); ++kk ) {
				Size const jj_2_kk = node_dist_( jj, kk );
				Size const jj_2_ii = node_dist_( jj, ii );
				Size const ii_2_kk = node_dist_( ii, kk );

				Size const jj_2_ii_2_kk = jj_2_ii + ii_2_kk;

				if ( jj_2_kk > jj_2_ii_2_kk ) {
					node_dist_( jj, kk ) =  jj_2_ii_2_kk;
					node_dist_( kk, jj ) =  jj_2_ii_2_kk;
				}
			}
		}
	}

	// produce some debug output
	//   if ( tr.Trace.visible() ) {
	//     tr.Trace << "jump_res distance table:\n";
	//     for ( Size ii = 1; ii <= jump_res_.size(); ++ ii ) {
	//       for ( Size jj = 1; jj <= jump_res_.size(); ++ jj ) {
	//     //tr.Trace << node_dist_( ii, jj ) << " ";
	//       }
	//       tr.Trace << std::endl;
	//     }
	//   }
}// compute_dist_map

/// @detail build table that gives for each residue the distance to
/// upstream and downstream jump-residues (if available)
///
/// format:
/// <edge_nr> <jump1> <dist1> <jump2> <dist2>
//
/// edge_nr is a unique number that identifes peptide-edges
/// jump1 and jump2 refers to our internally-numbered jump-residues (entries in node_dist_ )
/// dist -- distance in sequence to the respective jump-residues
void
ShortestPathInFoldTree::build_peptide_table( core::kinematics::FoldTree const& f ) {
	using namespace core::kinematics;
	using namespace  ObjexxFCL::format;
	res2jumps_.dimension( f.nres(), 5, -1 );   //5 entries per residue.
	utility::vector1< int > leaves;
	// go thru edges and fill res2jump array accordingly
	Size edge_nr = 1;
	for ( auto it=f.begin(), eit=f.end();
			it!=eit;
			++it, ++edge_nr ) {
		int start_jump = get_jump( it->start() ); // returns -1 if node is not a jump residue
		int stop_jump = get_jump( it->stop() );
		if ( start_jump < 0 ) leaves.push_back( it->start() );
		if ( stop_jump < 0 ) leaves.push_back( it->stop() );
		if ( !it->is_jump() ) { // a peptide edge
			for ( int seqpos = std::min(it->start(),it->stop()); seqpos<=std::max( it->stop(),it->start() ); seqpos++ ) {
				res2jumps_( seqpos, 1 ) = edge_nr;
				res2jumps_( seqpos, 2 ) = start_jump;
				res2jumps_( seqpos, 3 ) = std::abs( (int) seqpos - (int) it->start() );
				res2jumps_( seqpos, 4 ) = stop_jump;
				res2jumps_( seqpos, 5 ) = std::abs( (int) seqpos - (int) it->stop() );
			} // for seqpos
			//    int edge_length = std::abs( (int) it->stop() - (int) it->start() );
			//    if (( start_jump > 0 && stop_jump < 0 ) && ( maxdist2leave( start_jump ) < edge_length ) ) maxdist2leave( start_jump ) = edge_length;
			//    if (( start_jump < 0 && stop_jump > 0 ) && ( maxdist2leave( stop_jump ) < edge_length ) ) maxdist2leave( stop_jump ) = edge_length;
		} else { // a jump edge
			res2jumps_( it->start(), 1 ) = edge_nr;
			res2jumps_( it->start(), 2 ) = start_jump;
			res2jumps_( it->start(), 3 ) = 0;
			res2jumps_( it->start(), 4 ) = stop_jump;
			res2jumps_( it->start(), 5 ) = 1;

			res2jumps_( it->stop(), 1 ) = edge_nr;
			res2jumps_( it->stop(), 2 ) = start_jump;
			res2jumps_( it->stop(), 3 ) = 1;
			res2jumps_( it->stop(), 4 ) = stop_jump;
			res2jumps_( it->stop(), 5 ) = 0;
		} //if jump
	} //for edges
	// find maximum distance:
	// it will be from leave to leave ... go throu all combinations and get distance...
	max_dist_ = 0;
	for ( Size i1 = 1; i1 <= leaves.size() ; i1++ ) {
		for ( Size i2 = i1; i2 <= leaves.size() ; i2++ ) {
			max_dist_ = std::max( max_dist_, dist( leaves[ i1 ],leaves[ i2 ] ) );
		}
	}

	// already no maximum distance to leaves in fold-tree
	// now go throu all combination of
	//  for ( int seqpos = 1; seqpos <= f.nres(); seqpos++ ) {
	//   if ( res2jumps_( seqpos, 1) == -1 ) { //haven't found this residue in any peptide edge
	//    res2jump2_( seqpos, 1 ) =
	//   }
	//  }
	//produce some debug output
	//   if ( tr.Trace.visible() ) {
	//     tr.Trace << " edge_nr   jump1    dist1   jump2    dist2 \n";
	//     for ( Size ii = 1; ii<=f.nres(); ii++ ) {
	//       for ( Size k = 1; k<=5; k++ ) {
	//  tr.Trace << I(3, res2jumps_( ii, k ) );
	//       }
	//       tr.Trace << std::endl;
	//     }

	//     tr.Trace << "\n total distance list \n";
	//     for ( Size ii = 1; ii<=f.nres(); ii++ ) {
	//       for ( Size jj = 1; jj<=f.nres(); jj++ ) {
	//     tr.Trace << "(" << ii << "," << jj << ") " << dist( ii, jj ) << std::endl;
	//       }
	//     }
	// }
}

/// @detail distance between two residues
/// with the help of our pre-computed data
/// this only requires comparison of 4 possible pathways:
/// go via upstream/downstream jump-residue for pos1/pos2
Size
ShortestPathInFoldTree::dist( Size pos1, Size pos2 ) const {
	if ( simple_fold_tree_ ) {
		return  std::abs( (int) pos1 - (int) pos2 );
	}

	// on same edge ?
	if ( res2jumps_( pos1, 1 /*asks for edge-nr*/ ) == res2jumps_( pos2, 1 /*asks for edge-nr*/ ) ) {
		return std::abs( (int) pos1 - (int) pos2 );
	};

	// compute for possibilities and take smallest
	int const inf( 12345678 ); //assumption: fewer than 12 million nodes in the graph.

	Size min_dist = inf;

	// check 2x2 possibilities of up-/downstream jump-residues for pos1/pos2
	for ( Size ii=1; ii<=2 ; ++ii ) { //choose jump-node for pos1
		if ( res2jumps_( pos1, 2*ii ) < 0 ) continue; // is not a jump-residue
		for ( Size jj=1; jj<=2; ++jj ) { //choose jump-node for pos2
			if ( res2jumps_( pos2, 2*jj ) < 0 ) continue; // is not a jump-residue
			//   tr.Trace << "dist1: " << res2jumps_( pos1, 2*ii+1) << " " << pos1 << " --> " << res2jumps_( pos1, 2*ii ) << "\n";
			//   tr.Trace << "dist2: " << res2jumps_( pos2, 2*jj+1) << " " << pos2 << " --> " << res2jumps_( pos2, 2*jj ) << "\n";
			Size dist = res2jumps_( pos1, 2*ii+1) + res2jumps_( pos2, 2*jj+1 )
				+ node_dist_(  res2jumps_( pos1, 2*ii ),  res2jumps_( pos2, 2*jj ) );
			if ( dist < min_dist ) min_dist = dist;
		}
	}
	debug_assert ( min_dist <= nres_ );
	return min_dist;
}


} //kinematics
} //core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::kinematics::ShortestPathInFoldTree::ShortestPathInFoldTree() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::ShortestPathInFoldTree::save( Archive & arc ) const {
	arc( CEREAL_NVP( jump_res_ ) ); // std::map<Size, Size>
	arc( CEREAL_NVP( node_dist_ ) ); // ObjexxFCL::FArray2D_int
	arc( CEREAL_NVP( res2jumps_ ) ); // ObjexxFCL::FArray2D_int
	arc( CEREAL_NVP( nres_ ) ); // Size
	arc( CEREAL_NVP( simple_fold_tree_ ) ); // _Bool
	arc( CEREAL_NVP( max_dist_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::ShortestPathInFoldTree::load( Archive & arc ) {
	arc( jump_res_ ); // std::map<Size, Size>
	arc( node_dist_ ); // ObjexxFCL::FArray2D_int
	arc( res2jumps_ ); // ObjexxFCL::FArray2D_int
	arc( nres_ ); // Size
	arc( simple_fold_tree_ ); // _Bool
	arc( max_dist_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::ShortestPathInFoldTree );
CEREAL_REGISTER_TYPE( core::kinematics::ShortestPathInFoldTree )

CEREAL_REGISTER_DYNAMIC_INIT( core_kinematics_ShortestPathInFoldTree )
#endif // SERIALIZATION
