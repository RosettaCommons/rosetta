// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/graph/find_neighbors.cc
/// @brief  Sets up the residue neighbor information
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
///
/// @remarks Thanks to Will Sheffler for his ideas on refining this and extending it to atom neighbors
/// @remarks Adapting libRosetta code for generalized neighbor detection

#ifndef INCLUDED_core_conformation_find_neighbors_hh
#define INCLUDED_core_conformation_find_neighbors_hh

// Package Headers
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/find_neighbors.fwd.hh>
#include <core/types.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>


// Numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzTriple.hh>
#include <numeric/xyzVector.hh>

#include <numeric/geometry/hashing/xyzStripeHashWithMeta.hh>

// ObjexxFCL headers
//#include <ObjexxFCL/KeyFArray1D.hh>
//#include <ObjexxFCL/KeyFArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// Utility headers
//#include <utility/pointer/access_ptr.hh>

// boost headers
#include <boost/unordered_map.hpp>

// C++ headers
#include <utility/assert.hh>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <vector>

#include <utility/vector1.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/prof.hh>


namespace core {
namespace conformation {


// add the following to top of file:
/// @brief uses default boost::hash combine to hash Cubes
struct DefaultCubeHash : std::unary_function< CubeKey, std::size_t > {
	// use std::size_t instead of core::Size just to be
	// consistent with boost::hash types

	/// @brief return hash value given CubeKey
	std::size_t operator()( CubeKey const & key ) const {
		std::size_t seed = 0;
		boost::hash_combine( seed, key.x() );
		boost::hash_combine( seed, key.y() );
		boost::hash_combine( seed, key.z() );
		return seed;
	}
};


// Find neighbors and place them in a graph
template <class Vertex, class Edge>
void
find_neighbors(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	Strategy strategy = AUTOMATIC
)
{
	// the naive and the octree strategies are equally fast
	// when detecting neighbors for 150 points.
	core::Size const N_POINTS_BREAK_EVEN = 150;

	// PROF_START( basic::TEST2 );

	if ( strategy == STRIPEHASH || ( strategy == AUTOMATIC && basic::options::option[ basic::options::OptionKeys::score::find_neighbors_stripehash ]() ) ) {
		find_neighbors_stripe( point_graph, neighbor_cutoff );
	} else if ( strategy == THREEDGRID || ( strategy == AUTOMATIC && basic::options::option[ basic::options::OptionKeys::score::find_neighbors_3dgrid ]() ) ) {
		find_neighbors_3dgrid( point_graph, neighbor_cutoff );
	} else if ( strategy == NAIVE || point_graph->num_vertices() < N_POINTS_BREAK_EVEN ) {
		find_neighbors_naive( point_graph, neighbor_cutoff );
	} else { // Use automatic or an octree strategy
		find_neighbors_octree( point_graph, neighbor_cutoff, strategy );
	}
	// PROF_STOP( basic::TEST2 );

}

template <class Vertex, class Edge>
void
find_neighbors_naive(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff
)
{
	// Constants
	core::Size const n_points( point_graph->num_vertices() );
	if ( n_points == 0 ) return;
	core::Real neighbor_cutoff_sq = neighbor_cutoff * neighbor_cutoff;

	// Exclusion checks
	if ( n_points <= 1 ) return; // Nothing to do

	// Naive method: O( R^2 ) for R residues but faster for small, compact conformations
	for ( core::Size ii = 1; ii <= n_points; ++ii ) {
		PointPosition const & ii_pos( point_graph->get_vertex(ii).data().xyz() );
		for ( core::Size jj = ii + 1; jj <= n_points; ++jj ) {
			core::Real const d_sq( ii_pos.distance_squared( point_graph->get_vertex(jj).data().xyz() ) ); // Using member version of distance_squared to work around GCC 3.4.5 ADL bug
			if ( d_sq <= neighbor_cutoff_sq ) {
				// Add neighbor link
				point_graph->add_edge( ii, jj, Edge( d_sq ) );
			}
		}
	}
}

template <class Vertex, class Edge>
struct AddEdgeVisitor{
	graph::UpperEdgeGraph<Vertex, Edge> & point_graph;
	AddEdgeVisitor(
		graph::UpperEdgeGraph<Vertex, Edge> & point_graph_in
	) : point_graph(point_graph_in) {}
	void visit(
		numeric::xyzVector<Real> const & /*v*/, Real  const & vm,
		numeric::xyzVector<Real> const & /*c*/, Real  const & cm, Real const & d_sq
	){
		if ( vm < cm ) point_graph.add_edge( Size(vm), Size(cm), Edge( d_sq ) );
	}
};

template <class Vertex, class Edge>
void
find_neighbors_stripe(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff
){
	core::Size const n_points( point_graph->num_vertices() );
	if ( n_points <= 1 ) return; // Nothing to do
	std::cout << n_points << std::endl;
	utility::vector1<PointPosition> pts(n_points);
	for ( core::Size ii = 1; ii <= n_points; ++ii ) pts[ii] = point_graph->get_vertex(ii).data().xyz();
	utility::vector1<core::Real> dummy;
	numeric::geometry::hashing::xyzStripeHashWithMeta<Real> hash(neighbor_cutoff,pts,dummy);
	AddEdgeVisitor<Vertex,Edge> visitor(*point_graph);
	for ( core::Size ii = 1; ii <= n_points; ++ii ) {
		hash.visit(pts[ii],ii,visitor);
	}
}

/// @brief Finds the residue neighbors efficiently using an octree-like spatial sort
/// @remarks
///  @li The "octree" algorithm isn't a real octree since we don't want/need the tree structure:
///       we are only interested in one distance cutoff criterion
///  @li The octree algorithm is O( R log R ) for R residues vs. O( R^2 ) for the naive algorithm
///  @li The octree algorithm seems to be faster for R >= 150 or so: this may come down as it gets refined
///  @li This is an initial implementation of a more scalable neighbor detection algorithm: it will be
///       further tuned as more tests cases are run.  This type of algorithm may be of greater benefit
///       for atom neighbor detection because the numbers involved are greater and the neighbor regions
///       encompass a smaller fraction of the other atoms in typical structures.
///  @li The spatial sorting is essentially an octree method but with the important distinction that
///       is doesn't have full a full tree structure: That would change the log R search complexity to
///       log C where C is the number of cubes in the whole bounding box, which could be as much O( R^3 )
///       ( O( log C ) == O( log R ) would still be true but the constant would be 3x worse).  Instead
///       we let the map build a tree of just the active cubes, keeping the depth minimal.  If we needed
///       to access the parent meta-cubes we would need the full octree, but we don't here.
///  @li The use of std::map to hold the octree should be compared against hash maps and other methods.
template <class Vertex, class Edge>
void
find_neighbors_octree(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	Strategy strategy
)
{
	//using namespace residue_neighbor_strategy;

	using numeric::min;
	using numeric::max;
	using numeric::square;


	// Types
	typedef  numeric::xyzTriple< core::Size >  CubeDim; // Cube dimensions
	//typedef  numeric::xyzTriple< core::Size >  CubeKey; // Cube index-triple key
	typedef  utility::vector1< PointPosition > Points;
	typedef  std::vector< core::Size >  PointIDs;
	//typedef  std::map< CubeKey, PointIDs, std::less< CubeKey >, boost::pool_allocator< std::pair< CubeKey const, PointIDs > >  >  Cubes;
	typedef  std::map< CubeKey, PointIDs >  Cubes;

	/// Andrew Ban's boost version
	// within find_neighbors() change the Cubes typedef to:
	///typedef boost::unordered_map< CubeKey, PointIDs, DefaultCubeHash > Cubes; // uncomment to use boost version

	// Constants
	core::Size const n_points( point_graph->num_vertices() );

	core::Real neighbor_cutoff_sq( neighbor_cutoff*neighbor_cutoff);

	//local copy
	Points points( n_points );
	for ( core::Size ii = 1; ii <= n_points; ++ii ) { points[ ii ] = point_graph->get_vertex( ii ).data().xyz(); }

	// Exclusion checks
	if ( n_points <= 1 ) return; // Nothing to do

	// Use automatic or an octree strategy

	// Bounding box of residue cutoff positions
	PointPosition bbl( points[ 1 ] ), bbu( bbl ); // Lower and upper corners of bounding box
	for ( core::Size ii = 2; ii <= n_points; ++ii ) {
		bbl.min( points[ ii ] );
		bbu.max( points[ ii ] );
	}

	core::Size const epsilon_multiplier( 10 ); // Increase this if assert failures hit in finding a point's cube
	core::Real const epsilon( epsilon_multiplier * std::numeric_limits< core::Real >::epsilon() );
	bbl -= epsilon; // Expand bounding box to assure all points get assigned cubes in it
	bbu += epsilon;

	// Set cube size and dimensions within bounding box
	core::Size const side_factor( 1 ); // 1 factor => Check <= 27 adjacent cubes // 2 factor => Check <= 8 adjacent cubes
	// Might gain some speed by replacing max_residue_pair_cutoff below with the max cutoff for pairs present
	core::Real const side( side_factor * neighbor_cutoff );
	debug_assert( side > core::Real( 0 ) );
	core::Real const side_inv( core::Real( 1 ) / side );
	CubeDim const cube_dim( // Cube dimensions
		core::Size( std::ceil( ( bbu.x() - bbl.x() ) * side_inv ) ),             // Test that ceil values == core::Size values
		core::Size( std::ceil( ( bbu.y() - bbl.y() ) * side_inv ) ),
		core::Size( std::ceil( ( bbu.z() - bbl.z() ) * side_inv ) )
	);
	// We rounded up the number of cubes in each dimension
	// We use cubes of exactly side x side x side dimensions
	// We treat the (0,0,0) cube as touching bbl at its low corner
	// The "highest" cube generally extends beyond bbu
	// We call this the expanded bounding box

	// Number of potential cubes in expanded bounding box (we don't create them all)
	core::Size const n_cube( cube_dim.x() * cube_dim.y() * cube_dim.z() );

	// Find upper Residue neighbors of each residue
	if ( ( n_cube < core::Size( 27 ) ) && ( strategy < OCTREE ) ) { // Naive strategy //! Tune the n_cube threshold based on more real-world trials

		// Naive method: O( R^2 ) for R residues but faster for small, compact conformations
		find_neighbors_naive<Vertex, Edge>( point_graph, neighbor_cutoff );

	} else { // Octree O( R log R ) strategy

		// Add residues to bounding box cube tree: Only cubes with residues are added
		Cubes cubes; /// STL MAP cubes

		// take a look at boost doc and implementation on more info to
		// set number of buckets/loadfactor, briefly e.g.:

		// Andrew Ban's Boost version
		// init at least 128 buckets, boost::unordered_map uses
		// prime number sequence for bucket growth
		//Cubes cubes( 128 ); // uncomment to use boost version
		// change max load factor to 4.0
		//cubes.max_load_factor( 4.0 ); // uncomment to use boost version

		for ( core::Size i = 1; i <= n_points; ++i ) {
			//AminoAcid & res( p[ i ] );
			PointPosition const pp( points[ i ]); //( res.neighbor_graph_position() );

			// Find the residue's cube: Cube coords are indexed from 0 to cube_dim -1
			CubeKey const cube_key(
				core::Size( ( pp.x() - bbl.x() ) * side_inv ),
				core::Size( ( pp.y() - bbl.y() ) * side_inv ),
				core::Size( ( pp.z() - bbl.z() ) * side_inv )
			);

			// Check that it is within the expanded bounding box
			debug_assert( cube_key.x() < cube_dim.x() );
			debug_assert( cube_key.y() < cube_dim.y() );
			debug_assert( cube_key.z() < cube_dim.z() );

			// Add the point's position to the cube's collection
			cubes[ cube_key ].push_back( i ); // Creates the cube if it doesn't exist yet
		}

		// Find upper neighbors
		core::Real const D_ZERO( 0 );
		for ( core::Size i = 1; i <= n_points; ++i ) {
			//AminoAcid & res( p[ i ] );
			PointPosition const pp( points[ i ]); //( res.neighbor_graph_position() );
			//AminoAcidKey const & res_key( res.cat_key() );
			//core::Size const res_number( res.number() );

			// Find the residue's cube indexes
			core::Size const icx( core::Size( ( pp.x() - bbl.x() ) * side_inv ) );
			core::Size const icy( core::Size( ( pp.y() - bbl.y() ) * side_inv ) );
			core::Size const icz( core::Size( ( pp.z() - bbl.z() ) * side_inv ) );

			// Get cube-relative position (for fast cube exclusion tests)
			core::Real const cx( pp.x() - ( bbl.x() + ( icx * side ) ) );
			core::Real const cy( pp.y() - ( bbl.y() + ( icy * side ) ) );
			core::Real const cz( pp.z() - ( bbl.z() + ( icz * side ) ) );

			// Check its cube and adjacent cubes (<= all 27 of them with side_factor==1)
			for ( core::Size ix = max( icx, core::Size( 1 ) ) - 1,  ixe = min( icx + 1, cube_dim.x() - 1 ); ix <= ixe; ++ix ) {
				for ( core::Size iy = max( icy, core::Size( 1 ) ) - 1, iye = min( icy + 1, cube_dim.y() - 1 ); iy <= iye; ++iy ) {
					for ( core::Size iz = max( icz, core::Size( 1 ) ) - 1, ize = min( icz + 1, cube_dim.z() - 1 ); iz <= ize; ++iz ) {
						Cubes::iterator const ic( cubes.find( CubeKey( ix,iy, iz ) ) );
						if ( ic != cubes.end() ) { // Cube exists
							if ( // This test gave a ~10% speedup in trials
									( ix != icx ? square( cx - ( ix > icx ? side : D_ZERO ) ) : D_ZERO ) +
									( iy != icy ? square( cy - ( iy > icy ? side : D_ZERO ) ) : D_ZERO ) +
									( iz != icz ? square( cz - ( iz > icz ? side : D_ZERO ) ) : D_ZERO )
									<= neighbor_cutoff_sq ) {
								// Max cutoff sphere intersects this cube so check each residue in it
								for ( PointIDs::iterator ia = ic->second.begin(), iae = ic->second.end(); ia != iae; ++ia ) {
									core::Size const j( *ia );
									if ( i < j ) { // It is an upper neighbor
										core::Real const d_sq( pp.distance_squared( points[ j ] ) );
										if ( d_sq <= neighbor_cutoff_sq ) {
											point_graph->add_edge( i, j, Edge( d_sq ) );
										}
										//if ( d_sq < residue_neighbor_count_cutoff_sq ) { // Add to neighbor counts
										// res.increment_n_neighbor();
										// resu.increment_n_neighbor();
										//}
									}
								}
							}
						}
					}
				}
			}

		}
	}
}

/// @brief Create a 3D grid of points.  O(N^3).  For "spherical" conformations, Theta(N).  Speeds neighbor detection
/// in abinitio by a factor of 2.  Definition: Spherical = span of x,y and z all O(N**1/3).  Note circularity.
/// Adendum: if the 3D grid used a list of point indices instead of a vector, then this would be Theta(N) for
/// spherical conformations; however, with a vector, this is O(NlgN).  With the additional assumption that
/// each cube contains O(1) points, then this implementation is O(N).  Such an assumption is unneccessary
/// in the list implementation.
/// @details Shameless code duplication below based on Stuart's
/// stl-map-based-neighbor-detection code.  Note that the FArray3D is an index-from-1
/// data structure.
template <class Vertex, class Edge>
void
find_neighbors_3dgrid(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff
)
{
	using numeric::min;
	using numeric::max;
	using numeric::square;


	// Types
	typedef  numeric::xyzTriple< core::Size >  CubeDim; // Cube dimensions
	typedef  numeric::xyzTriple< core::Size >  CubeKey; // Cube index-triple key
	typedef  utility::vector1< PointPosition > Points;
	typedef  std::vector< core::Size >  PointIDs;
	typedef ObjexxFCL::FArray3D< PointIDs > Cubes; // The 3D array that will be indexed into.  Indexed from 1, not 0.

	// Constants
	core::Size const n_points( point_graph->num_vertices() );

	core::Real neighbor_cutoff_sq( neighbor_cutoff*neighbor_cutoff);

	//local copy
	Points points( n_points );
	for ( core::Size ii = 1; ii <= n_points; ++ii ) { points[ ii ] = point_graph->get_vertex( ii ).data().xyz(); }

	// Exclusion checks
	if ( n_points <= 1 ) return; // Nothing to do

	// Use automatic or an octree strategy

	// Bounding box of residue cutoff positions
	PointPosition bbl( points[ 1 ] ), bbu( bbl ); // Lower and upper corners of bounding box
	for ( core::Size ii = 2; ii <= n_points; ++ii ) {
		bbl.min( points[ ii ] );
		bbu.max( points[ ii ] );
	}

	core::Size const epsilon_multiplier( 10 ); // Increase this if assert failures hit in finding a point's cube
	core::Real const epsilon( epsilon_multiplier * std::numeric_limits< core::Real >::epsilon() );
	bbl -= epsilon; // Expand bounding box to assure all points get assigned cubes in it
	bbu += epsilon;

	// Set cube size and dimensions within bounding box
	core::Size const side_factor( 1 ); // 1 factor => Check <= 27 adjacent cubes // 2 factor => Check <= 8 adjacent cubes
	// Might gain some speed by replacing max_residue_pair_cutoff below with the max cutoff for pairs present
	core::Real const side( side_factor * neighbor_cutoff );
	debug_assert( side > core::Real( 0 ) );
	core::Real const side_inv( core::Real( 1 ) / side );
	CubeDim const cube_dim( // Cube dimensions
		core::Size( std::ceil( ( bbu.x() - bbl.x() ) * side_inv ) ),             // Test that ceil values == core::Size values
		core::Size( std::ceil( ( bbu.y() - bbl.y() ) * side_inv ) ),
		core::Size( std::ceil( ( bbu.z() - bbl.z() ) * side_inv ) )
	);
	// We rounded up the number of cubes in each dimension
	// We use cubes of exactly side x side x side dimensions
	// We treat the (1,1,1) cube as touching bbl at its low corner
	// The "highest" cube generally extends beyond bbu
	// We call this the expanded bounding box

	// Add residues to bounding box cube tree: Only cubes with residues are added

	/// NOT THREAD SAFE -- Static variable below would avoid allocation and deallocation costs of the 3D array.
	/// Does not seem to offer any appreciable speed advantages.
	//static Cubes cubes;
	///cubes.dimension( cube_dim.x(), cube_dim.y(), cube_dim.z() );

	/// keep track of the non-empty voxels so we can delete them later.  The cubes array must be empty
	/// at the beginning of neighbor detection.
	//utility::vector1< core::Size > nonempty_cube_indices;

	/// Thread safe version; potentially more expensive than the non-thread-safe version,
	/// but has not proven so in experimentation.
	Cubes cubes( cube_dim.x(), cube_dim.y(), cube_dim.z() );

	for ( core::Size i = 1; i <= n_points; ++i ) {
		PointPosition const pp( points[ i ]);

		// Find the residue's cube: Cube coords are indexed from 1 to cube_dim.
		CubeKey const cube_key(
			core::Size( ( pp.x() - bbl.x() ) * side_inv ) + 1,
			core::Size( ( pp.y() - bbl.y() ) * side_inv ) + 1,
			core::Size( ( pp.z() - bbl.z() ) * side_inv ) + 1
		);

		// Check that it is within the expanded bounding box
		debug_assert( cube_key.x() <= cube_dim.x() );
		debug_assert( cube_key.y() <= cube_dim.y() );
		debug_assert( cube_key.z() <= cube_dim.z() );

		// Add the point's position to the cube's collection
		//cubes[ cube_key ].push_back( i ); // Creates the cube if it doesn't exist yet
		core::Size i_index = cubes.index( cube_key.x(), cube_key.y(), cube_key.z() );
		if ( cubes[ i_index ].size() == 0 ) {
			/// In the statically-allocated version, the cubes object must be emptied
			/// at the conclusion of neighbor detection; keep track of those cubes which
			/// have some entry to avoid the expense of traversing the whole cubes object
			/// later.
			//nonempty_cube_indices.push_back( i_index );

			/// In the thread-safe version, guess that any cube with 1 point contained inside
			/// it will likely contain several.  Allocate a bit of space now.  O(NlgN) if the points
			/// are not well distributed.
			cubes[ i_index ].reserve( 10 );
		}
		cubes[ i_index ].push_back( i );
		///std::cout << "Cube " << i_index << " for residue " << i << " at coordinate: (" << pp.x() << "," << pp.y() <<"," << pp.z() << ")" << std::endl;
	}

	// Find upper neighbors
	//core::Real const D_ZERO( 0 );
	for ( core::Size i = 1; i <= n_points; ++i ) {
		//AminoAcid & res( p[ i ] );
		PointPosition const pp( points[ i ]);

		// Find the residue's cube indexes
		core::Size const icx( core::Size( ( pp.x() - bbl.x() ) * side_inv ) + 1 );
		core::Size const icy( core::Size( ( pp.y() - bbl.y() ) * side_inv ) + 1 );
		core::Size const icz( core::Size( ( pp.z() - bbl.z() ) * side_inv ) + 1 );

		// Check its cube and adjacent cubes (<= all 27 of them with side_factor==1)
		for ( core::Size ix = max( icx, core::Size( 2 ) ) - 1,  ixe = min( icx + 1, cube_dim.x() ); ix <= ixe; ++ix ) {
			for ( core::Size iy = max( icy, core::Size( 2 ) ) - 1, iye = min( icy + 1, cube_dim.y() ); iy <= iye; ++iy ) {
				for ( core::Size iz = max( icz, core::Size( 2 ) ) - 1, ize = min( icz + 1, cube_dim.z() ); iz <= ize; ++iz ) {

					//Cubes::iterator const ic( cubes.find( CubeKey( ix,iy, iz ) ) );
					core::Size cube_index = cubes.index( ix, iy, iz );

					///std::cout << "Searching for neighbors of point " << i << " in cube [" << ix << "," << iy << "," << iz << ") index: " << cube_index << std::endl;

					if ( cubes[ cube_index ].size() != 0 ) { // Cube exists
						for ( PointIDs::iterator ia = cubes[ cube_index ].begin(), iae = cubes[ cube_index ].end(); ia != iae; ++ia ) {
							core::Size const j( *ia );
							///std::cout << "point " << j << " found " << std::endl;
							if ( i < j ) { // It is an upper neighbor
								core::Real const d_sq( pp.distance_squared( points[ j ] ) );
								if ( d_sq <= neighbor_cutoff_sq ) {
									point_graph->add_edge( i, j, Edge( d_sq ) );
								}
								//if ( d_sq < residue_neighbor_count_cutoff_sq ) { // Add to neighbor counts
								// res.increment_n_neighbor();
								// resu.increment_n_neighbor();
								//}
							}
						}
					}
				}
			}
		}

	}

	/// Only necessary in the non-thread-safe version
	/// before returning, empty the cubes array so it's ready for the next round
	//for ( core::Size ii = 1; ii <= nonempty_cube_indices.size(); ++ii ) {
	// cubes[ nonempty_cube_indices[ ii ] ].clear();
	//}
}

template <class Vertex, class Edge>
void
find_neighbors_restricted(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const & residue_selection,
	Strategy strategy = AUTOMATIC
)
{
	// the naive and the octree strategies are equally fast
	// when detecting neighbors for 150 points.
	core::Size const N_POINTS_BREAK_EVEN = 150;

	// PROF_START( basic::TEST2 );

	if ( strategy == THREEDGRID || ( strategy == AUTOMATIC && basic::options::option[ basic::options::OptionKeys::score::find_neighbors_3dgrid ] ) ) {
		find_neighbors_3dgrid_restricted<Vertex,Edge>( point_graph, neighbor_cutoff, residue_selection );
	} else if ( strategy == NAIVE || point_graph->num_vertices() < N_POINTS_BREAK_EVEN ) {
		find_neighbors_naive_restricted<Vertex,Edge>( point_graph, neighbor_cutoff, residue_selection );
	} else { // Use automatic or an octree strategy
		find_neighbors_octree_restricted<Vertex,Edge>( point_graph, neighbor_cutoff, residue_selection, strategy );
	}
	// PROF_STOP( basic::TEST2 );
}

template <class Vertex, class Edge>
void
find_neighbors_naive_restricted(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const & residue_selection
)
{
	// Constants
	core::Size const n_points( point_graph->num_vertices() );
	if ( n_points == 0 ) return;
	core::Real neighbor_cutoff_sq = neighbor_cutoff * neighbor_cutoff;

	// Exclusion checks
	if ( n_points <= 1 ) return; // Nothing to do

	// Naive method: O( R^2 ) for R residues but faster for small, compact conformations
	for ( core::Size ii = 1; ii <= n_points; ++ii ) {
		if ( !residue_selection[ ii ] ) continue;
		PointPosition const & ii_pos( point_graph->get_vertex(ii).data().xyz() );
		for ( core::Size jj = 1; jj <= n_points; ++jj ) {
			if ( jj <= ii && residue_selection[ jj ] ) continue;
			core::Real const d_sq( ii_pos.distance_squared( point_graph->get_vertex(jj).data().xyz() ) ); // Using member version of distance_squared to work around GCC 3.4.5 ADL bug
			if ( d_sq <= neighbor_cutoff_sq ) {
				// Add neighbor link
				point_graph->add_edge( ii, jj, Edge( d_sq ) );
			}
		}
	}
}

/// @brief Finds the residue neighbors efficiently using an octree-like spatial sort
template <class Vertex, class Edge>
void
find_neighbors_octree_restricted(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const & residue_selection,
	Strategy strategy
)
{
	//using namespace residue_neighbor_strategy;

	using numeric::min;
	using numeric::max;
	using numeric::square;


	// Types
	typedef  numeric::xyzTriple< core::Size >  CubeDim; // Cube dimensions
	//typedef  numeric::xyzTriple< core::Size >  CubeKey; // Cube index-triple key
	typedef  utility::vector1< PointPosition > Points;
	typedef  std::vector< core::Size >  PointIDs;
	//typedef  std::map< CubeKey, PointIDs, std::less< CubeKey >, boost::pool_allocator< std::pair< CubeKey const, PointIDs > >  >  Cubes;
	typedef  std::map< CubeKey, PointIDs >  Cubes;

	/// Andrew Ban's boost version
	// within find_neighbors() change the Cubes typedef to:
	///typedef boost::unordered_map< CubeKey, PointIDs, DefaultCubeHash > Cubes; // uncomment to use boost version

	// Constants
	core::Size const n_points( point_graph->num_vertices() );

	core::Real neighbor_cutoff_sq( neighbor_cutoff*neighbor_cutoff);

	//local copy
	Points points( n_points );
	for ( core::Size ii = 1; ii <= n_points; ++ii ) { points[ ii ] = point_graph->get_vertex( ii ).data().xyz(); }

	// Exclusion checks
	if ( n_points <= 1 ) return; // Nothing to do

	// Use automatic or an octree strategy

	// Bounding box of residue cutoff positions
	PointPosition bbl( points[ 1 ] ), bbu( bbl ); // Lower and upper corners of bounding box
	for ( core::Size ii = 2; ii <= n_points; ++ii ) {
		bbl.min( points[ ii ] );
		bbu.max( points[ ii ] );
	}

	core::Size const epsilon_multiplier( 10 ); // Increase this if assert failures hit in finding a point's cube
	core::Real const epsilon( epsilon_multiplier * std::numeric_limits< core::Real >::epsilon() );
	bbl -= epsilon; // Expand bounding box to assure all points get assigned cubes in it
	bbu += epsilon;

	// Set cube size and dimensions within bounding box
	core::Size const side_factor( 1 ); // 1 factor => Check <= 27 adjacent cubes // 2 factor => Check <= 8 adjacent cubes
	// Might gain some speed by replacing max_residue_pair_cutoff below with the max cutoff for pairs present
	core::Real const side( side_factor * neighbor_cutoff );
	debug_assert( side > core::Real( 0 ) );
	core::Real const side_inv( core::Real( 1 ) / side );
	CubeDim const cube_dim( // Cube dimensions
		core::Size( std::ceil( ( bbu.x() - bbl.x() ) * side_inv ) ),             // Test that ceil values == core::Size values
		core::Size( std::ceil( ( bbu.y() - bbl.y() ) * side_inv ) ),
		core::Size( std::ceil( ( bbu.z() - bbl.z() ) * side_inv ) )
	);
	// We rounded up the number of cubes in each dimension
	// We use cubes of exactly side x side x side dimensions
	// We treat the (0,0,0) cube as touching bbl at its low corner
	// The "highest" cube generally extends beyond bbu
	// We call this the expanded bounding box

	// Number of potential cubes in expanded bounding box (we don't create them all)
	core::Size const n_cube( cube_dim.x() * cube_dim.y() * cube_dim.z() );

	// Find upper Residue neighbors of each residue
	if ( ( n_cube < core::Size( 27 ) ) && ( strategy < OCTREE ) ) { // Naive strategy //! Tune the n_cube threshold based on more real-world trials

		// Naive method: O( R^2 ) for R residues but faster for small, compact conformations
		find_neighbors_naive_restricted<Vertex, Edge>( point_graph, neighbor_cutoff,residue_selection );

	} else { // Octree O( R log R ) strategy

		// Add residues to bounding box cube tree: Only cubes with residues are added
		Cubes cubes; /// STL MAP cubes

		// take a look at boost doc and implementation on more info to
		// set number of buckets/loadfactor, briefly e.g.:

		// Andrew Ban's Boost version
		// init at least 128 buckets, boost::unordered_map uses
		// prime number sequence for bucket growth
		//Cubes cubes( 128 ); // uncomment to use boost version
		// change max load factor to 4.0
		//cubes.max_load_factor( 4.0 ); // uncomment to use boost version

		for ( core::Size i = 1; i <= n_points; ++i ) {
			//AminoAcid & res( p[ i ] );
			PointPosition const pp( points[ i ]); //( res.neighbor_graph_position() );

			// Find the residue's cube: Cube coords are indexed from 0 to cube_dim -1
			CubeKey const cube_key(
				core::Size( ( pp.x() - bbl.x() ) * side_inv ),
				core::Size( ( pp.y() - bbl.y() ) * side_inv ),
				core::Size( ( pp.z() - bbl.z() ) * side_inv )
			);

			// Check that it is within the expanded bounding box
			debug_assert( cube_key.x() < cube_dim.x() );
			debug_assert( cube_key.y() < cube_dim.y() );
			debug_assert( cube_key.z() < cube_dim.z() );

			// Add the point's position to the cube's collection
			cubes[ cube_key ].push_back( i ); // Creates the cube if it doesn't exist yet
		}

		// Find upper neighbors
		core::Real const D_ZERO( 0 );
		for ( core::Size i = 1; i <= n_points; ++i ) {
			if ( !residue_selection[ i ] ) continue;
			//AminoAcid & res( p[ i ] );
			PointPosition const pp( points[ i ]); //( res.neighbor_graph_position() );
			//AminoAcidKey const & res_key( res.cat_key() );
			//core::Size const res_number( res.number() );

			// Find the residue's cube indexes
			core::Size const icx( core::Size( ( pp.x() - bbl.x() ) * side_inv ) );
			core::Size const icy( core::Size( ( pp.y() - bbl.y() ) * side_inv ) );
			core::Size const icz( core::Size( ( pp.z() - bbl.z() ) * side_inv ) );

			// Get cube-relative position (for fast cube exclusion tests)
			core::Real const cx( pp.x() - ( bbl.x() + ( icx * side ) ) );
			core::Real const cy( pp.y() - ( bbl.y() + ( icy * side ) ) );
			core::Real const cz( pp.z() - ( bbl.z() + ( icz * side ) ) );

			// Check its cube and adjacent cubes (<= all 27 of them with side_factor==1)
			for ( core::Size ix = max( icx, core::Size( 1 ) ) - 1,  ixe = min( icx + 1, cube_dim.x() - 1 ); ix <= ixe; ++ix ) {
				for ( core::Size iy = max( icy, core::Size( 1 ) ) - 1, iye = min( icy + 1, cube_dim.y() - 1 ); iy <= iye; ++iy ) {
					for ( core::Size iz = max( icz, core::Size( 1 ) ) - 1, ize = min( icz + 1, cube_dim.z() - 1 ); iz <= ize; ++iz ) {
						Cubes::iterator const ic( cubes.find( CubeKey( ix,iy, iz ) ) );
						if ( ic != cubes.end() ) { // Cube exists
							if ( // This test gave a ~10% speedup in trials
									( ix != icx ? square( cx - ( ix > icx ? side : D_ZERO ) ) : D_ZERO ) +
									( iy != icy ? square( cy - ( iy > icy ? side : D_ZERO ) ) : D_ZERO ) +
									( iz != icz ? square( cz - ( iz > icz ? side : D_ZERO ) ) : D_ZERO )
									<= neighbor_cutoff_sq ) {
								// Max cutoff sphere intersects this cube so check each residue in it
								for ( PointIDs::iterator ia = ic->second.begin(), iae = ic->second.end(); ia != iae; ++ia ) {
									core::Size const j( *ia );
									if ( i < j || !residue_selection[ j ] ) { // It is an upper neighbor
										core::Real const d_sq( pp.distance_squared( points[ j ] ) );
										if ( d_sq <= neighbor_cutoff_sq ) {
											point_graph->add_edge( i, j, Edge( d_sq ) );
										}
										//if ( d_sq < residue_neighbor_count_cutoff_sq ) { // Add to neighbor counts
										// res.increment_n_neighbor();
										// resu.increment_n_neighbor();
										//}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

/// @brief Create a 3D grid of points.  O(N^3).  For "spherical" conformations, Theta(N).  Speeds neighbor detection
/// in abinitio by a factor of 2.  Definition: Spherical = span of x,y and z all O(N**1/3).  Note circularity.
/// Adendum: if the 3D grid used a list of point indices instead of a vector, then this would be Theta(N) for
/// spherical conformations; however, with a vector, this is O(NlgN).  With the additional assumption that
/// each cube contains O(1) points, then this implementation is O(N).  Such an assumption is unneccessary
/// in the list implementation.
template <class Vertex, class Edge>
void
find_neighbors_3dgrid_restricted(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const &  residue_selection
)
{
	using numeric::min;
	using numeric::max;
	using numeric::square;


	// Types
	typedef  numeric::xyzTriple< core::Size >  CubeDim; // Cube dimensions
	typedef  numeric::xyzTriple< core::Size >  CubeKey; // Cube index-triple key
	typedef  utility::vector1< PointPosition > Points;
	typedef  std::vector< core::Size >  PointIDs;
	typedef ObjexxFCL::FArray3D< PointIDs > Cubes; // The 3D array that will be indexed into.  Indexed from 1, not 0.

	// Constants
	core::Size const n_points( point_graph->num_vertices() );

	core::Real neighbor_cutoff_sq( neighbor_cutoff*neighbor_cutoff);

	//local copy
	Points points( n_points );
	for ( core::Size ii = 1; ii <= n_points; ++ii ) { points[ ii ] = point_graph->get_vertex( ii ).data().xyz(); }

	// Exclusion checks
	if ( n_points <= 1 ) return; // Nothing to do

	// Use automatic or an octree strategy

	// Bounding box of residue cutoff positions
	PointPosition bbl( points[ 1 ] ), bbu( bbl ); // Lower and upper corners of bounding box
	for ( core::Size ii = 2; ii <= n_points; ++ii ) {
		bbl.min( points[ ii ] );
		bbu.max( points[ ii ] );
	}

	core::Size const epsilon_multiplier( 10 ); // Increase this if assert failures hit in finding a point's cube
	core::Real const epsilon( epsilon_multiplier * std::numeric_limits< core::Real >::epsilon() );
	bbl -= epsilon; // Expand bounding box to assure all points get assigned cubes in it
	bbu += epsilon;

	// Set cube size and dimensions within bounding box
	core::Size const side_factor( 1 ); // 1 factor => Check <= 27 adjacent cubes // 2 factor => Check <= 8 adjacent cubes
	// Might gain some speed by replacing max_residue_pair_cutoff below with the max cutoff for pairs present
	core::Real const side( side_factor * neighbor_cutoff );
	debug_assert( side > core::Real( 0 ) );
	core::Real const side_inv( core::Real( 1 ) / side );
	CubeDim const cube_dim( // Cube dimensions
		core::Size( std::ceil( ( bbu.x() - bbl.x() ) * side_inv ) ),             // Test that ceil values == core::Size values
		core::Size( std::ceil( ( bbu.y() - bbl.y() ) * side_inv ) ),
		core::Size( std::ceil( ( bbu.z() - bbl.z() ) * side_inv ) )
	);
	// We rounded up the number of cubes in each dimension
	// We use cubes of exactly side x side x side dimensions
	// We treat the (1,1,1) cube as touching bbl at its low corner
	// The "highest" cube generally extends beyond bbu
	// We call this the expanded bounding box

	// Add residues to bounding box cube tree: Only cubes with residues are added

	/// NOT THREAD SAFE -- Static variable below would avoid allocation and deallocation costs of the 3D array.
	/// Does not seem to offer any appreciable speed advantages.
	//static Cubes cubes;
	///cubes.dimension( cube_dim.x(), cube_dim.y(), cube_dim.z() );

	/// keep track of the non-empty voxels so we can delete them later.  The cubes array must be empty
	/// at the beginning of neighbor detection.
	//utility::vector1< core::Size > nonempty_cube_indices;

	/// Thread safe version; potentially more expensive than the non-thread-safe version,
	/// but has not proven so in experimentation.
	Cubes cubes( cube_dim.x(), cube_dim.y(), cube_dim.z() );

	for ( core::Size i = 1; i <= n_points; ++i ) {
		PointPosition const pp( points[ i ]);

		// Find the residue's cube: Cube coords are indexed from 1 to cube_dim.
		CubeKey const cube_key(
			core::Size( ( pp.x() - bbl.x() ) * side_inv ) + 1,
			core::Size( ( pp.y() - bbl.y() ) * side_inv ) + 1,
			core::Size( ( pp.z() - bbl.z() ) * side_inv ) + 1
		);

		// Check that it is within the expanded bounding box
		debug_assert( cube_key.x() <= cube_dim.x() );
		debug_assert( cube_key.y() <= cube_dim.y() );
		debug_assert( cube_key.z() <= cube_dim.z() );

		// Add the point's position to the cube's collection
		//cubes[ cube_key ].push_back( i ); // Creates the cube if it doesn't exist yet
		core::Size i_index = cubes.index( cube_key.x(), cube_key.y(), cube_key.z() );
		if ( cubes[ i_index ].size() == 0 ) {
			/// In the statically-allocated version, the cubes object must be emptied
			/// at the conclusion of neighbor detection; keep track of those cubes which
			/// have some entry to avoid the expense of traversing the whole cubes object
			/// later.
			//nonempty_cube_indices.push_back( i_index );

			/// In the thread-safe version, guess that any cube with 1 point contained inside
			/// it will likely contain several.  Allocate a bit of space now.  O(NlgN) if the points
			/// are not well distributed.
			cubes[ i_index ].reserve( 10 );
		}
		cubes[ i_index ].push_back( i );
		///std::cout << "Cube " << i_index << " for residue " << i << " at coordinate: (" << pp.x() << "," << pp.y() <<"," << pp.z() << ")" << std::endl;
	}

	// Find upper neighbors
	//core::Real const D_ZERO( 0 );
	for ( core::Size i = 1; i <= n_points; ++i ) {
		//AminoAcid & res( p[ i ] );
		if ( !residue_selection[ i ] ) continue;
		PointPosition const pp( points[ i ]);

		// Find the residue's cube indexes
		core::Size const icx( core::Size( ( pp.x() - bbl.x() ) * side_inv ) + 1 );
		core::Size const icy( core::Size( ( pp.y() - bbl.y() ) * side_inv ) + 1 );
		core::Size const icz( core::Size( ( pp.z() - bbl.z() ) * side_inv ) + 1 );

		// Check its cube and adjacent cubes (<= all 27 of them with side_factor==1)
		for ( core::Size ix = max( icx, core::Size( 2 ) ) - 1,  ixe = min( icx + 1, cube_dim.x() ); ix <= ixe; ++ix ) {
			for ( core::Size iy = max( icy, core::Size( 2 ) ) - 1, iye = min( icy + 1, cube_dim.y() ); iy <= iye; ++iy ) {
				for ( core::Size iz = max( icz, core::Size( 2 ) ) - 1, ize = min( icz + 1, cube_dim.z() ); iz <= ize; ++iz ) {

					//Cubes::iterator const ic( cubes.find( CubeKey( ix,iy, iz ) ) );
					core::Size cube_index = cubes.index( ix, iy, iz );

					///std::cout << "Searching for neighbors of point " << i << " in cube [" << ix << "," << iy << "," << iz << ") index: " << cube_index << std::endl;

					if ( cubes[ cube_index ].size() != 0 ) { // Cube exists
						for ( PointIDs::iterator ia = cubes[ cube_index ].begin(), iae = cubes[ cube_index ].end(); ia != iae; ++ia ) {
							core::Size const j( *ia );
							///std::cout << "point " << j << " found " << std::endl;
							if ( i < j || !residue_selection[ j ] ) { // It is an upper neighbor
								core::Real const d_sq( pp.distance_squared( points[ j ] ) );
								if ( d_sq <= neighbor_cutoff_sq ) {
									point_graph->add_edge( i, j, Edge( d_sq ) );
								}
								//if ( d_sq < residue_neighbor_count_cutoff_sq ) { // Add to neighbor counts
								// res.increment_n_neighbor();
								// resu.increment_n_neighbor();
								//}
							}
						}
					}
				}
			}
		}
	}

	/// Only necessary in the non-thread-safe version
	/// before returning, empty the cubes array so it's ready for the next round
	//for ( core::Size ii = 1; ii <= nonempty_cube_indices.size(); ++ii ) {
	// cubes[ nonempty_cube_indices[ ii ] ].clear();
	//}

}

template <class Vertex, class Edge>
core::Size
get_nearest_neighbor(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Size node_id,
	core::Real neighbor_cutoff,
	Strategy strategy
)
{
	find_neighbors<Vertex,Edge>(point_graph,neighbor_cutoff,strategy);
	graph::UEVertex<Vertex,Edge> query_vertex =  point_graph->get_vertex(node_id);
	typename graph::UEVertex<Vertex,Edge>::UpperEdgeListIter query_it; //sometimes I really don't like C++...
	core::Real min_d_squared = 99999999.9;
	core::Size min_node_index = node_id;
	for ( query_it = query_vertex.upper_edge_list_begin(); query_it != query_vertex.upper_edge_list_end(); ++query_it ) {
		core::Real d_squared = query_it->data().dsq();
		if ( d_squared < min_d_squared ) {
			min_d_squared = d_squared;
			min_node_index = query_it->upper_vertex();
		}
	}
	return min_node_index;
}

template <class Vertex, class Edge>
void
find_neighbors_naive_surface(
	utility::pointer::shared_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< std::pair< Size, Size > > const & non_surface_ranges,
	utility::vector1< bool > const & is_surface
)
{
	//std::cout<<"finding neighbors...../n";
	// Constants
	core::Size const n_points( point_graph->num_vertices() );
	if ( n_points == 0 ) return;
	core::Real neighbor_cutoff_sq = neighbor_cutoff * neighbor_cutoff;

	// Exclusion checks
	if ( n_points <= 1 ) return; // Nothing to do

	// Naive method: O( R^2 ) for R residues but faster for small, compact conformations
	//std::cout<<"Protein Start Location:"<<conformation.num_jump();
	for ( Size ii = 1; ii <= non_surface_ranges.size(); ++ii ) {
		for ( Size jj = non_surface_ranges[ ii ].first, jjend = non_surface_ranges[ ii ].second; jj <= jjend; ++jj ) {
			PointPosition const & jj_pos( point_graph->get_vertex(jj).data().xyz() );
			for ( Size kk = 1; kk <= n_points; ++kk ) {
				if ( kk <= jj && ! is_surface[ kk ] ) continue;
				Real const d_sq( jj_pos.distance_squared( point_graph->get_vertex( kk ).data().xyz() ) );
				if ( d_sq <= neighbor_cutoff_sq ) {
					Size lower = kk < jj ? kk : jj;
					Size upper = kk < jj ? jj : kk;
					point_graph->add_edge( lower, upper, Edge( d_sq ) );
					//point_graph->add_edge( jj, kk, Edge( d_sq ) );
				}
			}
		}
	}
}


} //conformation
} //core

#endif
