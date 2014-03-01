// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/methods/fold_tree_functions.cc
/// @brief methods for manipulating FoldTrees
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/methods/fold_tree_functions.hh>

// package headers
#include <protocols/forge/methods/util.hh>

// project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/VariantType.hh>

#include <core/graph/DisjointSets.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/exit.hh>

// C++ headers
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <core/pose/util.hh>
#include <utility/vector1.hh>



namespace protocols {
namespace forge {
namespace methods {


static numeric::random::RandomGenerator RG( 2557316 ); // magic number, don't change


// Tracer instance for this file
// Named after the original location of this code
static basic::Tracer TR( "protocols.forge.methods.fold_tree_functions" );


/// @brief enforce Edge has start <= stop (swap if necessary)
/// @return true if start and stop were swapped
bool order( core::kinematics::Edge & e ) {
	if ( e.start() > e.stop() ) {
		int const tmp = e.stop();
		e.stop() = e.start();
		e.start() = tmp;

		if ( !e.start_atom().empty() || !e.stop_atom().empty() ) {
			std::string tmp_s = e.stop_atom();
			e.stop_atom() = e.start_atom();
			e.start_atom() = tmp_s;
		}

		return true;
	}

	return false;
}


/// @brief find the k'th closest larger peptide vertex of given vertex in fold tree
/// @param[in] v the given vertex
/// @param[in] ft fold tree to search
/// @param[in] k find the k'th closest, must be > 0
/// @return 0 if no such vertex
core::Size
closest_larger_peptide_vertex(
	core::Size const v,
	core::kinematics::FoldTree const & ft,
	core::Size const k
)
{
	using core::Size;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	assert( k > 0 );

	std::vector< Size > peptide_vertices;

	for ( FoldTree::const_iterator e = ft.begin(), ee = ft.end(); e != ee; ++e ) {
		if ( e->label() == Edge::PEPTIDE ) {
			peptide_vertices.push_back( e->start() );
			peptide_vertices.push_back( e->stop() );
		}
	}

	// run through in sorted order
	Size count = 0;
	std::sort( peptide_vertices.begin(), peptide_vertices.end() );
	std::unique( peptide_vertices.begin(), peptide_vertices.end() );
	for ( std::vector< Size >::const_iterator i = peptide_vertices.begin(), ie = peptide_vertices.end(); i != ie; ++i ) {
		if ( *i > v ) {
			++count;
			if ( count == k ) {
				return *i;
			}
		}
	}

	return 0;
}


/// @brief find the k'th closest smaller peptide vertex of given vertex in fold tree
/// @param[in] v the given vertex
/// @param[in] ft fold tree to search
/// @param[in] k find the k'th closest, must be > 0
/// @return 0 if no such vertex
core::Size
closest_smaller_peptide_vertex(
	core::Size const v,
	core::kinematics::FoldTree const & ft,
	core::Size const k
)
{
	using core::Size;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	assert( k > 0 );

	std::vector< Size > peptide_vertices;

	for ( FoldTree::const_iterator e = ft.begin(), ee = ft.end(); e != ee; ++e ) {
		if ( e->label() == Edge::PEPTIDE ) {
			peptide_vertices.push_back( e->start() );
			peptide_vertices.push_back( e->stop() );
		}
	}

	// run through in backwards sorted order
	Size count = 0;
	std::sort( peptide_vertices.begin(), peptide_vertices.end() );
	std::unique( peptide_vertices.begin(), peptide_vertices.end() );
	for ( std::vector< Size >::const_reverse_iterator i = peptide_vertices.rbegin(), ie = peptide_vertices.rend(); i != ie; ++i ) {
		if ( *i < v ) {
			++count;
			if ( count == k ) {
				return *i;
			}
		}
	}

	return 0;
}


/// @brief query if vertex already exists in fold tree
bool
vertex_exists(
	core::Size const v,
	core::kinematics::FoldTree const & ft
)
{
	using core::Size;
	using core::kinematics::FoldTree;

	for ( FoldTree::const_iterator e = ft.begin(), ee = ft.end(); e != ee; ++e ) {
		if ( static_cast< Size >( e->start() ) == v || static_cast< Size >( e->stop() ) == v ) {
			return true;
		}
	}

	return false;
}


/// @brief find all jump edges whose start or stop lands on the given position
/// @param[in] pos The position.
/// @param[in] ft The fold tree.
/// @return a list of jump edges
utility::vector1< core::kinematics::Edge > jumps_connected_to_position(
	core::Size const pos,
	core::kinematics::FoldTree const & ft
)
{
	using core::Size;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	typedef utility::vector1< Edge > Edges;

	Edges jump_edges;

	for ( FoldTree::const_iterator e = ft.begin(), ee = ft.end(); e != ee; ++e ) {
		if ( e->label() > 0 && ( static_cast< Size >( e->start() ) == pos || static_cast< Size >( e->stop() ) == pos ) ) {
			jump_edges.push_back( *e );
		}
	}

	return jump_edges;
}


/// @brief find the jump connecting two continuous segments of a fold tree
/// @param[in] u Any vertex on the first segment.
/// @param[in] v Any vertex on the second segment.
/// @param[in] ft The fold tree to query.
/// @return the jump number connecting the two segments, 0 if no such jump
core::Size find_connecting_jump(
	core::Size const u,
	core::Size const v,
	core::kinematics::FoldTree const & ft
)
{
	using core::Size;
	using core::graph::DisjointSets;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	// add all regular edges sans jumps to union-find to construct continuous segments
	DisjointSets uf( ft.nres() );
	for ( FoldTree::const_iterator e = ft.begin(), ee = ft.end(); e != ee; ++e ) {
		Edge edge( *e );
		order( edge );

		if ( e->label() < 0 ) {
			union_interval( edge.start(), edge.start(), edge.stop(), uf );
		}
	}

	// representatives of the segments of u and v
	Size const root_u = uf.ds_find( u );
	Size const root_v = uf.ds_find( v );

	for ( Size i = 1, ie = ft.num_jump(); i <= ie; ++i ) {
		Edge const & e = ft.jump_edge( i );

		Size const root_start = uf.ds_find( e.start() );
		Size const root_stop = uf.ds_find( e.stop() );
		if ( ( root_start == root_u && root_stop == root_v ) || ( root_start == root_v && root_stop == root_u ) ) {
			return i;
		}
	}

	return 0; // couldn't find a jump edge
}


/// @brief remove a cutpoint, merging the two sections that are adjacent to the cut
/// @return true if cutpoint removed, false if not a cutpoint
bool remove_cutpoint(
	core::Size const v,
	core::kinematics::FoldTree & ft
)
{
	using core::Size;
	using core::graph::DisjointSets;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	typedef utility::vector1< Edge > EdgeList;

	if ( !ft.is_cutpoint( v ) ) {
		return false;
	}

	// first:
	// - add all regular edges sans jumps to union-find to construct continuous
	//   segments
	// - collect all jump edges
	DisjointSets uf( ft.nres() );
	EdgeList jump_edges;

	for ( FoldTree::const_iterator e = static_cast< FoldTree const & >( ft ).begin(), ee = static_cast< FoldTree const & >( ft ).end(); e != ee; ++e ) {
		Edge edge( *e );
		order( edge );

		if ( e->label() > 0 ) {
			jump_edges.push_back( edge );
		} else {
			union_interval( edge.start(), edge.start(), edge.stop(), uf );
		}
	}

	// find the representatives of section containing cutpoint and cutpoint+1
	Size const left_root = uf.ds_find( v );
	Size const right_root = uf.ds_find( v + 1 );

	// Run through jump edges to find the jump connecting the left & right sets.
	// Jump edges were ordered when collecting them, so the start() must correspond
	// to the left set and the stop() must correspond to the right set.
	EdgeList::const_iterator j = jump_edges.begin(), je = jump_edges.end();
	while(
		j != je &&
		uf.ds_find( j->start() ) != left_root &&
		uf.ds_find( j->stop() ) != right_root
	)
	{
		++j;
	}

	// remove the jump and subsequently the cutpoint; this call should
	// automatically reorder the tree
	ft.delete_jump_and_intervening_cutpoint( j->start(), j->stop() );

	return true;
}


/// @brief seal a fold tree by removing all specified cutpoints
/// @param[in] cutpoints Cutpoints to remove.
/// @param[in,out] ft The input tree.
/// @return A new tree with cutpoints removed and all other topology kept
///  constant.
void remove_cutpoints(
	utility::vector1< core::Size > const & cutpoints,
	core::kinematics::FoldTree & ft
)
{
	using core::Size;
	using core::graph::DisjointSets;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	typedef utility::vector1< Edge > EdgeList;

	// first:
	// - add all regular edges sans jumps to union-find to construct continuous
	//   segments
	// - collect all jump edges
	DisjointSets uf( ft.nres() );
	EdgeList jump_edges;

	for ( FoldTree::const_iterator e = static_cast< FoldTree const & >( ft ).begin(), ee = static_cast< FoldTree const & >( ft ).end(); e != ee; ++e ) {
		Edge edge( *e );
		order( edge );

		if ( e->label() > 0 ) {
			jump_edges.push_back( edge );
		} else {
			union_interval( edge.start(), edge.start(), edge.stop(), uf );
		}
	}

	// run through all cutpoints
	for ( utility::vector1< Size >::const_iterator i = cutpoints.begin(), ie = cutpoints.end(); i != ie; ++i ) {
		// find the representatives of section containing cutpoint and cutpoint+1
		Size const v = *i;
		Size const left_root = uf.ds_find( v );
		Size const right_root = uf.ds_find( v + 1 );

		// Run through jump edges to find the jump connecting the left & right sets.
		// Jump edges were ordered when collecting them, so the start() must correspond
		// to the left set and the stop() must correspond to the right set.
		EdgeList::const_iterator j = jump_edges.begin(), je = jump_edges.end();
		while(
			j != je &&
			uf.ds_find( j->start() ) != left_root &&
			uf.ds_find( j->stop() ) != right_root
		)
		{
			++j;
		}

		// remove the jump and subsequently the cutpoint; this call should
		// automatically reorder the tree
		ft.delete_jump_and_intervening_cutpoint( j->start(), j->stop() );

		// track the merge between the two segments
		uf.ds_union( left_root, right_root );
	}
}


/// @brief attempt to shift jumps in a fold tree based on fixed positions in a MoveMap
/// @param[in] ft FoldTree to alter
/// @param[in] movemap MoveMap to try and honor
/// @return FoldTree with jumps shifted, if possible
/// @remarks Procedure will shift a jump position on a continuous segment to a
///  randomly selected fixed position on that segment specified by the MoveMap.
///  No changes will be made to a jump point if it is contained within a segment
///  lacking valid fixed positions.  This procedure has the side effect of
///  collapsing all jump points to a single point within segments that have
///  fixed positions.
core::kinematics::FoldTree
shift_jumps(
	core::kinematics::FoldTree const & ft,
	core::kinematics::MoveMap const & mm
)
{
	using core::Size;
	using core::graph::DisjointSets;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	typedef utility::vector1< Edge > EdgeList;
	typedef utility::vector1< Size > NodeList;
	typedef std::map< Size, Size > Root2Jump;
	typedef std::map< Size, NodeList > Root2Nodes;

	FoldTree aft; // altered FoldTree, return this
	DisjointSets uf( ft.nres() );
	EdgeList regular_edges; // non-jump edges
	EdgeList jump_edges;

	// separate edges into categories; construct continuous segments using union-find
	for ( FoldTree::const_iterator e = ft.begin(), ee = ft.end(); e != ee; ++e ) {
		if ( e->label() > 0 ) { // jump edge

			jump_edges.push_back( *e );
			order( jump_edges.back() );

		} else { // regular edge

			regular_edges.push_back( *e );
			order( regular_edges.back() );
			union_interval( regular_edges.back().start(), regular_edges.back().start(), regular_edges.back().stop(), uf );

		}
	}

	// pick a single jump point per segment
	Root2Jump r2j;
	Root2Nodes r2n = uf.sets();
	for ( Root2Nodes::iterator i = r2n.begin(), ie = r2n.end(); i != ie; ++i ) {
		NodeList & nodes = i->second;
		NodeList fixed;

		// find all useable fixed positions for jump
		for ( NodeList::const_iterator j = nodes.begin(), je = nodes.end(); j != je; ++j ) {
			if ( !mm.get_bb( *j ) ) {
				fixed.push_back( *j );
			}
		}

		// pick fixed position for jump
		if ( !fixed.empty() ) {
			r2j[ i->first ] = fixed[ RG.random_range( 1, fixed.size() ) ];
		}
	}

	// swap all jump points, split any necessary regular edges
	for ( EdgeList::iterator e = jump_edges.begin(), ee = jump_edges.end(); e != ee; ++e ) {
		Root2Jump::const_iterator new_start = r2j.find( uf.ds_find( e->start() ) );
		Root2Jump::const_iterator new_stop = r2j.find( uf.ds_find( e->stop() ) );

		// swap start
		if ( new_start != r2j.end() ) {
			e->start() = new_start->second;
		}

		// swap stop
		if ( new_stop != r2j.end() ) {
			e->stop() = new_stop->second;
		}

		// split facets if necessary
		add_vertex( e->start(), regular_edges );
		add_vertex( e->stop(), regular_edges );
	}

	// add all regular edges
	for ( EdgeList::const_iterator e = regular_edges.begin(), ee = regular_edges.end(); e != ee; ++e ) {
		aft.add_edge( *e );
	}

	// add all jump edges
	for ( EdgeList::const_iterator e = jump_edges.begin(), ee = jump_edges.end(); e != ee; ++e ) {
		aft.add_edge( *e );
	}

	// finalize
	aft.reorder( ft.root() );

	return aft;
}


/// @brief construct a fold tree from Pose wrt chain endings and residue types
/// @remarks Determines edge types from residues (polymer vs non-polymer).
///  Each chain will end up as one edge in the fold tree.  Jumps will be made
///  from the new root to a random fixed backbone residue as specified by the
///  movemap.  If all residues in a chain are moveable, will choose any random
///  residue.
/// @param[in] pose Pose.
/// @param[in] ft_root Root of the new fold tree.
/// @param[in] mm MoveMap used to select jump positions.
/// @return fold tree wrt chain endings and residue types
core::kinematics::FoldTree
fold_tree_from_pose(
	core::pose::Pose const & pose,
	core::Size const ft_root,
	core::kinematics::MoveMap const & mm
)
{
	using core::Size;
	using core::conformation::Residue;
	using core::graph::DisjointSets;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	//typedef utility::vector1< Edge > Edges;
	typedef utility::vector1< Size > Nodes;
	typedef std::map< Size, Nodes > Root2Nodes;

	FoldTree ft;

	// NOTE: This procedure could be significantly simpler but the current fold
	// tree implementation does not allow it; hence what is essentially a series
	// of workarounds below.

	// use union-find to track connected components
	DisjointSets uf( pose.n_residue() );

	// first connect sequential polymer vertices and do connections
	// between chemical vertices
	for ( Size i = 1, ie = pose.conformation().num_chains(); i <= ie; ++i ) {
		Size const chain_begin = pose.conformation().chain_begin( i );
		Size const chain_end = pose.conformation().chain_end( i );

		Size prior_polymer_r = pose.residue( chain_begin ).is_polymer() ? chain_begin : 0;
		for ( Size r = chain_begin; r <= chain_end; ++r ) {
			Residue const & res = pose.residue( r );

			if ( res.is_polymer() ) {

				// if polymer and sequential, union with prior residue;
				// cannot add polymer edge yet (also, some FoldTree internal
				// topology checks barf when too many [i,i+1] vertices)
				if ( prior_polymer_r > 0 && r == prior_polymer_r + 1 ) { // safe for r == 1
					uf.ds_union( r, prior_polymer_r );
				}

				prior_polymer_r = r;

			} else { // chemical connection

				for ( Size c = 1, ce = res.n_residue_connections(); c <= ce; ++c ) {
					Size const cr = res.connected_residue_at_resconn( c );
					Residue const & cres = pose.residue( cr );

					assert( static_cast< Size >( res.chain() ) == i );

					if ( r < cr && res.chain() == cres.chain() ) {

						// go ahead and add the edge to fold tree;
						// uses implicit chemical edge constructor
						ft.add_edge(
							Edge(
								r, cr,
								res.atom_name( res.connect_atom( cres ) ),
								cres.atom_name( cres.connect_atom( res ) ) )
						);

						uf.ds_union( cr, r );
					}
				} // foreach residue connection

			}
		} // foreach chain
	}

	// grab components up to this point
	Root2Nodes r2n = uf.sets();

	// run through each component, adding polymer edges and assigning jumps
	Size jump_count = 0;
	for ( Root2Nodes::iterator i = r2n.begin(), ie = r2n.end(); i != ie; ++i ) {
		Nodes & nodes = i->second;
		std::sort( nodes.begin(), nodes.end() );
		bool const component_contains_ft_root = std::find( nodes.begin(), nodes.end(), ft_root ) != nodes.end();

		if ( nodes.size() > 1 ) {

			// Grab only polymer vertices.  Need to do this to separate any
			// mixed vertex scenarios.
			Nodes polymer_vertices;
			for ( Nodes::const_iterator v = nodes.begin(), ve = nodes.end(); v != ve; ++v ) {
				if ( pose.residue( *v ).is_polymer() ) {
					polymer_vertices.push_back( *v );
				}
			}

			if ( polymer_vertices.empty() ) {

				// everything is chemical, pick a jump point and be done
				if ( !component_contains_ft_root ) {
					Size const jump_point = RG.random_range( *nodes.begin(), *nodes.rbegin() );
					ft.add_edge( Edge( ft_root, jump_point, ++jump_count ) );
					uf.ds_union( jump_point, ft_root );
				}

			} else if ( polymer_vertices.size() == 1 ) {

				// single vertex, jump point is the single polymer vertex
				if ( !component_contains_ft_root ) {
					ft.add_edge( Edge( ft_root, *polymer_vertices.begin(), ++jump_count ) );
					uf.ds_union( *polymer_vertices.begin(), ft_root );
				}

			} else {

				// construct re-indexed uf to track connected polymer subsections
				DisjointSets indexed_uf( polymer_vertices.size() );
				for ( Size j = 1, je = polymer_vertices.size(); j < je; ++j ) {
					if ( polymer_vertices[ j ] + 1 == polymer_vertices[ j+1 ] ) {
						// for simplicity do the union here instead of when the
						// polymer edge is added
						indexed_uf.ds_union( j+1, j );
					}
				}

				Root2Nodes indexed_r2n = indexed_uf.sets();
				for ( Root2Nodes::iterator j = indexed_r2n.begin(), je = indexed_r2n.end(); j != je; ++j ) {
					Nodes subsection = j->second;
					std::sort( subsection.begin(), subsection.end() );
					Size const ss_begin = polymer_vertices[ *subsection.begin() ];
					Size const ss_end = polymer_vertices[ *subsection.rbegin() ];
					bool const subsection_contains_ft_root = ss_begin <= ft_root && ft_root <= ss_end;

					// find appropriate fixed bb positions for jump
					utility::vector1< Size > fixed;
					for ( Size k = ss_begin; k <= ss_end; ++k ) {
						if ( !mm.get_bb( k ) ) {
							fixed.push_back( k );
						}
					}

					// pick a jump point
					Size jump_point = ft_root;
					if ( !subsection_contains_ft_root ) {
						if ( fixed.empty() ) {
							jump_point = RG.random_range( ss_begin, ss_end );
						} else {
							jump_point = fixed[ RG.random_range( 1, fixed.size() ) ];
						}
					}

					// order necessary vertices, remove duplicates
					std::set< Size > vertices;
					vertices.insert( ss_begin );
					vertices.insert( ss_end );
					vertices.insert( jump_point );

					// add polymer edges, there are only three cases
					switch ( vertices.size() ) {
						case 3:
							ft.add_edge( Edge( *vertices.begin(), *( ++vertices.begin() ), Edge::PEPTIDE ) );
							ft.add_edge( Edge( *( ++vertices.begin() ), *( ++( ++vertices.begin() ) ), Edge::PEPTIDE ) );
							break;
						case 2:
							ft.add_edge( Edge( *vertices.begin(), *( ++vertices.begin() ), Edge::PEPTIDE ) );
							break;
						case 1:
							// do nothing
							break;
						default:
							TR.Fatal << "FATAL: fold_tree_from_pose() : vertices.size() not in [1, 3]" << std::endl;
							utility_exit_with_message( "should not be here" );
							break;
					}

					// add jump edge
					if ( !subsection_contains_ft_root ) {
						ft.add_edge( Edge( ft_root, jump_point, ++jump_count ) );
						uf.ds_union( jump_point, ft_root );
					}

				}

			}

		} else { // nodes.size() == 1
			// with one vertex we only need to connect a jump
			if ( !component_contains_ft_root ) {
				ft.add_edge( Edge( ft_root, *nodes.begin(), ++jump_count ) );
				uf.ds_union( *nodes.begin(), ft_root );
			}
		}
	}

	assert( uf.n_disjoint_sets() == 1 ); // everything should now be connected

	// finalize
	bool const success = ft.reorder( ft_root );
	if ( !success ) {
		TR.Fatal << "FATAL: fold_tree_from_pose(): " << ft << std::endl;
		utility_exit_with_message( "bad fold tree topology" );
	}

	return ft;
}


/// @brief merge two fold trees by jump between their roots
/// @param[in] left_tree
/// @param[in] right_tree
/// @return Merged FoldTree with all vertices of right_tree placed to the
///  right of vertices in left_tree.  Jump labels in right_tree will be
///  renumbered += left_tree.num_jump().  New tree is rooted in the same
///  place as left_tree.
core::kinematics::FoldTree
merge(
	core::kinematics::FoldTree const & left_tree,
	core::kinematics::FoldTree const & right_tree
)
{
	return merge(
		left_tree, left_tree.root(), std::string(),
		right_tree, right_tree.root(), std::string(),
		false
	);
}


/// @brief merge two fold trees connecting by jump via specified positions
/// @param[in] left_tree
/// @param[in] left_position position on left_tree to connect
/// @param[in] right_tree
/// @param[in] right_position position on right_tree to connect
/// @return Merged FoldTree with all vertices of right_tree placed to the
///  right of vertices in left_tree.  Jump labels in right_tree will be
///  renumbered += left_tree.num_jump().  New tree is rooted in the same
///  place as left_tree.
core::kinematics::FoldTree
merge(
	core::kinematics::FoldTree const & left_tree,
	core::Size const left_position,
	core::kinematics::FoldTree const & right_tree,
	core::Size const right_position
)
{
	return merge(
		left_tree, left_position, std::string(),
		right_tree, right_position, std::string(),
		false
	);
}


/// @brief merge two fold trees connecting by jump via specified positions
/// @param[in] left_tree
/// @param[in] left_position position on left_tree to connect
/// @param[in] left_jump_atom Use this atom for the left side of the jump.
///  Use empty string to indicate default setting.
/// @param[in] right_tree
/// @param[in] right_position position on right_tree to connect
/// @param[in] right_jump_atom Use this atom for the right set of the jump.
///  Use empty string to indicate default setting.
/// @param[in] keep_stub_in_residue Attempt to keep generated stubs of the
///  jump within their respective residues.  default False
/// @return Merged FoldTree with all vertices of right_tree placed to the
///  right of vertices in left_tree.  Jump labels in right_tree will be
///  renumbered += left_tree.num_jump().  New tree is rooted in the same
///  place as left_tree.
core::kinematics::FoldTree
merge(
	core::kinematics::FoldTree const & left_tree,
	core::Size const left_position,
	std::string const & left_jump_atom,
	core::kinematics::FoldTree const & right_tree,
	core::Size const right_position,
	std::string const & right_jump_atom,
	bool const keep_stub_in_residue
)
{
	using core::Size;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;

	assert( left_position <= left_tree.nres() );
	assert( right_position <= right_tree.nres() );

	FoldTree new_ft = left_tree;

	// add all shifted edges from right_tree
	for ( FoldTree::const_iterator r = right_tree.begin(), re = right_tree.end(); r != re; ++r ) {
		Edge new_edge = *r;
		new_edge.start() += left_tree.nres();
		new_edge.stop() += left_tree.nres();
		if ( new_edge.label() > 0 ) {
			new_edge.label() += left_tree.num_jump();
		}
		new_ft.add_edge( new_edge );
	}

	// compute jump positions
	Size const jleft = left_position;
	Size const jright = right_position + left_tree.nres();

	// split edge originating from left tree
	if ( !vertex_exists( jleft, new_ft ) ) {
		FoldTree::const_iterator f = facet_containing_position(
			jleft,
			static_cast< FoldTree const & >( new_ft ).begin(),
			static_cast< FoldTree const & >( new_ft ).end(),
			Edge::PEPTIDE
		);
		assert( f != static_cast< FoldTree const & >( new_ft ).end() ); // paranoia

		Edge const base = *f;
		new_ft.delete_edge( base );
		new_ft.add_edge( base.start(), jleft, base.label() );
		new_ft.add_edge( jleft, base.stop(), base.label() );
	}

	// split edge originating from right tree
	if ( !vertex_exists( jright, new_ft ) ) {
		FoldTree::const_iterator f = facet_containing_position(
			jright,
			static_cast< FoldTree const & >( new_ft ).begin(),
			static_cast< FoldTree const & >( new_ft ).end(),
			Edge::PEPTIDE
		);
		assert( f != static_cast< FoldTree const & >( new_ft ).end() ); // paranoia

		Edge const base = *f;
		new_ft.delete_edge( base );
		new_ft.add_edge( base.start(), jright, base.label() );
		new_ft.add_edge( jright, base.stop(), base.label() );
	}

	// add new jump
	new_ft.add_edge(
		Edge(
			jleft, jright,
			left_tree.num_jump() + right_tree.num_jump() + 1,
			left_jump_atom, right_jump_atom,
			keep_stub_in_residue
		)
	);

	// re-root
	new_ft.reorder( left_tree.root() );

	return new_ft;
}


/// @brief replace a section of one fold tree with another fold tree, connecting by
///  jump between their roots
/// @param[in] original_tree
/// @param[in] replace_begin residue starting the section of original_tree to replace
/// @param[in] replace_end residue ending the section of original_tree to replace
/// @param[in] movemap MoveMap whose fixed backbone positions dictates where new jumps
///  may be placed.
/// @param[in] replacement_tree
/// @return FoldTree with section replaced.
/// @remarks The procedure will attempt to honor the MoveMap as much as it can.  The caveat
///  is that sequences of calls to some FoldTree routines may shift the jumps internally in
///  a way that is not easily predictable.  If the procedure cannot find an allowed
///  residue for a jump, it will make a jump to the median residue in the disconnected
///  fold tree interval.
core::kinematics::FoldTree
replace(
	core::kinematics::FoldTree const & original_tree,
	int const replace_begin,
	int const replace_end,
	core::kinematics::MoveMap const & movemap,
	core::kinematics::FoldTree const & replacement_tree
)
{
	using core::Size;
	using core::graph::DisjointSets;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;

	//typedef utility::vector1< Edge > EdgeList;
	typedef utility::vector1< Size > NodeList;
	typedef std::map< Size, NodeList > Root2Nodes;

	assert( replace_begin <= replace_end );
	assert( original_tree.root() < replace_begin || original_tree.root() > replace_end );

	Size const replace_length = replace_end - replace_begin + 1;
	Size const final_nres = original_tree.nres() - replace_length + replacement_tree.nres();
	Size const final_ft_root = original_tree.root() < replace_begin ?
	                           original_tree.root() :
	                           original_tree.root() - replace_length + replacement_tree.nres();

	// create altered fold tree and delete the section to be replaced
	FoldTree aft = original_tree;
	for ( int r = replace_begin; r <= replace_end; ++r ) {
		aft.delete_seqpos( replace_begin );
	}

	TR.Debug << "original_tree: " << original_tree << std::endl;
	TR.Debug << "aft: " << aft << std::endl;

	DisjointSets uf( final_nres ); // union-find
	FoldTree tracking_ft = aft;
	FoldTree new_ft;

	// left section
	for ( FoldTree::const_iterator e = static_cast< FoldTree const & >( aft ).begin(),
	      ee = static_cast< FoldTree const & >( aft ).end(); e != ee; ++e )
	{
		if ( e->label() < 0 && e->start() < replace_begin && e->stop() < replace_begin ) {
			Edge new_edge = *e;
			order( new_edge ); // for union_interval

			new_ft.add_edge( new_edge );
			tracking_ft.delete_edge( *e );

			union_interval( new_edge.start(), new_edge.start(), new_edge.stop(), uf );
		}
	}

	TR.Debug << "new_ft after left: " << new_ft << std::endl;

	// right section, shifted to take into account replacement tree
	for ( FoldTree::const_iterator e = static_cast< FoldTree const & >( aft ).begin(),
	      ee = static_cast< FoldTree const & >( aft ).end(); e != ee; ++e )
	{
		if ( e->label() < 0 && e->start() >= replace_begin && e->stop() >= replace_begin ) {
			Edge new_edge = *e;
			order( new_edge ); // for union_interval
			new_edge.start() += replacement_tree.nres();
			new_edge.stop() += replacement_tree.nres();

			new_ft.add_edge( new_edge );
			tracking_ft.delete_edge( *e );

			union_interval( new_edge.start(), new_edge.start(), new_edge.stop(), uf );
		}
	}

	TR.Debug << "new_ft after right: " << new_ft << std::endl;

	// jumps
	for ( FoldTree::const_iterator e = static_cast< FoldTree const & >( aft ).begin(),
	      ee = static_cast< FoldTree const & >( aft ).end(); e != ee; ++e )
	{
		if ( e->label() > 0 ) {
			Edge new_edge = *e;

			if ( new_edge.start() >= replace_begin ) {
				new_edge.start() += replacement_tree.nres();
			}

			if ( new_edge.stop() >= replace_begin ) {
				new_edge.stop() += replacement_tree.nres();
			}

			new_ft.add_edge( new_edge );
			tracking_ft.delete_edge( *e );

			uf.ds_union( new_edge.start(), new_edge.stop() );
		}
	}

	assert( tracking_ft.size() <= 2 );

	TR.Debug << "new_ft after jump: " << new_ft << std::endl;

	// run through remaining edges, splitting if necessary
	for ( FoldTree::const_iterator e = static_cast< FoldTree const & >( tracking_ft ).begin(),
	      ee = static_cast< FoldTree const & >( tracking_ft ).end(); e != ee; ++e )
	{
		Edge edge = *e;
		order( edge ); // for union_interval

		if ( edge.stop() == replace_begin ) { // on right vertex
			Edge new_edge = edge;
			--new_edge.stop();

			new_ft.add_edge( new_edge );

			union_interval( new_edge.start(), new_edge.start(), new_edge.stop(), uf );

		} else if ( edge.start() == replace_begin ) { // on left vertex
			Edge new_edge = edge;
			new_edge.start() += replacement_tree.nres();
			new_edge.stop() += replacement_tree.nres();

			new_ft.add_edge( new_edge );

			union_interval( new_edge.start(), new_edge.start(), new_edge.stop(), uf );

		} else if ( edge.start() < replace_begin ) { // split
			Edge left_edge = edge;
			left_edge.stop() = replace_begin - 1;

			Edge right_edge = edge;
			right_edge.start() = replace_begin + replacement_tree.nres();
			right_edge.stop() += replacement_tree.nres();

			new_ft.add_edge( left_edge );
			new_ft.add_edge( right_edge );

			union_interval( left_edge.start(), left_edge.start(), left_edge.stop(), uf );
			union_interval( right_edge.start(), right_edge.start(), right_edge.stop(), uf );
		}
	}

	// remove any single residue edges that might exist after splitting
	new_ft.delete_self_edges();

	TR.Debug << "new_ft after remaining: " << new_ft << std::endl;

	// replacement tree indexing shift
	Size const rshift = replace_begin - 1;
	Size const jshift = new_ft.num_jump();

	// add edges from replacement tree
	for ( FoldTree::const_iterator e = replacement_tree.begin(), ee = replacement_tree.end(); e != ee; ++e ) {
		Edge new_edge = *e;
		new_edge.start() +=  rshift;
		new_edge.stop() += rshift;
		if ( new_edge.label() > 0 ) {
			new_edge.label() += jshift;
		}

		new_ft.add_edge( new_edge );

		if ( new_edge.label() > 0 ) {
			uf.ds_union( new_edge.start(), new_edge.stop() );
		} else {
			order( new_edge ); // for union_interval
			union_interval( new_edge.start(), new_edge.start(), new_edge.stop(), uf );
		}
	}

	// make new jump from final_ft_root to root of replacement tree
	new_ft.add_edge( final_ft_root, replacement_tree.root() + rshift, new_ft.num_jump() + 1 );
	uf.ds_union( final_ft_root, replacement_tree.root() + rshift );

	TR.Debug << "new_ft after replacement: " << new_ft << std::endl;

	// construct re-indexed movemap
	MoveMap mm;
	mm.set_bb_true_range( 1, final_nres );
	for ( Size i = 1, ie = original_tree.nres(); i <= ie; ++i ) {
		if ( i < static_cast< Size >( replace_begin ) ) {
			mm.set_bb( i, movemap.get_bb( i ) );
		} else if ( static_cast< Size >( replace_end ) < i ) {
			mm.set_bb( i - replace_length + replacement_tree.nres(), movemap.get_bb( i ) );
		}
	}

	// find disconnected components and make jumps from final_root to an
	// allowed residue governed by the MoveMap
	Root2Nodes r2n = uf.sets();
	for ( Root2Nodes::iterator i = r2n.begin(), ie = r2n.end(); i != ie; ++i ) {
		if ( i->first != uf.ds_find( final_ft_root ) ) {
			NodeList & nodes = i->second;
			NodeList fixed;

			// find all useable fixed positions for jump
			for ( NodeList::const_iterator j = nodes.begin(), je = nodes.end(); j != je; ++j ) {
				if ( !mm.get_bb( *j ) ) {
					fixed.push_back( *j );
				}
			}

			// new edge
			Edge new_edge;
			new_edge.start() = final_ft_root;
			new_edge.label() = new_ft.num_jump() + 1;

			// add jump
			if ( fixed.empty() ) { // problematic case

				// we have no fixed positions, so just add the jump to the
				// (lower) median node
				std::sort( nodes.begin(), nodes.end() );
				if ( nodes.size() % 2 == 1 ) { // median
					new_edge.stop() = nodes[ ( nodes.size() + 1 ) / 2 ];
				} else { // lower median
					new_edge.stop() = nodes[ nodes.size() / 2 ];
				}

			} else { // normal case
				new_edge.stop() = fixed[ RG.random_range( 1, fixed.size() ) ];
			}

			// split existing edge if necessary
			FoldTree::const_iterator f_connecting = facet_containing_position(
				new_edge.stop(),
				static_cast< FoldTree const & >( new_ft ).begin(),
				static_cast< FoldTree const & >( new_ft ).end(),
				Edge::PEPTIDE
			);
			if ( f_connecting != static_cast< FoldTree const & >( new_ft ).end() ) {
				Edge const connecting_edge = *f_connecting;
				new_ft.delete_edge( connecting_edge );
				new_ft.add_edge( connecting_edge.start(), new_edge.stop(), connecting_edge.label() );
				new_ft.add_edge( connecting_edge.stop(), new_edge.stop(), connecting_edge.label() );
			}

			// add to fold tree
			new_ft.add_edge( new_edge );

			// add to uf
			uf.ds_union( new_edge.start(), new_edge.stop() );

		}
	}

	assert( uf.n_disjoint_sets() == 1 );

	TR.Debug << "new_ft after disconnected: " << new_ft << std::endl;

	// finalize
	bool const success = new_ft.reorder( final_ft_root );
	runtime_assert( success );

	return new_ft;
}

// set up star foldtree for special manipulation
void make_star_foldtree(
		core::pose::Pose & pose,
		//core::kinematics::MoveMap & mm,
		protocols::loops::Loops loops ) {
	using namespace core::chemical;
	using namespace core::kinematics;

	core::Size nres = pose.total_residue()-1;
	core::kinematics::FoldTree newF;

	core::Size prev_cut = 0, this_cut, out_midpt;
	int njump =1 ;
	for ( core::Size i=1; i <= loops.size(); ++i ) {
		core::Size loop_start = loops[i].start() ;
		core::Size loop_end = loops[i].stop() -1 ; // special case to treat loop objects as definition of a cut and not a real range. only use this for two chain tree!

		bool start_is_cut = pose.fold_tree().is_cutpoint( loop_start-1 );
		bool end_is_cut = pose.fold_tree().is_cutpoint( loop_end );

		if ( loop_start == 1) continue;
		if ( loop_end == nres) continue;

		/// some really weird cases
		if ( start_is_cut && !end_is_cut && prev_cut == loop_start-1 )
			continue;
		if ( start_is_cut && end_is_cut && prev_cut != loop_start-1 ) {
			// need to add two cuts for this loop
			this_cut   = loop_start-1;
			//out_midpt  = (prev_cut + this_cut+1)/4;
			out_midpt  = prev_cut +1; // try to root the jump at beginning
			newF.add_edge( nres+1, out_midpt, njump );
			if (out_midpt != prev_cut+1)
				newF.add_edge( out_midpt, prev_cut+1, Edge::PEPTIDE );
			if (out_midpt != this_cut)
				newF.add_edge( out_midpt, this_cut  , Edge::PEPTIDE );
			prev_cut = this_cut;

			njump++;
		}

		if ( start_is_cut && !end_is_cut && prev_cut != loop_start-1 )
			this_cut = loop_start-1;
		else if ( end_is_cut )
			this_cut = loop_end;
		else
			this_cut = loops[i].cut();

		//out_midpt  = (prev_cut + this_cut+1)/2;
		out_midpt = prev_cut + 1; //try to root the jump at beginning.

		if ( !start_is_cut && !end_is_cut ) {
			core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, this_cut   );
			core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, this_cut+1 );
		}
		newF.add_edge( nres+1, out_midpt, njump );
		if (out_midpt != prev_cut+1)
			newF.add_edge( out_midpt, prev_cut+1, Edge::PEPTIDE );
		if (out_midpt != this_cut)
			newF.add_edge( out_midpt, this_cut  , Edge::PEPTIDE );
		TR << "add edge " << njump << " : " << prev_cut+1 << "..." << out_midpt << "..." << this_cut << std::endl;

		njump++;
		prev_cut = this_cut;
	}

	// c term
	if (prev_cut != nres) {
		out_midpt  = (prev_cut + nres+1)/2;
		newF.add_edge( prev_cut+1, nres+1, njump );
		newF.add_edge( prev_cut+1, nres  , Edge::PEPTIDE );
		TR << "add edge " << njump << " : " << prev_cut+1 << "..." << out_midpt << "..." << nres << std::endl;
	}

	newF.reorder( nres+1 );  // root the tree on the VRT res
	TR << newF << std::endl;
	pose.fold_tree( newF );
}

void jumps_and_cuts_from_pose( core::pose::Pose & pose, utility::vector1< std::pair<core::Size, core::Size > > & jumps, utility::vector1< core::Size > & cuts){

	core::kinematics::FoldTree f_orig = pose.fold_tree();

	for ( core::Size i = 1; i<= f_orig.num_jump(); ++i ) {
    core::Size down ( f_orig.downstream_jump_residue(i) );
    core::Size up ( f_orig.upstream_jump_residue(i) );
        jumps.push_back( std::pair<int,int>( down, up ) );
	}
 cuts =  f_orig.cutpoints();
}


} // methods
} // forge
} // protocols
