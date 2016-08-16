// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/methods/fold_tree_functions.hh
/// @brief methods for manipulating FoldTrees
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_methods_fold_tree_functions_hh
#define INCLUDED_protocols_forge_methods_fold_tree_functions_hh

// type headers
#include <core/types.hh>

// project headers
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// utility headers

// C++ headers
#include <string>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace methods {


/// @brief enforce Edge has start <= stop (swap if necessary)
/// @return true if start and stop were swapped
bool order( core::kinematics::Edge & e );


/// @brief add a set of Edges to a FoldTree
/// @param[in] begin iterator pointing to the first edge to add
/// @param[in] end iterator pointing just past the last edge to add
/// @param[out ft the FoldTree to modify
template< typename EdgeIterator >
void
add_edges(
	EdgeIterator begin,
	EdgeIterator end,
	core::kinematics::FoldTree & ft
)
{
	for ( EdgeIterator e = begin; e != end; ++e ) {
		ft.add_edge( *e );
	}
}


/// @brief find regular facet (non-jump edge sans face vertices) that contains given position
/// @param[in] pos position
/// @param[in] begin iterator pointing to the first edge to search
/// @param[in] end iterator pointing just past the last edge to search
/// @return iterator pointing to the Edge defining the facet that contains the
///  position; if vertex lies on the face of an edge or vertex not found, returns
///  'end'
template< typename EdgeIterator >
EdgeIterator
regular_facet_containing_position(
	core::Size const pos,
	EdgeIterator const begin,
	EdgeIterator const end
)
{
	using core::Size;
	using core::kinematics::Edge;

	for ( EdgeIterator e = begin; e != end; ++e ) {
		if ( e->label() < 0 ) {

			if ( ( static_cast< Size >( e->start() ) < pos && pos < static_cast< Size >( e->stop() ) ) ||
					( static_cast< Size >( e->start() ) > pos && pos > static_cast< Size >( e->stop() ) ) ) {
				return e;
			}

		}
	}

	return end;
}


/// @brief add a vertex, splitting edges if necessary
/// @param[in] v the vertex to add
/// @param[in,out] edges the list of edges to modify
/// @return true if vertex was added and an edge was split, false otherwise
template< typename EdgeList >
bool
add_vertex(
	core::Size const v,
	EdgeList & edges
)
{
	using core::kinematics::Edge;

	// attempt to find the facet
	typename EdgeList::iterator f = regular_facet_containing_position(
		v,
		edges.begin(),
		edges.end()
	);

	// if the facet exists, it needs to be split
	if ( f != edges.end() ) {
		Edge left = *f;
		Edge right = *f;

		order( left );
		order( right );

		left.stop() = v;
		right.start() = v;

		edges.erase( f );

		edges.push_back( left );
		edges.push_back( right );

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
	core::Size const k = 1
);


/// @brief find the k'th closest smaller peptide vertex of given vertex in fold tree
/// @param[in] v the given vertex
/// @param[in] ft fold tree to search
/// @param[in] k find the k'th closest, must be > 0
/// @return 0 if no such vertex
core::Size
closest_smaller_peptide_vertex(
	core::Size const v,
	core::kinematics::FoldTree const & ft,
	core::Size const k = 1
);


/// @brief query if vertex already exists in fold tree
bool
vertex_exists(
	core::Size const v,
	core::kinematics::FoldTree const & ft
);


/// @brief find facet (edge sans face vertices) that contains given position
/// @param[in] pos position
/// @param[in] begin iterator pointing to the first edge to search
/// @param[in] end iterator pointing just past the last edge to search
/// @param[in] label the types of Edge to search, default PEPTIDE
/// @return iterator pointing to the Edge defining the facet that contains the
///  position; if vertex lies on the face of an edge or vertex not found, returns
///  'end'
template< typename EdgeIterator >
EdgeIterator
facet_containing_position(
	core::Size const pos,
	EdgeIterator const begin,
	EdgeIterator const end,
	int const label = core::kinematics::Edge::PEPTIDE
)
{
	using core::Size;
	using core::kinematics::Edge;

	for ( EdgeIterator e = begin; e != end; ++e ) {
		if ( e->label() == label ) {

			if ( ( static_cast< Size >( e->start() ) < pos && pos < static_cast< Size >( e->stop() ) ) ||
					( static_cast< Size >( e->start() ) > pos && pos > static_cast< Size >( e->stop() ) ) ) {
				return e;
			}

		}
	}

	return end;
}


/// @brief find regular (non-jump) edges contained within the interval [left, right]
template< typename EdgeIterator >
utility::vector1< core::kinematics::Edge >
regular_edges_within_interval(
	core::Size const left,
	core::Size const right,
	EdgeIterator const begin,
	EdgeIterator const end
)
{
	using core::Size;
	using core::kinematics::Edge;

	utility::vector1< Edge > regular_edges;
	for ( EdgeIterator e = begin; e != end; ++e ) {
		if ( e->label() < 0 ) {

			Edge tmp( *e );
			order( tmp );

			if ( !( right < static_cast< Size >( tmp.start() ) || static_cast< Size >( tmp.stop() ) < left ) ) {
				regular_edges.push_back( *e );
			}

		}
	}

	return regular_edges;
}


/// @brief find specific type of edges contained within the interval [left, right]
template< typename EdgeIterator >
utility::vector1< core::kinematics::Edge >
edges_within_interval(
	core::Size const left,
	core::Size const right,
	int const edge_type,
	EdgeIterator const begin,
	EdgeIterator const end
)
{
	using core::Size;
	using core::kinematics::Edge;

	utility::vector1< Edge > collected_edges;
	for ( EdgeIterator e = begin; e != end; ++e ) {
		if ( e->label() == edge_type ) {

			Edge tmp( *e );
			order( tmp );

			if ( !( right < static_cast< Size >( tmp.start() ) || static_cast< Size >( tmp.stop() ) < left ) ) {
				collected_edges.push_back( *e );
			}

		}
	}

	return collected_edges;
}


/// @brief find all jump edges whose start or stop lands on the given position
/// @param[in] pos The position.
/// @param[in] ft The fold tree.
/// @return a list of jump edges
/// @remarks jump edges start/stop will be ordered such that 'pos' will always
///  be at 'start'.
utility::vector1< core::kinematics::Edge > jumps_connected_to_position(
	core::Size const pos,
	core::kinematics::FoldTree const & ft
);


/// @brief find the jump connecting two continuous segments of a fold tree
/// @param[in] u Any vertex on the first segment.
/// @param[in] v Any vertex on the second segment.
/// @param[in] ft The fold tree to query.
/// @return the jump number connecting the two segments, 0 if no such jump
core::Size find_connecting_jump(
	core::Size const u,
	core::Size const v,
	core::kinematics::FoldTree const & ft
);


/// @brief remove a cutpoint, merging the two sections that are adjacent to the cut
/// @return true if cutpoint removed, false if not a cutpoint
bool remove_cutpoint(
	core::Size const v,
	core::kinematics::FoldTree & ft
);


/// @brief seal a fold tree by removing all specified cutpoints
/// @param[in] cutpoints Cutpoints to remove.
/// @param[in,out] ft The input tree.
/// @return A new tree with cutpoints removed and all other topology kept
///  constant.
void remove_cutpoints(
	utility::vector1< core::Size > const & cutpoints,
	core::kinematics::FoldTree & ft
);


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
);


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
);


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
);


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
);


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
	bool const keep_stub_in_residue = false
);


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
);

void make_star_foldtree( core::pose::Pose & pose, protocols::loops::Loops loops);

void jumps_and_cuts_from_pose( core::pose::Pose & pose, utility::vector1< std::pair< core::Size, core::Size > > & jumps, utility::vector1< core::Size > & cuts);


} // methods
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_methods_fold_tree_functions_HH */
