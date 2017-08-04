// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/FoldTree.cc
/// @brief  Fold tree class
/// @author Phil Bradley

// Unit headers
#include <core/kinematics/FoldTree.hh>

// Rosetta Headers
#include <core/kinematics/util.hh>
#include <core/id/SequenceMapping.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/py/PyAssert.hh>

// ObjexxFCL formating
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ Headers
#include <algorithm>
#include <list>
#include <string>
#include <iostream>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>
#include <utility/serialization/ObjexxFCL/FArray1D.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#endif // SERIALIZATION

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

static THREAD_LOCAL basic::Tracer TR( "core.kinematics.FoldTree" );

namespace core {
namespace kinematics {

// @brief Auto-generated virtual destructor
FoldTree::~FoldTree() {}


///////////////////////////////////////////////////////////////////////
// ensure that stored data depending on the topology are up-to-date
// modifications by add_edge, delete_edge, reorder, etc, etc
// are indicated by setting new_topology and/or new_order to true

bool
operator==(
	FoldTree const & a,
	FoldTree const & b
)
{
	if ( a.edge_list_.size() != b.edge_list_.size() ) return false;
	for ( auto
			edge_a( a.edge_list_.begin() ), a_end( a.edge_list_.end() ),
			edge_b( b.edge_list_.begin() );
			edge_a != a_end; ++edge_a, ++edge_b ) {
		if ( *edge_a != *edge_b ) return false;
	}
	return true;
}

bool operator!=(const FoldTree& a, const FoldTree& b) {
	return !(a == b);
}

bool
FoldTree::is_equivalent(
	FoldTree const & b
) const {
	if ( this->root() != b.root() ) { return false; }
	if ( this->edge_list_.size() != b.edge_list_.size() ) { return false; }

	// make copy, as we're modifying them with sorting
	EdgeList mine( this->edge_list_ );
	EdgeList theirs( b.edge_list_ );

	// Edge has operator < defined, so we can sort properly.
	std::sort( mine.begin(), mine.end() );
	std::sort( theirs.begin(), theirs.end() );

	for ( core::Size ii(0); ii < mine.size(); ++ii ) {
		if ( mine[ii] != theirs[ii] ) {
			return false;
		}
	}

	return true;
}

/// @details Computes a fixed-length, hash-based identifier for this FoldTree,
/// permitting efficient comparison of a pair of FoldTrees. The need for this
/// functionality arose from a desire to reuse an object that was unaware of
/// changes to the FoldTree. Rather than perform a costly deep comparison by
/// evaluating edge lists, we wanted a simple method for quickly testing
/// whether the naive object should be reinstantiated. This method is most
/// useful in situations where there are many edges in the FoldTree.
size_t FoldTree::hash_value() const {
	// compute hash(repr)
	return hasher( to_string() );
}

/// @details easy output of string
std::string
FoldTree::to_string() const{
	// obtain a string representation of the FoldTree
	std::string repr;
	std::stringstream ss;
	ss << *this;
	repr = ss.str();
	return repr;
}

/// @details  Delete self-edges in the foldtree, allowing the edge 1->1 for a single residue tree
void
FoldTree::delete_self_edges()
{
	if ( edge_list_.size() <= 1 ) return;

	utility::vector1< Edge > delete_me;

	for ( auto & it : edge_list_ ) {
		if ( it.start() == it.stop() ) delete_me.push_back( it );
	}

	if ( !delete_me.empty() ) {
		for ( Size i=1; i<= delete_me.size(); ++i ) {
			delete_edge( delete_me[i] );
		}
		new_topology = true;
	}
}


/// @details  Delete a sequence position from a foldtree. If the residue is a jump point or the root of the tree,
/// we will have to rearrange the topology.
void
FoldTree::delete_seqpos( int const seqpos )
{
	PyAssert( (seqpos>0) && (seqpos<=static_cast<int>(nres())), "FoldTree::delete_seqpos( int const seqpos ): input variable seqpos has a meaningless value");
	if ( is_jump_point( seqpos ) || is_root( seqpos ) ) {
		delete_jump_seqpos( seqpos );
	} else {
		delete_seqpos_simple( seqpos );
	}
}

/// @details  Useful for removing a loop modeling jump+cut
/// @note  This will trigger a renumbering of the jumps if jump_number < num_jump()
void
FoldTree::delete_jump_and_intervening_cutpoint( int const jump_number )
{
	delete_jump_and_intervening_cutpoint( jump_point(1,jump_number), jump_point(2,jump_number) );
}

/// @details  Useful for removing a loop modeling jump+cut
/// @note  This will trigger a renumbering of the jumps if jump number < num_jump()
void
FoldTree::delete_jump_and_intervening_cutpoint( int jump_begin, int jump_end, Size cut/*= 0*/ )
{
	PyAssert( is_jump_point(jump_begin),
		"FoldTree::delete_jump_and_intervening_cutpoint( int jump_begin , int jump_end ): "
		"input variable jump_begin has a meaningless value");
	PyAssert( is_jump_point(jump_end),
		"FoldTree::delete_jump_and_intervening_cutpoint( int jump_begin , int jump_end ): "
		"input variable jump_end has a meaningless value");

	debug_assert( is_jump_point( jump_begin ) && is_jump_point( jump_end ) );

	if ( jump_begin > jump_end ) {
		int tmp( jump_begin );
		jump_begin = jump_end;
		jump_end = tmp;
	}

	Size const old_num_jump( num_jump() );
	Size const old_root( root() );
	Size jump_number(0);

	for ( Size i=1; i<= old_num_jump; ++i ) {
		if ( jump_point(1,i) == jump_begin && jump_point(2,i) == jump_end ) jump_number = i;
	}
	debug_assert( jump_number ); // will fail if there's no jump between jump_begin and jump_end

	// look for a cutpoint between if unspecified
	if ( !cut ) {
		//Size cut(0);
		for ( int i=jump_begin; i<= jump_end-1; ++i ) {
			if ( is_cutpoint( i ) ) {
				if ( cut ) utility_exit_with_message( "multiple cutpoints between jump positions!" );
				cut = i;
			}
		}
	}

	debug_assert( cut );

	add_edge( cut, cut+1, Edge::PEPTIDE );
	delete_unordered_edge( jump_begin, jump_end, jump_number );

	if ( jump_number < old_num_jump ) {
		TR << "delete_jump_and_intervening_cutpoint: renumbering the remaining jumps!" << std::endl;
		for ( auto & it : *this ) {
			if ( it.is_jump() && (Size)(it.label()) > jump_number ) it.label() = it.label() - 1;
		}
	}

	delete_extra_vertices();

	reorder( old_root );
	debug_assert( check_fold_tree() );
}


/// @details  Slide a cutpoint from one position to another.
void
FoldTree::slide_cutpoint( Size const current_cut, Size const target_cut )
{
	debug_assert( is_cutpoint( current_cut ) && !is_cutpoint( target_cut ) );
	PyAssert( (is_cutpoint( current_cut )) && (!is_cutpoint( target_cut )), "FoldTree::slide_cutpoint( Size const current_cut , Size const target_cut ): input variable current_cut or target_cut has a meaningless value");
	Size const current_root( root() );

	TR.Trace << "slide_cutpoint: current= "<< current_cut << " target= " << target_cut << ' ' << *this;
	add_vertex( target_cut );
	add_vertex( target_cut+1 );
	add_edge( current_cut, current_cut+1, Edge::PEPTIDE );
	delete_unordered_edge( target_cut, target_cut+1, Edge::PEPTIDE );
	reorder( current_root );
	delete_extra_vertices();
	delete_self_edges();
	debug_assert( check_fold_tree() );
	TR.Trace << "slide_cutpoint: final tree= " << *this;
}

/// @details  Switches a given jump to connect two different residues.
void
FoldTree::slide_jump( Size const jump_number, Size const new_res1, Size const new_res2 )
{
	TR.Debug << "slide_jump: starting tree= " << *this << std::endl;
	Edge const old_jump_edge( jump_edge( jump_number ) );
	Size const pos1( std::min( new_res1, new_res2 ) );
	Size const pos2( std::max( new_res1, new_res2 ) );
	utility::vector1< Edge > new_edges, remove_edges;
	Size const original_root( root() );
	for ( auto & it : *this ) {
		core::Size const start( (core::Size)it.start() );
		core::Size const stop( (core::Size)it.stop() );
		if ( it.label() != Edge::PEPTIDE ) continue;
		if ( (start <= pos1 && stop >= pos1) || (stop <= pos1 && start >= pos1) ) { // edges not always in sequential order (eg - jump in middle of chain)
			//TR.Debug << "start-pos1-stop " << start <<" " << pos1 << " " << stop << std::endl;
			new_edges.push_back( Edge( start, pos1, Edge::PEPTIDE ) );
			new_edges.push_back( Edge( pos1, stop, Edge::PEPTIDE ) );
			remove_edges.push_back( it );
		}
		if ( (start <= pos2 && stop >= pos2) || (stop <= pos2 && start >= pos2) ) {
			//TR.Debug << "start-pos2-stop " << start <<" " << pos2 << " " << stop << std::endl;
			new_edges.push_back( Edge( start, pos2, Edge::PEPTIDE ) );
			new_edges.push_back( Edge( pos2, stop, Edge::PEPTIDE ) );
			remove_edges.push_back( it );
		}
	}
	for ( auto & remove_edge : remove_edges ) {
		delete_edge( remove_edge );
	}
	for ( auto & new_edge : new_edges ) {
		if ( (core::Size)new_edge.start() == original_root ) {
			prepend_edge( new_edge ); // preserve root!
		} else {
			add_edge( new_edge );
		}
	}
	runtime_assert( (core::Size)root() == original_root );

	Edge const new_jump_edge( pos1, pos2, jump_number );
	delete_edge( old_jump_edge );
	if ( (core::Size)new_jump_edge.start() == original_root ) {
		prepend_edge( new_jump_edge ); // preserve root!
	} else {
		add_edge( new_jump_edge );
	}
	delete_extra_vertices();
	delete_self_edges();
	runtime_assert( (core::Size)root() == original_root );
	TR.Debug << "slide_jump: final tree= " << *this << std::endl;
}

/// @details  Delete a position from the foldtree that is either a jump_point or the root. Will require some
/// rearranging of the topology.
/// LOGIC: note that there are 0 or 1 incoming edges to this vertex
/// need a replacement vertex for edges involving this guy: (note that this number will have to be adjusted)
/// I.   if seqpos is polymer residue (ie contained in a polymer edge)
///      1. if there's an incoming polymer segment, choose the previous rsd in this segment
///      2. choose the next residue in the first outgoing polymer edge
/// II.  if there's an incoming edge, choose the start of this edge
/// III. (non-polymer root residue) choose stop of first edge in foldtree
void
FoldTree::delete_jump_seqpos( int const seqpos )
{

	TR.Trace << "delete_jump_seqpos: " << seqpos << ' ' << *this << std::endl;

	Size const old_nres( nres() );
	debug_assert( is_jump_point( seqpos ) || is_root( seqpos ) );
	PyAssert( ( is_jump_point( seqpos )) || ( is_root( seqpos ) ), "FoldTree::delete_jump_seqpos( int const seqpos): input variable seqpos has a meaningless value");

	/// this could confuse us
	delete_self_edges();

	/// 1st order of business: determine the replacement vertex for this guy
	Edge incoming_edge;
	utility::vector1< Edge > outgoing_edges, outgoing_polymer_edges;

	for ( const_iterator it = begin(), ite = end(); it != ite; ++it ) {
		if ( it->stop() == seqpos ) incoming_edge = *it;
		else if ( it->start() == seqpos ) {
			outgoing_edges.push_back( *it );
			if ( it->is_polymer() ) {
				outgoing_polymer_edges.push_back( *it );
			}
		}
	}

	bool const is_root_rsd( is_root( seqpos ) );
	debug_assert( is_root_rsd == !incoming_edge.valid() );
	bool const is_polymer_rsd( ( !is_root_rsd && incoming_edge.is_polymer() ) || !outgoing_polymer_edges.empty());

	Size new_seqpos( 0 ), deleted_jump_number( 0 );
	if ( is_polymer_rsd ) {
		if ( !is_root_rsd && incoming_edge.is_polymer() ) {
			// case I.1 (see comments at start)
			debug_assert( incoming_edge.stop() == seqpos );
			new_seqpos = seqpos - incoming_edge.polymer_direction();
		} else {
			// case I.2
			Edge const & e( outgoing_polymer_edges.front() );
			new_seqpos = seqpos + e.polymer_direction();
		}
	} else {
		if ( is_root_rsd ) {
			// case III
			new_seqpos = begin()->stop();
			if ( begin()->is_jump() ) deleted_jump_number = begin()->label();
		} else {
			// case II
			new_seqpos = incoming_edge.start();
			if ( incoming_edge.is_jump() ) deleted_jump_number = incoming_edge.label();
		}
	}
	debug_assert( new_seqpos );

	// NO MODIFICATIONS TO THE TOPOLOGY OR ORDER OF THE TREE UP TO THIS POINT (except for delete_self_edges at start) ///
	Size const old_num_jump( num_jump() );
	if ( deleted_jump_number && deleted_jump_number != old_num_jump ) {
		// have to relabel
		TR.Info << "delete_jump_seqpos: deleting jump " << deleted_jump_number << std::endl;
		TR.Info << "delete_jump_seqpos: renumbering jump " << old_num_jump << " to " << deleted_jump_number << std::endl;
		jump_edge( old_num_jump        ).label() = deleted_jump_number;
		jump_edge( deleted_jump_number ).label() = old_num_jump;
		new_topology = true;
	}

	// now remap the edges that contain seqpos
	for ( auto & it : *this ) {
		if ( it.stop() == seqpos ) it.stop() = new_seqpos;
		else if ( it.start() == seqpos ) it.start() = new_seqpos;
	}
	new_topology = true;

	// delete the edges with start == stop after removing seqpos
	delete_self_edges();

	// now adjust the sequence numbering of all edges to reflect deletion of seqpos
	id::SequenceMapping old2new( id::SequenceMapping::identity( old_nres ) );
	old2new.delete_target_residue( seqpos );

	apply_sequence_mapping( old2new );

	// special case
	if ( edge_list_.size() == 1 && begin()->is_jump() &&
			begin()->start() == begin()->stop() ) {
		edge_list_.clear();
		edge_list_.push_back( Edge( 1, 1, Edge::PEPTIDE ) );
		new_topology = true;
	}

	debug_assert( check_fold_tree() );

	TR.Trace << "delete_jump_seqpos: after " << seqpos << ' ' << new_seqpos << *this << std::endl;
}

/// @details  Get the number of the jump that builds (connects to) a given residue
int
FoldTree::get_jump_that_builds_residue( int const seqpos ) const
{
	check_topology();
	Edge const & edge( get_residue_edge( seqpos ) );
	if ( !edge.is_jump() ) utility_exit_with_message( "Error in core::kinematics::FoldTree::get_jump_that_builds_residue(): This residue is not the child of (built by) a jump!" );
	return edge.label();
}

/// @brief  Get the residue that is immediately upstream of this residue (and tell us whether connection is jump or bond).
int
FoldTree::get_parent_residue( int const seqpos, bool & connected_by_jump ) const {
	// (1) Root
	if ( seqpos == root() ) {
		connected_by_jump = true;
		return 0;
	}

	// (2) Jump.
	kinematics::Edge const & edge = get_residue_edge( seqpos );
	if ( edge.is_jump() ) {
		connected_by_jump = true;
		runtime_assert( edge.stop() == seqpos );
		return edge.start();
	}

	// (3) Covalent connection
	int parent_res( 0 );
	if ( edge.start() < seqpos ) {
		runtime_assert( edge.start() < edge.stop() );
		parent_res = seqpos - 1;
		//  if ( parent_res < edge.start() ) return 0;
	} else {
		runtime_assert( edge.start() > edge.stop() );
		parent_res = seqpos + 1;
		//  if ( parent_res > edge.start() ) return 0;
	}
	return parent_res;
}

/// @brief  Get the residue that is immediately upstream of this residue.
int
FoldTree::get_parent_residue( int const seqpos ) const {
	bool connected_by_jump( true );
	return get_parent_residue( seqpos, connected_by_jump );
}

/// @details  Delete a sequence position from a foldtree. This will not work
/// at positions that are jump points, ie start or stop vertices for jump edges. (or the root of the tree!)
/// So basically only works for polymer residues.
void
FoldTree::delete_seqpos_simple( int const seqpos )
{
	TR.Trace << "delete_seqpos_simple: before " << seqpos << ' ' << *this << std::endl;
	debug_assert( !is_jump_point( seqpos ) && !is_root( seqpos ) );

	Size const old_nres( nres() );

	// always a good idea for safety
	delete_self_edges();

	// first remap the edge (if it exists) that contains seqpos as a vertex
	// do this before renumbering everything
	for ( auto it = begin(), ite = end(); it != ite; ++it ) {
		debug_assert( it->start() != seqpos );
		if ( it->stop() == seqpos ) {
			debug_assert( it->is_polymer() );
			it->stop() = it->stop() - it->polymer_direction();
			new_topology = true;
			break; // there should be only one incoming edge
		}
	}

	// now adjust the sequence numbering of all edges to reflect deletion of seqpos
	id::SequenceMapping old2new( id::SequenceMapping::identity( old_nres ) );
	old2new.delete_target_residue( seqpos );debug_assert( !old2new[ seqpos ] );

	apply_sequence_mapping( old2new );

	// they may have been introduced
	delete_self_edges();

	// sanity
	debug_assert( check_fold_tree() );

	TR.Trace << "delete_seqpos_simple: after " << seqpos << ' ' << *this << std::endl;
}


/// @details  Renumber all vertices according to an input sequence mapping
void
FoldTree::apply_sequence_mapping( id::SequenceMapping const & old2new )
{

	for ( auto & it : *this ) {
		it.start() = old2new[ it.start() ];
		it.stop () = old2new[ it.stop () ];
	}

	new_topology = true;
}


/////////////////////////////////////////////////////////////////////////////
//
// the residue at position seqpos moves to position seqpos+1
//
// vertices remapped, only question is cutpoint at seqpos-1, should it move to seqpos?
//
/// @details (ie between current rsds seqpos-1 and seqpos, so that the sequence position of the new residue is seqpos)
/// if seqpos-1 is a cutpoint in the current fold_tree -- we have a choice about how to connect the new
/// residue: it could be joined to the preceding peptide segment (join to seqpos-1) or to the following
/// segment (joined to the residue currently at seqpos). join_upper and join_lower control the behavior in
/// this case.
///
/// @note  seqpos may be greater than current nres, ie we may be "inserting" at end
void
FoldTree::insert_polymer_residue(
	int const seqpos,
	bool const join_lower,
	bool const join_upper
)
{
	int const old_size( nres() );
	utility::vector1< int > old2new( old_size, 0 );
	for ( int i=1; i<= old_size; ++i ) {
		if ( i<seqpos ) old2new[i] = i;
		else old2new[i] = i+1;
	}

	// this is the only case where we actually need to know join_lower and join_upper,
	// since in this case we have a choice of how to attach the new residue
	bool const special_case( is_cutpoint( seqpos -1 ) );

	// in either of these cases we actually need to add a new edge to the tree, since we are gluing
	// to a jump_point at a cutpoint, which means that there's no polymer edge extending in our
	// direction yet
	bool const seqpos_minus_1_is_vertex( seqpos > 1         && (is_jump_point( seqpos-1 ) || is_root(seqpos-1)));
	bool const seqpos_is_vertex        ( seqpos <= old_size && (is_jump_point( seqpos   ) || is_root(seqpos)));
	bool add_edge_lower( special_case && join_lower && seqpos_minus_1_is_vertex );
	bool add_edge_upper( special_case && join_upper && seqpos_is_vertex         );

	for ( auto & it : edge_list_ ) {
		it.start() = old2new[ it.start() ];

		if ( special_case ) {
			if ( join_lower && it.is_polymer() && it.stop() == seqpos - 1 && !seqpos_minus_1_is_vertex ) {
				it.stop() = seqpos; // extend this segment to include seqpos
				debug_assert( !add_edge_lower );
				continue;
			}
			if ( join_upper && it.is_polymer() && it.stop() == seqpos     && !seqpos_is_vertex ) {
				// dont remap seqpos to seqpos+1 as we would otherwise do, thereby extending this peptide edge to include seqpos
				debug_assert( !add_edge_upper );
				continue;
			}
		}

		it.stop () = old2new[ it.stop () ];
	}

	new_topology = true;

	if ( add_edge_lower ) {
		add_edge( seqpos-1, seqpos, Edge::PEPTIDE );
	} else if ( add_edge_upper ) {
		add_edge( seqpos+1, seqpos, Edge::PEPTIDE );
	}

	debug_assert( join_lower == !is_cutpoint( seqpos-1 ) && join_upper == !is_cutpoint( seqpos ) );
}


/////////////////////////////////////////////////////////////////////////////
//
// the residue at position seqpos moves to position seqpos+1
//
// vertices remapped, only question is cutpoint at seqpos-1, should it move to seqpos?
//
/// @details (ie between current rsds seqpos-1 and seqpos, so that the sequence position of the new residue is seqpos)
/// if seqpos-1 is a cutpoint in the current fold_tree -- we have a choice about how to connect the new
/// residue: it could be joined to the preceding peptide segment (join to seqpos-1) or to the following
/// segment (joined to the residue currently at seqpos). join_upper and join_lower control the behavior in
/// this case.
///
/// @note  seqpos may be greater than current nres, ie we may be "inserting" at end
void
FoldTree::insert_residue_by_chemical_bond(
	int const seqpos,
	int const anchor_residue,
	std::string const& anchor_atom,
	std::string const& root_atom )
{
	debug_assert( is_cutpoint( seqpos - 1 ) );
	int const old_size( nres() );
	//int const new_jump_number( num_jump() + 1 );

	utility::vector1< int > old2new( old_size, 0 );
	for ( int i=1; i<= old_size; ++i ) {
		if ( i<seqpos ) old2new[i] = i;
		else old2new[i] = i+1;
	}
	int anchor_pos = old2new[ anchor_residue ];

	for ( auto & it : edge_list_ ) {
		it.start() = old2new[ it.start() ];
		it.stop () = old2new[ it.stop () ];
	}

	add_edge( Edge( anchor_pos, seqpos, anchor_atom, root_atom ) ); // different from add_edge call in insert_residue_by_jump, here Edge constructor assumes it's a chemical edge
	debug_assert( check_fold_tree() );
	new_topology = true;
}

/////////////////////////////////////////////////////////////////////////////
/// @details  Insert a new residue into the tree at position seqpos, anchoring it
/// to the rest of the tree by a jump
///
/// the residue at position seqpos moves to position seqpos+1
///
/// vertices remapped, only question is cutpoint at seqpos-1, should it move to seqpos?
void
FoldTree::insert_residue_by_jump(
	int const seqpos,
	int anchor_pos, // in the current numbering system
	std::string const& anchor_atom, // = 0
	std::string const& root_atom // = 0
)
{
	debug_assert( is_cutpoint( seqpos - 1 ) );
	int const old_size( nres() );
	int const new_jump_number( num_jump() + 1 );

	utility::vector1< int > old2new( old_size, 0 );
	for ( int i=1; i<= old_size; ++i ) {
		if ( i<seqpos ) old2new[i] = i;
		else old2new[i] = i+1;
	}
	anchor_pos = old2new[ anchor_pos ];

	for ( auto & it : edge_list_ ) {
		it.start() = old2new[ it.start() ];
		it.stop () = old2new[ it.stop () ];
	}

	add_edge( anchor_pos, seqpos, new_jump_number );
	debug_assert( check_fold_tree() );
	new_topology = true;
	if ( anchor_atom.size() ) {
		debug_assert( root_atom.size() );
		set_jump_atoms( new_jump_number, anchor_atom, root_atom );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Insert another fold_tree as a subtree. Residues are inserted as a contiguous block beginning at
/// insert_seqpos. Jumps are inserted as a contiguous block beginning at insert_jumppos.
/// Note that insert_seqpos could be equal to nres()+1, ie subtree is being appended at the end.
/// The jump anchoring subtree runs from the residue currently numbered "anchor_pos" to the residue
/// insert_seqpos + subtree.root() - 1, and has label/number anchor_jump_number
void
FoldTree::insert_fold_tree_by_jump(
	FoldTree const & subtree,
	int const insert_seqpos,         // rsd 1 in subtree goes here
	int const insert_jumppos,        // jump 1 in subtree goes here
	int const anchor_pos,            // in the old numbering system
	int anchor_jump_number,          // in the new jump numbering system, default=0
	std::string const & anchor_atom, // could be ""
	std::string const & root_atom )  // ditto
{
	if ( !is_cutpoint( anchor_pos -1 ) && !is_cutpoint( anchor_pos ) ) {
		TR.Warning << "insert_fold_tree_by_jump: anchor_pos is not a vertex of the tree!! anchor_pos= " << anchor_pos <<
			std::endl;
		TR.Warning << "*this= " << *this << std::endl;
		TR.Warning << "Adding anchor_pos as a new vertex!" << std::endl;
		add_vertex( anchor_pos );
	}

	int const old_nres( nres() );
	int const old_njump( num_jump() );
	int const insert_nres( subtree.nres() );
	int const insert_njump( subtree.num_jump() );
	int const new_njump( old_njump + insert_njump + 1 );
	if ( !anchor_jump_number ) anchor_jump_number = new_njump; // put the anchor_jump at the end

	TR.Trace << "insert_fold_tree_by_jump: insert_seqpos= " << insert_seqpos << " insert_jumppos= " << insert_jumppos <<
		" anchor_pos= " << anchor_pos << " anchor_jump_number= " << anchor_jump_number <<
		" old_nres= " << old_nres << " old_njump= " << old_njump <<
		" insert_nres= " << insert_nres << " insert_njump= " << insert_njump << std::endl;

	TR.Trace << "insert_fold_tree_by_jump: old_fold_tree: " << *this << std::endl;
	TR.Trace << "insert_fold_tree_by_jump: subtree: " << subtree <<std::endl;

	utility::vector1< int > old2new_res ( old_nres , 0 );
	utility::vector1< int > old2new_jump( old_njump, 0 );
	for ( int i=1; i<= old_nres; ++i ) {
		if ( i < insert_seqpos ) old2new_res[ i ] = i;
		else old2new_res[ i ] = i+insert_nres;
	}
	if ( anchor_jump_number >= insert_jumppos && anchor_jump_number < insert_jumppos + insert_njump ) {
		utility_exit_with_message("anchor_jump_number cant fall within the range of inserted jumps!");
	}
	/// for debugging
	utility::vector1< bool > labeled( new_njump, false );
	labeled[ anchor_jump_number ] = true;
	for ( int i= insert_jumppos; i< insert_jumppos+insert_njump; ++i ) labeled[i] = true;
	for ( int i=1; i<= old_njump; ++i ) {
		int o2n( i < insert_jumppos ? i : i+insert_njump );
		if ( o2n >= anchor_jump_number ) ++o2n;
		if ( o2n >= insert_jumppos && o2n < insert_jumppos+insert_njump ) o2n = insert_jumppos + insert_njump;
		debug_assert( !labeled[ o2n ] );
		labeled[ o2n ] = true;
		old2new_jump[ i ] = o2n;
	}

	// remap the numbers of our edges
	for ( auto & it : edge_list_ ) {
		it.start() = old2new_res[ it.start() ];
		it.stop () = old2new_res[ it.stop () ];
		if ( it.is_jump() ) {
			it.label() = old2new_jump[ it.label() ];
		}
	}

	// add the anchoring jump
	int const root_pos( subtree.root() + insert_seqpos - 1 );
	edge_list_.push_back( Edge( old2new_res[ anchor_pos ], root_pos, anchor_jump_number ) );

	// now append the edges from subtree
	for ( auto new_edge : subtree ) {
		new_edge.start() = new_edge.start() + insert_seqpos-1;
		new_edge.stop () = new_edge.stop () + insert_seqpos-1;
		if ( new_edge.is_jump() ) new_edge.label() = new_edge.label() + insert_jumppos-1;

		if ( new_edge.start() != new_edge.stop() ) {
			// fixes the case where a single-residue FT creates an invalid single-residue peptide edge
			// led to by a jump.
			edge_list_.push_back( new_edge );
		}
	}

	new_topology = true;

	if ( anchor_atom.size() ) {
		debug_assert( root_atom.size() );
		set_jump_atoms( anchor_jump_number, anchor_atom, root_atom );
	}

	TR.Trace << "insert_fold_tree_by_jump: new fold_tree: " << *this << std::endl;
	debug_assert( check_fold_tree() );
}


/////////////////////////////////////////////////////////////////////////////
/// @details  Add a new residue -- either polymer or jump -- to the end of the tree
/// if the new residue is appended by polymer connection, add it at the end of polymer.
/// if the new residue is appended by jump connection, the cutpoint is the original polymer end
void
FoldTree::append_residue(
	bool const attach_by_jump, // = false,
	int const jump_anchor_residue, // = 0,
	std::string const& jump_upstream_atom, // = "",
	std::string const& jump_downstream_atom // = ""
)
{
	// NOTE -- all public method calls, no worries about data being
	// up to date

	int const old_nres( nres() );

	// add a new edge (maybe temporary)
	add_edge( old_nres, old_nres+1, Edge::PEPTIDE );

	if ( old_nres == 1 ) delete_self_edges();

	if ( attach_by_jump ) {
		new_jump( jump_anchor_residue, old_nres + 1, old_nres );
		if ( jump_upstream_atom.size() ) {
			set_jump_atoms( num_jump(), jump_upstream_atom, jump_downstream_atom );
		}
	} else {
		delete_extra_vertices();
	}
}

//////////////////////////////////////////////////////////////////////
void
FoldTree::append_residue_by_chemical_bond(
	int const anchor_residue,
	std::string const& anchor_atom,
	std::string const& root_atom
)
{
	int const old_nres( nres() );

	// chemical edge constructor
	add_edge( Edge( anchor_residue, old_nres+1, anchor_atom, root_atom ) );

	new_topology = true;

	debug_assert( is_cutpoint( old_nres ) );
}


/////////////////////////////////////////////////////////////////////////////
/// @details
/// do not allow self edge and do not order the edge list here
void
FoldTree::add_edge(
	int const start,
	int const stop,
	int const label
)
{
	// jump out if self-edge, exception: a one residue pose (happens at setup of pdb-read)
	if ( start == stop && !( start == 1 && edge_list_.empty() ) ) return;
	runtime_assert( start > 0 );
	runtime_assert( stop  > 0 );
	new_topology = true; // book-keeping
	edge_list_.push_back( Edge( start, stop, label ) );
}

/// @details This 'add_edge' calls the edge constructor with the same args and is
/// used to form chemical edges.
void
FoldTree::add_edge(
	int const start,
	int const stop,
	std::string const & start_atom,
	std::string const & stop_atom
)
{
	// jump out if self-edge, exception: a one residue pose (happens at setup of pdb-read)
	if ( start == stop && !( start == 1 && edge_list_.empty() ) ) return;
	runtime_assert( start > 0 );
	runtime_assert( stop  > 0 );
	new_topology = true; // book-keeping
	edge_list_.push_back( Edge( start, stop, start_atom, stop_atom ) );
}


/////////////////////////////////////////////////////////////////////////////
/// @details  Does not ensure proper folding order
void
FoldTree::add_edge( Edge const & new_edge )
{
	new_topology = true; // book-keeping
	edge_list_.push_back( new_edge );
}

/////////////////////////////////////////////////////////////////////////////
/// @brief Prepend the edge  <new_edge>. Useful alternative to add_edge for setting root.
/// @details  Does not ensure proper folding order
void
FoldTree::prepend_edge(
	Edge const & new_edge
)
{
	new_topology = true; // book-keeping
	edge_list_.insert( edge_list_.begin(), new_edge );
}


/// @details  This function is used primarily to restore CHEMICAL Edges that
/// have been replaced by JUMPs after other FoldTree manipulations.
/// @author   Labonte <JWLabonte@jhu.edu>
void
FoldTree::replace_edge( Edge const & old_edge, Edge const & replacement_edge  )
{
	new_topology = true;
	bool found( false );
	for ( auto & edge : edge_list_ ) {
		if ( edge == old_edge ) {
			edge = replacement_edge;
			found = true;
			break;
		}
	}
	if ( ! found ) {
		TR.Fatal << "FoldTree::replace_edge(...): edge not in tree." << std::endl;
		utility_exit();
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details
/// die if the iterator is out of range
void
FoldTree::delete_edge( FoldTree::iterator edge )
{
	new_topology = true; // book-keeping
	if ( edge < edge_list_.begin() || edge >= edge_list_.end() ) {
		TR.Fatal << "FoldTree::delete_edge(...) edge not contained in edge_list_." << std::endl;
		utility_exit();
	}
	edge_list_.erase( edge );
}

/////////////////////////////////////////////////////////////////////////////
/// @details
void
FoldTree::delete_edge( Edge const & edge )
{
	delete_unordered_edge( edge.start(), edge.stop(), edge.label() );
}

/////////////////////////////////////////////////////////////////////////////
/// @details
/// needs to match start residue number, end residue number and label index.
/// abort if the edge is not found.
void
FoldTree::delete_unordered_edge(
	int const start,
	int const stop,
	int const label
)
{
	new_topology = true;
	bool found(false);
	for ( auto it=edge_list_.begin(), it_end=edge_list_.end();
			it != it_end; ++it ) {
		if ( it->label() == label &&
				( ( it->start() == start && it->stop() == stop ) ||
				( it->start() == stop  && it->stop() == start ) ) ) {
			edge_list_.erase( it );
			found = true;
			break;
		}
	}
	if ( !found ) {
		TR.Fatal << "FoldTree::delete_unordered_edge(...) edge not in tree: " <<
			' ' << start << ' ' << stop << ' ' << label << std::endl;
		TR.Fatal << *this;
		utility_exit();
	}
}


//////////////////////////////////////////////////////////////////
/// @details
/// it assumes that the segment is completely contained in a single edge
/// of the tree. Only edge_list is updated and new topology is set true.
/// No derived data is updated.
void
FoldTree::delete_segment(
	int const seg_begin,
	int const seg_end
)
{
	int const n2c(1), c2n(-1);

	TR.Info << "FoldTree::delete_segment: " << seg_begin << ' ' <<
		seg_end << std::endl;

	int const size( seg_end - seg_begin + 1 );

	FArray1D_int mapping( nres_, -1 );
	for ( int i=1; i<= nres_; ++i ) {
		if ( i < seg_begin ) {
			mapping(i) = i;
		} else if ( i > seg_end ) {
			mapping(i) = i - size;
		}
	}

	std::vector< Edge > new_edge_list_;

	for ( auto it = edge_list_.begin(), it_end = edge_list_.end();
			it != it_end; ++it ) {
		int pos1( mapping( it->start() ) );
		int pos2( mapping( it->stop()  ) );
		int const dir( it->start() < it->stop() ? n2c :c2n );
		if ( pos1 == -1 || pos2 == -1 ) debug_assert( it->is_polymer() );

		if ( pos1 == -1 ) {
			debug_assert( (dir == n2c && it->start() == seg_begin && it->stop() > seg_end) ||
				(dir == c2n && it->start() == seg_end   && it->stop() < seg_begin) );
			if ( dir == n2c ) {
				pos1 = seg_begin; // n2c
			} else {
				pos1 = seg_begin - 1; // c2n
			}
		} else if ( pos2 == -1 ) {
			debug_assert( (dir == c2n && it->stop() == seg_begin && it->start() > seg_end) ||
				(dir == n2c && it->stop() == seg_end   && it->start() < seg_begin) );
			if ( dir == c2n ) {
				pos2 = seg_begin; // c2n
			} else {
				pos2 = seg_begin - 1; // n2c
			}
		}
		if ( pos1 != pos2 ) {
			new_edge_list_.push_back( Edge( pos1, pos2, it->label() ) );
			//std::cout << "remap edge: " << *it << ' ' << pos1 << ' ' << pos2 <<
			//std::endl;
		}
	}

	new_topology = true;
	edge_list_ = new_edge_list_;
}

// should this set new_topology TRUE???
//
//////////////////////////////////////////////////////////////////
/// @details  this is an internal function, used for testing if an edge is separating
void
FoldTree::update_edge_label(
	int const start,
	int const stop,
	int const old_label,
	int const new_label
)
{
	bool found(false);
	for ( auto & it : edge_list_ ) {
		if ( it.label() == old_label &&
				( ( it.start() == start && it.stop() == stop ) ||
				( it.start() == stop  && it.stop() == start ) ) ) {
			it.label() = new_label;
			found = true;
			break;
		}
	}
	if ( !found ) {
		TR.Fatal << "FoldTree::update_edge_label(...) edge not in tree: " <<
			' ' << start << ' ' << stop << ' ' << old_label << std::endl;
		TR.Fatal << *this;
		utility_exit();
	}
}


// should this set new_topology TRUE???
//
//////////////////////////////////////////////////////////////////
/// @details  this is an internal function, used for testing if an edge is separating
int
FoldTree::edge_label(
	int const start,
	int const stop
)
{
	bool found(false);
	int label(-1000);
	for ( auto & it : edge_list_ ) {
		if ( ( it.start() == start && it.stop() == stop ) ||
				( it.start() == stop  && it.stop() == start ) ) {
			label = it.label();
			found = true;
			break;
		}
	}
	if ( !found ) {
		TR.Fatal << "FoldTree::edge_label(...) edge not in tree: " <<
			' ' << start << ' ' << stop << ' '  << std::endl;
		TR.Fatal << *this;
		utility_exit();
	}
	return label;
}

void FoldTree::split_existing_edge_at_residue( int const resNo )
{
	get_residue_edge( resNo ); // ensures that an edge containing `resNo` exists
	add_vertex( resNo );
}

///////////////////////////////////////////////////////
/// @details
/// after deleting a jump, there may be vertices of the
/// tree which are neither jumps nor cutpoints. So delete them!
/// this will combine two adjacent short edges into a long one
void
FoldTree::delete_extra_vertices()
{
	// first get rid of any self-edges, eg left over from the tree when it was a single-residue tree
	delete_self_edges();

	// want to use is_jump_point_, is_cutpoint_
	// so ensure that this data is up-to-date
	check_topology();

	int const _root = root(); //need to keep that locally, since after killing the first edge the root() might have changed
	//
	while ( true ) {
		int kill_vertex( 0 );
		for ( auto it=edge_list_.begin(),it_end = edge_list_.end();
				it != it_end && !kill_vertex; ++it ) {
			// are either of these vertices extraneous ( neither jump-point nor cutpoint)
			//      TR.Trace << "superfluous? check Edge " << *it << std::endl;
			//    TR.Trace << "isjump: " << is_jump_point_( it->start() ) << " " << is_jump_point_( it->stop() ) << std::endl;
			//    TR.Trace << "iscut: "
			//         << is_cutpoint_( it->start() ) << " " << is_cutpoint_( it->start()-1  ) << " "
			//         << is_cutpoint_( it->stop() ) << " " << is_cutpoint_( it->stop()-1 ) << std::endl;
			//    TR.Trace << "isroot: " << ( _root == it->start() ) << " " << ( _root == it->stop() ) << std::endl;
			if ( ( it->start() != _root ) &&
					!is_jump_point_[ it->start() ] &&
					!is_cutpoint_  ( it->start() ) &&
					!is_cutpoint_  ( it->start()-1 ) ) kill_vertex = it->start();
			if ( ( it->stop() != _root ) &&
					!is_jump_point_[ it->stop() ] &&
					!is_cutpoint_  ( it->stop() ) &&
					!is_cutpoint_  ( it->stop()-1  ) ) kill_vertex = it->stop();
			if ( kill_vertex ) break;
		}
		if ( !kill_vertex ) break;
		//  TR.Trace << "delete vertex: " << kill_vertex << std::endl;
		int nbounds(0);
		FArray1D_int bounds(2,0);
		for ( auto it=edge_list_.begin(),it_end = edge_list_.end();
				it != it_end && nbounds<2; ++it ) {
			if ( it->start() == kill_vertex ) {
				nbounds++;
				bounds( nbounds ) = it->stop();
			} else if ( it->stop() == kill_vertex ) {
				nbounds++;
				bounds( nbounds ) = it->start();
			}
		}
		//std::cout << "bounds: " << bounds(1) << ' ' << bounds(2) << std::endl << *this;
		delete_unordered_edge( bounds(1), kill_vertex, Edge::PEPTIDE );
		delete_unordered_edge( bounds(2), kill_vertex, Edge::PEPTIDE );
		add_edge( bounds(1), bounds(2), Edge::PEPTIDE );
	}

	// need to call reorder since directionality of new edges may be incorrect
	reorder( _root );
}


/////////////////////////////////////////////////////////////////////////////
/// @details  Get the sequence position of the downstream vertex of the jump indicated by
/// the jump_number argument. Downstream means that if we traverse the tree starting at the root
/// then we hit that vertex second.
///
/// return 0 if failed
int
FoldTree::downstream_jump_residue( int const jump_number ) const
{
	check_order();
	debug_assert( jump_number >= 1 && jump_number <= num_jump_ );
	for ( const auto & it : edge_list_ ) {
		if ( it.label() == jump_number ) return it.stop();
	}
	return 0;
}

/////////////////////////////////////////////////////////////////////////////
/// @details  Get the sequence position of the upstream vertex of the jump indicated by
/// the jump_number argument. Upstream means that if we traverse the tree starting at the root
/// then we hit that vertex first.
///
/// return 0 if failed
int
FoldTree::upstream_jump_residue( int const jump_number ) const
{
	check_order();
	debug_assert( jump_number >= 1 && jump_number <= num_jump_ );
	for ( const auto & it : edge_list_ ) {
		if ( it.label() == jump_number ) return it.start();
	}
	return 0;
}

// will this be too slow? creates a new edge_list_ object,
// then at the end copies it into edge_list_; To
///////////////////////////////////////////////////////
/// @details  Reorder the tree so that start_residue is the new root.
/// returns false if no re-ordering allowed! To reorder
/// successfully, start_residue needs to be a vertex in the
/// original fold tree.
bool
FoldTree::reorder( int const start_residue, bool const verbose_if_fail /* = true */ )
{
	if ( new_topology ) update_nres(); // need nres

	// Scratch FArray1D -- only re-allocate upon changes of nres
	if ( nres_ != int( linked_.size1() ) ) {
		linked_.dimension( nres_ );
	}

	EdgeList new_edge_list_;

	// keep track of which residues have been added to the new list
	linked_ = false;

	linked_( start_residue) = true;

	bool new_member (true);

	while ( new_member ) {
		new_member = false;
		for ( const_iterator it = edge_list_.begin(),
				it_end = edge_list_.end(); it != it_end; ++it ) {
			Edge const& old_edge( *it );
			Edge edge( old_edge ); //makes sure everything in Edge gets copied !!!
			if ( linked_( edge.start() ) && !linked_( edge.stop() ) ) {
				new_edge_list_.push_back( edge );
				if ( !edge.is_polymer() ) {
					linked_( edge.stop() ) = true;
				} else {
					for ( core::Size i=std::min(edge.start(),edge.stop()), i_end=std::max(edge.start(),edge.stop()); i<=i_end; ++i ) linked_( i ) = true;
				}
				new_member = true;
			} else if ( linked_( edge.stop() ) && !linked_( edge.start() ) ) {
				//switch start / stop information
				edge.start() = old_edge.stop();
				edge.stop() = old_edge.start();
				edge.start_atom() = old_edge.stop_atom();
				edge.stop_atom() = old_edge.start_atom();
				new_edge_list_.push_back( edge );
				if ( !edge.is_polymer() ) {
					linked_( edge.stop() ) = true;
				} else {
					for ( core::Size i=std::min(edge.start(),edge.stop()), i_end=std::max(edge.start(),edge.stop()); i<=i_end; ++i ) linked_( i ) = true;
				}
				new_member = true;
			}
		}
	}// while ( new_member )

	if ( new_edge_list_.size() != edge_list_.size() ) {
		if ( verbose_if_fail && nres() > 1 /* nres = 1 is innocuous edge case*/ ) {
			TR.Error << "FoldTree::reorder( " << start_residue << " ) failed, new/old edge_list_ size mismatch" << std::endl;
			TR.Error << "old_edge_list.size() " << edge_list_.size() << "  new_edge_list.size()" << new_edge_list_.size() << std::endl;
			TR.Error << *this << std::endl;

			//  TR.Error << "show old edge list " << std::endl;
			//  for( FoldTree::const_iterator it(edge_list_.begin()), end(edge_list_.end()); it!=end; ++it){
			//   TR.Error << *it << std::endl;
			//  }

			//  TR.Error << "show new edge list " << std::endl;
			//  for( FoldTree::const_iterator it(new_edge_list_.begin()), end(new_edge_list_.end()); it!=end; ++it){
			//   TR.Error << *it << std::endl;
			//  }
		}

		return false;
	}

	// slow: clear edge_list_ then copy from new_edge_list_...
	edge_list_ = new_edge_list_;
	reassign_atoms_for_intra_residue_stubs();
	new_order = true;
	return true; // success
} // FoldTree::reorder(...)

///////////////////////////////////////////////////////////////////////////////
/// @brief Make a simple, 1->total_residue tree.
void
FoldTree::simple_tree( int const nres_in )
{
	PyAssert( (nres_in>0), "FoldTree::simple_tree( int const nres_in ): input variable nres_in has a meaningless value");
	new_topology = true; // ensure that derived data are re-calculated
	edge_list_.clear();
	add_edge(1, nres_in, Edge::PEPTIDE); // from 1->total_residue
}

//////////////////////////////////////////////////////////////////////////////
/// @brief Returns true if this tree is a simple 1->total_residue FoldTree,
/// returns false otherwise.
bool
FoldTree::is_simple_tree() const {
	bool is_simple( false );
	if ( edge_list_.size() == 1 ) {
		auto e = begin();
		if ( e->start() == 1 && (Size) e->stop() == nres() && e->is_polymer() ) {
			is_simple = true;
		}
	}

	return is_simple;
}


//////////////////////////////////////////////////////////////////////////////
/// @details  After this call you're guaranteed that v is a vertex of the tree, ie not contained in the
/// interior of a peptide edge
void
FoldTree::add_vertex( int const v )
{
	// find the edge that contains new_cutpoint and cut it
	for ( auto it= edge_list_.begin(), ite= edge_list_.end(); it != ite; ++it ) {
		if ( it->is_polymer() &&
				( ( it->start() < v && it->stop () > v ) ||
				( it->stop () < v && it->start() > v ) ) ) {
			int const start( std::min( it->start(), it->stop()) );
			int const stop ( std::max( it->start(), it->stop()) );
			delete_edge( it );
			if ( start < v ) add_edge( start, v, Edge::PEPTIDE );
			if ( stop  > v ) add_edge( v,  stop, Edge::PEPTIDE );
			break;
		}
	}
	new_topology = true;
}


/////////////////////////////////////////////////////////////////////////////
/// @details  Add a new jump to an existing fold tree, returns the jump_number of the new jump.
int
FoldTree::new_jump(
	int const jump_pos1,
	int const jump_pos2,
	int const new_cutpoint
)
{
	debug_assert( !is_cutpoint( new_cutpoint ) );
	// otherwise causes an error with delete_unordered_edge
	PyAssert( !is_cutpoint( new_cutpoint ), "FoldTree::new_jump( int const new_cutpoint ): new_cutpoint is already a cutpoint!" );
	// Commenting out because it get in the way when building symmetry Pose
	//PyAssert( (jump_pos1-new_cutpoint)*(jump_pos2-new_cutpoint)<0, "FoldTree::new_jump( int const new_cutpoint ): new_cutpoint is not between the jump points!" );
	PyAssert( (jump_pos1>0) && (jump_pos1<=nres_), "FoldTree::new_jump( int const jump_pos1 ): jump_pos1 is out of range!" );
	PyAssert( (jump_pos2>0) && (jump_pos2<=nres_), "FoldTree::new_jump( int const jump_pos2 ): jump_pos2 is out of range!" );

	int const root( edge_list_.begin()->start() );
	int const new_jump_number( num_jump() + 1 );

	add_vertex( jump_pos1 );
	add_vertex( jump_pos2 );
	add_vertex( new_cutpoint   );
	add_vertex( new_cutpoint+1 );

	add_edge( jump_pos1, jump_pos2, new_jump_number );
	delete_unordered_edge( new_cutpoint, new_cutpoint+1, Edge::PEPTIDE );

	reorder( root );
	new_topology = true;

	debug_assert( is_cutpoint( new_cutpoint ) && is_jump_point( jump_pos1 ) && is_jump_point( jump_pos2 ) );

	return new_jump_number;
}

/// @details  Add a new jump to an existing fold tree, returns the jump_number of the new jump.
void
FoldTree::new_chemical_bond(
	int const anchor_pos,
	int const root_pos,
	std::string const & anchor_atom,
	std::string const & root_atom,
	int const new_cutpoint
)
{
	debug_assert( !is_cutpoint( new_cutpoint ) );

	int const root( edge_list_.begin()->start() );

	add_vertex( anchor_pos );
	add_vertex( root_pos   );
	add_vertex( new_cutpoint   );
	add_vertex( new_cutpoint+1 );

	Edge new_edge( anchor_pos, root_pos, anchor_atom, root_atom ); // not obvious that this is a chemical bond c-tor
	debug_assert( new_edge.is_chemical_bond() );
	add_edge( new_edge );
	delete_unordered_edge( new_cutpoint, new_cutpoint+1, Edge::PEPTIDE );

	reorder( root );
	new_topology = true;

	debug_assert( is_cutpoint( new_cutpoint ) && is_jump_point( anchor_pos ) && is_jump_point( root_pos ) );
}

/////////////////////////////////////////////////////////////////////////////
/// @details  Construct a new tree (self) from a set of jump points and cutpoints
/// this assumes that we can make a tree, ie that the number of cuts
/// is equal to the number of jumps.
///
/// @note The root vertex of new tree is 1.
bool
FoldTree::tree_from_jumps_and_cuts(
	int const nres_in,
	int const num_jump_in,
	FArray2D_int const & jump_point_in, /* 2 x njump */
	FArray1D_int const & cuts,
	int const root_in,
	bool const verbose /* = true */
)
{
	// tells routines that need derived data to re-update,eg connected() needs
	// nres
	new_topology = true;

	if ( num_jump_in == 0 ) {
		// make a simple tree. this could also have been done by simple_tree()
		edge_list_.clear();
		if ( root_in != 1 ) {
			//this is for re-rooting the tree
			add_edge( 1, root_in, Edge::PEPTIDE );
			add_edge( root_in, nres_in, Edge::PEPTIDE );
			reorder(root_in);
		} else {
			add_edge( 1, nres_in, Edge::PEPTIDE );
		}
		return true; // success
	}

	//jjh Report jumps for reporting purposes
	if ( verbose ) {
		for ( int i = 1, iend = num_jump_in ; i <= iend ; ++i ) {
			TR.Debug << "Jump #" << i << " from " << jump_point_in( 1, i ) <<
				" to " << jump_point_in( 2, i ) << std::endl;
		}
	}

	// make a list of the unique jump_points in increasing order:
	// so we can construct the peptide edges
	typedef std::list< int > Int_list;
	Int_list vertex_list;
	// vertex_list.push_back( 1 );
	FArray1D_bool is_cut( nres_in, false ); // keep track of cuts
	for ( int i = 1; i <= num_jump_in; ++i ) {
		for ( int j = 1; j <= 2; ++j ) {
			int const pos ( jump_point_in(j,i) );
			//  debug_assert( pos >= 1 && pos <= nres_in );
			//   if ( jump_point_in( j, i ) == 1 || jump_point_in( j, i ) == nres_in ) {
			//    if ( verbose ) TR.Warning << "attempt to create jump with residue 1 or NRES: not supported.. returning invalid tree" << std::endl;
			//    return false;
			//   }
			vertex_list.push_back( pos );
		}
		runtime_assert( jump_point_in(1,i) < jump_point_in(2,i) );
		int const cut( cuts(i) );
		debug_assert( cut >= 1 && cut < nres_in );
		is_cut( cut ) = true;
		vertex_list.push_back( cut );
		vertex_list.push_back( cut+1 );
	}
	// vertex_list.push_back( nres_in );

	vertex_list.sort();
	vertex_list.unique(); // something like that...

	// start building the tree, add peptide edges
	edge_list_.clear();

	int const jump_stop( *( vertex_list.begin() ) );
	if ( jump_stop > 1 ) add_edge( 1, jump_stop, Edge::PEPTIDE );

	for ( auto it = vertex_list.begin(),
			it_end = vertex_list.end(); it != it_end; ++it ) {
		auto it_next (it);
		++it_next;
		if ( it_next == it_end ) break;

		int const start ( *it );
		int const stop ( *it_next );
		debug_assert( start >= 1 && start < stop && stop <= nres_in );
		if ( !is_cut(start) ) {
			add_edge( start, stop, Edge::PEPTIDE );
		} else {
			debug_assert( stop == start + 1 );
		}
	}

	// Add final edge.
	auto last_jump_it = vertex_list.end();
	int const jump_start( *(--last_jump_it) );
	if ( jump_start < nres_in ) add_edge( jump_start, nres_in, Edge::PEPTIDE );

	// now add the edges corresponding to jumps
	for ( int i=1; i<= num_jump_in; ++i ) {
		add_edge( jump_point_in(1,i), jump_point_in(2,i), i );
	}

	bool reorder_success = reorder(root_in, verbose);
	if ( !reorder_success ) return false;

	return check_fold_tree();
}

//////////////////////////////////////////////////////////////////////////////
/// @details
/// random_tree_from_jumps builds a tree from a list of jump_points and a
/// bias array telling us where we should prefer to cut:
///         P(cut at i ) ~ cut_bias(i).
/// It returns false if choosing fails (incompatible jumps? cycles?....).
/// as noted in update_edge_labels(), without corresponding cutpoint for each
/// jump added, the fold tree will be cyclic ( e.g., some of the edges are non-
/// separating with label==-2). This function keeps finding those edges and choose
/// cutpoints based on the biased probability.
/// note that during the process of building out tree, we cant use
/// any of the derived data, since conceptually the derived data comes
/// after the tree.
///
/// Note (from rhiju): updated routine to allow jumps at 1 or NRES, but for now
///  keep this behavior inactive by default --  other routines used by
///  protein jumpers (get_residue_edge?) appear to stumble if jumps start
///  at the root residue.
bool FoldTree::random_tree_from_jump_points(
	int const nres_in,
	int const num_jump_in,
	FArray2D_int const & jump_point_in, /* 2 x njump */
	FArray1D_float const & cut_bias,
	int const root_in /* = 1 */,
	bool const allow_jump_at_1_or_NRES /* = false */)
{
	std::vector< int > obligate_cut_points; //empty.
	return random_tree_from_jump_points( nres_in, num_jump_in, jump_point_in, obligate_cut_points, cut_bias, root_in, allow_jump_at_1_or_NRES );
}


bool
FoldTree::random_tree_from_jump_points(
	int const nres_in,
	int const num_jump_in,
	FArray2D_int const & jump_point_in,
	std::vector< int > const & obligate_cut_points,
	FArray1D_float const & cut_bias,
	int const root_in /* = 1 */,
	bool const allow_jump_at_1_or_NRES /* = false */)
{
	// tells routines that need derived data to re-update,eg connected() needs
	// nres
	new_topology = true;

	if ( num_jump_in == 0 ) {
		// make a simple tree. this could also have been done by simple_tree()
		edge_list_.clear();
		add_edge( 1, nres_in, Edge::PEPTIDE );
		return true; // success
	}

	// this array is passed in to the edge-cutting routines
	FArray1D_float cut_bias_sum( DRange(0,nres_in) );
	cut_bias_sum(0) = 0.0;
	for ( int i=1; i<= nres_in; ++i ) {
		cut_bias_sum(i) = cut_bias_sum( i-1 ) + cut_bias(i);
	}


	// make a list of the unique jump_points in increasing order:
	// so we can construct the peptide edges
	typedef std::list< int > Int_list;
	Int_list jump_list;
	// jump_list.push_back( 1 );
	for ( int i = 1; i <= num_jump_in; ++i ) {
		for ( int j = 1; j <= 2; ++j ) {
			if ( !allow_jump_at_1_or_NRES && (jump_point_in( j, i ) == 1 || jump_point_in( j, i ) == nres_in )  ) {
				TR.Warning << "attempt to create jump with residue 1 or NRES: not supported.. returning invalid tree" << std::endl;
				return false;
			}
			jump_list.push_back( jump_point_in(j,i) );

		}
	}
	// jump_list.push_back( nres_in );

	jump_list.sort();
	jump_list.unique(); // something like that...

	// start building the tree, add peptide edges
	edge_list_.clear();

	//Add beginning edge.
	int const jump_stop( *jump_list.begin() );
	if ( jump_stop > 1 ) add_edge( 1, jump_stop, Edge::PEPTIDE );

	// Add intervening segments.
	for ( auto it = jump_list.begin(),
			it_end = jump_list.end(); it != it_end; ++it ) {
		auto it_next (it);
		++it_next;
		if ( it_next == it_end ) break;

		int const start ( *it );
		int const stop ( *it_next );
		int const label ( -2 );//(start == 1 || stop == nres_in) ? Edge::PEPTIDE : -2 );
		debug_assert( start >= 1 && start < stop && stop <= nres_in );
		add_edge( start, stop, label );
	}

	//Add final edge.
	auto last_jump_it = jump_list.end();
	int const jump_start( *(--last_jump_it) );
	if ( jump_start < nres_in ) add_edge(  jump_start, nres_in, Edge::PEPTIDE );

	debug_assert( edge_list_[0].start() == 1 &&
		edge_list_[ edge_list_.size()-1 ].stop() == nres_in );

	// now add the edges corresponding to jumps
	for ( int i=1; i<= num_jump_in; ++i ) {
		add_edge( jump_point_in(1,i), jump_point_in(2,i), i );
	}

	//Add cuts that the user has given ... there could be none.
	int const num_user_cut_points = obligate_cut_points.size();
	for ( int i = 0; i < num_user_cut_points; i++ ) {
		int const cut_point = obligate_cut_points[i];
		bool const success = cut_edge( cut_point );

		if ( !success ) {
			return false;
			// this is a problem!
			//   TR.Fatal << "Problem with user-defined cutpoint: " << cut_point << " " << (*this) << std::endl;
			//   utility_exit();
		}

	}

	// keep cutting edges until we get a tree
	bool is_a_tree ( false );
	while ( ! is_a_tree ) {
		update_edge_labels();
		is_a_tree = true;
		for ( auto & it : edge_list_ ) {
			if ( it.label() == -2 ) {
				is_a_tree = false;
				break;
			}
		}

		if ( ! is_a_tree ) {
			if ( ! cut_random_edge( cut_bias_sum, nres_in ) ) {
				// failed to cut
				TR.Error << "failure in FoldTree::choose_random_cuts(...)." << std::endl;
				return false; // signal failure
			}
		}
	} // while ( ! is_a_tree )

	reorder(root_in);

	/// wrong behaviour:
	//debug_assert( check_fold_tree() ); // tree should be valid now
	// return true;
	//This assertion is out of place:
	// the interface details that function returns false if it fails...
	return check_fold_tree();
} // random_tree_from_jump_points(...)


///////////////////////////////////////////////////////////
bool FoldTree::cut_edge( int const cut_point ) {

	for ( auto it = edge_list_.begin(), it_end = edge_list_.end();
			it != it_end; ++it ) {
		if ( it ->label() < 0 &&
				( ( it->start() <= cut_point && it->stop()  >= cut_point+1 ) ||
				( it->stop()  <= cut_point && it->start() >= cut_point+1 ) ) ) {
			if ( it->label() == Edge::PEPTIDE ) {
				// you cant cut at this cutpoint
				break;
			} else {
				debug_assert( it->label() == -2 );
				TR.Debug << "cutting at " << cut_point << std::endl;

				int const start( std::min( it->start(), it->stop()) );
				int const stop ( std::max( it->start(), it->stop()) );
				delete_edge( it );
				if ( start  < cut_point ) add_edge( start, cut_point, Edge::PEPTIDE); // separating
				if ( cut_point+1 < stop ) add_edge( cut_point+1, stop, Edge::PEPTIDE); // separating
				return true; // successfully cut the tree
			}
		}
	} // loop over edge_list
	return false;
}

/////////////////////////////////////////////////////////////////////////////
/// @details  Fill a vector of cutpoints.
/// @note  Slow: just for convenience routines...
utility::vector1< int >
FoldTree::cutpoints() const
{
	check_topology();
	utility::vector1< int > cuts;
	for ( int i=1; i<= num_cutpoint_; ++i ) {
		cuts.push_back( cutpoint_[i] );
	}
	return cuts;
}

///////////////////////////////////////////////////////////////////////////////
/// @details
///
/// this routine assigns labels to the edges of a graph based
/// on whether or not those edges are separating -- ie whether
/// they can be cut without disconnecting the graph.
/// we know we have a tree when all the edges are separating
///
/// edge labels:
/// - +N means that the edge corresponds to jump #N (=> uncuttable)
/// - -2 means cuttable
/// - -1 means separating == PEPTIDE
/// - 0 means cut
///
/// we assume that the only possible change in edge labeling that we
/// need to make is from a -2 to a -1
/// ie, the -1's are still correct
/// also, there shouldn't be any 0's before this is called
/// 0's are for communicating between this function and
/// the logical function connected_graph(g)
void
FoldTree::update_edge_labels()
{
	debug_assert( Edge::PEPTIDE != -2 );
	for ( auto & it : edge_list_ ) {
		if ( it.label() == -2 ) { // labeled as not separating
			it.label() = 0;
			if ( ! connected() ) {
				it.label() = Edge::PEPTIDE; // now it's separating
			} else {
				it.label() = -2; // still not separating
			}
		} else if ( it.label() == 0 ) { // debug
			TR.Fatal << "zero edge label in graph:" << std::endl;
			utility_exit();
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details  Cuts a random edge chosen with per-rsd frequency given by cut_bias_sum.
/// private: returns true if success, false if failure.
/// operates on the edge_list_.
/// doesnt assume any derived data is up-to-date.
/// only an non-separating edge (label==-2) can be cut, otherwise fold tree will
/// not be valid.
bool
FoldTree::cut_random_edge(
	FArray1D_float const & cut_bias_sum,
	int const nres_in
)
{
	int tries(0);

	while ( tries < 100000 ) {
		++tries;

		int const cut_point ( pick_loopy_cutpoint( nres_in, cut_bias_sum ) );
		bool success = cut_edge( cut_point );
		TR.Debug << "Trying cut_point: " << cut_point << " " << success << std::endl;
		if ( success ) return true;
	} // keep trying

	return false; // too many tries
}

/////////////////////////////////////////////////////////////////////////////
/// @details  Assign new numbers to the jumps.
/// after we delete a jump, we may want to re-number the others.
/// note of course this will invalidate the jump_transform array of
/// any pose with this fold_tree, so be sure to call jumps_from_positions
/// or something
void
FoldTree::renumber_jumps()
{
	int counter(0);
	new_topology = true;
	for ( auto & it : edge_list_ ) {
		if ( it.is_jump() ) {
			++counter;
			TR.Debug << "renumber jumps:: from,to " << it.label() << ' ' <<
				counter << std::endl;
			it.label() = counter;
		}
	}
}

//@details Assigns new numbers to the jumps while maintaining the relative
// order of the labels.
void
FoldTree::renumber_jumps_ordered()
{
	int counter(0);
	new_topology = true;
	int highest_label(0);
	core::Size edgecounter = 0;
	utility::vector1< int > jump_labels;
	for ( auto & it : edge_list_ ) {
		edgecounter++;
		if ( it.is_jump() ) {
			jump_labels.push_back(it.label());
			++counter;
			if ( it.label() > highest_label ) highest_label = it.label();
		}
	}

	for ( auto & it : edge_list_ ) {
		if ( it.is_jump() ) {
			//determine the number of jumps with lower labels to account for missing
			int lower_count = 0;
			for ( Size i=1; i<=jump_labels.size(); i++ ) {
				if ( jump_labels[i] <= it.label() ) lower_count +=1;
			}
			//TR << "renumber jumps:: from,to " << it->label() << ' ' <<
			// it->label()-(it->label()-lower_count) << std::endl;
			it.label() = it.label()-(it.label()-lower_count);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details Is the tree connected?
/// returns true if fold_tree is connected
/// doesn't assume that the fold_tree is in valid folding order, or even a tree
bool
FoldTree::connected() const
{
	if ( new_topology ) update_nres();
	if ( int( linked_.size1() ) != nres_ ) linked_.dimension( nres_ );

	if ( edge_list_.size() <1 ) return true;

	linked_ = false;

	// grow out a single connected component from the start vertex:
	const_iterator const it_begin ( edge_list_.begin() );
	const_iterator const it_end ( edge_list_.end() );

	// mark start vertex as linked:
	linked_( it_begin->start() ) = true;

	bool new_member ( true );

	while ( new_member ) {     // keep adding new members
		new_member = false;
		for ( const_iterator it = it_begin; it != it_end; ++it ) {
			if ( it->label() == 0 ) continue;
			if ( linked_( it->start() ) && ! linked_( it->stop()) ) {
				linked_(it->stop()) = true;
				new_member = true;
			} else if ( linked_( it->stop() ) && ! linked_( it->start() ) ) {
				linked_( it->start() ) = true;
				new_member = true;
			}
		}
	}

	for ( const_iterator it = it_begin; it != it_end; ++it ) {
		if ( ! linked_( it->start() ) || ! linked_( it->stop() ) ) {
			// it's disconnected!
			return false;
		}
	}
	return true;
} // FoldTree::connected()


/////////////////////////////////////////////////////////////////////////////
/// @details  Create two new foldtrees f1 and f2 by splitting myself at jump jump_number
/// Uses the following routine to figure out which vertices should be in each tree.
///
/// @note The N-terminal vertex of jump "jump_number" goes to tree f1
void
FoldTree::partition_by_jump(
	int const jump_number,
	FoldTree & f1, //contains N-terminal vertex in jump, like partner1 in partition_by_jump below
	FoldTree & f2
) const
{
	check_topology(); // update derived data if necessary

	// partition the residues
	FArray1D_bool partner1( nres_, false );
	partition_by_jump( jump_number, partner1 );
	f1.clear();
	f2.clear();

	utility::vector1< int > seqpos1( nres_, 0 );
	utility::vector1< int > seqpos2( nres_, 0 );
	Size nres1(0), nres2(0);

	for ( Size i=1; i<= Size(nres_); ++i ) {
		if ( partner1(i) ) {
			++nres1;
			seqpos1[i] = nres1;
		} else {
			++nres2;
			seqpos2[i] = nres2;
		}
	}

	int jump1(0);
	int jump2(0);
	for ( auto edge : edge_list_ ) {
		if ( edge.is_jump() && edge.label() == jump_number ) continue; // the split jump itself

		if ( partner1( edge.start() ) ) {
			debug_assert( partner1( edge.stop() ) );

			edge.start() = seqpos1[ edge.start() ];
			edge.stop () = seqpos1[ edge.stop () ];
			debug_assert( edge.start() && edge.stop() );
			if ( edge.is_jump() ) {
				++jump1;
				edge.label() = jump1;
			}
			f1.add_edge( edge );

		} else {
			debug_assert( !partner1( edge.stop() ) );

			edge.start() = seqpos2[ edge.start() ];
			edge.stop () = seqpos2[ edge.stop () ];
			debug_assert( edge.start() && edge.stop() );
			if ( edge.is_jump() ) {
				++jump2;
				edge.label() = jump2;
			}
			f2.add_edge( edge );
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details
/// when a jump edge is removed, the fold tree is separated into two parts. This
/// fucntion is to find all residues connecting to the jump starting residue and flag
/// them in the partner1(n_res) array as true. The residues on the other side are
/// flagged as false. Useful to distinguish two docking partners when fold tree is
/// properly set up.
void
FoldTree::partition_by_jump(
	int const jump_number,
	FArray1D_bool & partner1
) const
{
	check_topology(); // update derived data if necessary

	debug_assert( jump_number <= num_jump_ );
	debug_assert( int(partner1.size1()) >= nres_ );

	// find n-terminal jump vertex
	int const pos1( jump_point_[jump_number ].first );

	// mark start vertex as linked:
	partner1 = false;
	//for ( int i=1; i<= nres_; ++i ) {
	// partner1(i) = false;
	//}
	partner1( pos1 ) = true;

	bool new_member ( true );

	// get pointers to the beginning and end of the edge_list_
	// const_iterator is a typedef in fold_tree.h:
	//
	// typedef std::vector< Edge > EdgeList;
	// typedef EdgeList::iterator iterator;
	// typedef EdgeList::const_iterator const_iterator;

	auto it_begin( edge_list_.begin() );
	auto it_end  ( edge_list_.end() );

	while ( new_member ) {     // keep adding new members
		new_member = false;
		for ( auto it = it_begin; it != it_end; ++it ) {
			if ( it->label() == jump_number ) continue; // skip jump

			int const start( std::min( it->start(), it->stop() ) );
			int const stop ( std::max( it->start(), it->stop() ) );
			if ( (partner1( start ) && !partner1( stop )) ||
					(partner1( stop ) && !partner1( start )) ) {
				new_member = true;
				if ( it->is_polymer() ) {
					// all the residues
					for ( int i=start; i<= stop; ++i ) {
						partner1( i ) = true;
					}
				} else {
					// just the vertices
					partner1( start ) = true;
					partner1( stop ) = true;
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
/// when a jump edge is removed, the fold tree is separated into two parts. This
/// fucntion is to find all residues connecting to the jump starting residue and flag
/// them in the partner1(n_res) array as true. The residues on the other side are
/// flagged as false. Useful to distinguish two docking partners when fold tree is
/// properly set up.
utility::vector1< bool >
FoldTree::partition_by_jump( Size const jump_nr ) const {
	ObjexxFCL::FArray1D<bool> partition_definition( nres_, false );
	partition_by_jump( jump_nr, partition_definition );

	//silly conversion. There may be a faster way to do this actually.
	utility::vector1< bool > partition_definition_vector1;
	for ( int n = 1; n <= nres_; n++ ) partition_definition_vector1.push_back( partition_definition(n) );

	return partition_definition_vector1;
}

/////////////////////////////////////////////////////////////////////////////
/// @details partition the fold tree into n parts based on specified jumps.
///  Just uses partition_by_jump over and over again.
utility::vector1< Size >
FoldTree::partition_coloring( utility::vector1< Size > const & jump_numbers ) const {

	utility::vector1< Size > partition_color( nres(), 0 );

	for ( Size n = 1; n <= jump_numbers.size(); n++ ) {
		Size const jump_nr = jump_numbers[ n ];
		Size const res1 = upstream_jump_residue( jump_nr );
		Size const res2 = downstream_jump_residue( jump_nr );
		runtime_assert( partition_color[ res1 ] == partition_color[ res2 ] ); // should not have been split yet.

		Size original_color = partition_color[ res1 ]; // or res2, they're the same color

		utility::vector1< bool > partition_definition = partition_by_jump( jump_nr );
		for ( Size k = 1; k <= nres(); k++ ) {
			if ( partition_definition[ k ] == partition_definition[ res2 ] &&
					partition_color[ k ] == original_color ) {
				partition_color[ k ] = n; // update coloring on downstream side of the partition.
			}
		}
		runtime_assert( partition_color[ res1 ] != partition_color[ res2 ] ); // better be split.
	}

	return partition_color;
}

/////////////////////////////////////////////////////////////////////////////
/// @details partition the fold tree in two parts if a cut would be introduced between seqpos and seqpos+1.
/// Function is an analog to partition_by_jump() -- its goal to find all residues
///  connecting to the jump starting residue and flag
/// them in the partner1(n_res) array as true. The residues on the other side are
/// flagged as false. Useful to distinguish two docking partners when fold tree is
/// properly set up.
void
FoldTree::partition_by_residue(
	int const seqpos,
	FArray1D_bool & partner1
) const
{
	check_topology(); // update derived data if necessary

	debug_assert( seqpos <= nres_ );
	debug_assert( int(partner1.size1()) >= nres_ );

	partner1 = false;

	// mark part of starting edge as linked.

	// First have to find edge. Standard get_residue_edge() has problem with root? Why?
	// Edge edge = get_residue_edge( seqpos );
	Edge edge;
	bool found_edge( false );
	partner1( seqpos ) = true;
	for ( const auto & it : *this ) {

		int const start( std::min( it.start(), it.stop() ) );
		int const stop ( std::max( it.start(), it.stop() ) );
		if ( it.is_polymer() && start <= seqpos && stop >= seqpos ) {
			edge = it;
			found_edge = true;

			// std::cout << "FOUND START EDGE: " << edge.start() << " " << edge.stop() << std::endl;
			int const seqpos_edge_start( std::min( edge.start(), edge.stop() ) );
			for ( int i = seqpos_edge_start; i< seqpos; ++i ) partner1( i ) = true;

		}
	}
	if ( !found_edge ) return;

	bool new_member ( true );

	// get pointers to the beginning and end of the edge_list_
	// const_iterator is a typedef in fold_tree.h:

	auto it_begin( edge_list_.begin() );
	auto it_end  ( edge_list_.end() );

	while ( new_member ) {     // keep adding new members
		new_member = false;
		for ( auto it = it_begin; it != it_end; ++it ) {

			int const start( std::min( it->start(), it->stop() ) );
			int const stop ( std::max( it->start(), it->stop() ) );

			if ( it->is_polymer() && start <= seqpos && stop >= seqpos ) continue; // skip input edge.

			if ( ( partner1( start ) && !partner1( stop ) ) ||
					( partner1( stop ) && !partner1( start ) ) ) {
				new_member = true;
				if ( it->is_polymer() ) {
					// all the residues
					for ( int i=start; i<= stop; ++i ) {
						partner1( i ) = true;
					}
				} else {
					// just the vertices
					partner1( start ) = true;
					partner1( stop ) = true;
				}
			}
		}
	}
}

int
FoldTree::jump_point(
	int const lower_higher, // = 1 or 2
	int const jump_number
) const
{
	PyAssert(((jump_number > 0) || (jump_number <= num_jump_)),
		"FoldTree::jump_point( int const lower_higher, int const jump_number ): Input variable jump_number is not a valid value.");
	check_topology();
	if ( lower_higher == 1 ) {
		return jump_point_[jump_number].first;
	} else if ( lower_higher == 2 ) {
		return jump_point_[jump_number].second;
	} else {
		std::cout << "FoldTree::jump_point() lower_higher needs to be 1 or 2" << std::endl;
		std::exit(9);
	}
}


////////////////////////////////////////////////////////////////////////////
/// @details
/// To keep the fold tree non-cyclic, for each jump added, there should be
/// a corresponding cutpoint. First call partition_by_jump() and the cutpoint
/// for this jump would be those two sequentially adjacent residues which are
/// not connected any more if the jump is disconnected.
///
/// @note  This cutpoint is not necessarily unique if the foldtree is sufficiently
/// complex. Chooses the first cutpoint with the desired property, starting at the N-terminus.
/// Will be unique eg if the jump is the intra-template jump used to support a single loop
/// region during loop modeling.
int
FoldTree::cutpoint_by_jump(
	int const jump_number
) const
{
	check_topology(); // update derived data if necessary

	FArray1D_bool partner1( nres_, false );
	partition_by_jump( jump_number, partner1 );
	int const i_min = std::min( upstream_jump_residue( jump_number ), downstream_jump_residue( jump_number ) );
	int const i_max = std::max( upstream_jump_residue( jump_number ), downstream_jump_residue( jump_number ) ) - 1;
	for ( int i = i_min; i <= i_max; i++ ) {
		if ( partner1(i) != partner1(i+1) ) {
			return i;
		}
	}

	TR.Fatal << " FoldTree::cutpoint_by_jump error: "
		<< "can not find the cutpoint! for jump_number: " << jump_number
		<< std::endl << (*this) << std::endl;
	utility_exit();
	return 0;

}


//////////////////////////////////////////////////////////////////////////
/// @details
/// - edge_count(i): delete the edge (i-1,i), how many residues are in the
///   component of the graph containing i-1?
/// - jump_edge_count(i): delete jump_number i. How many residues are in
///   the component of the graph containing the jump point on the
///   N-terminal side of the jump.
///
/// @note edge_count(cutpoint+1) doesn't really make sense
///         currently set to 0 but routines should avoid looking at this
///         value (see eg refold_reorder(...) )
/// @note Not checked out for chemical links
void
FoldTree::setup_edge_counts() const
{
	// redimension?
	if ( static_cast<int>(edge_count.size()) != nres_ ) {
		edge_count = utility::vector1<int>(nres_, 0);
	}
	if ( (int)jump_edge_count.size() != num_jump_ ) {
		jump_edge_count = utility::vector1<int>(num_jump_, -1);
	}
	if ( nres_ != static_cast<int>( linked_.size1() ) ) {
		linked_.dimension( nres_ );
	}

	std::fill( edge_count.begin(), edge_count.end(), 0 );
	std::fill( jump_edge_count.begin(), jump_edge_count.end(), -1 );

	min_edge_count = nres_;

	const_iterator const it_begin ( edge_list_.begin() );
	const_iterator const it_end   ( edge_list_.end()   );

	for ( const_iterator it = it_begin; it != it_end; ++it ) {
		linked_ = false;
		int const begin_res ( std::min( it->start(), it->stop()) );
		int link_count (0);
		linked_( begin_res ) = true;

		// find all the residues linked to begin_res, when we arent allowed
		// to traverse the edge *it
		bool new_member = true;
		while ( new_member ) {
			new_member = false;

			for ( const_iterator it2 = it_begin; it2 != it_end; ++it2 ) {
				if ( it2 == it ) continue;
				int const start ( std::min( it2->start(), it2->stop() ) );
				int const stop  ( std::max( it2->start(), it2->stop() ) );
				int new_link(0);
				if ( linked_(start) && !linked_(stop) ) {
					new_link = stop;
				} else if ( linked_(stop) && !linked_(start) ) {
					new_link = start;
				}
				if ( new_link > 0 ) {
					new_member = true;
					linked_(new_link) = true;
					// how many new residues does this edge add?
					link_count +=
						( it2->is_jump() ) ? 1 : stop - start;
				}
			}
		} // while ( new_member )

		if ( it->is_jump() ) {
			// jump
			int const jump_number ( it->label() );
			jump_edge_count[ jump_number ] = link_count + 1;
			min_edge_count = std::min( min_edge_count, std::max(
				jump_edge_count[ jump_number ],
				nres_ - jump_edge_count[ jump_number ] ) );
		} else {
			// peptide edge
			int const end_res  ( std::max(it->start(), it->stop()) );

			for ( int i= begin_res+1; i<= end_res; ++i ) {
				edge_count[i] = link_count + i - begin_res;
				min_edge_count = std::min( min_edge_count, std::max(
					edge_count[i], nres_ - edge_count[i] ) );
			}
		}
	}

	for ( int i=1; i<= num_jump_; ++i ) debug_assert( jump_edge_count[ i ] >= 1 );

} // FoldTree::setup_edge_counts(...)


Edge const &
FoldTree::jump_edge( int const jump_number ) const
{
	PyAssert(((jump_number > 0) || (jump_number <= num_jump_)),
		"FoldTree::jump_edge( int const jump_number ): Input variable jump_number is not a valid value.");
	check_order();
	return edge_list_[ jump_edge_[ jump_number ] ];
}


Edge &
FoldTree::jump_edge( int const jump_number )
{
	PyAssert(((jump_number > 0) || (jump_number <= num_jump_)),
		"FoldTree::jump_edge( int const jump_number ): Input variable jump_number is not a valid value.");
	check_order();
	return edge_list_[ jump_edge_[ jump_number ] ];
}


Edge const &
FoldTree::get_residue_edge( int const seqpos ) const
{
	if ( seqpos == root() ) utility_exit_with_message( "FoldTree:: residue_edge is undefined for root vertex" );

	for ( const auto & it : *this ) {
		if ( seqpos == it.stop() ||
				( it.is_peptide() &&
				( ( seqpos > it.start() && seqpos <= it.stop() ) ||
				( seqpos < it.start() && seqpos >= it.stop() ) ) ) ) {
			return it;
		}
	}
	utility_exit_with_message( "no edge found that contains seqpos!" );
	return *begin();
}

/// @details  Returns all edges that build a residue directly off of seqpos
utility::vector1< Edge >
FoldTree::get_outgoing_edges( int const seqpos ) const
{
	utility::vector1< Edge > outgoing;
	for ( const auto & it : *this ) {
		//SML 9/30/08 it->stop() != seqpos protects one-residue poses from having outgoing edges
		if ( (it.start() == seqpos && it.stop() != seqpos) ||
				( it.is_polymer() && ( ( seqpos >= it.start() && seqpos < it.stop() ) ||
				( seqpos <= it.start() && seqpos > it.stop() ) ) ) ) {
			outgoing.push_back( it );
		}
	}
	return outgoing;
}


/// @author  Labonte <JWLabonte@jhu.edu>
utility::vector1< Edge >
FoldTree::get_jump_edges( ) const
{
	check_order();
	utility::vector1< Edge > edges;
	for ( const auto & it : *this ) {
		if ( it.is_jump() ) {
			edges.push_back( it );
		}
	}
	return edges;
}

/// @author  Morgan Nance
utility::vector1< Edge >
FoldTree::get_chemical_edges( ) const
{
	check_order();
	utility::vector1< Edge > edges;
	for ( const auto & it : *this ) {
		if ( it.is_chemical_bond() ) {
			edges.push_back( it );
		}
	}
	return edges;
}

///////////////////////////////////////////////////////////////////////
/// @details  Returns the folding direction of a given polymer (peptide) residue. If the residue
/// is in a peptide edge this is the direction in which that edge is traveled if
/// we traverse the tree starting at the root.
/// Will die if residue is root or if residue is built by jump or chemical bond.
int
FoldTree::get_polymer_residue_direction( int const seqpos ) const
{
	for ( const auto & it : *this ) {
		if ( it.is_peptide() ) {
			// peptide edge
			if ( seqpos > it.start() && seqpos <= it.stop() ) {
				return 1; // forward
			} else if ( seqpos < it.start() && seqpos >= it.stop() ) {
				return -1; // backward
			}
		}
	}
	utility_exit_with_message( "no peptide edge found that contains (builds) seqpos!" );
	return 0;
}

///////////////////////////////////////////////////////////////////////
/// @details  Internal routine for updating data that is derived from the edge list (which is the only primary data).
void
FoldTree::update_cutpoints() const
{

	// first: is_cutpoint_
	is_cutpoint_.dimension( DRange(0,nres_) );
	is_cutpoint_ = true;

	// loop through the peptide edges, each implies a range of NON-cutpoints:
	for ( const auto & it : edge_list_ ) {
		if ( it.is_polymer() ) {
			for ( int j = std::min( it.start(), it.stop() ),
					j_end = std::max( it.start(), it.stop() ); j < j_end; ++j ) {
				is_cutpoint_(j) = false;
			}
		}
	}
	debug_assert( is_cutpoint_( 0 ) && is_cutpoint_( nres_ ) );

	// count the cutpoints. note that 0,total_residue dont count as cutpoints
	// for num_cutpoint_
	num_cutpoint_ = 0;

	for ( int i = 1; i < nres_; ++i ) {
		if ( is_cutpoint_(i) ) {
			++num_cutpoint_;
		}
	}

	// now cutpoint_ and cutpoint_map_ (which are inverses)
	cutpoint_ = utility::vector1<int>( num_cutpoint_, 0);
	cutpoint_map_ = utility::vector1<int>( nres_, 0);

	for ( int i = 1, cut=0; i < nres_; ++i ) {
		if ( is_cutpoint_(i) ) {
			cutpoint_[ ++cut ] = i;
			cutpoint_map_[ i ] = cut;
		}
		debug_assert( i<nres_-1 || cut == num_cutpoint_ );
	}

}


///////////////////////////////////////////////////////////////////////////////
/// @details  Internal routine for updating data that is derived from the edge list (which is the only primary data).
void
FoldTree::update_nres() const
{
	int tmp_nres (0);
	for ( const auto & it : edge_list_ ) {
		tmp_nres = std::max( tmp_nres, std::max( it.start(), it.stop()) );
	}
	if ( tmp_nres != nres_ ) {
		//std::cout << "FoldTree::update_nres: nres has changed from: " << nres_
		//     << " to: " << tmp_nres << std::endl;
		nres_ = tmp_nres;
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details  Internal routine for updating data that is derived from the edge list (which is the only primary data).
void
FoldTree::update_num_jump() const
{
	int tmp_num_jump (0);
	int biggest_label (0); // for debugging
	for ( const auto & it : edge_list_ ) {
		if ( it.is_jump() ) {
			++tmp_num_jump;
			biggest_label = std::max( biggest_label, it.label() );
		}
	}

	if ( biggest_label != tmp_num_jump ) {
		TR.Fatal << "problem with the fold_tree: biggest_label != num_jump " <<
			biggest_label << ' ' << tmp_num_jump << std::endl;
		TR.Fatal << *this;
		utility_exit();
	}

	if ( tmp_num_jump != num_jump_ ) {
		//std::cout << "FoldTree::update_num_jump: num_jump has changed from: "
		//     << num_jump_ << " to: " << tmp_num_jump << std::endl;
		num_jump_ = tmp_num_jump;
	}
} // FoldTree::update_num_jump


///////////////////////////////////////////////////////////////////////////////
/// @details  Internal routine for updating data that is derived from the edge list (which is the only primary data).
/// fills is_jump_point, jump_point
void
FoldTree::update_jump_points() const
{
	// re-dimension?
	if ( (int)is_jump_point_.size() != nres_ ) {
		is_jump_point_ = utility::vector1<bool>(nres_, false);
	} else {
		std::fill(is_jump_point_.begin(), is_jump_point_.end(), false);
	}
	if ( (int)jump_point_.size() != num_jump_ ) {
		jump_point_ = utility::vector1<std::pair<int,int> >(num_jump_);
	}

	for ( const auto & it : edge_list_ ) {
		if ( it.is_jump() ) {
			int const jump_number ( it.label() );
			debug_assert( jump_number <= num_jump_ );

			is_jump_point_[ it.start() ] = true;
			is_jump_point_[ it.stop () ] = true;

			jump_point_[jump_number].first = std::min( it.start(), it.stop());
			jump_point_[jump_number].second = std::max( it.start(), it.stop());
		} else if ( it.is_chemical_bond() ) {
			// this is a little hacky -- calling chemical bond connections jump_points
			// need this for internal use, eg operations like insert_polymer_residue and delete_extra_vertices
			is_jump_point_[ it.start() ] = true;
			is_jump_point_[ it.stop () ] = true;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details fills jump_edge
/// routines that use jump_edge should call check_order to ensure that
/// its up to date
///
/// @note that this is sensitive to the order of the tree
/// so it's updated in check_order()
///
/// Internal routine for updating data that is derived from the edge list (which is the only primary data).
void
FoldTree::update_jump_edge() const
{
	if ( (int)jump_edge_.size() != num_jump_ ) {
		jump_edge_ = utility::vector1<int>(num_jump_, 0);
	}
	int jump_index(0);
	for ( auto it = edge_list_.begin(), it_end = edge_list_.end();
			it != it_end; ++it, ++jump_index ) {
		if ( it->is_jump() ) {
			int const jump_number ( it->label() );
			debug_assert( jump_number <= num_jump_ );
			debug_assert( edge_list_[ jump_index ] == *it );

			jump_edge_[ jump_number ] = jump_index;
		}
	}
}

void
FoldTree::show(std::ostream & out) const
{
	out << "   Edge   \t   Jump     Jump #\n";
	for ( const auto & it : *this ) {
		// how do I reference vector member?
		if ( it.is_jump() ) {
			out << "          \t" << I(4,4,it.start()) << "--" << I(4,4,it.stop()) << "  " << I(3,3,it.label()) << '\n';
		} else {
			out << I(4,4,it.start()) << "--" << I(4,4,it.stop()) << '\n';
		}
	}
	out << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
/// @details  Foldtree output to stream
std::ostream &
operator <<( std::ostream & os, FoldTree const & t )
{
	os << "FOLD_TREE ";
	for ( const auto & it : t ) {
		os << it;
	}
	return os;
}


/////////////////////////////////////////////////////////////////////////////
/// @details  Foldtree input from stream
std::istream &
operator >>( std::istream & is, FoldTree & t )
{
	t.new_topology = true;
	t.edge_list_.clear();

	std::string tag;
	is >> tag;
	if ( !is.fail() && tag == "FOLD_TREE" ) {
		while ( !is.fail() ) {
			Edge e;
			is >> e;
			if ( is.fail() ) break;
			t.edge_list_.push_back( e );
		}
		is.clear();
	}

	if ( t.edge_list_.size() == 0 ) {
		is.setstate( std::ios_base::failbit );
		TR.Error << "no fold_tree info in this stream." << std::endl;
	} else {
		if ( ! t.check_fold_tree() ) {
			TR.Error << "bad fold_tree, reordering." << std::endl;
			t.reorder( t.edge_list_.begin()->start() );
			if ( ! t.check_fold_tree() ) {
				TR.Error << "bad fold_tree still bad" << std::endl;
				TR.Error << t;
			}
		}
	}
	return is;
}

/////////////////////////////////////////////////////////////////////////////
/// @details  Check to see if a foldtree is in valid folding order.
/// To be valid, a fold tree needs to be connected, but not cyclic. So the tree
/// is traversed from the root residue and if any residue has not been visited or
/// has been visited multiple times, the fold tree is bad.
bool
FoldTree::check_fold_tree() const
{
	if ( edge_list_.size() <= 0 ) return false;
	check_topology();
	//if ( new_topology ) update_nres(); // largest vertex
	if ( int( seen_.size1() ) != nres_ ) seen_.dimension( nres_ );

	seen_ = false;
	seen_( root() ) = true; // Root is first to build
	for ( Edge const & edge: edge_list_ ) {
		int const start( edge.start() );
		int const stop ( edge.stop() );
		if ( ! seen_( start ) ) {
			TR.Error << "Bad fold tree at edge " << edge << std::endl;
			TR.Error << "Start residue " << start << " not built yet. " << std::endl;
			TR.Error << *this << std::endl;
			return false;
		}
		if ( start == stop ) continue;
		if ( seen_( stop ) ) {
			TR.Error << "Bad fold tree at edge " << edge << std::endl;
			TR.Error << "Stop residue " << stop << " has already been built. " << std::endl;
			TR.Error << *this << std::endl;
			return false;
		}
		if ( ! edge.is_polymer() ) {
			seen_( stop ) = true;
		} else {
			int const dir( start < stop ? 1 : -1 );
			for ( int i=start + dir; i!= stop + dir; i+= dir ) {
				if ( seen_( i ) ) {
					TR.Error << "Bad fold tree at edge " << edge << std::endl;
					TR.Error << "Residue " << i << " is in a polymeric edge, but has already been built from another edge." << std::endl;
					TR.Error << *this << std::endl;
					return false;
				}
				// for debugging purposes, do not uncomment unless you want it to
				// print out very often!
				//std::cout << "i=" << i << std::endl;
				seen_( i ) = true;
			}
		}
	}
	for ( int i=1; i<= nres_; ++i ) {
		if ( !seen_(i) ) {
			TR.Error << "Bad fold tree!" << std::endl;
			TR.Error << "Residue " << i << " is not built by any edge." << std::endl;
			TR.Error << *this << std::endl;
			return false;
		}
	}
	return true;
} // check_fold_tree()

bool
FoldTree::check_edges_for_atom_info() const
{
	// if ( edge_list_.size() <= 0 ) return false;
	// if ( new_topology ) update_nres(); // largest vertex
	auto it ( edge_list_.begin() );
	for ( auto it_end = edge_list_.end(); it != it_end; ++it ) {
		if ( it->label()== -2 && ! it->has_atom_info() ) {
			TR<< "bad chemical edge from"<< it->start() << " to "<< it->stop();
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
/// @details  Set connection atoms for a jump. This is not used by the foldtree, only to communicate to the
/// AtomTree during construction of an atomtree from a foldtree.
void
FoldTree::set_jump_atoms(
	int const jump_number,
	std::string const&  upstream_atom,
	std::string const&  downstream_atom,
	bool bKeepStubInResidue /* false */
)
{
	Edge & edge( jump_edge( jump_number ) );
	edge.upstream_atom() = upstream_atom;
	edge.downstream_atom() = downstream_atom;
	edge.keep_stub_in_residue() = bKeepStubInResidue;
	// either set both or none
	debug_assert( ( upstream_atom.size() && downstream_atom.size() )
		|| ( !upstream_atom.size() && !downstream_atom.size() ) );
}


//version of above but makes it permutation safe!
void
FoldTree::set_jump_atoms(
	int const jump_number,
	core::Size res1,
	std::string const&  atom1,
	core::Size res2,
	std::string const&  atom2,
	bool bKeepStubInResidue /* false */
)
{
	runtime_assert( res1 != res2 );
	Edge & edge( jump_edge( jump_number ) );
	if ( Size(edge.start()) == res1 ) {
		edge.upstream_atom() = atom1;
	} else {
		runtime_assert( Size(edge.stop()) == res1 );
		edge.downstream_atom() = atom1;
	}

	if ( Size(edge.start()) == res2 ) {
		edge.upstream_atom() = atom2;
	} else {
		runtime_assert( Size(edge.stop()) == res2 );
		edge.downstream_atom() = atom2;
	}

	edge.keep_stub_in_residue() = bKeepStubInResidue;
	// either set both or none
	debug_assert( ( atom1.size() && atom2.size() )
		|| ( !atom1.size() && !atom2.size() ) );
}


/////////////////////////////////////////////////////////////////////////////
/// @details  Get the upstream connection resid (connection atom # at the "start" vertex)
/// If it hasn't been set return 0.
/// Also see set_jump_atoms, which sets this data.
std::string
FoldTree::upstream_atom( int const jump_number ) const
{
	Edge const & edge( jump_edge( jump_number ) );
	if ( edge.has_atom_info() ) {
		return edge.upstream_atom();
	} else {
		return "";
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details  Get the downstream connection atomno (connection atom # at the "stop" vertex)
/// If it hasn't been set return 0.
/// Also see set_jump_atoms, which sets this data.
std::string
FoldTree::downstream_atom( int const jump_number ) const
{
	Edge const & edge( jump_edge( jump_number ) );
	if ( edge.has_atom_info() ) {
		return edge.downstream_atom();
	} else {
		return "";
	}
}


/////////////////////////////////////////////////////////////////////////////
// this assumes no fragment insertions across chainbreaks which is
// guaranteed by the settings of the insert_size map
//
// borrowed some code from refold_reorder
// returns the size of the largest single fixed region after
// a fragment is inserted from begin_res to begin_res+size-1
Size
FoldTree::count_fixed_residues(
	Size const begin_res,
	Size const size,
	Size & min_edge_count_out
) const
{
	check_topology();

	//moving the call to setup_edge_counts() from check_topology() to here
	//after thorough testing, this seems to be the only function that requires up-to-date edge
	//count data, so moving it here should cut down on unecessary edge count calculations
	setup_edge_counts();
	// pass out the value for min_edge_count
	// this is a measure of the magnitude of the largest single-residue or
	// single-jump
	// move we could possibly make. ie, its the number of fixed residues for the
	// move with the smallest number of fixed residues
	min_edge_count_out = min_edge_count;
	debug_assert( size > 0 );
	Size const end_res ( begin_res + size - 1);
	debug_assert( begin_res >= 1 && end_res <= static_cast<Size> ( nres_ ) );

	int best = 0;
	if ( ! is_cutpoint_( begin_res-1 ) ) {
		int const n_fixed ( edge_count[ begin_res ] );
		if ( n_fixed > best ) {
			best = n_fixed;
		}
	}

	if ( ! is_cutpoint_( end_res ) ) {
		int const c_fixed ( nres_ - edge_count[ end_res + 1] );
		if ( c_fixed > best ) {
			best = c_fixed;
		}
	}

	// how to test this stuff?
	for ( int i = 1; i<= num_jump_; ++i ) {
		for ( int j = 1; j <= 2; ++j ) {
			Size const pos = j == 1 ? jump_point_[i].first : jump_point_[i].second;
			if ( begin_res <= pos && pos <= end_res ) {
				int const fixed
					( j==1 ? nres_ - jump_edge_count[ i ] : jump_edge_count[ i ] );
				if ( fixed > best ) {
					best = fixed;
				}
			}
		}
	}
	return best;
}

void FoldTree::reassign_atoms_for_intra_residue_stubs() {
	debug_assert( check_fold_tree() ); // necessary?

	for ( Size jump_nr = 1; jump_nr <= num_jump(); ++jump_nr ) {
		if ( !jump_edge( jump_nr ).keep_stub_in_residue() ) continue; // do nothing

		std::string anchor = "";

		if ( jump_edge( jump_nr ).start() == root() ) {
			anchor = "N";
		} else {

			Edge anchor_edge = get_residue_edge( jump_edge( jump_nr ).start() );
			//work out upstream Jump Atom from upstream folding direction

			if ( !anchor_edge.is_jump() ) {
				bool bN2C = anchor_edge.start() < anchor_edge.stop();
				if ( bN2C ) {
					anchor = "C";
				} else {
					anchor = "N";
				}
			} else {
				//if it is a jump it will be an N now or later when we get to it.
				anchor = "N";
			}
		}

		// choosing the root to be N and setting keep_Stub_in_resiude makes N-CA-C the stub
		// C-->CA-->N
		std::string root = "N";
		TR.Debug << "set anchor and root atom for jump " << jump_nr << " to " << anchor << " and " << root << std::endl;

		std::string const upstream_atom_name = ObjexxFCL::strip_whitespace( jump_edge( jump_nr ).upstream_atom() );
		if ( upstream_atom_name != "" && upstream_atom_name != "N" && upstream_atom_name != "C" && upstream_atom_name != "CA"  ) {
			TR.Debug << "UPSTREAM_ATOM_NAME" <<  upstream_atom_name << std::endl;
			anchor = upstream_atom_name;
		}
		std::string const downstream_atom_name = ObjexxFCL::strip_whitespace( jump_edge( jump_nr ).downstream_atom() );
		if ( downstream_atom_name != "" && downstream_atom_name != "N" && downstream_atom_name != "C" && downstream_atom_name != "CA"  ) {
			TR.Debug << "DOWNSTREAM_ATOM_NAME" <<  downstream_atom_name << std::endl;
			root = downstream_atom_name;
		}

		set_jump_atoms( jump_nr, anchor, root, true /* keep_stub_in_residue */ );
	} // for loop
}

void FoldTree::put_jump_stubs_intra_residue() {
	//  TR.Trace << (*this) << std::endl;
	debug_assert( check_fold_tree() ); // necessary?
	for ( Size jump_nr = 1; jump_nr <= num_jump(); ++jump_nr ) {
		if ( ! jump_edge( jump_nr ).has_atom_info() ) {
			jump_edge( jump_nr ).keep_stub_in_residue() = true;
		} // if there is atom_info already, do not reassign
	} // for loop
	reassign_atoms_for_intra_residue_stubs();
}

} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::FoldTree::save( Archive & arc ) const {
	arc( CEREAL_NVP( edge_list_ ) ); // EdgeList
	arc( CEREAL_NVP( new_topology ) ); // _Bool
	arc( CEREAL_NVP( new_order ) ); // _Bool
	arc( CEREAL_NVP( nres_ ) ); // int
	arc( CEREAL_NVP( num_jump_ ) ); // int
	arc( CEREAL_NVP( num_cutpoint_ ) ); // int
	arc( CEREAL_NVP( jump_point_ ) ); // utility::vector1<std::pair<int, int> >
	arc( CEREAL_NVP( is_jump_point_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( cutpoint_ ) ); // utility::vector1<int>
	arc( CEREAL_NVP( cutpoint_map_ ) ); // utility::vector1<int>
	arc( CEREAL_NVP( is_cutpoint_ ) ); // ObjexxFCL::FArray1D_bool
	arc( CEREAL_NVP( jump_edge_ ) ); // utility::vector1<int>
	arc( CEREAL_NVP( edge_count ) ); // utility::vector1<int>
	arc( CEREAL_NVP( min_edge_count ) ); // int
	arc( CEREAL_NVP( jump_edge_count ) ); // utility::vector1<int>
	arc( CEREAL_NVP( linked_ ) );
	arc( CEREAL_NVP( seen_ ) );
	// turns out the "hasher" is never written to, it just hashes
	// strings as needed; it probably doesn't need to be a member
	// variable!
	// arc( CEREAL_NVP( hasher ) ); // boost::hash<std::string>
	// EXEMPT hasher
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::FoldTree::load( Archive & arc ) {
	arc( edge_list_ ); // EdgeList
	arc( new_topology ); // _Bool
	arc( new_order ); // _Bool
	arc( nres_ ); // int
	arc( num_jump_ ); // int
	arc( num_cutpoint_ ); // int
	arc( jump_point_ ); // utility::vector1<std::pair<int, int> >
	arc( is_jump_point_ ); // utility::vector1<_Bool>
	arc( cutpoint_ ); // utility::vector1<int>
	arc( cutpoint_map_ ); // utility::vector1<int>
	arc( is_cutpoint_ ); // ObjexxFCL::FArray1D_bool
	arc( jump_edge_ ); // utility::vector1<int>
	arc( edge_count ); // utility::vector1<int>
	arc( min_edge_count ); // int
	arc( jump_edge_count ); // utility::vector1<int>
	arc( linked_ );
	arc( seen_ );
	// hasher is not serialized/deserialized; see above note
	// arc( hasher ); // boost::hash<std::string>
	// EXEMPT hasher
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::FoldTree );
CEREAL_REGISTER_TYPE( core::kinematics::FoldTree )

CEREAL_REGISTER_DYNAMIC_INIT( core_kinematics_FoldTree )
#endif // SERIALIZATION
