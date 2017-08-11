// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/FoldTree.hh
/// @brief  Fold tree class
/// @author Phil Bradley

#ifndef INCLUDED_core_kinematics_FoldTree_HH
#define INCLUDED_core_kinematics_FoldTree_HH

// Unit headers
#include <core/kinematics/FoldTree.fwd.hh>

// Package Headers
#include <core/kinematics/Edge.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/SequenceMapping.fwd.hh>
#ifdef PYROSETTA
#include <core/id/SequenceMapping.hh>
#endif

// utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// External Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>

#include <boost/functional/hash/hash.hpp>         // for hash

// C++ Headers
#include <vector>
#include <iostream>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {

/// @brief The FoldTree is a residue-based tree-like representation of a molecule
///
/// @note all the derived data is "mutable", so that we can
/// update them as needed on the fly inside "const" member functions
/// this emphasizes that the core data is just the edge_list_
///
/// @note see @ref foldtree_overview "FoldTree overview and concepts" for details
///
/// Common Methods:
///     FoldTree.check_fold_tree
///     FoldTree.clear
///     FoldTree.new_jump
///     FoldTree.nres
///     FoldTree.num_jump
///     FoldTree.simple_tree
///     Foldtree.size
class FoldTree : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~FoldTree() override;
	// types
	typedef std::vector< Edge > EdgeList;
	typedef EdgeList::iterator iterator;
	typedef EdgeList::const_iterator const_iterator;

	/// @brief constructor
	FoldTree():
		ReferenceCount(),
		new_topology (true), // new_topology == true ==> new_order == true
		new_order(true),
		nres_(0),
		num_jump_(0),
		num_cutpoint_(0),
		min_edge_count(0)
	{}

	/// @brief Constructs a simple FoldTree
	///
	/// ft = FoldTree( nres_in )
	///
	/// int  nres_in   /the number of residues in this simple tree
	FoldTree( int const nres_in ):
		ReferenceCount(),
		new_topology (true), // new_topology == true ==> new_order == true
		new_order(true),
		nres_(0),
		num_jump_(0),
		num_cutpoint_(0),
		min_edge_count(0)
	{
		simple_tree( nres_in );
	}

	/// @brief operator=
	/// @note this version doesn't copy any of the derived data!
	/// will this be too slow? it will trigger re-calculating of everything
	/// every time we reject a move....
	FoldTree & operator =( FoldTree const & src )
	{
		edge_list_ = src.edge_list_;
		new_topology = true;
		return *this;
	}

	FoldTree clone()
	{
		return *this;
	}

	// non-modifying access /////////////////////////////////////////////////////

	/// @brief Return the starting residue of the first kinematic Edge to which res belongs.
	Size boundary_left( Size res ) const
	{
		return get_residue_edge( int( res ) ).start();
	}

	/// @brief Return the ending residue of the first kinematic Edge to which res belongs.
	Size boundary_right( Size res ) const
	{
		return get_residue_edge( int( res ) ).stop();
	}


	/// @brief Returns the number of edges in the FoldTree
	///
	/// example(s):
	///     ft.size()
	/// See also:
	///     FoldTree
	///     FoldTree.add_edge
	///     FoldTree.check_fold_tree
	///     FoldTree.delete_edge
	///     FoldTree.jump_edge
	///     FoldTree.nres
	///     FoldTree.num_jump
	inline
	int
	size() const
	{
		return int( edge_list_.size() );
	}

	/// @brief begin iterator of the edge_list
	inline
	const_iterator
	begin() const
	{
		return edge_list_.begin();
	}

	/// @brief end iterator of the edge_list
	inline
	const_iterator
	end() const
	{
		return edge_list_.end();
	}

	/// @brief Displays the FoldTree information
	///
	/// example(s):
	///    ft.show()
	/// See Also:
	///    Pose
	void
	show( std::ostream & out ) const;

	// modifiers //////////////////////////////////////////////////////////
	// add,delete edges:

	/// @brief Adds an edge from  <start>  to  <stop>
	void
	add_edge(
		core::Size const start,
		core::Size const stop,
		int const label
	);

	/// @brief Especially useful version of add_edge for chemical edge construction
	void
	add_edge(
		Size const start,
		Size const stop,
		std::string const & start_atom,
		std::string const & stop_atom
	);

	/// @brief Add the edge  <new_edge>
	///
	/// example(s):
	///     ft.add_edge(edge1)
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.delete_edge
	///     FoldTree.jump_edge
	///     FoldTree.nres
	///     FoldTree.num_jump
	void add_edge( Edge const & new_edge );

	/// @brief Prepend the edge  <new_edge>. Useful alternative to add_edge for setting root.
	void prepend_edge( Edge const & new_edge );

	/// @brief  Find and replace an Edge in the FoldTree.
	void replace_edge( Edge const & old_edge, Edge const & replacement_edge  );

	/// @brief Deletes the edge  <edge>  in the FoldTree by iterator
	void delete_edge( iterator edge );

	/// @brief Delete the edge  <edge>  in the fold tree by example edge
	///
	/// example(s):
	///     ft.delete_edge(edge1)
	/// See also:
	///     FoldTree
	///     FoldTree.add_edge
	///     FoldTree.check_fold_tree
	///     FoldTree.jump_edge
	///     FoldTree.nres
	///     FoldTree.num_jump
	void delete_edge( Edge const & edge );

	/// @brief Find an edge in fold tree and delete it
	void delete_unordered_edge( Size const start, Size const stop,int const label);

	/// @brief Changes the label of an edge in fold tree
	void update_edge_label( Size const start, Size const stop, int const old_label, int const new_label );

	/// @brief Returns the edge label of the edge from  <start>  to  <stop>
	int edge_label( Size const start, Size const stop );

	/// @brief Splits an edge into two at a specified position.
	void split_existing_edge_at_residue( Size const resNo );

	/// @brief Deletes all edge in the FoldTree
	///
	/// example(s):
	///     ft.clear()
	/// See also:
	///     FoldTree
	///     FoldTree.add_edge
	///     FoldTree.check_fold_tree
	///     FoldTree.delete_edge
	///     FoldTree.new_jump
	///     FoldTree.nres
	///     FoldTree.num_jump
	///     FoldTree.simple_tree
	///     FoldTree.size
	void clear() {edge_list_.clear(); new_topology=true;}

	/// @brief Renumbers the jump edges in the FoldTree
	/// How?
	void renumber_jumps();

	///@brief Renumbers the jump edges in the FoldTree while
	// maintaining relative labels.
	void renumber_jumps_ordered();

	/// @brief Deletes edges with start==stop
	/// allowable 1->1 edge for single residue FoldTree
	void delete_self_edges();

	/// @brief Delete vertices that are no longer necessary any more
	/// How is this determined?
	void delete_extra_vertices();

	/// @brief Deletes a continuous segment from  <seq_begin>  to  <seq_end>
	///
	/// example(s):
	///     ft.delete_segment(13,37)
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.delete_edge
	///     FoldTree.new_jump
	///     FoldTree.nres
	///     FoldTree.simple_tree
	void delete_segment( Size const seg_begin, Size const seg_end );

	/// @brief Deletes the residue  <seqpos>  from the FoldTree.
	/// Will rearrange topology if necessary.
	///
	/// example(s):
	///     ft.delete_seqpos(3)
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.clear
	///     FoldTree.new_jump
	///     FoldTree.nres
	///     FoldTree.num_jump
	///     FoldTree.simple_tree
	void
	delete_seqpos( Size const seqpos );

	/// @brief Inserts a polymer residue at position  <seqpos>
	/// How?
	void
	insert_polymer_residue(
		Size const seqpos,
		bool const join_lower,
		bool const join_upper
	);

	/// @brief Inserts a bonded residue at position  <seqpos>
	void
	insert_residue_by_chemical_bond(
		Size const seqpos,
		Size const anchor_residue,
		std::string const& anchor_atom,
		std::string const& root_atom
	);

	/// @brief Inserts a residue attached only by a jump. precondition is that seqpos-1 is a cutpoint
	/// @note that anchor_pos is wrt the current numbering system (ie before insertion)
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	void
	insert_residue_by_jump(
		Size const seqpos,
		Size anchor_pos,
		std::string const& anchor_atom = "",
		std::string const& root_atomno = ""
	);

	/// @brief Inserts a fold_tree as a subtree
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	void
	insert_fold_tree_by_jump(
		FoldTree const & subtree,
		Size const insert_seqpos,               // rsd 1 in subtree goes here
		Size const insert_jumppos,              // jump 1 in subtree goes here
		Size const anchor_pos,                  // in the old numbering system
		Size anchor_jump_number = 0,            // in the new jump numbering system, default=0
		std::string const & anchor_atom = "",  // could be ""
		std::string const & root_atom   = ""   // ditto
	);

	/// @brief Renumber all vertices according to an input sequence mapping
	void
	apply_sequence_mapping( id::SequenceMapping const & old2new );

	/// @brief Adds a new jump edge from  <pos1>  to  <pos2>  with cutpoint  <cutpoint>
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	core::Size new_jump( Size const jump_pos1, Size const jump_pos2, Size const cutpoint );

	void
	new_chemical_bond(
		Size const anchor_pos,
		Size const root_pos,
		std::string const & anchor_atom,
		std::string const & root_atom,
		Size const new_cutpoint
	);

	// ways to build an entire tree:

	/// @brief Produces a 1-edge FoldTree that is  <nres_in>  long
	/// No jumps or extraneous features
	///
	/// example(s):
	///     ft.simple_tree()
	/// See also:
	///     FoldTree
	///     FoldTree.add_edge
	///     FoldTree.check_fold_tree
	///     FoldTree.clear
	///     FoldTree.delete_edge
	///     FoldTree.is_simple_tree
	///     FoldTree.new_jump
	///     FoldTree.nres
	///     FoldTree.num_jump
	///     FoldTree.size
	void simple_tree( Size const nres_in ); // makes a 1-->total_residue tree

	/// @brief Returns true if the FoldTree has 1-edge (non-jump)
	///
	/// example(s):
	///     ft.is_simple_tree()
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.num_jump
	///     FoldTree.simple_tree
	bool is_simple_tree() const;

	/// @brief Builds a FoldTree from a list of  <jump_points>  and random cut
	/// points based on some biased probability
	/// Returns bool about success
	bool random_tree_from_jump_points(
		Size const nres_in, Size const num_jump_in,
		ObjexxFCL::FArray2D< Size > const & jump_point,
		ObjexxFCL::FArray1D_float const & cut_bias,
		Size const root_in=1,
		bool const allow_jump_at_1_or_NRES = false
	);

	/// @brief Builds a FoldTree from a list of  <jump_points>  and random cut
	/// points based on some biased probability and any user-defined obligate cutpoints
	/// Returns bool about success
	bool random_tree_from_jump_points(
		Size const nres_in, Size const num_jump_in,
		ObjexxFCL::FArray2D< Size > const & jump_point,
		std::vector< Size > const & obligate_cut_points,
		ObjexxFCL::FArray1D_float const & cut_bias,
		Size const root_in = 1,
		bool const allow_jump_at_1_or_NRES = false );

	/// @brief Constructs a FoldTree from listed  <jump point>  and  <cuts>
	/// Returns bool about success
	bool tree_from_jumps_and_cuts(
		Size const nres_in, Size const num_jump_in,
		ObjexxFCL::FArray2D< Size > const & jump_point,
		ObjexxFCL::FArray1D< Size > const & cuts,
		Size const root_in = 1,
		bool const verbose = false
	);

	/// @brief Appends a new residue to the tree, either by a jump or as a continuation of a peptide segment
	void
	append_residue(
		bool const attach_by_jump = false,
		core::Size const jump_anchor_residue = 0,
		std::string const& jump_upstream_atom = "",
		std::string const&  jump_downstream_atom = ""
	);

	/// @brief Appends a new residue to the tree using a chemical (APL-style) connection
	void
	append_residue_by_chemical_bond(
		core::Size const anchor_residue,
		std::string const& anchor_atom,
		std::string const& root_atom
	);


	/// @brief Reorders the FoldTree to start at residue  <start_residue>
	bool reorder( Size const start_residue, bool const verbose_if_fail = true );


	// check status ////////////////////////////////////////////////////////
	/// @brief Returns true if this is a valid FoldTree
	///
	/// example(s):
	///     ft.check_fold_tree()
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.clear
	///     FoldTree.is_simple_tree
	///     FoldTree.new_jump
	///     FoldTree.nres
	///     FoldTree.num_jump
	///     FoldTree.simple_tree
	bool check_fold_tree() const; // foldable?

	/// @brief Returns true if the FoldTree is connected
	bool connected() const; // connected

	/// @brief chemical edges should have atom info
	bool check_edges_for_atom_info() const;

	/// @brief the starting residue for this jump
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	Size upstream_jump_residue( Size const jump_number ) const;

	/// @brief the stopping residue for this jump
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	Size downstream_jump_residue( Size const jump_number ) const;

	/// @brief partition into two foldtrees by cutting at jump= jump_number
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	void
	partition_by_jump(
		Size const jump_number,
		FoldTree & f1, //contains N-terminal vertex in jump, like partner1 in partition_by_jump below
		FoldTree & f2
	) const;

	/// @brief partition the fold tree in two parts if the jump is disconnected.
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	void
	partition_by_jump( Size const jump_number, ObjexxFCL::FArray1D_bool & partner1 ) const ;

	/// @brief partition the fold tree in two parts if the jump is disconnected.
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	utility::vector1< bool >
	partition_by_jump( Size const jump_nr ) const;

	/// @brief partition the fold tree into n parts based on specified jumps.
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	utility::vector1< Size >
	partition_coloring( utility::vector1< Size > const & jump_numbers ) const;

	/// @brief partition the fold tree in two parts if a cut would be introduced between seqpos and seqpos+1
	void
	partition_by_residue(
		Size const seqpos,
		ObjexxFCL::FArray1D_bool & partner1
	) const ;

	/// @brief Returns the corresponding cutpoint position for jump  <jump_number>
	/// WARNING: if you look for all cutpoints by cycling thru jump_numbers you may be dissapointed
	/// you will get most likely the same cutpoint for several different jump_numbers
	/// however: the method cutpoint( nr ) will give you the number you are looking for
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	Size
	cutpoint_by_jump( Size const jump_number ) const;

	// these routines are for storing extra information about the jumps
	// you can specify which atoms should be the up/downstream atoms for
	// each jump. Note that this assumes that the topology and order of the
	// fold_tree are fixed, at least until this data is needed.
	//
	// the data are cleared upon re-ordering/wiring of the tree.

	/// @brief the jump atom on the staring side
	std::string
	upstream_atom( Size const jump_number ) const;

	/// @brief the jump atom on the stopping side
	std::string
	downstream_atom( Size const jump_number ) const;

	/// @brief Returns the direction (n2c, c2n) in which the given (peptide) residue is built during folding.
	int
	get_polymer_residue_direction( Size const seqpos ) const;

	/// @brief Returns the edge that builds the residue  <seqpos>
	/// Does not work for root atom (will fail)
	Edge const &
	get_residue_edge( Size const seqpos ) const;

	/// @brief Returns all edges that build a residue directly off of  <seqpos>
	utility::vector1< Edge >
	get_outgoing_edges( Size const seqpos ) const;

	/// @brief Return all jump Edges from the FoldTree
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	utility::vector1< Edge >
	get_jump_edges() const;

	/// @brief Returns all chemical edges from fold tree
	utility::vector1< Edge >
	get_chemical_edges() const;

	/// @brief  Get the number of the jump that builds (connects to) a given residue
	/// It's an error if the residue isn't built directly by a jump.
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	core::Size
	get_jump_that_builds_residue( Size const seqpos ) const;

	/// @brief  Get the residue that is immediately upstream of this residue (and tell us whether connection is jump or bond).
	core::Size
	get_parent_residue( Size const seqpos, bool & connected_by_jump ) const;

	/// @brief  Get the residue that is immediately upstream of this residue.
	core::Size
	get_parent_residue( Size const seqpos ) const;

	/// @brief define the specific atoms that should be connected by this jump
	/**
	This information can then be used in setting up the AtomTree from the
	FoldTree. Data is stored in the Edge corresponding to this Jump.
	If not specified, residue-specific defaults will be used.
	*/
	void
	set_jump_atoms(
		Size const jump_number,
		std::string const& upstream_atom,
		std::string const& downstream_atom,
		bool bKeepStubInResidue = false
	);

	//version of above but makes it permutation safe!
	void
	set_jump_atoms(
		Size const jump_number,
		core::Size res1,
		std::string const&  atom1,
		core::Size res2,
		std::string const&  atom2,
		bool bKeepStubInResidue = false
	);

	/// @brief this reorganizes upstream/downstream atoms of all jumps such that
	/// stubs are N-CA-C
	void put_jump_stubs_intra_residue();

	/// @brief this reorganizes upstream/downstream atoms of jumps that have flag
	/// keep_stub_in_resiue = true such that
	/// stubs are N-CA-C
	void reassign_atoms_for_intra_residue_stubs();


	// stream I/O //////////////////////////////////////////////////////////
	// these two should be inverses!
	// depend on I/O for the class: Edge

	/// @brief input operator
	friend std::istream & operator >>(std::istream & is, FoldTree & t);

	/// @brief output operator
	friend std::ostream & operator <<(std::ostream & os, FoldTree const & t);

	///////////////////////////// derived data //////////////////////
	// routines for retrieving derived data, typically fast

	/// @brief computes a fixed-length, hash-based identifier for this FoldTree,
	/// permitting efficient comparison between a pair of FoldTrees
	std::size_t hash_value() const;

	/// @brief easy output of string
	std::string to_string() const;

	/// @brief Returns the number of residues in the FoldTree
	///
	/// example(s):
	///     ft.nres()
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.num_jump
	///     FoldTree.simple_tree
	///     FoldTree.size
	inline Size nres() const;

	/// @brief Returns the number of jumps in the FoldTree
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	/// example(s):
	///     ft.num_jump()
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.jump_edge
	///     FoldTree.new_jump
	///     FoldTree.nres
	inline Size num_jump() const;

	/// @brief starting or stopping residue of a jump edge
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	Size jump_point( Size const lower_higher, Size const jump_number ) const;

	/// @brief Returns true if  <seqpos>  is a starting or stopping residue of a jump edge
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	/// example(s):
	///     ft.is_jump_point()
	/// See also:
	///     FoldTree
	///     FoldTree.is_cutpoint
	///     FoldTree.new_jump
	///     FoldTree.num_jump
	inline bool is_jump_point( Size const seqpos ) const;

	/// @brief Returns the cutpoint position of jump number  <cut>
	///
	/// example(s):
	///     ft.cutpoint(1)
	/// See also:
	///     FoldTree
	///     FoldTree.jump_edge
	///     FoldTree.is_cutpoint
	///     FoldTree.is_jump_point
	///     FoldTree.num_cutpoint
	///     FoldTree.num_jump
	inline Size cutpoint( Size const cut ) const;

	/// @brief Returns the number of cutpoints in the FoldTree
	///
	/// example(s):
	///     ft.num_cutpoint()
	/// See also:
	///     FoldTree
	///     FoldTree.cutpoint
	///     FoldTree.is_cutpoint
	///     FoldTree.nres
	///     FoldTree.num_jump
	inline Size num_cutpoint() const;

	/// @brief Returns true is position  <seqpos>  is a cutpoint
	///
	/// example(s):
	///     ft.is_cutpoint(37)
	/// See also:
	///     FoldTree
	///     FoldTree.cutpoint
	///     FoldTree.new_jump
	///     FoldTree.nres
	///     FoldTree.num_cutpoint
	///     FoldTree.num_jump
	inline bool is_cutpoint( Size const seqpos ) const;

	/// @brief cutpoint number for this residue
	inline Size cutpoint_map( Size const seqpos ) const;

	/// @brief Returns the jump edge with jump number  <jump_number>  (const)
	///
	/// example(s):
	///     ft.jump_edge(1)
	/// See also:
	///     FoldTree
	///     FoldTree.new_jump
	///     FoldTree.num_jump
	Edge const & jump_edge( Size const jump_number ) const;

	/// @brief Returns the jump edge with jump number  <jump_number>  (non-const)
	Edge & jump_edge( Size const jump_number );

	/// @brief get the jump_nr connected to jump upstream->downstream, returns 0 if not found
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	inline core::Size jump_nr( core::Size upstream_res, core::Size downstream_res ) const;


	/// @brief Returns true if the FoldTree is empty
	///
	/// example(s):
	///     ft.empty()
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.clear
	///     FoldTree.nres
	///     FoldTree.simple_tree
	///     FoldTree.size
	inline
	bool
	empty() const
	{
		return edge_list_.empty();
	}


	/// @brief Returns true if  <seqpos>  the the root
	///
	/// example(s):
	///     ft.empty()
	/// See also:
	///     FoldTree
	///     FoldTree.is_cutpoint
	///     FoldTree.is_jump_point
	///     FoldTree.nres
	inline
	bool
	is_root( Size const seqpos ) const
	{
		return ( !edge_list_.empty() && seqpos == begin()->start() );
	}


	/// @brief Returns true if  <seqpos>  is the root
	///
	inline
	bool
	possible_root( Size const & seqpos ) const{
		if ( seqpos == 1 ) return true;
		if ( seqpos == nres() ) return true;
		if ( is_cutpoint( seqpos ) ) return true;
		if ( is_cutpoint( seqpos - 1 ) ) return true;
		return false;
	}


	/// @brief Returns the root vertex position of the FoldTree
	///
	/// example(s):
	///     ft.empty()
	/// See also:
	///     FoldTree
	///     FoldTree.is_root
	///     FoldTree.clear
	///     FoldTree.simple_tree
	inline
	Size
	root() const
	{
		debug_assert( !empty() );
		return edge_list_.begin()->start();
	}


	/// @brief Return true if a jump exists between  <pos1>  and  <pos2>
	///
	/// example(s):
	///     ft.empty()
	/// See also:
	///     FoldTree
	///     FoldTree.check_fold_tree
	///     FoldTree.jump_edge
	///     FoldTree.new_jump
	///     FoldTree.nres
	inline
	bool
	jump_exists( Size const pos1, Size const pos2 ) const
	{
		for ( int i=1; i<= static_cast<int>( num_jump()); ++i ) {
			if ( ( ( pos1 == jump_point(1,i) && pos2 == jump_point(2,i) ) ||
					( pos1 == jump_point(2,i) && pos2 == jump_point(1,i) ) ) ) return true;
		}
		return false;
	}


	//  // old stuff -- still used? yeah abinitio mode- endbias check
	//inline ObjexxFCL::FArray1D_int const & get_jump_edge_count() const;
	Size count_fixed_residues(
		Size const begin_res,
		Size const size,
		Size & min_edge_count_out
	) const;

	// more routines for getting derived data
	// slow: creates vector each time, just for convenience
	/// @brief get all cutpoints
	utility::vector1< Size >
	cutpoints() const;

	/// @brief equal to operator
	/// @details  Checks that edges are in the same order and are equal
	/// That is, the order of the edges in the edge list is important to equality here
	friend
	bool
	operator==( FoldTree const & a, FoldTree const & b );

	/// @brief Not equal to operator
	friend bool operator!=( FoldTree const & a, FoldTree const & b );

	/// @brief Check if the two FoldTrees build in a similar fashion
	/// That is, if their roots are the same and if all edges correspond to each other
	/// -- Does not necessarily check that the order of the edges are the same.
	bool
	is_equivalent( FoldTree const & b ) const;

	/// @brief  Slide a polymer cutpoint from one location to another
	void
	slide_cutpoint( Size const current_cut, Size const target_cut );

	/// @brief change an existing jump to start and end in new positions
	void
	slide_jump( Size const jump_number, Size const new_res1, Size const new_res2 );

	/// @brief  Useful for removing a loop modeling jump+cut
	void
	delete_jump_and_intervening_cutpoint(
		Size jump_begin,
		Size jump_end,
		Size cut = 0
	);

	/// @brief  Useful for removing a loop modeling jump+cut
	void
	delete_jump_and_intervening_cutpoint(
		Size const jump_number
	);

	// AMW: I am elevating this to public - this is a concise and
	// helpful transformation when you don't want to destroy
	// your current FT state
	/// @brief helper function to try cutting an edge in a tree.
	bool cut_edge( Size const cut_point );


	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// private methods
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////

private:

	/// @brief delete a polymer residue (non-jump,non-root)
	void
	delete_seqpos_simple( Size const seqpos );

	/// @brief delete a root/jump_point residue
	void
	delete_jump_seqpos( Size const seqpos );


	/// @brief non-const begin iterator of edge_list
	inline
	iterator
	begin()
	{
		return edge_list_.begin();
	}

	/// @brief non-const end iterator of edge_list
	inline
	iterator
	end()
	{
		return edge_list_.end();
	}

	/// @brief  Helper function
	void
	add_vertex( Size const v );

	///////////////////////////// internal book-keeping /////////////
	// maintain things that depend on both topology and/or order
	/// @brief update fold tree topology (when edges are changed) if necessary
	inline bool check_topology() const;

	/// @brief update fold tree order (when edges are same but the order in the edge_list is changed) if necessary
	inline bool check_order() const;

	// private methods for updating derived data
	/// @brief update total number residues in the fold tree
	void update_nres() const;

	/// @brief update number of jumps in the fold tree
	void update_num_jump() const;

	/// @brief update jump residues list
	void update_jump_points() const;

	/// @brief update the index of jump edges in the edge list
	void update_jump_edge() const;

	/// @brief update cutpoints info in the fold tree
	void update_cutpoints() const;

	/// @brief update edge counts info
	void setup_edge_counts() const;


	// private edge_list_ modifiers
	/// @brief update edge labels based on whether edges are separating or not.
	void update_edge_labels();

	/// @brief cut an edge randomly based on probability without disconnecting fold tree
	bool
	cut_random_edge(
		ObjexxFCL::FArray1D_float const & cut_bias_sum,
		Size const nres_in
	);

	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////


	// the list of edges.
	/// @note vector for fast traversal, but re-ordering, deleting are slow.
	EdgeList edge_list_;

	// book-keeping, so we know when to update derived data
	/// @brief edges in the edge_list_ have been changed.
	mutable bool new_topology;

	/// @brief edges in the edge_list_ have been reordered.
	mutable bool new_order;

	// some derived data ////////////////////////////////////////
	// you don't need to worry about setting any derived data, they are all
	// calculated as needed from the fold_tree.
	// accessed by get_XXX where XXX is the data name
	// note that these are MUTABLE so they can be synced with the
	// edge_list_ on the fly inside "const" access methods

	/// @brief just the largest vertex in edge_list_
	mutable Size nres_;

	/// @brief number of jump edges (edges in edge_list_ with label>0)
	mutable Size num_jump_;

	/// @brief number of cutpoints in the fold tree.
	/// @note number_cutpoint_ == num_jump_ (if a connected tree)
	mutable Size num_cutpoint_;

	/// @brief jump number to jump residue number. dimensioned as (2,num_jump_)
	mutable utility::vector1< std::pair< Size, Size > > jump_point_;

	/// @brief whehter a residue is a jump_point, dimensioned as nres_
	mutable utility::vector1<bool> is_jump_point_;

	/// @brief cutpoint number to cutpoint residue number, dimesioned as num_cutpoint_.
	mutable utility::vector1<Size> cutpoint_;

	/// @brief residue number of cutpoint number, 0 if it is not a cutpoint. dimensioned as nres_.
	mutable utility::vector1<Size> cutpoint_map_;

	/// @brief whether a residue is a cutpoint, dimensioned as nres_
	mutable ObjexxFCL::FArray1D_bool is_cutpoint_; // nres

	/// @brief jump number to edge index number in the edge_list_, dimensioned as num_jump_.
	mutable utility::vector1<Size> jump_edge_;

	/// @brief dimensioned as nres_, see setup_edge_counts for more info
	mutable utility::vector1<Size> edge_count;

	/// @brief the minimum number in edge_count and jump_edge_count.
	mutable Size min_edge_count;

	/// @brief dimensioned as num_jump, see setup_edge_counts for more info
	mutable utility::vector1<Size> jump_edge_count;

	/// @brief Used as scratch space by a couple of functions.
	mutable ObjexxFCL::FArray1D_bool linked_, seen_;

	/// @brief computes fixed-size identifier for a string input
	boost::hash<std::string> hasher;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // FoldTree


///////////////////////////////////////////////////////////////////////////////
/// @brief check_topology and check_order handle the updating of data that depends on
/// tree topology and/or tree order.
/// @details any routine that depends on the stored derived data (eg any of the access
/// methods ) should call check_topology() or check_order() at the beginning.
/// @note also, any method function that changes the tree topology or order should set
/// the private data members new_topology and/or new_topology to true.
inline
bool
FoldTree::check_topology() const
{
	if ( !new_topology ) return false;
	new_topology = false;
	new_order = true;

	// insert updating calls here: ///////////////////////////
	update_nres();
	update_num_jump();
	update_jump_points();
	update_cutpoints();
	//setup_edge_counts(); moving this call to count_fixed_residues() since it is the only
	//function that requires up-to-date edge count data
	//////////////////////////////////////////////////////////
	return true;
}

/// @brief returns true if order has changed
/// see details for check_topology
inline
bool
FoldTree::check_order() const
{
	check_topology();
	if ( !new_order ) return false;
	new_order = false;

	// insert updating calls here: ///////////////////////////
	//////////////////////////////////////////////////////////
	update_jump_edge();

	//////////////////////////////////////////////////////////
	return true;
}


///////////////////////////////////////////////////////////////////////////////
/// @brief  routines for retrieving the derived data
/// will call check_topology and/or check_order first
inline
Size
FoldTree::nres() const
{
	check_topology();
	return nres_;
}

///////////////////////////////////////////////////////////////////////////////
inline
Size
FoldTree::num_jump() const
{
	check_topology();
	return num_jump_;
}


///////////////////////////////////////////////////////////////////////////////
/// @brief is seqpos a jump-point?
inline
bool
FoldTree::is_jump_point( Size const seqpos ) const
{
	check_topology();
	return is_jump_point_[seqpos];
}


///////////////////////////////////////////////////////////////////////////////
inline
Size
FoldTree::cutpoint( Size const cut ) const
{
	check_topology();
	return cutpoint_[cut];
}


///////////////////////////////////////////////////////////////////////////////
// number of cutpoints
inline
Size
FoldTree::num_cutpoint() const
{
	check_topology();
	return num_cutpoint_;
}

/// @brief whether a jump exists between these residues
inline
core::Size
FoldTree::jump_nr( core::Size const pos1, core::Size const pos2 ) const
{
	for ( Size i=1; i<= num_jump(); ++i ) {
		if ( ( ( pos1 == jump_point(1,i) && pos2 == jump_point(2,i) ) ||
				( pos1 == jump_point(2,i) && pos2 == jump_point(1,i) ) ) ) return i;
	}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
/// @brief cutpoint_map is the inverse of cutpoint_, ie it assigns each
/// sequence position that is a cutpoint to the cutpoint number
/// associated with that cutpoint (cutpoints are numbered in increasing
/// residue number from the beginning of the chain)
inline
Size
FoldTree::cutpoint_map( Size const seqpos ) const
{
	check_topology();
	// sanity
	debug_assert( (  is_cutpoint_(seqpos) && cutpoint_map_[seqpos]  > 0 &&
		cutpoint_[cutpoint_map_[seqpos]] == seqpos ) ||
		( !is_cutpoint_(seqpos) && cutpoint_map_[seqpos] == 0 ) );
	return cutpoint_map_[seqpos];
}


///////////////////////////////////////////////////////////////////////////////
inline
bool
FoldTree::is_cutpoint( Size const seqpos ) const
{
	check_topology();

	if ( seqpos > 0 && seqpos < nres_ ) {
		// 0 and nres count as cutpoints for is_cutpoint but they aren't internal cutpoints
		// in the foldtree so they dont appear in cutpoint_map_
		debug_assert( (  is_cutpoint_(seqpos) && cutpoint_map_[seqpos]  > 0 &&
			cutpoint_[cutpoint_map_[seqpos]] == seqpos ) ||
			( !is_cutpoint_(seqpos) && cutpoint_map_[seqpos] == 0 ) );
	}
	return is_cutpoint_(seqpos);
}

} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_kinematics_FoldTree )
#endif // SERIALIZATION


#endif // INCLUDED_core_kinematics_FoldTree_HH
