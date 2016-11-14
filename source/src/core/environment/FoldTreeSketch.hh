// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/environment/FoldTreeSketch.hh
/// @brief A temporary "sketch" of a fold tree that keeps track of everything as a graph
///        and doesn't mind invalid states (cycles, detached bits, etc.)
///
/// @author Justin Porter

#ifndef INCLUDED_core_environment_FoldTreeSketch_hh
#define INCLUDED_core_environment_FoldTreeSketch_hh

// Unit Headers
#include <core/environment/FoldTreeSketch.fwd.hh>

// Package headers

// Project headers
#include <core/kinematics/FoldTree.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <core/types.hh>

// C++ Headers
#include <stack>
#include <set>

// ObjexxFCL Headers

namespace core {
namespace environment {

class EXCN_FTSketchGraph : public utility::excn::EXCN_Msg_Exception {
	typedef utility::excn::EXCN_Msg_Exception Parent;

public:
	EXCN_FTSketchGraph( Size by, Size on, std::string const& action, std::string const& reason );

	EXCN_FTSketchGraph( std::string const& message );
};

class FoldTreeSketch : public utility::pointer::ReferenceCount {

	class Node;
	typedef utility::pointer::shared_ptr< Node > NodeOP;
	typedef utility::pointer::shared_ptr< Node const > NodeCOP;
	typedef utility::pointer::weak_ptr< Node > NodeAP;
	typedef utility::pointer::weak_ptr< Node const > NodeCAP;

public:

	FoldTreeSketch();

	FoldTreeSketch( core::Size const length );

	FoldTreeSketch( core::kinematics::FoldTree const& );

	/// @brief copy the FoldTreeSketch
	/// @warning This is not a full deep copy. It doesn't copy the Node pointers deeply, which means
	//           you can accidentally modify the parent FTS through this one.
	FoldTreeSketch( FoldTreeSketch const& );

	/// @brief insert a jump between positions p1 and p2 in the graph
	void insert_jump( Size const p1,
		Size const p2 );

	/// @brief insert a cut after position l (i.e. between position l and l+1)
	void insert_cut( Size const l );

	/// @brief add a new (graph-detached) peptide-bonded strech of length l at the end of the sequence
	void append_peptide( Size const l );

	void append_residue() { return append_peptide( 1 ); }

	bool has_cut( Size const ) const;

	bool has_jump( Size const p1, Size const p2 ) const;

	utility::vector1< Size > const cycle( Size const start = 1 ) const;

	std::set< Size > remove_cycles();

	std::set< Size > remove_cycles( utility::vector1< Real > const& bias );

	/// @brief wrapper for the other render()
	core::kinematics::FoldTreeOP render() const;

	/// @brief build this graph into the FoldTree ft. The state of the graph must be valid here.
	void render( core::kinematics::FoldTree& ft ) const;

	Size nres() const;

	Size num_jumps() const;

	int num_cuts() const;

	/// @brief randomly insert a cut using the bias passed in.
	core::Size insert_cut( utility::vector1< Real > bias );

private:

	bool cuttable( utility::vector1< Size > const& cycle ) const;

	void range_check( Size const seqpos ) const;

	void visit( NodeCAP n, NodeCAP parent ) const;

	utility::vector1< NodeOP > nodes_;

	core::Size n_jumps_;
	int n_cuts_;

	class Node : public ReferenceCount {
	public:
		typedef std::set< NodeCAP, utility::pointer::owner_less< NodeCAP > > EdgeList;

	private:
		Node( core::Size i );

	public:
		static NodeOP newNode( core::Size i);

		void add_peptide_neighbor( NodeOP n );

		void add_jump_neighbor( NodeAP n );

		void rm_peptide_neighbor( NodeAP n );

		void rm_jump_neighbor( NodeAP n );

		bool has_peptide_neighbor( NodeCAP n ) const;

		bool has_jump_neighbor( NodeCAP n ) const;

		bool has_neighbor( NodeCAP n ) const;

		void collect_jumps( std::set< std::pair< Size, Size > >& v ) const;

		bool has_cycle( std::stack< NodeCAP >& path, NodeCAP start ) const;

		core::Size seqid() const;

		void unvisit() const { parent_.reset(); }

		bool visited() const { return !parent_.expired(); }

	private:

		void rm_neighbor( FoldTreeSketch::Node::EdgeList& list, Node* n );

		core::Size seqid_;

		NodeAP pep_prev_;
		NodeAP pep_next_;

		// Used during cycle search.
		mutable NodeCAP parent_;

		EdgeList jump_neighbors_;

		NodeAP this_weak_ptr_;

	};

}; // end FoldTreeSketch base class

} // environment
} // core

#endif //INCLUDED_core_environment_FoldTreeSketch_HH
