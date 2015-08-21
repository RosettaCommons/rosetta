// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/FoldGraph.hh
/// @brief FoldGraph - a fold-tree which takes two-way dependencies into account
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_FoldGraph_hh
#define INCLUDED_protocols_denovo_design_components_FoldGraph_hh

// Unit headers
#include <protocols/denovo_design/components/FoldGraph.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Core headers
#include <core/graph/Graph.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Numeric/Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <set>
#include <stack>

namespace protocols {
namespace denovo_design {
namespace components {

// Fold graph types
typedef std::set< core::Size > NodeSet;
typedef std::set< std::string > NamedNodeSet;
typedef utility::vector1< NodeSet > Solution;
typedef utility::vector1< NamedNodeSet > NamedSolution;

class FoldGraph : public utility::pointer::ReferenceCount {
public:
	FoldGraph( StructureData const & perm, core::pose::PoseCOP pose );

	/// @brief gives a fold tree based on the segments in the given permutation
	core::kinematics::FoldTree
	fold_tree( StructureData const & perm, utility::vector1< std::string > const & start_segment ) const;

	/// @brief given a segment name, returns the associated graph node index
	core::Size
	nodeidx( std::string const & segment ) const;

	/// @brief given a node index, returns the associated segment
	std::string const &
	segment( core::Size const nodeidx ) const;

	/// @brief adds an non-peptide edge to the foldgraph to indicate a non-covalent interaction
	void add_edge( std::string const & seg1, std::string const & seg2 );

	/// @brief adds a peptide edge to the foldgraph, indicating a direct covalent interaction
	void add_peptide_edge( std::string const & seg1, std::string const & seg2 );

	/// @brief generates a solution of movable segments to be used in remodel, based on the foldgraph
	Solution compute_best_solution(
		StructureData const & perm,
		utility::vector1< std::string > const & staple_loops,
		utility::vector1< std::string > const & cut_loops ) const;

	/// @brief takes a solution based on nodeidx and uses names instead
	NamedSolution named_solution( Solution const & solution ) const;

	/// @brief generates a loops object to be used in remodel, based on the foldgraph
	protocols::loops::LoopsOP create_loops(
		StructureData const & perm,
		utility::vector1< std::string > const & staple_loops,
		utility::vector1< std::string > const & cut_loops ) const;

	protocols::loops::LoopsOP create_loops(
		StructureData const & perm,
		Solution const & solution,
		utility::vector1< std::string > const & cut_loops ) const;

	/// @brief returns true if there is a peptide edge between the two given segments
	bool has_peptide_edge( std::string const & seg1, std::string const & seg2 ) const;

	/// @brief returns true if there is a non-peptide edge between the two given segments
	bool has_edge( std::string const & seg1, std::string const & seg2 ) const;

private:
	/// @brief recursive function to traverse graphs and build fold tree
	/// @param[parent_direction] -1 if the previous edge is a peptide edge going backward, 1 if the previous edge is going forward, and 0 if the previous edge is a jump
	void fold_tree_rec(
		core::kinematics::FoldTree & ft,
		NodeSet & visited,
		std::stack< core::Size > & node_stack,
		StructureData const & perm,
		std::string const & segment_name,
		int const parent_direction,
		bool const polymer_only ) const;

	/// @brief recursive inner function that traverses the foldgraph to create loops objects
	void create_loops_dfs(
		Solution & solutions,
		NodeSet & visited,
		core::Size const current_node,
		NodeSet const & cut_loop_nodes,
		StructureData const & perm ) const;

	/// @brief returns true if there is a peptide edge between the two given nodes
	bool has_peptide_edge( core::Size const n1, core::Size const n2 ) const;

	/// @brief returns true if there is a non-peptide edge between the two given nodes
	bool has_edge( core::Size const n1, core::Size const n2 ) const;

	/// @brief checks to see whether solutions found by DFS search can be combined.
	void add_combined_solutions( utility::vector1< Solution > & solutions ) const;

private:
	/// @brief checks a solution to ensure that covalently bound segments are not found both inside and outside of loops
	bool check_solution( StructureData const & perm, Solution const & solution ) const;

	std::map< std::string, core::Size > seg2node_;
	std::map< core::Size, std::string > node2seg_;
	core::graph::Graph g_;
	core::graph::Graph gpeptide_;
};

} // components
} // denovo_design
} // protocols

#endif
