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
/// @details
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_FoldGraph_hh
#define INCLUDED_protocols_denovo_design_components_FoldGraph_hh

// Unit headers
#include <protocols/denovo_design/components/FoldGraph.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/loops/Loops.fwd.hh>

// Core headers
#include <core/graph/Graph.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Numeric/Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <map>
#include <set>
#include <stack>

namespace protocols {
namespace denovo_design {
namespace components {

// Fold graph types
class FoldGraph : public utility::pointer::ReferenceCount {
public:
	typedef std::stack< core::Size > NodeStack;
	typedef std::set< core::Size > NodeSet;
	typedef utility::vector1< NodeSet > Solution;
	typedef utility::vector1< Solution > Solutions;
	typedef std::map< core::Size, std::string > NodeToSegmentMap;
	typedef std::map< std::string, core::Size > SegmentToNodeMap;
	typedef std::map< core::Size, SegmentNameList > MovableGroupToSegmentNameListMap;

public:
	FoldGraph( StructureData const & perm );
	virtual ~FoldGraph();

	/// @brief gives a fold tree based on the segments in the given permutation
	core::kinematics::FoldTree
	fold_tree( SegmentNames const & start_segment ) const;

	/// @brief given a segment name, returns the associated graph node index
	core::Size
	nodeidx( std::string const & segment ) const;

	/// @brief given a node index, returns the associated segment
	std::string const &
	segment( core::Size const nodeidx ) const;

	/// @brief generates a solution of movable segments to be used in remodel, based on the foldgraph
	Solution
	compute_best_solution( SegmentNames const & staple_loops ) const;

	/// @brief generates a loops object to be used in remodel, based on the foldgraph
	protocols::loops::LoopsOP
	create_loops( SegmentNames const & staple_loops ) const;

	protocols::loops::LoopsOP
	create_loops( Solution const & solution ) const;

	/// @brief const access to StructureData
	StructureData const &
	sd() const;

private:
	/// @brief adds an non-peptide edge to the foldgraph to indicate a non-covalent interaction
	void
	add_edge( std::string const & seg1, std::string const & seg2 );

	/// @brief adds a peptide edge to the foldgraph, indicating a direct covalent interaction
	void
	add_peptide_edge( std::string const & seg1, std::string const & seg2 );

	/// @brief returns true if there is a peptide edge between the two given segments
	bool
	has_peptide_edge( std::string const & seg1, std::string const & seg2 ) const;

	/// @brief returns true if there is a non-peptide edge between the two given segments
	bool
	has_edge( std::string const & seg1, std::string const & seg2 ) const;

	/// @brief populates segment to node number map using the StructureData object
	void
	build_seg2node();

	/// @brief populates node number to segment map using a segment to node number map
	void
	build_node2seg( SegmentToNodeMap const & seg2node );

	/// @brief propagaes the internal set of cutpoint-containing segments
	void
	find_cutpoints();

	/// @brief recursive function to traverse graphs and build fold tree
	/// @param[parent_direction] -1 if the previous edge is a peptide edge going backward, 1 if the previous edge is going forward, and 0 if the previous edge is a jump
	/// @param[parent_resid] parent residue
	void
	fold_tree_rec(
		core::kinematics::FoldTree & ft,
		NodeSet & visited,
		NodeStack & node_stack,
		std::string const & segment_name,
		core::Size const parent_resid,
		int const parent_direction,
		bool const polymer_only ) const;

	/// @brief recursive inner function that traverses the foldgraph to create loops objects
	void
	create_loops_dfs(
		Solution & solutions,
		NodeSet & visited,
		core::Size const current_node ) const;

	/// @brief returns true if there is a peptide edge between the two given nodes
	bool
	has_peptide_edge( core::Size const n1, core::Size const n2 ) const;

	/// @brief returns true if there is a non-peptide edge between the two given nodes
	bool
	has_edge( core::Size const n1, core::Size const n2 ) const;

	/// @brief checks to see whether solutions found by DFS search can be combined.
	void
	add_combined_solutions( Solutions & solutions ) const;

	void
	add_non_polymeric_connections( Solutions & solutions ) const;

	Solution const &
	select_best_solution( Solutions const & solutions ) const;

	/// @brief sorts solutions such that the ones maximizing segments w/ same MG are fixed.
	void
	sort_solutions( Solutions & solutions ) const;

	/// @brief checks a solution to ensure that covalently bound segments are not found both inside and outside of loops
	bool
	check_solution( Solution const & solution ) const;

	int
	find_cutpoint_in_range( int const min_res, int const max_res ) const;

private:
	StructureDataCOP sd_;
	SegmentToNodeMap seg2node_;
	NodeToSegmentMap node2seg_;
	NodeSet cutpoints_;
	core::graph::Graph g_;
	core::graph::Graph gpeptide_;
};

// helper types for outputting sets of nodes
typedef std::set< std::string > NamedNodeSet;

class NamedSolution : public utility::vector1< NamedNodeSet > {
public:
	NamedSolution( FoldGraph const & fg, FoldGraph::Solution const & solution );

	friend std::ostream &
	operator<<( std::ostream & os, NamedSolution const & solution );
};

class NamedSolutions : public utility::vector1< NamedSolution > {
public:
	NamedSolutions( FoldGraph const & fg, FoldGraph::Solutions const & solutions );

	friend std::ostream &
	operator<<( std::ostream & os, NamedSolutions const & solutions );
};


} // components
} // denovo_design
} // protocols

#endif
