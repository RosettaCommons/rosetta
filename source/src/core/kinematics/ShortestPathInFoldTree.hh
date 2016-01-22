// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/ShortestPathInFoldTree.hh
/// @brief  Fold tree helper class
/// @author Oliver Lange


#ifndef INCLUDED_core_kinematics_ShortestPathInFoldTree_hh
#define INCLUDED_core_kinematics_ShortestPathInFoldTree_hh


// Unit headers
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>

// Package Headers
#include <core/kinematics/FoldTree.fwd.hh>

// Rosetta Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

// C++ Headers
#include <map>

#include <core/kinematics/Edge.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {

/// @brief A helper class that can tell the shortest distance of residues in a given Foldtree
/// The fold-tree is parsed only in the beginning such that individual queries are fast
class ShortestPathInFoldTree;
typedef utility::pointer::shared_ptr< ShortestPathInFoldTree > ShortestPathInFoldTreeOP;

class ShortestPathInFoldTree : public utility::pointer::ReferenceCount {
public:
	// types
	typedef utility::vector1< Edge > EdgeList;

public:
	/// @brief cs-tor, parses fold-tree and caches important distances:
	///        memory N^2+M*5    N: 2 x number of jumps     M: number of residues
	ShortestPathInFoldTree( core::kinematics::FoldTree const& f );

	/// @brief copy constructor
	ShortestPathInFoldTree( ShortestPathInFoldTree const & src );

	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~ShortestPathInFoldTree();


	/// @brief returns the shortest distance of two residues going along Fold-Tree edges.
	Size dist( Size pos1, Size pos2 ) const;

	/// @brief returns the shortest distance for the two residues that are furthest apart
	Size max_dist() const {
		return max_dist_;
	}

private:

	/// @brief build table that gives for each residue the distance to the next jump(s) (1 or 2)
	void build_peptide_table( core::kinematics::FoldTree const& f);

	/// @brief build N^2 table of distances between jump-residues
	void build_jumpres_distmap( core::kinematics::FoldTree const& f);

	/// @brief if seqpos is a jump_res resturn its internal number
	int get_jump( Size seqpos ) const {
		std::map< Size, Size >::const_iterator fit = jump_res_.find ( seqpos );
		if ( fit != jump_res_.end() ) {
			return fit->second;
		} else {
			return -1;
		}
	}

	/// @brief helper of build_jumpres_distmap
	void compute_dist_map( core::kinematics::FoldTree const&);
	void init_dist_map( EdgeList const& );

	/// private data ===============================

	/// @brief map from pose-numbering to internal numbering --> jump_res are counted from 1...N
	/// first: pose-numbering, second internal number
	std::map< Size, Size > jump_res_;

	/// @brief 2D array for distances between each pair of jump_residues
	ObjexxFCL::FArray2D_int node_dist_;

	/// @brief 1D array for distances of residue j to
	/// each of max 2 jump_residues that are on the same peptide edge in the fold-tree
	///  --> for each residue 1+2x2 entries.
	///  < edgenr >  <jump1> <dist1>  <jump2> <dist2>
	///  jump1 upstream jump_residue (internal numbering) jump2 downstream jump-residues
	///  no up/downstream jumps --> -1
	ObjexxFCL::FArray2D_int res2jumps_;

	/// @brief for error checking what is nres_ in the fold-tree
	Size nres_;

	/// @brief if fold-tree is simple (no jumps) don't compute anything.
	bool simple_fold_tree_;

	/// @brief the furthest distance a query to dist() can return
	mutable Size max_dist_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ShortestPathInFoldTree();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // kinematics
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_kinematics_ShortestPathInFoldTree )
#endif // SERIALIZATION


#endif
