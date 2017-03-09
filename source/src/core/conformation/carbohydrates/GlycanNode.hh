// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core.pose.carbohydrates.GlycanNode.hh
/// @brief Class to store info a node (residue) within a glycan tree
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author raemisch (raemisch@scripps.edu)


#ifndef INCLUDED_core_conformation_carbohydrates_GlycanNode_hh
#define INCLUDED_core_conformation_carbohydrates_GlycanNode_hh

#include <core/conformation/carbohydrates/GlycanNode.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

// Project headers
#include <core/types.hh>

// C++ headers
#include <vector>
#include <utility>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace core {
namespace conformation {
namespace carbohydrates {

///@brief Class to store info a node (residue) within a glycan tree
class GlycanNode : public utility::pointer::ReferenceCount {

	///   Typedefs for storing connections
	typedef Size upstream_atom;
	typedef std::pair< Size, Size > downstream_atom; // residue number, atom number
	typedef std::pair<upstream_atom,downstream_atom> connection; // this is an edge :)


public:

	GlycanNode();
	//GlycanNode(GlycanNode const & src);

	/// @brief Standard constructor
	///  tree_start_pos is the first residue of the tree in which this node is located in.
	///  pos is the glycan residue we are setting up.
	GlycanNode(conformation::Conformation const & conf, Size const tree_start_pos, Size const pos);

	virtual ~GlycanNode();

	GlycanNode( GlycanNode const & src );

	GlycanNodeOP
	clone() const;


public:

	///@brief Get the residue corresponding to this node.
	core::Size
	get_resnum() const {
		return node_residue_;
	}
	
	///@brief Get the residue distance from this glycan residue to the start of the glycan.
	/// Used for Layer-based glycan sampling.
	core::Size
	get_distance_to_start() const {
		return distance_to_start_;
	}



	///@brief Get parent residue number
	core::Size
	get_parent() const {
		return parent_residue_;
	}

	/// @brief Linkage number on the parent residue.
	/// @details an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment point on the
	/// parent monosaccharide residue; e.g., 4 specifies O4; n = 0 specifies that the residue at <seqpos> is a lower
	/// terminus or connected to a non-sugar.
	core::Size
	get_linkage_position() const {
		return linkage_position_;
	}

	/// @brief Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic
	/// carbon.
	bool
	has_exocyclic_linkage() const {
		return has_exocyclic_linkage_;
	}

public:

	///@brief Get all connected downstream residues
	/// (Downstream = this->end of glycan)
	utility::vector1< Size > const &
	get_children() const {
		return children_;
	}

	//core::Size
	//get_child_at_linkage_position() const;

	///@brief Get the downstream residue connecting to this residue that is part of the mainchain.
	/// If this has NO mainchain connection, we return zero.
	/// (Downstream = this->end of glycan)
	core::Size
	get_mainchain_child() const;

	///@brief Get all downstream connections
	/// (Downstream = this->end of glycan
	utility::vector1< connection > const &
	get_downstream_connections() {
		return downstream_connections_;
	}

public:

	/// @brief Setup all data for this glycan residue.
	///
	///  tree_start_pos is the first residue of the tree in which this node is located in.
	///  pos is the glycan residue we are setting up.
	void
	setup_info( conformation::Conformation const & conf, Size const tree_start_pos, Size const pos);

	///@brief Update connectivity following a length-change event.
	void
	update_connectivity_data( conformation::Conformation const & conf );

	///@brief Remap residue following length change.
	void
	remap_residue( Size new_start, Size new_node );

private:
	///   Member variables

	///@brief The node's pose-internal residue number
	Size node_residue_;

	///@brief The node's parent pose-internal residue number
	Size parent_residue_;

	///@brief The starting residue of this Node's Tree.  Instead of holding a whole AP.
	Size tree_start_residue_; ///

	///@brief directly connected residues (downstream residues)
	/// (Downstream = this->end of glycan
	utility::vector1< Size > children_;

	///@brief All connections of this node, i.e. connections to children:
	/// (Downstream = this->end of glycan
	utility::vector1< connection > downstream_connections_;

	///@brief Position representation of the glycan residue.  From 1 to N, N being the size of the glyan tree.
	Size glycan_position_;

	///@brief Residue distance to the glycan root.  Used for Layer-Based Tree sampling.
	Size distance_to_start_;
	
	/// @brief Whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic
	/// carbon.
	bool has_exocyclic_linkage_;

	/// @brief Linkage number on the parent residue.
	/// @details an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment point on the
	/// parent monosaccharide residue; e.g., 4 specifies O4; n = 0 specifies that the residue at <seqpos> is a lower
	/// terminus or connected to a non-sugar.
	Size linkage_position_;

	Size mainchain_child_;


};


} //core
} //chemical
} //carbohydrates



#endif //INCLUDED_core_pose_carbohydrates_GlycanNode_hh





