// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core.pose.carbohydrates.GlycanTreeSet.hh
/// @brief Class to store info on all glycan trees of a pose
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author raemisch (raemisch@scripps.edu)

#ifndef INCLUDED_core_conformation_carbohydrates_GlycanTreeSet_hh
#define INCLUDED_core_conformation_carbohydrates_GlycanTreeSet_hh

//unit headers
#include <core/conformation/carbohydrates/GlycanTreeSet.fwd.hh>


// Core headers
#include <core/conformation/carbohydrates/GlycanTree.fwd.hh>
#include <core/conformation/carbohydrates/GlycanNode.fwd.hh>
#include <core/types.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include<map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {
namespace carbohydrates {


///@brief Class to store info on all glycan trees of a pose
///
/// This class contains a map, where the key is the first glycan residue
/// If the key would be the root residue (Asn), we couldn't work with pure glycans.
///
class GlycanTreeSet : public utility::VirtualBase  {

public:


	GlycanTreeSet();
	//GlycanTreeSet(GlycanTreeSet const & src);

	/// @brief Construct which populates the glycan trees.
	GlycanTreeSet( conformation::Conformation const & conf );

	~GlycanTreeSet() override;

	GlycanTreeSet( GlycanTreeSet const & src );

	GlycanTreeSetOP
	clone() const;

public:


	bool
	has_tree( core::Size const glycan_start_position ) const;

	/// @brief  Is this residue a part of a glycan tree?
	bool is_residue_in_tree( const core::uint glycan_residue ) const;

	///@brief Get a glycan tree corresponding to a particular starting residue.
	GlycanTreeCOP
	get_tree( core::Size const glycan_start_position ) const;

	///@brief Convenience function to the tree of a particular residue
	GlycanTreeCOP
	get_tree_containing_residue( core::Size const glycan_residue) const;


	///@brief Get a map of the tree start to the glycan tree.
	std::map< Size, GlycanTreeOP> const &
	get_tree_map() const;

	///@brief Get a list of all the glycan trees.
	utility::vector1< GlycanTreeCOP > const
	get_all_trees() const;

	utility::vector1< Size >
	get_start_points() const;

	///@brief Setup the glycan trees.  Done by other classes if contained in Pose.
	/// Unless you know what you are doing, you should not need to use this function.
	void
	setup_glycan_trees( conformation::Conformation const & pose);

public:

	///@brief Get the number of glycan trees
	core::Size
	n_trees() const;

	///@brief Get the number of glycan trees
	core::Size
	get_size() const;

	///@brief Get the number of glycan trees
	core::Size
	size() const;

	core::Size
	get_largest_glycan_tree_length() const;


	///@brief Get the largest glycan tree layer (from 0 to N).
	core::Size
	get_largest_glycan_tree_layer() const;

	///@brief Get the largest glycan tree layer (from 0 to N). of the glycan node residues passed in
	core::Size
	get_largest_glycan_tree_layer( utility::vector1< bool > const & subset ) const;


	///@brief Get the smallest glycan tree layer (from 0 to N).
	core::Size
	get_smallest_glycan_tree_layer() const;

	///@brief Get the smallest glycan tree layer (from 0 to N).
	core::Size
	get_smallest_glycan_tree_layer( utility::vector1< bool > const & subset ) const;

public:

	///////////////////////////////////////////////////////////////////////
	///                     ///
	///        Convenience functions that access stored data.           ///
	///                    ///
	///////////////////////////////////////////////////////////////////////

	///@brief Does the current glycan tree set have the node?
	bool
	has_node( core::Size glycan_residue ) const;

	///@brief Convenience function to get the node of a particular residue
	GlycanNodeCOP
	get_node( core::Size glycan_residue ) const;

	///@brief Convenience function to get the parent residue number from a GlycanNode
	core::Size
	get_parent( core::Size glycan_residue ) const;

	/// @brief Convenience function to get the main-chain child residue number from a GlycanNode.
	core::uint
	get_mainchain_child( core::Size glycan_residue ) const;

	///@brief Convenience function to get the starting position of a particular tree which contains the residue.
	/// This is the first Glycan Residue of the tree.
	/// Accessed from a stored GlycanTree class.
	core::Size
	get_tree_start_of_glycan_residue( core::Size resnum) const;

	///@brief Convenience function to get the root of a particular tree which contains the residue.
	/// This is the residue that 'roots' the glycan.  It should be protein.
	///  If this is 0, it means the glycan is a 'free' glycan.
	/// Accessed from a stored GlycanTree class.
	core::Size
	get_tree_root_of_glycan_residue( core::Size resnum ) const;

	///@brief Get the residue distance from this glycan residue to the root of the glycan.
	/// Used for Layer-based glycan sampling.
	core::Size
	get_distance_to_start( core::Size glycan_residue ) const;

	/// @brief Linkage number on the parent residue.
	/// @details an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment point on the
	/// parent monosaccharide residue; e.g., 4 specifies O4; n = 0 specifies that the residue at <seqpos> is a lower
	/// terminus or connected to a non-sugar.
	core::uint
	get_linkage_position( core::Size const resnum ) const;

	/// @brief Convenience function to get whether the glycosidic linkage between the
	///  residue and previous residue (parent residue) has an exocyclic
	///  carbon.
	bool
	has_exocyclic_glycosidic_linkage( core::Size resnum ) const;

public:

	///@brief Respond to a length change event.  Do not use this manually.  It is used by the GlcyanTreeObserver.
	///  The observer is attached to your pose and will respond.
	void
	on_length_change( core::conformation::signals::LengthEvent const & event );

private:

	std::map< Size, GlycanTreeOP > glycan_tree_set_; // The actual tree set - Maps first glycan residue to the tree
	std::map< Size, GlycanTreeOP > glycan_res_to_tree_; //Maps glycan residues to their glycan tree.

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //core
} //chemical
} //carbohydrates

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_carbohydrates_GlycanTreeSet )
#endif // SERIALIZATION

#endif //INCLUDED_core_pose_carbohydrates_GlycanTreeSet_hh





