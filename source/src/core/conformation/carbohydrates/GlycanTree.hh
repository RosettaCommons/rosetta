// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core.pose.carbohydrates.GlycanTree.hh
/// @brief Class to store info a glycan tree
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Sebastian Rämisch (raemisch@scripps.edu)

#ifndef INCLUDED_core_conformation_carbohydrates_GlycanTree_hh
#define INCLUDED_core_conformation_carbohydrates_GlycanTree_hh

#include <core/conformation/carbohydrates/GlycanTree.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <core/conformation/carbohydrates/GlycanNode.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <core/id/SequenceMapping.hh>

// C++ headers
#include <map>
#include <utility>
#include <vector>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

/////////////////
// This class gathers and stores information on entire glycan trees. Typical glycan trees
// contain multiple chains but with consective residue numbers. Gathering all the information about
// connectivity of multiple chains invokes quite a number of functions that do funky things to find out which
// residue belongs to which glycan tree and where are they connected to the protein. This container is meant for
// storing this information instead of running all those functions (see carbohydrates/util.cc) repeatedly.

namespace core {
namespace conformation {
namespace carbohydrates {



///@brief Class to store info a glycan tree
class GlycanTree : public utility::pointer::ReferenceCount {

public:



	GlycanTree();
	//GlycanTree(GlycanTree const & src);

	virtual ~GlycanTree();

	GlycanTree( GlycanTree const & src );

	GlycanTreeOP
	clone() const;

	// standard constructor
	GlycanTree( conformation::Conformation const & conf, Size const start_pos);

	// constructor that still needs updating by pose
	GlycanTree(Size const start_pos);


public:

	///@brief Is this tree connected to a [protein] or is it free-standing?
	bool
	is_connected() const;

	///@brief Does this Tree have a particular residue?
	bool
	has_glycan_residue( core::Size resnum ) const;

	///@brief Get the length of this tree
	core::Size
	get_size() const;

	///@brief Get the length of this tree
	core::Size
	size() const;

	///@brief Get the first glycan residue number of this tree.
	/// This is also used to identify the tree.
	/// This is different than the root!
	core::Size
	get_start() const;

	///@brief Get the root of this tree. The connecting residue.
	/// If it is zero, it means this tree is a free-glycan and not connected!
	core::Size
	get_root() const;

public:

	///@brief Get the GlycanNode of a particular glycan residue
	/// GlycanNode has lots of information on the particular residue as part of this glycan tree.
	///
	GlycanNodeCOP
	get_node( core::Size resnum ) const;

public:

	///@brief Get the tips (ending residues) of each foliage end.
	utility::vector1< core::Size > const &
	get_tips() const;

	///@brief Get all the residues of this glycan
	utility::vector1< core::Size >
	get_residues() const;

public:

	///@brief Populate the glycan nodes of this GlycanTree
	void
	setup_glycan_nodes(conformation::Conformation const & conf, core::Size const start_pos);

	///@brief Update the nodes connectivity and the tree root after a length-change event.
	void
	update_on_length_change( core::conformation::signals::LengthEvent const & event);

	///@brief Update the starting position if the start has been deleted.
	void
	update_start_position( core::Size const start_pos );


private:
	// The glycan tree istself is just a map< 'residue number', GlycanNode instance >
	std::map< Size, GlycanNodeOP > tree_;

	utility::vector1< Size > branch_tips_;

	core::Size start_pos_;
	Size root_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_carbohydrates_GlycanTree )
#endif // SERIALIZATION


#endif //INCLUDED_core_pose_carbohydrates_GlycanTree_hh





