// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core.pose.carbohydrates.GlycanTree.cc
/// @brief Class to store info a glycan tree
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author raemisch (raemisch@scripps.edu)

#include <core/conformation/carbohydrates/GlycanTree.hh>
#include <core/types.hh>
#include <core/conformation/carbohydrates/GlycanNode.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/LengthEvent.hh>


// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>


// C++ headers
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <basic/Tracer.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>

#endif // SERIALIZATION


static THREAD_LOCAL basic::Tracer TR( "core.conformation.carbohydrates.GlycanTree" );


namespace core {
namespace conformation {
namespace carbohydrates {

GlycanTree::GlycanTree():
	utility::pointer::ReferenceCount()
{

}


// @brief standard constructor
// -> Get all nodes of a glycan tree
GlycanTree::GlycanTree( conformation::Conformation const & conf, Size const start_pos ):
	utility::pointer::ReferenceCount()

{
	setup_glycan_nodes( conf, start_pos );

}// constructor

GlycanTree::GlycanTree( Size const start_pos ):
	utility::pointer::ReferenceCount(),
	start_pos_(start_pos)
{

}

GlycanTree::GlycanTree( GlycanTree const & src ):
	utility::pointer::ReferenceCount(),
	branch_tips_(src.branch_tips_),
	start_pos_(src.start_pos_),
	root_(src.root_)
{
	tree_.clear();
	for ( auto const & kv : src.tree_ ) {
		tree_[kv.first] = GlycanNodeOP( new GlycanNode( *kv.second));
	}

}




GlycanTreeOP
GlycanTree::clone() const {
	return GlycanTreeOP( new GlycanTree( *this ) );
}


GlycanTree::~GlycanTree(){}

void
GlycanTree::setup_glycan_nodes(conformation::Conformation const & conf, Size const start_pos) {
	tree_.clear();
	start_pos_ = start_pos;


	// Find all residues belonging to the tree starting at start_pos
	std::pair<utility::vector1< Size >, utility::vector1< Size>> branch_and_tips = get_carbohydrate_residues_and_tips_of_branch(conf, start_pos, true /*include_start_position*/);

	utility::vector1< Size > branch_residues = branch_and_tips.first;

	branch_tips_ = branch_and_tips.second;

	if ( *(branch_residues.begin()) != start_pos ) {
		std::string curr_pos = utility::Real2string(start_pos,0);
		std::string msg = "ERROR: Finding tree residues for start position " + curr_pos + " failed!";
		std::string msg2 = utility::to_string(branch_residues);

		utility_exit_with_message(msg+" "+msg2);
	}

	for ( Size pos : branch_residues ) {
		// Create a new GlycanNode instance for this residue
		// The node will contain info about all direct downstream connections
		// Add the node to the tree
		tree_[pos] = GlycanNodeOP( new GlycanNode(conf, start_pos, pos));
	}

	root_ = find_seqpos_of_saccharides_parent_residue( conf.residue(start_pos) );
}

bool
GlycanTree::is_connected() const {
	if ( root_ == 0 ) {
		return false;
	} else {
		return true;
	}
}

bool
GlycanTree::has_glycan_residue(core::Size resnum) const {
	bool result = tree_.count( resnum ) != 0;
	return result;
}

core::Size
GlycanTree::get_size() const {
	return tree_.size();
}

core::Size
GlycanTree::size() const {
	return tree_.size();
}

core::Size
GlycanTree::get_root() const {
	return root_;
}

core::Size
GlycanTree::get_start() const {
	return start_pos_;
}

GlycanNodeCOP
GlycanTree::get_node(core::Size resnum) const{
	return tree_.at( resnum );
}

utility::vector1< Size > const &
GlycanTree::get_tips() const {
	return branch_tips_;
}

utility::vector1< Size >
GlycanTree::get_residues() const {
	utility::vector1< Size > keys;

	for ( auto const & kv : tree_ ) {
		keys.push_back(kv.first);
	}
	return keys;
}

void
GlycanTree::update_on_length_change(core::conformation::signals::LengthEvent const & event){

	using namespace core::conformation::signals;

	core::conformation::Conformation const & conf = *event.conformation;
	core::id::SequenceMapping smap( event );

	// Takes just as long as repopulating.
	//New GlycanTree that came from adding a new free-standing tree OR
	// We deleted a branch point, which created two glycans!
	if ( tree_.size() == 0 ) {
		setup_glycan_nodes(conf, start_pos_);
	} else {

		std::map< Size, GlycanNodeOP > new_trees;

		//Update root
		root_ = find_seqpos_of_saccharides_parent_residue( conf.residue(start_pos_) );

		//Update branch residues
		std::pair<utility::vector1< Size >, utility::vector1< Size>> branch_and_tips = get_carbohydrate_residues_and_tips_of_branch( conf, start_pos_, true );

		utility::vector1< Size > branch_residues = branch_and_tips.first;
		branch_tips_ = branch_and_tips.second;

		if ( *(branch_residues.begin()) != start_pos_ ) {
			std::string curr_pos = utility::Real2string(start_pos_,0);
			std::string msg = "ERROR: Finding tree residues for start position " + curr_pos + " failed!";
			utility_exit_with_message(msg);
		}

		//Update numbering and connectivity
		for ( auto & kv : tree_ ) {

			GlycanNodeOP node = kv.second;
			Size old_resnum = kv.first;


			if ( (event.tag == LengthEvent::RESIDUE_PREPEND || event.tag == LengthEvent::RESIDUE_APPEND) && old_resnum == event.residue->seqpos() ) {
				node->update_connectivity_data( conf );
				new_trees[ old_resnum ] = node;
			} else if ( event.tag == LengthEvent::RESIDUE_DELETE && old_resnum == event.position ) {
				continue;
			} else {
				Size new_resnum = smap.get_corresponding_residue_in_current(old_resnum);
				if ( new_resnum != 0 ) {
					node->remap_residue(start_pos_, new_resnum);
					node->update_connectivity_data( conf );
				}
				new_trees[ new_resnum ] = node;
			}




		}
		//Add any new nodes from append or prepend residue operations.
		for ( Size node_residue : branch_residues ) {

			//New residue.
			if ( ! new_trees.count( node_residue ) ) {
				GlycanNodeOP new_tree = GlycanNodeOP( new GlycanNode( conf, start_pos_, node_residue ) );
				new_trees[node_residue] = new_tree;

			}

		}
		tree_ = new_trees;
	}

}

void
GlycanTree::update_start_position( core::Size const start_pos){
	start_pos_ = start_pos;

}





} //core
} //chemical
} //carbohydrates

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::carbohydrates::GlycanTree::save( Archive & arc ) const {
	arc( CEREAL_NVP( tree_ ) ); // std::map< Size, GlycanNodeOP >
	arc( CEREAL_NVP( branch_tips_ ) ); // utility::vector< core::Size >
	arc( CEREAL_NVP( start_pos_ ) ); // core::Size
	arc( CEREAL_NVP( root_ ) ); //core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::carbohydrates::GlycanTree::load( Archive & arc ) {
	arc( tree_ );
	arc( branch_tips_ );
	arc( start_pos_ );
	arc( root_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::carbohydrates::GlycanTree );
CEREAL_REGISTER_TYPE( core::conformation::carbohydrates::GlycanTree )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_carbohydrates_GlycanTree )
#endif // SERIALIZATION




