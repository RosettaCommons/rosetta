// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core.pose.carbohydrates.GlycanNode.cc
/// @brief Class to store info a node (residue) within a glycan tree
/// @author Jared Adolf-Bryfogle(jadolfbr@gmail.com)
/// @author raemisch (raemisch@scripps.edu)

#include <core/conformation/carbohydrates/GlycanNode.hh>
#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResConnID.hh>


// C++ headers
#include <utility>

static THREAD_LOCAL basic::Tracer TR( "core.pose.carbohydrates.GlycanNode" );


namespace core {
namespace conformation {
namespace carbohydrates {

GlycanNode::GlycanNode():
	utility::pointer::ReferenceCount()
{

}

GlycanNode::~GlycanNode(){}

// Copy constructor
//GlycanNode::GlycanNode( GlycanNode const & ) {
//}





// @brief standard constructor
GlycanNode::GlycanNode( conformation::Conformation const & conf, Size const tree_start_pos, Size const pos ) {
	setup_info( conf, tree_start_pos, pos);
}// construtor}

GlycanNode::GlycanNode( GlycanNode const & src ):
	utility::pointer::ReferenceCount(),
	node_residue_( src.node_residue_ ),
	parent_residue_( src.parent_residue_ ),
	tree_start_residue_( src.tree_start_residue_ ),
	children_( src.children_ ),
	downstream_connections_( src.downstream_connections_ ),
	glycan_position_( src.glycan_position_ ),
	distance_to_root_( src.distance_to_root_ ),
	has_exocyclic_linkage_( src.has_exocyclic_linkage() ),
	linkage_position_( src.linkage_position_ ),
	mainchain_child_( src.mainchain_child_ )
{

}

GlycanNodeOP
GlycanNode::clone() const {
	return GlycanNodeOP( new GlycanNode( *this ) );
}

void
GlycanNode::setup_info( conformation::Conformation const & conf, Size const tree_start_pos, Size const pos){


	//////////////////////////////////
	// Initialize member variables  //
	//////////////////////////////////

	node_residue_ = pos;
	tree_start_residue_ = tree_start_pos;
	update_connectivity_data(conf);

	downstream_connections_.clear();
}

void
GlycanNode::remap_residue( Size new_start, Size new_node){
	tree_start_residue_ = new_start;
	node_residue_ = new_node;
}


void
GlycanNode::update_connectivity_data(conformation::Conformation const & conf){
	typedef std::pair< Size, Size > Downstream; // residue number, connection number (index)
	typedef std::pair<Size,Downstream> Connection; // Connection ID,dowstream connection


	conformation::Residue const & this_residue  = conf.residue(node_residue_);

	parent_residue_ = find_seqpos_of_saccharides_parent_residue( this_residue );

	// Find children, if any, and populate connections
	//chemical::ResidueType const & residue_type = conf.residue_type(node_residue_);
	Size connections = this_residue.n_possible_residue_connections(); //Want the index to match here.

	children_.clear();
	for ( core::Size conn = 1; conn <= connections; ++conn ) {

		Size child = conf.residue( node_residue_ ).connected_residue_at_resconn( conn );

		if ( child != 0 && child != parent_residue_ && conf.residue( child ).is_carbohydrate() ) {

			children_.push_back( child );

			Size upstream_conn_id = conn;
			Size downstream_conn_id = conn;

			if ( this_residue.connect_map(conn).resid() == child ) {
				downstream_conn_id = this_residue.connect_map(conn).connid();
			}
			Downstream downstream_conn = std::make_pair(child,downstream_conn_id);
			Connection child_connection = std::make_pair(upstream_conn_id, downstream_conn);
			// Store this connction in the downstream_connection member variable
			downstream_connections_.push_back(child_connection);

		}
	}

	glycan_position_ = get_glycan_position_from_resnum( conf, tree_start_residue_, node_residue_ );
	distance_to_root_ = carbohydrates::get_distance_to_root( conf, node_residue_ );
	mainchain_child_ = find_seqpos_of_saccharides_mainchain_child( conf.residue( node_residue_ ) );
	has_exocyclic_linkage_ = has_exocyclic_glycosidic_linkage(conf, node_residue_ );

	if ( parent_residue_ != 0 ) {
		linkage_position_ = get_linkage_position_of_saccharide_residue( conf.residue( node_residue_), conf.residue(parent_residue_));
	} else {
		linkage_position_ = 0;
	}

}

Size
GlycanNode::get_mainchain_child() const {
	return mainchain_child_;
}

} //core
} //chemical
} //carbohydrates






