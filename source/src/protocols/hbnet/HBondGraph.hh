// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/hbnet/HBondGraph.hh
/// @brief class headers for HBondGraph, HBondNode, and HBondEdge
/// @details This class is used to store and traverse data used for HBNet's monte carlo branching protocol. Most (if not all) of the information held in this graph is from HBNet's RotamerSets and InteractionGraph. Nodes in this graph represent rotamers from the rotamer set (node id == global rotamer id) and edges respresent hydrogen bonds.
/// @author Jack Maguire, jack@med.unc.edu

#ifndef INCLUDED_protocols_hbnet_HBondGraph_HH
#define INCLUDED_protocols_hbnet_HBondGraph_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/hbnet/HBondGraph.fwd.hh>

#include <utility/graph/Graph.hh>
#include <core/types.hh>

#include <utility/graph/unordered_object_pool.hpp>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <set>

namespace protocols {
namespace hbnet {

///@brief Each HBondNode represents a rotamer from the RotamerSets object
class HBondNode : public utility::graph::Node {

public:

	//constructor
	HBondNode( utility::graph::Graph*, core::Size node_id );

	HBondNode( utility::graph::Graph*, core::Size node_id, core::Size mres_id, core::Size rotamer_id );

	//destructor
	~HBondNode();

public:

	void copy_from( utility::graph::Node const * source ) override;

	void print() const override;

	///@brief get molten residue id for this rotamer
	inline unsigned int moltenres() const{
		return mres_id_;
	}

	///@brief set molten residue id for this rotamer.
	inline void set_moltenres( core::Size mres ){
		mres_id_ = ( unsigned int ) ( mres );
	}

	///@brief get local rotamer id (local to the residue position)
	///@details this is equivalent to core::pack::rotamer_set::RotamerSetsBase::rotid_on_moltenresidue( this->get_node_index() )
	inline unsigned int local_rotamer_id() const{
		return rotamer_id_;
	}

	///@brief set local rotamer id (local to the residue position).
	inline void set_local_rotamer_id( core::Size rot_id ){
		rotamer_id_ = ( unsigned int ) ( rot_id );
	}

	///@brief duplicate interface for getting the global rotamer id. Identical to this->get_node_index()
	///details this is equivalent to core::pack::rotamer_set::RotamerSetsBase::nrotamer_offset_for_moltenres( this->moltenres() ) + this->local_rotamer_id()
	inline core::Size global_rotamer_id() const{
		return get_node_index();
	}

	///@brief keep track of another node (rotamer) that this clashes with. You do not need to call this for all of the rotamers that share a residue position.
	inline void register_clash( core::Size node_id ){
		ids_of_clashing_nodes_.insert(
			std::upper_bound( ids_of_clashing_nodes_.begin(), ids_of_clashing_nodes_.end(), node_id ),
			node_id
		);
	}

	///@brief does this node clash with another node (at another residue position)?
	inline bool clashes( core::Size node_id ) const{
		return std::binary_search( ids_of_clashing_nodes_.begin(), ids_of_clashing_nodes_.end(), node_id );
	}

	core::Size count_static_memory() const override;
	core::Size count_dynamic_memory() const override;

private:

	//There are going to be a lot of nodes so I am only using 32 bits here. If you are reading this in a dystopian future where a Pose can have more than 4294967295 molten residues or a single residue position could have more than 4294967295 rotamers, please update these.
	unsigned int mres_id_;
	unsigned int rotamer_id_;

	utility::vector1< unsigned int > ids_of_clashing_nodes_;
};

///@brief Each HBondEdge represents a hydrogen bond
class HBondEdge : public utility::graph::Edge {

public:

	//constructor
	HBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind );

	HBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind, core::Real energy );

	//destructor
	~HBondEdge();

public:

	void copy_from( utility::graph::Edge const * source ) override;

	core::Size count_static_memory() const override;
	core::Size count_dynamic_memory() const override;

	///@brief this is intended to be the raw energy from the interaction graph between the rotamers represented by this->get_first_node_ind() and this->get_second_node_ind()
	inline float energy() const{
		return energy_;
	}

	inline void set_energy( core::Real energy ){
		energy_ = energy;
	}

	///@brief redundant interface for energy getter and setter. I find myself forgetting if the method is called "score" or "energy" so this way both are right
	inline float score() const{
		return energy_;
	}

	inline void set_score( core::Real energy ){
		energy_ = energy;
	}


private:

	float energy_;

};


class HBondGraph : public utility::graph::Graph {

public:

	//constructor
	HBondGraph();
	HBondGraph( core::Size num_nodes );

	//destructor
	~HBondGraph() override;

protected:

	utility::graph::Node * create_new_node( platform::Size node_index ) override;

	utility::graph::Edge * create_new_edge( core::Size index1, core::Size index2 ) override;

	utility::graph::Edge * create_new_edge( utility::graph::Edge const * example_edge ) override;

public:
	void delete_edge( utility::graph::Edge * edge ) override;

	core::Size count_static_memory() const override;
	core::Size count_dynamic_memory() const override;

private:

	boost::unordered_object_pool< HBondEdge > * hbond_edge_pool_;

};

} //hbnet
} //protocols

#endif
