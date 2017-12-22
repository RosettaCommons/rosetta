// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/graph/HBondGraph.hh
/// @brief class headers for HBondGraph, HBondNode, and HBondEdge
/// @details This class is used to store and traverse data used for HBNet's monte carlo branching protocol. Most (if not all) of the information held in this graph is from HBNet's RotamerSets and InteractionGraph. Nodes in this graph represent rotamers from the rotamer set (node id == global rotamer id) and edges respresent hydrogen bonds.
/// @details See core/pack/hbonds/HBondGraph_util.hh for helper functions. Here is the inteded way to use an HBondGraph:
/// (1) Call ctor
/// (2) Immediately call core::pack::hbonds::init_node_info()
/// (3) Use MCHBNetInteractionGraph and RotamerSets::compute_energies() to populate edges into this graph
/// (4) If you are using an AtomLevelHBondGraph, call core::pack::hbonds::determine_atom_level_edge_info_for_all_edges() and core::pack::hbonds::determine_atom_level_node_info_for_all_nodes()
/// (5) Optional: If you are using an AtomLevelHBondGraph and you only care to analyze unsatisfied atoms, call find_satisfying_interactions_with_background()
/// @author Jack Maguire, jack@med.unc.edu

#ifndef INCLUDED_core_scoring_hbonds_graph_HBondGraph_HH
#define INCLUDED_core_scoring_hbonds_graph_HBondGraph_HH

#include <core/scoring/hbonds/graph/HBondGraph.fwd.hh>
#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/graph/unordered_object_pool.hpp>
#include <utility/graph/Graph.hh>

#include <boost/pool/pool.hpp>
#include <set>

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

///@brief Each HBondNode represents a rotamer from the RotamerSets object
class HBondNode : public utility::graph::Node {

public:
	//Please do not use these. We need this to exist only so that we can reserve space in the vector< HBondNode >
	HBondNode();
	HBondNode( const HBondNode& );
	///////////////////////////////////////

	//constructor
	HBondNode( utility::graph::Graph*, core::Size node_id );

	HBondNode( utility::graph::Graph*, core::Size node_id, core::Size mres_id, core::Size rotamer_id );

	//destructor
	~HBondNode() override;

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

	unsigned int mres_id_;
	unsigned int rotamer_id_;

	utility::vector1< unsigned int > ids_of_clashing_nodes_;

public://please do not use this. It is required for pyrosetta compilation
	HBondNode & operator = ( HBondNode const & src ){
		mres_id_ = src.mres_id_;
		rotamer_id_ = src.rotamer_id_;
		ids_of_clashing_nodes_ = src.ids_of_clashing_nodes_;
		return *this;
	}
};

///@brief Each HBondEdge represents a hydrogen bond
class HBondEdge : public utility::graph::Edge {

public:

	//constructor
	HBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind );

	HBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind, core::Real energy );

	//destructor
	~HBondEdge() override;

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

public://please do not use this. It is required for pyrosetta compilation
	HBondEdge & operator = ( HBondEdge const & src ){
		energy_ = src.energy_;
		return *this;
	}

};

class AbstractHBondGraph : public utility::graph::Graph {
public:
	//constructor
	//AbstractHBondGraph();
	//AbstractHBondGraph( core::Size num_nodes );

	//destructor
	//~AbstractHBondGraph() override;

	virtual HBondNode const * get_hbondnode( platform::Size index ) const = 0;
	virtual HBondNode * get_hbondnode( platform::Size index ) = 0;

	virtual HBondEdge * find_hbondedge( platform::Size node1, platform::Size node2 ) = 0;
	virtual HBondEdge const * find_hbondedge( platform::Size node1, platform::Size node2 ) const = 0;

	inline HBondEdge * register_hbond( core::Size rotamerA, core::Size rotamerB, core::Real score ){
		HBondEdge * new_edge = static_cast< HBondEdge * >( add_edge( rotamerA, rotamerB ) );
		new_edge->set_energy( score );
		return new_edge;
	}

protected:
	///@brief This is only called by Graph::delete_everything() as it deletes its ptrs to nodes.
	/// Since we will be deleting all nodes upon our destruction, we do not need to do anything here.
	/// The absence of the override results in a double-free runtime error.
	void delete_node( utility::graph::Node * ) override {}

};

class HBondGraph : public AbstractHBondGraph {

public:

	//constructor
	HBondGraph();
	HBondGraph( core::Size num_nodes );

	//destructor
	~HBondGraph() override;

	void set_num_nodes( platform::Size num_nodes ) override;

protected:

	utility::graph::Node * create_new_node( platform::Size node_index ) override;

	utility::graph::Edge * create_new_edge( core::Size index1, core::Size index2 ) override;
	utility::graph::Edge * create_new_edge( utility::graph::Edge const * example_edge ) override;

public:
	void delete_edge( utility::graph::Edge * edge ) override;

	core::Size count_static_memory() const override;
	core::Size count_dynamic_memory() const override;

public: //inline access methods
	HBondNode const * get_hbondnode( platform::Size index ) const override
	{
		return &all_nodes_[ index ];
	}

	HBondNode * get_hbondnode( platform::Size index ) override
	{
		return &all_nodes_[ index ];
	}

	HBondEdge * find_hbondedge( platform::Size node1, platform::Size node2 ) override
	{
		return static_cast< HBondEdge * >( find_edge( node1, node2 ) );
	}

	HBondEdge const * find_hbondedge( platform::Size node1, platform::Size node2 ) const override
	{
		return static_cast< HBondEdge const * >( find_edge( node1, node2 ) );
	}


private:

	boost::unordered_object_pool< HBondEdge > * hbond_edge_pool_;
	utility::vector1< HBondNode > all_nodes_;

};

} //graph
} //hbonds
} //scoring
} //core

#endif
