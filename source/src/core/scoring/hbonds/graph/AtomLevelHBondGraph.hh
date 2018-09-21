// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/graph/AtomLevelHBondGraph.hh
/// @brief class headers for AtomLevelHBondGraph, AtomLevelHBondNode, and AtomLevelHBondEdge
/// @details This graph is different from a HBondGraph because it stores atom information for hbonds and residues.
/// See core/pack/hbonds/HBondGraph_util.hh (especially create_init_and_create_edges_for_atom_level_hbond_graph() ) and HBondGraph.hh for instructions on how to best use this data structure. Here is the inteded way to use an HBondGraph:
/// (1) Call ctor
/// (2) Immediately call core::pack::hbonds::init_node_info()
/// (3) Use MCHBNetInteractionGraph and RotamerSets::compute_energies() to populate edges into this graph
/// (4) If you are using an AtomLevelHBondGraph, call core::pack::hbonds::determine_atom_level_edge_info_for_all_edges() and core::pack::hbonds::determine_atom_level_node_info_for_all_nodes()
/// (5) Optional: If you are using an AtomLevelHBondGraph and you only care to analyze unsatisfied atoms, call find_satisfying_interactions_with_background()
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_scoring_hbonds_graph_AtomLevelHBondGraph_HH
#define INCLUDED_core_scoring_hbonds_graph_AtomLevelHBondGraph_HH

#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.fwd.hh>
#include <core/scoring/hbonds/graph/AtomInfo.hh>
#include <core/scoring/hbonds/graph/HBondInfo.hh>

#include <numeric/xyzVector.hh>

#include <utility/graph/unordered_object_pool.hpp>
#include <utility/graph/LowMemGraph.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>

#include <set>
#include <algorithm>

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {


///@brief Each AtomLevelHBondNode represents a rotamer from the RotamerSets object
class AtomLevelHBondNode : public utility::graph::LowMemNode {

public:

	//Please do not use these. We need this to exist only so that we can reserve space in the vector< AtomLevelHBondNode >
	AtomLevelHBondNode();
	AtomLevelHBondNode( const AtomLevelHBondNode& );
	/////////////////////////////////

	//constructor
	AtomLevelHBondNode( Size node_id );

	AtomLevelHBondNode( Size node_id, Size mres_id, Size rotamer_id );

	//destructor
	~AtomLevelHBondNode();

public:

	// void copy_from( utility::graph::Node const * source ) override;

	void print() const /* override */;

	Size count_static_memory() const /*override*/;
	Size count_dynamic_memory() const /*override*/;

public://getters, setters, accessors

	///@brief get molten residue id for this rotamer
	unsigned int moltenres() const{
		return mres_id_;
	}

	///@brief set molten residue id for this rotamer.
	void set_moltenres( Size mres ){
		mres_id_ = ( unsigned int ) ( mres );
	}

	///@brief get local rotamer id (local to the residue position)
	///@details this is equivalent to pack_basic::RotamerSetsBase::rotid_on_moltenresidue( this->get_node_index() )
	unsigned int local_rotamer_id() const{
		return rotamer_id_;
	}

	///@brief set local rotamer id (local to the residue position).
	void set_local_rotamer_id( Size rot_id ){
		rotamer_id_ = ( unsigned int ) ( rot_id );
	}

	///@brief duplicate interface for getting the global rotamer id. Identical to this->get_node_index()
	///details this is equivalent to pack_basic::RotamerSetsBase::nrotamer_offset_for_moltenres( this->moltenres() ) + this->local_rotamer_id()
	Size global_rotamer_id() const{
		return get_node_index();
	}

	///@brief keep track of another node (rotamer) that this clashes with. You do not need to call this for all of the rotamers that share a residue position.
	void register_clash( Size node_id ){
		ids_of_clashing_nodes_.insert(
			std::upper_bound( ids_of_clashing_nodes_.begin(), ids_of_clashing_nodes_.end(), node_id ),
			node_id
		);
	}

	///@brief does this node clash with another node (at another residue position)?
	bool clashes( Size node_id ) const{
		return std::binary_search( ids_of_clashing_nodes_.begin(), ids_of_clashing_nodes_.end(), node_id );
	}

	utility::vector1< AtomInfo > & polar_sc_atoms_not_satisfied_by_background() {
		return polar_sc_atoms_not_satisfied_by_background_;
	}

	utility::vector1< AtomInfo > const & polar_sc_atoms_not_satisfied_by_background() const {
		return polar_sc_atoms_not_satisfied_by_background_;
	}

	void add_polar_atom(
		unsigned short int local_atom_id,
		numeric::xyzVector< float > const & atom_position,
		bool is_hydrogen,
		bool is_donor,
		bool is_acceptor,
		bool is_hydroxyl,
		bool is_backbone
	){
		polar_sc_atoms_not_satisfied_by_background_.emplace_back( local_atom_id, atom_position, is_hydrogen, is_donor, is_acceptor, is_hydroxyl, is_backbone );
	}

	/// @brief The polar atoms are sorted, so we have to cleverly insert and maintain the order
	/// @details Used by merge
	void add_polar_atom_if_doesnt_exist( AtomInfo const & info ) {
		auto iter = std::lower_bound( polar_sc_atoms_not_satisfied_by_background_.begin(), polar_sc_atoms_not_satisfied_by_background_.end(), info );
		if ( iter->local_atom_id() == info.local_atom_id() ) return;
		polar_sc_atoms_not_satisfied_by_background_.insert( iter, info );
	}

	//Unstable means that the order is not preserved. I think this should not be used because it is handy to have all of the hydrogens at the end.
	//for example, seeing if there are any unsatisfied heavy atoms is an O(1) lookup of ( ! atom_vec.size() || atom_vec.front().is_hydrogen() )
	//These vectors are so short that I imagine we can pay the cost of down-shifting the higher elements on every removal

	/*static bool remove_atom_info_from_vec_unstable( utility::vector1< AtomInfo > & atom_vec, unsigned short int local_atom_id ){
	unsigned int const size = atom_vec.size();
	for( unsigned int ii = 1; ii <= size; ++ii ){
	if( atom_vec[ ii ].local_atom_id() == local_atom_id ){
	std::swap( atom_vec[ ii ], atom_vec[ size ] );
	debug_assert( atom_vec[ size ].local_atom_id() == local_atom_id );
	atom_vec.pop_back();
	return true;
	}
	}
	return false;
	}*/

	///@brief returns true if the atom was present, false is absent. False does not indicate failure!
	static bool remove_atom_info_from_vec_stable( utility::vector1< AtomInfo > & atom_vec, unsigned short int local_atom_id ){
		unsigned int const size = atom_vec.size();
		for ( unsigned int ii = 0; ii < size; ++ii ) {//THIS IS 0 ON PURPOSE
			if ( atom_vec[ ii + 1 ].local_atom_id() == local_atom_id ) {
				atom_vec.erase( atom_vec.begin() + ii );
				return true;
			}
		}
		return false;
	}

	void remove_atom_info_stable( unsigned short int local_atom_id ){
		remove_atom_info_from_vec_stable( polar_sc_atoms_not_satisfied_by_background_, local_atom_id );
	}

	/// @brief Only merges new data added by this class. Does not merge edges!!!
	/// @detail There are two ways to merge the atoms.
	///         OR_logic - Take the union of the atoms at each node
	///         AND_logic - Take the intersection of the atoms at each node
	///
	///   other_node_to_my_node is a map of node indices from the other graph
	///         to node indices in this graph.
	///   other_node_to_my_node.size() == other_graph.num_nodes()
	void merge_data(
		AtomLevelHBondNode const & other,
		utility::vector1< Size > const & other_node_to_my_node,
		bool merge_with_OR_logic
	);

private:
	unsigned int mres_id_;
	unsigned int rotamer_id_;

	utility::vector1< unsigned int > ids_of_clashing_nodes_;

	utility::vector1< AtomInfo > polar_sc_atoms_not_satisfied_by_background_;

public://please do not use this. It is required for pyrosetta compilation
	AtomLevelHBondNode & operator = ( AtomLevelHBondNode const & src ) {
		mres_id_ = src.mres_id_;
		rotamer_id_ = src.rotamer_id_;
		ids_of_clashing_nodes_ = src.ids_of_clashing_nodes_;
		polar_sc_atoms_not_satisfied_by_background_ = src.polar_sc_atoms_not_satisfied_by_background_;
		runtime_assert( false );
		return *this;
	}

};


///@brief Each AtomLevelHBondEdge represents a hydrogen bond
class AtomLevelHBondEdge : public utility::graph::LowMemEdge {

public:

	//constructor
	AtomLevelHBondEdge( Size first_node_ind, Size second_node_ind );

	AtomLevelHBondEdge( Size first_node_ind, Size second_node_ind, Real energy );

	//destructor
	~AtomLevelHBondEdge();

public:

	// void copy_from( utility::graph::Edge const * source ) override;

	Size count_static_memory() const /*override*/;
	Size count_dynamic_memory() const /*override*/;

	void register_hbond(
		bool first_node_is_donor,
		unsigned short int local_atom_id_A,
		unsigned short int local_atom_id_D,
		unsigned short int local_atom_id_H
	){
		hbonds_.emplace_back( first_node_is_donor, local_atom_id_A, local_atom_id_D, local_atom_id_H );
	}

	///@brief this is intended to be the raw energy from the interaction graph between the rotamers represented by this->get_first_node_ind() and this->get_second_node_ind()
	float energy() const {
		return energy_;
	}

	void set_energy( Real energy ){
		energy_ = energy;
	}

	///@brief redundant interface for energy getter and setter. I find myself forgetting if the method is called "score" or "energy" so this way both are right
	float score() const{
		return energy_;
	}

	void set_score( Real energy ){
		energy_ = energy;
	}

public://getters and setters

	std::vector< HBondInfo > & hbonds() {
		return hbonds_;
	}

	std::vector< HBondInfo > const & hbonds() const {
		return hbonds_;
	}

	/// @brief Only merges new data added by this class. Flips first_node_is_donor if necessary
	/// @detail other_node_to_my_node is a map of node indices from the other graph
	///         to node indices in this graph.
	///   other_node_to_my_node.size() == other_graph.num_nodes()
	void merge_data( AtomLevelHBondEdge const & other, utility::vector1< Size > const & other_node_to_my_node );

private:
	float energy_;
	std::vector< HBondInfo > hbonds_; // sizeof(std::vector) == 24 ; sizeof(utility::vector1) == 32
	//  If memory is a bottleneck, it's because there are too many edges

public://please do not use this. It is required for pyrosetta compilation
	AtomLevelHBondEdge & operator = ( AtomLevelHBondEdge const & src ) {
		energy_ = src.energy_;
		hbonds_ = src.hbonds_;
		runtime_assert( false );
		return *this;
	}

};

class AtomLevelHBondGraph : public utility::graph::LowMemGraph<AtomLevelHBondNode,AtomLevelHBondEdge> {

public:
	typedef utility::graph::LowMemGraph<AtomLevelHBondNode,AtomLevelHBondEdge> PARENT;

	//constructor
	AtomLevelHBondGraph();
	AtomLevelHBondGraph( Size num_nodes );

	//destructor
	~AtomLevelHBondGraph() override;

	Size count_static_memory() const override;
	Size count_dynamic_memory() const override;


	/// @brief Merges all info from other into this graph. All nodes from other must
	///        be present in this graph (by rotamer_id and mres_id).
	/// @detail There are two ways to merge the atoms of a node.
	///         OR_logic - Take the union of the atoms at each node
	///         AND_logic - Take the intersection of the atoms at each node
	void merge( AtomLevelHBondGraph const & other, bool merge_nodes_with_OR_logic );

public://old methods & methods used for unit tests
	AtomLevelHBondEdge * register_hbond( Size rotamerA, Size rotamerB, Real score );

};


} //graph
} //hbonds
} //scoring
} //core

#endif
