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
/// See core/pack/hbonds/HBondGraph_util.hh (especially create_init_and_create_edges_for_atom_level_hbond_graph() ) and HBondGraph.hh for instructions on how to best use this data structure
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_scoring_hbonds_graph_AtomLevelHBondGraph_HH
#define INCLUDED_core_scoring_hbonds_graph_AtomLevelHBondGraph_HH

#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.fwd.hh>
#include <core/scoring/hbonds/graph/HBondGraph.hh>
#include <core/scoring/hbonds/graph/AtomInfo.hh>
#include <core/scoring/hbonds/graph/HBondInfo.hh>

#include <numeric/xyzVector.hh>

#include <utility/graph/unordered_object_pool.hpp>
#include <utility/graph/Graph.hh>
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
class AtomLevelHBondNode : public HBondNode {

public:

	//Please do not use these. We need this to exist only so that we can reserve space in the vector< AtomLevelHBondNode >
	AtomLevelHBondNode();
	AtomLevelHBondNode( const AtomLevelHBondNode& );
	/////////////////////////////////

	//constructor
	AtomLevelHBondNode( utility::graph::Graph*, core::Size node_id );

	AtomLevelHBondNode( utility::graph::Graph*, core::Size node_id, core::Size mres_id, core::Size rotamer_id );

	//destructor
	~AtomLevelHBondNode() override;

public:

	void copy_from( utility::graph::Node const * source ) override;

	//void print() const override;

	core::Size count_static_memory() const override;
	core::Size count_dynamic_memory() const override;

public://getters, setters, accessors

	inline utility::vector1< AtomInfo > & polar_sc_atoms_not_satisfied_by_background() {
		return polar_sc_atoms_not_satisfied_by_background_;
	}

	inline utility::vector1< AtomInfo > const & polar_sc_atoms_not_satisfied_by_background() const {
		return polar_sc_atoms_not_satisfied_by_background_;
	}

	inline void add_polar_atom(
		unsigned short int local_atom_id,
		numeric::xyzVector< float > const & atom_position,
		bool is_hydrogen,
		bool is_donor,
		bool is_acceptor,
		bool is_hydroxyl
	){
		polar_sc_atoms_not_satisfied_by_background_.emplace_back( local_atom_id, atom_position, is_hydrogen, is_donor, is_acceptor, is_hydroxyl );
	}

	//Unstable means that the order is not preserved. I think this should not be used because it is handy to have all of the hydrogens at the end.
	//for example, seeing if there are any unsatisfied heavy atoms is an O(1) lookup of ( ! atom_vec.size() || atom_vec.front().is_hydrogen() )
	//These vectors are so short that I imagine we can pay the cost of down-shifting the higher elements on every removal

	/*inline static bool remove_atom_info_from_vec_unstable( utility::vector1< AtomInfo > & atom_vec, unsigned short int local_atom_id ){
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
	inline static bool remove_atom_info_from_vec_stable( utility::vector1< AtomInfo > & atom_vec, unsigned short int local_atom_id ){
		unsigned int const size = atom_vec.size();
		for ( unsigned int ii = 0; ii < size; ++ii ) {//THIS IS 0 ON PURPOSE
			if ( atom_vec[ ii + 1 ].local_atom_id() == local_atom_id ) {
				atom_vec.erase( atom_vec.begin() + ii );
				return true;
			}
		}
		return false;
	}

	inline void remove_atom_info_stable( unsigned short int local_atom_id ){
		remove_atom_info_from_vec_stable( polar_sc_atoms_not_satisfied_by_background_, local_atom_id );
	}

	/*virtual void remove_atom_info_unstable( unsigned short int local_atom_id ){
	remove_atom_info_unstable( polar_sc_atoms_not_satisfied_by_background_, local_atom_id );
	}*/

private:
	//utility::vector1< AtomInfo > all_polar_sc_atoms_;
	utility::vector1< AtomInfo > polar_sc_atoms_not_satisfied_by_background_;

public://please do not use this. It is required for pyrosetta compilation
	AtomLevelHBondNode & operator = ( AtomLevelHBondNode const & src ){
		HBondNode::operator = ( src );
		polar_sc_atoms_not_satisfied_by_background_ = src.polar_sc_atoms_not_satisfied_by_background_;
		return *this;
	}

};


///@brief Each AtomLevelHBondEdge represents a hydrogen bond
class AtomLevelHBondEdge : public HBondEdge {

public:

	//constructor
	AtomLevelHBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind );

	AtomLevelHBondEdge( utility::graph::Graph* owner, core::Size first_node_ind, core::Size second_node_ind, core::Real energy );

	//destructor
	~AtomLevelHBondEdge();

public:

	void copy_from( utility::graph::Edge const * source ) override;

	core::Size count_static_memory() const override;
	core::Size count_dynamic_memory() const override;

	inline void register_hbond( bool first_node_is_donor, unsigned short int local_atom_id_A, unsigned short int local_atom_id_D, unsigned short int local_atom_id_H ){
		hbonds_.emplace_back( first_node_is_donor, local_atom_id_A, local_atom_id_D, local_atom_id_H );
	}

	inline utility::vector1< HBondInfo > & hbonds() {
		return hbonds_;
	}

	inline utility::vector1< HBondInfo > const & hbonds() const {
		return hbonds_;
	}

private:
	utility::vector1< HBondInfo > hbonds_;

public://please do not use this. It is required for pyrosetta compilation
	AtomLevelHBondEdge & operator = ( AtomLevelHBondEdge const & src ){
		HBondEdge::operator = ( src );
		hbonds_ = src.hbonds_;
		return *this;
	}

};


///@brief AtomLevelHBondGraph does not derive directly from the HBondGraph because HBondGraph would create a ton of HBondNodes that we will never use.
class AtomLevelHBondGraph : public AbstractHBondGraph {

public:

	//constructor
	AtomLevelHBondGraph();
	AtomLevelHBondGraph( core::Size num_nodes );

	//destructor
	~AtomLevelHBondGraph() override;

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

	inline AtomLevelHBondNode const * get_atomlevel_hbondnode( platform::Size index ) const
	{
		return &all_nodes_[ index ];
	}

	inline AtomLevelHBondNode * get_atomlevel_hbondnode( platform::Size index )
	{
		return &all_nodes_[ index ];
	}

	inline AtomLevelHBondEdge * find_atomlevel_hbondedge( platform::Size node1, platform::Size node2 )
	{
		return static_cast< AtomLevelHBondEdge * >( find_edge( node1, node2 ) );
	}

	inline AtomLevelHBondEdge const * find_atomlevel_hbondedge( platform::Size node1, platform::Size node2 ) const
	{
		return static_cast< AtomLevelHBondEdge const * >( find_edge( node1, node2 ) );
	}

public: //virtual access methods
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
		return static_cast< AtomLevelHBondEdge * >( find_edge( node1, node2 ) );
	}

	HBondEdge const * find_hbondedge( platform::Size node1, platform::Size node2 ) const override
	{
		return static_cast< AtomLevelHBondEdge const * >( find_edge( node1, node2 ) );
	}

private:

	boost::unordered_object_pool< AtomLevelHBondEdge > * hbond_edge_pool_;
	utility::vector1< AtomLevelHBondNode > all_nodes_;

};


} //graph
} //hbonds
} //scoring
} //core

#endif
