// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hbnet/NetworkState.hh
/// @brief specialized structs for MC-HBNet to store hbond network info
/// @author Jack Maguire, jack@med.unc.edu

#ifndef INCLUDED_protocols_hbnet_NetworkState_hh
#define INCLUDED_protocols_hbnet_NetworkState_hh

#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/hbonds/graph/HBondGraph.hh>
#include <core/scoring/hbonds/graph/HBondInfo.hh>
#include <core/scoring/hbonds/graph/AtomInfo.hh>
#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>
//#include <core/scoring/hbonds/graph/LKAtomLevelHBondGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <numeric/xyzVector.hh>

namespace protocols {
namespace hbnet {

struct polar_atom {
	polar_atom( unsigned int sequence_position, numeric::xyzVector< float > const & atom_position, bool hydroxyl, bool satisfied=false /*, bool buried */ ){
		seqpos = sequence_position;
		xyz = atom_position;
		is_hydroxyl = hydroxyl; // worth it to store here, because so many special cases revolve around OH's
		is_satisfied = satisfied;
	}

	bool operator < ( polar_atom const & rhs ) const {
		//x is the primary metric for sorting. The y and z parts are just to make sure that two atoms with the same x value are still included
		if ( xyz.x() != rhs.xyz.x() ) return xyz.x() < rhs.xyz.x();
		if ( xyz.y() != rhs.xyz.y() ) return xyz.y() < rhs.xyz.y();
		return xyz.z() < rhs.xyz.z();
	}

	unsigned int seqpos;
	numeric::xyzVector< float > xyz;
	bool is_hydroxyl;
	bool is_satisfied;
};

struct NetworkState{

	NetworkState(
		core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed_in,
		core::scoring::hbonds::graph::AbstractHBondGraphOP hbond_graph,
		core::pack::rotamer_set::RotamerSetsOP rotsets
	) :
		monte_carlo_seed( monte_carlo_seed_in ),
		full_twobody_energy( monte_carlo_seed->energy() ),
		score ( 0 )
	{
		using namespace core::scoring::hbonds::graph;

		edges.push_back( monte_carlo_seed );

		nodes.clear();
		core::scoring::hbonds::graph::HBondNode const * const first_node = hbond_graph->get_hbondnode( monte_carlo_seed->get_first_node_ind() );
		nodes.push_back( first_node );
		add_polar_atoms( nodes.back(), rotsets );

		//in the case of symmetry, this might evaluate to false:
		if ( monte_carlo_seed->get_first_node_ind() != monte_carlo_seed->get_second_node_ind() ) {
			core::scoring::hbonds::graph::HBondNode const * const second_node = hbond_graph->get_hbondnode( monte_carlo_seed->get_second_node_ind() );
			nodes.push_back( second_node );
			add_polar_atoms( nodes.back(), rotsets );
		}

		//find satisfied atoms
		AtomLevelHBondEdge const * al_edge = static_cast< AtomLevelHBondEdge const * >( monte_carlo_seed );
		for ( HBondInfo const & hbond : al_edge->hbonds() ) {
			bool const first_node_is_donor = hbond.first_node_is_donor();
			unsigned short int const local_atom_id_A = hbond.local_atom_id_A();//Acceptor
			unsigned short int const local_atom_id_D = hbond.local_atom_id_D();//Donor
			unsigned short int const local_atom_id_H = hbond.local_atom_id_H();//Hydrogen

			//AtomLevelHBondNode const * al_node1 = static_cast< AtomLevelHBondEdge const * >( al_edge->get_other_node() )

			unsigned int const mres1 = first_node->moltenres();
			unsigned int const mres2 = hbond_graph->get_hbondnode( monte_carlo_seed->get_second_node_ind() )->moltenres();

			if ( first_node_is_donor ) {
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( unsatisfied_sc_atoms[ mres2 ], local_atom_id_A );

				utility::vector1< AtomInfo > & donor_atom_vec = unsatisfied_sc_atoms[ mres1 ];
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_D );
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_H );
			} else {
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( unsatisfied_sc_atoms[ mres1 ], local_atom_id_A );

				utility::vector1< AtomInfo > & donor_atom_vec = unsatisfied_sc_atoms[ mres2 ];
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_D );
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_H );
			}
		}
	}

	bool operator < ( NetworkState const & rhs) const {
		return full_twobody_energy < rhs.full_twobody_energy;
	}

	void add_polar_atoms(
		core::scoring::hbonds::graph::HBondNode const * node,
		core::pack::rotamer_set::RotamerSetsOP rotsets
	){
		core::Size const mres = rotsets->moltenres_for_rotamer( node->get_node_index() );
		if ( unsatisfied_sc_atoms.find( mres ) != unsatisfied_sc_atoms.end() ) return;
		unsatisfied_sc_atoms[ mres ] = static_cast< core::scoring::hbonds::graph::AtomLevelHBondNode const * >( node )->polar_sc_atoms_not_satisfied_by_background(); //copy ctor
	}


	utility::vector1< core::scoring::hbonds::graph::HBondNode const * > nodes;
	utility::vector1< core::scoring::hbonds::graph::HBondEdge const * > edges;

	core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed;//"Seed" hbond to branch off of
	core::Real full_twobody_energy;//Sum of hbond score + clash score for all residue pairs in "residues" data object
	core::Real score;//This holds whatever metric is used for sorting

	//This is a map of all unsatisfied atoms ( Heavy and H both included ) for every moltenres in the network:
	std::map < unsigned int /*mres*/, utility::vector1< core::scoring::hbonds::graph::AtomInfo > > unsatisfied_sc_atoms;

};

struct NetworkStateScoreComparator{
	static bool compare( NetworkState const & a, NetworkState const & b ){
		return a.score < b.score;
	}
};

}//hbnet
}//protocols
#endif
