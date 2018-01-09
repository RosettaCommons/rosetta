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
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <numeric/xyzVector.hh>

namespace protocols {
namespace hbnet {

using mres_unsat_pair=std::pair< unsigned int /*mres*/, utility::vector1< core::scoring::hbonds::graph::AtomInfo > >;

struct polar_atom {
	polar_atom(
		unsigned int sequence_position,
		numeric::xyzVector< float > const & atom_position,
		bool hydroxyl,
		bool satisfied = false
	){
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

struct compare_mres_unsat_pair : public std::binary_function< mres_unsat_pair, mres_unsat_pair, bool >{
	bool operator()( mres_unsat_pair const & a, mres_unsat_pair const & b ) const {
		return a.first < b.first;
	}

	bool operator()( mres_unsat_pair const & a, unsigned int b ) const {
		return a.first < b;
	}

	bool operator()( unsigned int a, mres_unsat_pair const & b ) const {
		return a < b.first;
	}

};


struct NetworkState{

	NetworkState(
		core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed_in,
		core::scoring::hbonds::graph::AbstractHBondGraphOP const & hbond_graph,
		core::pack::rotamer_set::RotamerSetsOP const & rotsets
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
				//AtomLevelHBondNode::remove_atom_info_from_vec_stable( unsatisfied_sc_atoms[ mres2 ], local_atom_id_A );
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( get_unsats_for_mres( mres2 )->second, local_atom_id_A );

				utility::vector1< AtomInfo > & donor_atom_vec = get_unsats_for_mres( mres1 )->second;
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_D );
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( donor_atom_vec, local_atom_id_H );
			} else {
				AtomLevelHBondNode::remove_atom_info_from_vec_stable( get_unsats_for_mres( mres1 )->second, local_atom_id_A );

				utility::vector1< AtomInfo > & donor_atom_vec = get_unsats_for_mres( mres2 )->second;
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
		core::pack::rotamer_set::RotamerSetsOP const & rotsets
	){
		core::Size const mres = rotsets->moltenres_for_rotamer( node->get_node_index() );
		auto iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
		auto const & atoms_to_copy = static_cast< core::scoring::hbonds::graph::AtomLevelHBondNode const * >( node )->polar_sc_atoms_not_satisfied_by_background();

		if ( iter == unsatisfied_sc_atoms_.end() || iter->first != mres ) {
			iter = unsatisfied_sc_atoms_.insert(
				iter,
				std::make_pair( mres, atoms_to_copy )
			);
		} else {
			iter->second = atoms_to_copy;
		}

	}

	utility::vector1< mres_unsat_pair >::iterator get_unsats_for_mres( unsigned int mres ) {
		auto iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
		runtime_assert( iter != unsatisfied_sc_atoms_.end() );
		runtime_assert( iter->first == mres );
		return iter;
	}

	utility::vector1< mres_unsat_pair >::const_iterator get_unsats_for_mres( unsigned int mres ) const {
		auto const iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
		runtime_assert( iter != unsatisfied_sc_atoms_.end() );
		runtime_assert( iter->first == mres );
		return iter;
	}

	inline bool mres_has_unsats( unsigned int mres ) const {
		auto iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
		if ( iter == unsatisfied_sc_atoms_.end() ) return false;
		if ( iter->first != mres ) return false;
		return ! iter->second.empty();
	}

	inline bool mres_has_heavy_unsats ( unsigned int mres ) const {
		auto iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
		if ( iter == unsatisfied_sc_atoms_.end() ) return false;
		if ( iter->first != mres ) return false;
		if ( iter->second.empty() ) return false;
		return iter->second[ 1 ].is_hydrogen();
	}

	inline utility::vector1< mres_unsat_pair > const & unsatisfied_sc_atoms_const() const {
		return unsatisfied_sc_atoms_;
	}

	utility::vector1< core::scoring::hbonds::graph::HBondNode const * > nodes;
	utility::vector1< core::scoring::hbonds::graph::HBondEdge const * > edges;

	core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed;//"Seed" hbond to branch off of
	core::Real full_twobody_energy;//Sum of hbond score + clash score for all residue pairs in "residues" data object
	core::Real score;//This holds whatever metric is used for sorting

private:
	utility::vector1< mres_unsat_pair > unsatisfied_sc_atoms_;
	compare_mres_unsat_pair sorter_;

	/*utility::vector1< mres_unsat_pair >::iterator get_unsats_for_mres( unsigned int mres ) {
	auto const iter = std::lower_bound( unsatisfied_sc_atoms.begin(), unsatisfied_sc_atoms.end(), mres );
	}*/

};

struct NetworkStateScoreComparator{
	static bool compare( NetworkState const & a, NetworkState const & b ){
		return a.score < b.score;
	}
};

}//hbnet
}//protocols
#endif
