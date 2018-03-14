// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/hbnet/NetworkState.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com

//Headers
#include <protocols/hbnet/NetworkState.hh>

namespace protocols {
namespace hbnet {

NetworkState::NetworkState(
	core::scoring::hbonds::graph::AtomLevelHBondEdge const * monte_carlo_seed,
	core::scoring::hbonds::graph::AtomLevelHBondGraphOP const & hbond_graph
) :
	nodes_ ( 0 ),
	edges_ ( 0 ),
	monte_carlo_seed_( monte_carlo_seed ),
	full_twobody_energy_( monte_carlo_seed->energy() ),
	score_ ( 0 ),
	unsatisfied_sc_atoms_( 0 ),
	sorter_()
{
	using namespace core::scoring::hbonds::graph;

	edges_.push_back( monte_carlo_seed_ );

	core::scoring::hbonds::graph::AtomLevelHBondNode const * const first_node =
		hbond_graph->get_hbondnode( monte_carlo_seed_->get_first_node_ind() );
	nodes_.push_back( first_node );
	add_polar_atoms( first_node );

	core::scoring::hbonds::graph::AtomLevelHBondNode const * const second_node =
		hbond_graph->get_hbondnode( monte_carlo_seed_->get_second_node_ind() );
	//in the case of symmetry, this might evaluate to false:
	if ( monte_carlo_seed_->get_first_node_ind() != monte_carlo_seed_->get_second_node_ind() ) {
		nodes_.push_back( second_node );
		add_polar_atoms( second_node );
	}

	unsigned int const mres1 = first_node->moltenres();
	unsigned int const mres2 = second_node->moltenres();

	//find satisfied atoms
	for ( HBondInfo const & hbond : monte_carlo_seed->hbonds() ) {
		bool const first_node_is_donor = hbond.first_node_is_donor();
		unsigned short int const local_atom_id_A = hbond.local_atom_id_A();//Acceptor
		unsigned short int const local_atom_id_D = hbond.local_atom_id_D();//Donor
		unsigned short int const local_atom_id_H = hbond.local_atom_id_H();//Hydrogen

		//Remove hbonding atoms from vector of unsats:
		if ( first_node_is_donor ) {
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

void NetworkState::add_polar_atoms(
	core::scoring::hbonds::graph::AtomLevelHBondNode const * node
){
	auto const & atoms_to_copy = node->polar_sc_atoms_not_satisfied_by_background();
	core::Size const mres = node->moltenres();

	auto iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
	if ( iter == unsatisfied_sc_atoms_.end() || iter->first != mres ) {
		//element does not exist yet, so inset it
		iter = unsatisfied_sc_atoms_.insert(
			iter, std::make_pair( mres, atoms_to_copy )
		);
	} else {
		iter->second = atoms_to_copy;
	}

}


} //hbnet
} //protocols
