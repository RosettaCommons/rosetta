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
	core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed,
	core::scoring::hbonds::graph::HBondGraphOP const & hbond_graph
) :
	nodes_ ( 0 ),
	edges_ ( 0 ),
	monte_carlo_seed_( monte_carlo_seed ),
	full_twobody_energy_( monte_carlo_seed->energy() ),
	score_ ( 0 ),
	unsatisfied_sc_atoms_()
{
	using namespace core::scoring::hbonds::graph;

	edges_.push_back( monte_carlo_seed_ );

	core::scoring::hbonds::graph::HBondNode const * const first_node =
		hbond_graph->get_node( monte_carlo_seed_->get_first_node_ind() );
	nodes_.push_back( first_node );
	add_polar_atoms( first_node );

	core::scoring::hbonds::graph::HBondNode const * const second_node =
		hbond_graph->get_node( monte_carlo_seed_->get_second_node_ind() );
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
			get_unsats_for_mres( mres2 )->remove( local_atom_id_A );

			core::scoring::hbonds::graph::AtomInfoSet * const donor_unsats =
				get_unsats_for_mres( mres1 );
			donor_unsats->remove( local_atom_id_D );
			donor_unsats->remove( local_atom_id_H );
		} else {
			get_unsats_for_mres( mres1 )->remove( local_atom_id_A );

			core::scoring::hbonds::graph::AtomInfoSet * const donor_unsats =
				get_unsats_for_mres( mres2 );
			donor_unsats->remove( local_atom_id_D );
			donor_unsats->remove( local_atom_id_H );
		}
	}
}

void NetworkState::add_polar_atoms(
	core::scoring::hbonds::graph::HBondNode const * node
){
	using namespace core::scoring::hbonds::graph;

	AtomInfoSet const & atoms_to_copy = node->polar_sc_atoms_not_satisfied_by_background();
	core::Size const mres = node->moltenres();

	AtomInfoSet * const unsat_vec = get_unsats_for_mres( mres );

	if ( unsat_vec == nullptr ) {
		//element does not exist yet, so inset it
		unsatisfied_sc_atoms_[ mres ] = atoms_to_copy;
	} else {
		(* unsat_vec) = atoms_to_copy;
	}

}


} //hbnet
} //protocols
