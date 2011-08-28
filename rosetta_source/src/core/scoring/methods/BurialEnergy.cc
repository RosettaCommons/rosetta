// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/BurialEnergy.cc
/// @brief  Radius of gyration energy function definition.
/// @author James Thompson

// Unit headers
#include <core/scoring/methods/BurialEnergy.hh>
#include <core/scoring/methods/BurialEnergyCreator.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/prof.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>

namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the BurialEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
BurialEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new BurialEnergy;
}

ScoreTypes
BurialEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rg );
	return sts;
}


/// c-tor
BurialEnergy::BurialEnergy() :
	parent( new BurialEnergyCreator )
{}


/// clone
EnergyMethodOP
BurialEnergy::clone() const
{
	return new BurialEnergy;
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
BurialEnergy::finalize_total_energy(
	core::pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( rsd.has_variant_type( "REPLONLY" ) ) continue;
		if ( ! rsd.is_protein() ) continue;

		core::conformation::Residue const & rsd;

		TwelveANeighborGraph const & graph ( pose.energies().twelveA_neighbor_graph() );
		Size const atomindex_i = rsd.atom_index( representative_atom_name( rsd.aa() ));

		core::conformation::Atom const & atom_i = rsd.atom(atomindex_i);

		// iterate across neighbors within 12 angstroms, count number <10A
		Size countN(0);
		for ( graph::Graph::EdgeListConstIter
				ir  = graph.get_node( rsd.seqpos() )->const_edge_list_begin(),
				ire = graph.get_node( rsd.seqpos() )->const_edge_list_end();
				ir != ire; ++ir ) {
			Size const j( (*ir)->get_other_ind( rsd.seqpos() ) );
			core::conformation::Residue const & rsd_j( pose.residue(j) );

			Size atomindex_j( rsd_j.nbr_atom() );
			core::conformation::Atom const & atom_j = rsd_j.atom(atomindex_j);
			Real sqdist = atom_i.xyz().distance_squared( atom_j.xyz() );
			if ( sqdist < 10 ) {
				countN++;
			}
		}

		// burial_prediction is negative for exposed predictions and positive for buried
		// magnitude of burial_prediction is related to prediction confidence
		Real const & burial_prediction( pred_burial_[ii] );
		Real const score( -1 * burial_prediction * countN );
		emap[ burial ] += score;
	}
}

core::Size
BurialEnergy::version() const
{
	return 1;
}

void
BurialEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ twelve_A_neighbor_graph ] = true;
}

void
BurialEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	pose.update_residue_neighbors();
}

void BurialEnergy::init_from_file() {
	utility::io::izstream input;
}

} // methods
} // scoring
} // core
