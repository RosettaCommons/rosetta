// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/HybridVDW_Energy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/HybridVDW_Energy.hh>
#include <core/scoring/methods/HybridVDW_EnergyCreator.hh>

// Package headers
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the HybridVDW_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
HybridVDW_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new HybridVDW_Energy );
}

ScoreTypes
HybridVDW_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( hybrid_vdw );
	return sts;
}


Real const vdw_scale_factor( 0.8 );

/// @details  C-TOR
HybridVDW_Energy::HybridVDW_Energy() :
	parent( methods::EnergyMethodCreatorOP( new HybridVDW_EnergyCreator ) ),
	atom_vdw_( ScoringManager::get_instance()->get_AtomVDW( chemical::HYBRID_FA_STANDARD_CENTROID ) )
{}


/// clone
EnergyMethodOP
HybridVDW_Energy::clone() const
{
	return EnergyMethodOP( new HybridVDW_Energy( *this ) );
}

/// @details  copy c-tor
HybridVDW_Energy::HybridVDW_Energy( HybridVDW_Energy const & /*src*/ ) = default;


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
HybridVDW_Energy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	// is this really necessary?
	pose.update_residue_neighbors();
}


void
HybridVDW_Energy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	// is this really necessary?
	pose.update_residue_neighbors();
}


void
HybridVDW_Energy::residue_pair_energy(
	conformation::Residue const & rsd1, //_in,
	conformation::Residue const & rsd2, //_in,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using namespace etable::count_pair;


	if ( ( rsd1.is_protein() && ! rsd2.is_protein() ) ||
			( rsd2.is_protein() && ! rsd1.is_protein() ) ) {
		//OK
	} else {
		return;
	}

	Real score(0.0);
	debug_assert( ! rsd1.is_bonded( rsd2 ) );

	// no countpair!
	for ( Size i = 1, i_end = rsd1.natoms(); i <= i_end; ++i ) {
		Vector const & i_xyz( rsd1.xyz(i) );
		Size const i_type( rsd1.atom_type_index(i) );
		Real const i_radius( atom_vdw_.approximate_vdw_radius( i_type ) );
		for ( Size j = 1, j_end = rsd2.natoms(); j <= j_end; ++j ) {
			Real const bump_dsq( numeric::square( i_radius + atom_vdw_.approximate_vdw_radius( rsd2.atom_type_index( j ) ) ) );
			Real const clash( bump_dsq - i_xyz.distance_squared( rsd2.xyz(j) ) );
			if ( clash > 0.0 ) {
				score += ( clash * clash ) / bump_dsq;
				//std::cout << "BUMP: " << I(4,rsd1.seqpos() ) << I(4,rsd2.seqpos() ) <<
				// ' ' << rsd1.atom_name(i) << ' ' << rsd2.atom_name(j) << ' ' << ( clash * clash ) / bump_dsq << std::endl;
			}
		}
	}
	debug_assert( rsd1.type().mode() == chemical::HYBRID_FA_STANDARD_CENTROID_t );
	debug_assert( rsd2.type().mode() == chemical::HYBRID_FA_STANDARD_CENTROID_t );
	emap[ hybrid_vdw ] += score * vdw_scale_factor;

}


void
HybridVDW_Energy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using namespace etable::count_pair;

	// what is my charge?
	Size const pos1( atom_id.rsd() );
	Size const i   ( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( pos1 ) );

	int const pos1_map( domain_map( pos1 ) );
	bool const pos1_fixed( pos1_map != 0 );

	Vector const & i_xyz( rsd1.xyz(i) );
	Size const i_type( rsd1.atom_type_index(i) );
	Real const i_radius( atom_vdw_.approximate_vdw_radius( i_type ) );

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// loop over *all* nbrs of rsd1 (not just upper or lower)
	for ( utility::graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			irue = energy_graph.get_node( pos1 )->const_edge_list_end();
			iru != irue; ++iru ) {
		Size const pos2( (*iru)->get_other_ind( pos1 ) );

		if ( pos1_fixed && pos1_map == domain_map( pos2 ) ) continue; // fixed wrt one another

		conformation::Residue const & rsd2( pose.residue( pos2 ) );

		if ( ! ( ( rsd1.is_protein() && ! rsd2.is_protein() ) ||
				( rsd2.is_protein() && ! rsd1.is_protein() ) ) ) continue;
		debug_assert( pos2 != pos1 );
		debug_assert( ! rsd1.is_bonded( rsd2 ) );
		// no countpair!
		for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {

			Vector const & j_xyz( rsd2.xyz(j) );
			Vector const f2( i_xyz - j_xyz );
			Real const dis2( f2.length_squared() );
			//Real const bump_dsq( i_atom_vdw[ rsd2.atom_type_index(j) ] );

			Real const bump_dsq( numeric::square( i_radius + atom_vdw_.approximate_vdw_radius( rsd2.atom_type_index( j ) ) ) );

			if ( dis2 < bump_dsq ) {
				// E += vdw_scale_factor_ * weights[vdw] * ( ( dis2 - bump_dsq ) **2 ) / bump_dsq
				Real const dE_dr_over_r = vdw_scale_factor * weights[ hybrid_vdw ] * 4.0 * ( dis2 - bump_dsq ) / bump_dsq;
				Vector const f1( i_xyz.cross( j_xyz ) );
				F1 += dE_dr_over_r * f1;
				F2 += dE_dr_over_r * f2;
			}
		}

	} // loop over nbrs of rsd1

}


/// @brief HybridVDW_Energy distance cutoff
Distance
HybridVDW_Energy::atomic_interaction_cutoff() const
{
	return 6.0; /// now subtracted off 3.0 from cutoffs in centroid params files
	//return 0.0; /// since all the cutoffs for centroid mode are rolled into the cendist check
}

/// @brief HybridVDW_Energy
void
HybridVDW_Energy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}
core::Size
HybridVDW_Energy::version() const
{
	return 1; // Initial versioning
}


}
}
}
