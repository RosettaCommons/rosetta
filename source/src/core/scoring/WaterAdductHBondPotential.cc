// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/WaterAdductHBondPotential.cc
/// @brief
/// @author Jim Havranek

// Project headers
#include <core/scoring/WaterAdductHBondPotential.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>


// // Project headers
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoreType.hh>

// Numeric
#include <numeric/conversions.hh>

// Utility headers
#include <utility/exit.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


namespace core {
namespace scoring {

WaterAdductHBondPotential::WaterAdductHBondPotential() :
	hbondoptions_( scoring::hbonds::HBondOptionsOP( new hbonds::HBondOptions ) ),
	hb_database_( scoring::hbonds::HBondDatabase::get_database())
{
	hbondoptions_->use_sp2_chi_penalty( false ); // do not use the chi penalty -- override command line
}

WaterAdductHBondPotential::~WaterAdductHBondPotential() {}

Real
WaterAdductHBondPotential::water_adduct_hbond_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	// Take turns assuming the waters are on either residue
	return h2o_hbond_score_1way( rsd1, rsd2 ) + h2o_hbond_score_1way( rsd2, rsd1 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A rough translation based on Phil's implementation in rosetta++

Real
WaterAdductHBondPotential::h2o_hbond_score_1way(
	conformation::Residue const & h2o_rsd,
	conformation::Residue const & other_rsd
) const
{


	Real h2o_hbond_score( 0.0 );
	hbonds::HBondDerivs deriv;

	static Real const HO_dist( 0.96 );
	static Real const theta( numeric::conversions::radians( 107.0 ) );
	static Real const HO_x_dist( HO_dist * std::cos( theta ) );
	static Real const HO_y_dist( HO_dist * std::sin( theta ) );

	// Loop over waters
	for ( Size i = 1, ei = h2o_rsd.natoms() ; i <= ei ; ++i ) {
		if ( !h2o_rsd.atom_type(i).is_h2o() ) continue;

		Size const h2o_index( i );

		Vector const h2o_coord( h2o_rsd.xyz( i ) );
		Vector const h2o_base_coord( h2o_rsd.xyz( h2o_rsd.atom_base(i) ) );
		Vector const h2o_base2_coord( h2o_rsd.xyz( h2o_rsd.abase2(i) ) );
		Vector const normal_from_base( Vector( h2o_coord - h2o_base_coord ).normalized() );

		// First check for interactions with donors
		for ( Size j = 1, ej = other_rsd.Hpos_polar().size() ; j <= ej ; ++j ) {
			Size const hyd_index( other_rsd.Hpos_polar()[j] );
			Size const donor_index( other_rsd.atom_base( hyd_index ) );
			Vector const hyd_coord( other_rsd.xyz( hyd_index ) );
			Vector const donor_coord( other_rsd.xyz( donor_index ) );

			//  std::cout << "other rsd type " << other_rsd.type().name() << " hyd name " << other_rsd.atom_name( hyd_index ) << " donor name " << other_rsd.atom_name( donor_index ) << std::endl;

			Vector const separation( hyd_coord - h2o_coord );
			if ( separation.length_squared() > scoring::hbonds::MAX_R2 ) continue;

			scoring::hbonds::HBEvalTuple const hbe_type( donor_index, other_rsd, h2o_index, h2o_rsd );

			Real h2o_hbond_energy( 0.0 );
			bool const eval_deriv( false );
			hb_energy_deriv( *hb_database_, *hbondoptions_,
				hbe_type, donor_coord, hyd_coord, h2o_coord, h2o_base_coord, h2o_base2_coord,
				h2o_hbond_energy, eval_deriv, deriv );

			if ( h2o_hbond_energy < 0.0 ) {
				h2o_hbond_score += h2o_hbond_energy;
			}

		}

		// Second check for interactions with acceptors
		for ( Size j = 1, ej = other_rsd.accpt_pos().size() ; j <= ej ; ++j ) {
			Size accpt_index( other_rsd.accpt_pos()[j] );
			Vector accpt_coord( other_rsd.xyz( accpt_index ) );

			Vector const separation( accpt_coord - h2o_coord );
			if ( separation.length_squared() > scoring::hbonds::MAX_R2 ) continue;

			Size accpt_base_index( other_rsd.atom_base( accpt_index ) );
			Vector accpt_base_coord( other_rsd.xyz( accpt_base_index ) );
			Size accpt_base2_index( other_rsd.abase2( accpt_index ) );
			Vector accpt_base2_coord( other_rsd.xyz( accpt_base2_index ) );
			scoring::hbonds::HBEvalTuple const hbe_type( h2o_index, h2o_rsd, accpt_index, other_rsd );

			// build a fictitious hydrogen for the water
			Vector const z( cross( normal_from_base, separation ).normalized() );
			Vector const y( cross( z, normal_from_base ) );
			Vector const faux_hyd_coord( h2o_coord + HO_x_dist * normal_from_base + HO_y_dist * y );

			Real h2o_hbond_energy( 0.0 );
			bool const eval_deriv( false );
			hb_energy_deriv(
				*hb_database_, *hbondoptions_,
				hbe_type, h2o_coord, faux_hyd_coord, accpt_coord, accpt_base_coord, accpt_base2_coord,
				h2o_hbond_energy, eval_deriv, deriv );

			if ( h2o_hbond_energy < 0.0 ) {
				h2o_hbond_score += h2o_hbond_energy;
			}
		}
	}

	return h2o_hbond_score;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
This routine fills an hbond-set with hbonds. All hbonds are included,
even ones which might be excluded later based on the backbone-hbond
exclusion.
**/

void
WaterAdductHBondPotential::fill_h2o_hbond_set(
	pose::Pose const & pose,
	hbonds::HBondSet & hbond_set
) const
{
	// clear old data
	hbond_set.clear();
	hbond_set.set_hbond_options( *hbondoptions_ );

	// need to know which residues are neighbors
	// and what the neighbor-numbers are for each residue since some of the
	// weights are environment-dependent.
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

	// loop over all nbr-pairs
	for ( Size res1 = 1; res1 <= pose.total_residue(); ++res1 ) {
		int const nb1 = tenA_neighbor_graph.get_node( res1 )->num_neighbors_counting_self();
		conformation::Residue const & rsd1( pose.residue( res1 ) );

		for ( graph::Graph::EdgeListConstIter
				iru = energy_graph.get_node(res1)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(res1)->const_upper_edge_list_end();
				iru != irue; ++iru ) {

			int const res2( (*iru)->get_second_node_ind() );

			conformation::Residue const & rsd2( pose.residue( res2 ) );

			int const nb2 = tenA_neighbor_graph.get_node( res2 )->num_neighbors_counting_self();

			// rsd1 as water, rsd2 as other
			get_residue_residue_h2o_hbonds_1way( rsd1, rsd2, nb1, nb2, hbond_set );

			// rsd2 as water, rsd1 as other
			get_residue_residue_h2o_hbonds_1way( rsd2, rsd1, nb2, nb1, hbond_set );

		} // nbrs of res1
	} // res1
}


/**
compiles list of hbonds
and evaluates the derivative
**/


void
WaterAdductHBondPotential::get_residue_residue_h2o_hbonds_1way(
	// input
	conformation::Residue const & h2o_rsd,
	conformation::Residue const & other_rsd,
	int const & /*h2o_nb*/,
	int const & /*other_nb*/,
	// output
	hbonds::HBondSet & hbond_set
) const
{
	debug_assert( h2o_rsd.seqpos() != other_rsd.seqpos() ); // otherwise include in allow

	// <f1,f2> -- derivative vectors
	hbonds::HBondDerivs deriv;

	static Real const HO_dist( 0.96 );
	static Real const theta( numeric::conversions::radians( 107.0 ) );
	static Real const HO_x_dist( HO_dist * std::cos( theta ) );
	static Real const HO_y_dist( HO_dist * std::sin( theta ) );

	// Loop over waters
	for ( Size i = 1, ei = h2o_rsd.natoms() ; i <= ei ; ++i ) {
		if ( !h2o_rsd.atom_type(i).is_h2o() ) continue;

		Size const h2o_index( i );

		Vector const h2o_coord( h2o_rsd.xyz( i ) );
		Vector const h2o_base_coord( h2o_rsd.xyz( h2o_rsd.atom_base(i) ) );
		Vector const h2o_base2_coord( h2o_rsd.xyz( h2o_rsd.abase2(i) ) );
		Vector const normal_from_base( Vector( h2o_coord - h2o_base_coord ).normalized() );

		// First check for interactions with donors
		for ( Size j = 1, ej = other_rsd.Hpos_polar().size() ; j <= ej ; ++j ) {
			Size const hyd_index( other_rsd.Hpos_polar()[j] );
			Size const donor_index( other_rsd.atom_base( hyd_index ) );
			Vector const hyd_coord( other_rsd.xyz( hyd_index ) );
			Vector const donor_coord( other_rsd.xyz( donor_index ) );

			//  std::cout << "other rsd type " << other_rsd.type().name() << " hyd name " << other_rsd.atom_name( hyd_index ) << " donor name " << other_rsd.atom_name( donor_index ) << std::endl;

			Vector const separation( hyd_coord - h2o_coord );
			if ( separation.length_squared() > scoring::hbonds::MAX_R2 ) continue;

			scoring::hbonds::HBEvalTuple const hbe_type( donor_index, other_rsd, h2o_index, h2o_rsd );

			Real h2o_hbond_energy( 0.0 );
			bool const eval_deriv( true );
			hb_energy_deriv( *hb_database_, *hbondoptions_,
				hbe_type, donor_coord, hyd_coord, h2o_coord, h2o_base_coord, h2o_base2_coord,
				h2o_hbond_energy, eval_deriv, deriv );

			if ( h2o_hbond_energy < 0.0 ) {
				//    Real const weight ( get_environment_dependent_weight( hbe_type, h2o_nb, other_nb, *hbond_set.options() ) );
				Real const weight ( 1.0 );
				hbond_set.append_hbond(
					hyd_index, other_rsd, h2o_index, h2o_rsd, hbe_type,
					h2o_hbond_energy, weight, deriv );
				//    std::cout << "Stashing h2o hbond, water is acceptor " << std::endl;
				//    std::cout << "Energy is " << h2o_hbond_energy << std::endl;
				//    std::cout << "F1 is " <<
				//    deriv.first[0] << " " <<
				//    deriv.first[1] << " " <<
				//    deriv.first[2] << " " <<
				//    std::endl;
				//    std::cout << "F2 is " <<
				//    deriv.second[0] << " " <<
				//    deriv.second[1] << " " <<
				//    deriv.second[2] << " " <<
				//    std::endl;
			}
		}

		// Second check for interactions with acceptors
		for ( Size j = 1, ej = other_rsd.accpt_pos().size() ; j <= ej ; ++j ) {
			Size accpt_index( other_rsd.accpt_pos()[j] );
			Vector accpt_coord( other_rsd.xyz( accpt_index ) );

			Vector const separation( accpt_coord - h2o_coord );
			if ( separation.length_squared() > scoring::hbonds::MAX_R2 ) continue;

			Size accpt_base_index( other_rsd.atom_base( accpt_index ) );
			Vector accpt_base_coord( other_rsd.xyz( accpt_base_index ) );
			Size accpt_base2_index( other_rsd.abase2( accpt_index ) );
			Vector accpt_base2_coord( other_rsd.xyz( accpt_base2_index ) );
			scoring::hbonds::HBEvalTuple const hbe_type( h2o_index, h2o_rsd, accpt_index, other_rsd );

			// build a fictitious hydrogen for the water
			Vector const z( cross( normal_from_base, separation ).normalized() );
			Vector const y( cross( z, normal_from_base ) );
			Vector const faux_hyd_coord( h2o_coord + HO_x_dist * normal_from_base + HO_y_dist * y );

			Real h2o_hbond_energy( 0.0 );
			bool const eval_deriv( true );
			hb_energy_deriv( *hb_database_, *hbondoptions_,
				hbe_type, h2o_coord, faux_hyd_coord,
				accpt_coord, accpt_base_coord, accpt_base2_coord,
				h2o_hbond_energy, eval_deriv, deriv );

			if ( h2o_hbond_energy < 0.0 ) {
				//    Real const weight ( get_environment_dependent_weight( hbe_type, h2o_nb, other_nb, *hbond_set.options() ) );
				Real const weight ( 1.0 );
				hbond_set.append_hbond( h2o_index, h2o_rsd, accpt_index, other_rsd, hbe_type,
					h2o_hbond_energy, weight, deriv );
				//    std::cout << "Stashing h2o hbond, water is donor " << std::endl;
				//    std::cout << "Energy is " << h2o_hbond_energy << std::endl;
				//    std::cout << "F1 is " <<
				//    deriv.first[0] << " " <<
				//    deriv.first[1] << " " <<
				//    deriv.first[2] << " " <<
				//    std::endl;
				//    std::cout << "F2 is " <<
				//    deriv.second[0] << " " <<
				//    deriv.second[1] << " " <<
				//    deriv.second[2] << " " <<
				//    std::endl;
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace scoring
} // namespace core
