// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/VDW_Energy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/VDW_Energy.hh>
#include <core/scoring/methods/VDW_EnergyCreator.hh>

// Package headers
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/prof.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the VDW_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
VDW_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new VDW_Energy( options );
}

ScoreTypes
VDW_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( vdw );
	return sts;
}


/// @details  C-TOR with method options object
VDW_Energy::VDW_Energy( EnergyMethodOptions const & options ):
	parent( new VDW_EnergyCreator ),
	atom_vdw_( ScoringManager::get_instance()->get_AtomVDW( options.atom_vdw_atom_type_set_name() ) ),
	atom_type_set_name_( options.atom_vdw_atom_type_set_name() ),
	vdw_scale_factor_( 0.8 ) // hack from rosetta++
{
}


/// clone
EnergyMethodOP
VDW_Energy::clone() const
{
	return new VDW_Energy( *this );
}

/// @details  copy c-tor
VDW_Energy::VDW_Energy( VDW_Energy const & src ):
	parent( src ),
	atom_vdw_( src.atom_vdw_ ),
	atom_type_set_name_( src.atom_type_set_name_ ),
	vdw_scale_factor_( src.vdw_scale_factor_ )
{}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
VDW_Energy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


///
void
VDW_Energy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


///
void
VDW_Energy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using namespace etable::count_pair;
	Real score(0.0);
	basic::ProfileThis doit( basic::VDW_ENERGY );
	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		// assuming only a single bond right now -- generalizing to arbitrary topology
		// also assuming crossover of 4, should be closest (?) to classic rosetta
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		for ( Size i = 1, i_end = rsd1.natoms(); i <= i_end; ++i ) {
			Vector const & i_xyz( rsd1.xyz(i) );
			Size const i_type( rsd1.atom_type_index(i) );
			utility::vector1< Real > const & i_atom_vdw( atom_vdw_( i_type ) );
			for ( Size j = 1, j_end = rsd2.natoms(); j <= j_end; ++j ) {
				Real weight(1.0);
				Size path_dist( 0 );
				if ( cpfxn->count( i, j, weight, path_dist ) ) {
					if ( weight < 0.99 ) continue; // don't count half-weight interxns in vdw_compute
					if ( rsd2.atom_type_index(j) <= i_atom_vdw.size() ){
						Real const bump_dsq( i_atom_vdw[ rsd2.atom_type_index(j) ] );
						Real const clash( bump_dsq - i_xyz.distance_squared( rsd2.xyz(j) ) );
						if ( clash > 0.0 ) {
							score += ( clash * clash ) / bump_dsq;
						}
					}else{
						std::cerr << "Etable: "  <<  rsd2.atom_type_index(j) << " " <<  i_atom_vdw.size() << "  " << i << " " <<  j << std::endl;
					  utility_exit_with_message( "Fatal Error in VDW_Energy" );
					}
				}
			}
		}
	} else {
		// no countpair!
		for ( Size i = 1, i_end = rsd1.natoms(); i <= i_end; ++i ) {
			Vector const & i_xyz( rsd1.xyz(i) );
			Size const i_type( rsd1.atom_type_index(i) );
			utility::vector1< Real > const & i_atom_vdw( atom_vdw_( i_type ) );
			for ( Size j = 1, j_end = rsd2.natoms(); j <= j_end; ++j ) {
				if ( rsd2.atom_type_index(j) <= i_atom_vdw.size() ){
					Real const bump_dsq( i_atom_vdw[ rsd2.atom_type_index(j) ] );
					Real const clash( bump_dsq - i_xyz.distance_squared( rsd2.xyz(j) ) );
					if ( clash > 0.0 ) {
						score += ( clash * clash ) / bump_dsq;
//  					std::cout << "BUMP: " << I(4,rsd1.seqpos() ) << I(4,rsd2.seqpos() )
// 										<<' ' << rsd1.atom_name(i) << ' ' << rsd2.atom_name(j) << ' '
// 										<< ( clash * clash ) / bump_dsq * vdw_scale_factor_
// 										<< ' ' << i_xyz.distance_squared( rsd2.xyz(j) ) <<  std::endl;
					}
				}else{
					std::cerr << "Etable: " <<  rsd2.atom_type_index(j) << " " <<  i_atom_vdw.size() << "  " << i << " " <<  j << std::endl;
				  utility_exit_with_message( "Fatal Error in VDW_Energy" );
				}
			}
		}
	}
	emap[ vdw ] += score * vdw_scale_factor_; // vdw prefactor!

// 	if ( score*0.8 > 0.001 && rsd1.seqpos() <= rsd2.seqpos() ) {
// 		using namespace ObjexxFCL::fmt;
// 		std::cout << "vdw_ij: " << I(4,rsd1.seqpos() ) << I(4,rsd2.seqpos() ) << F(9,3,score*0.8) << std::endl;
// 	}
}


void
VDW_Energy::eval_atom_derivative(
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
	basic::ProfileThis doit( basic::VDW_ENERGY );
	// what is my charge?
	Size const pos1( atom_id.rsd() );
	Size const i   ( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( pos1 ) );

	int const pos1_map( domain_map( pos1 ) );
	bool const pos1_fixed( pos1_map != 0 );

	Vector const & i_xyz( rsd1.xyz(i) );
	Size const i_type( rsd1.atom_type_index(i) );
	utility::vector1< Real > const & i_atom_vdw( atom_vdw_( i_type ) );

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// loop over *all* nbrs of rsd1 (not just upper or lower)
	for ( graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			irue = energy_graph.get_node( pos1 )->const_edge_list_end();
			iru != irue; ++iru ) {
		Size const pos2( (*iru)->get_other_ind( pos1 ) );

		if ( pos1_fixed && pos1_map == domain_map( pos2 ) ) continue; // fixed wrt one another

		conformation::Residue const & rsd2( pose.residue( pos2 ) );

		assert( pos2 != pos1 );

		if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
			// generalizing to arbitrary topology
			// also assuming crossover of 4, should be closest (?) to classic rosetta
			CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

			for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {

				Real cp_weight(1.0);
				Size path_dist( 0 );
				if ( cpfxn->count( i, j, cp_weight, path_dist ) ) {
					if ( cp_weight < 0.99 ) continue; // dont count half-weight interxns in vdw_compute
					Vector const & j_xyz( rsd2.xyz(j) );
					Vector const f2( i_xyz - j_xyz );
					Real const dis2( f2.length_squared() );
					Real const bump_dsq( i_atom_vdw[ rsd2.atom_type_index(j) ] );
					if ( dis2 < bump_dsq ) {
						// E += vdw_scale_factor_ * weights[vdw] * cp_weight * ( ( dis2 - bump_dsq ) **2 ) / bump_dsq
						Real const dE_dr_over_r = vdw_scale_factor_ * weights[ vdw ] * cp_weight * 4.0 * ( dis2 - bump_dsq ) / bump_dsq;
						Vector const f1( i_xyz.cross( j_xyz ) );
						F1 += dE_dr_over_r * f1;
						F2 += dE_dr_over_r * f2;
					}
				}
			}
		} else {
			// no countpair!
			for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {

				Vector const & j_xyz( rsd2.xyz(j) );
				Vector const f2( i_xyz - j_xyz );
				Real const dis2( f2.length_squared() );
				Real const bump_dsq( i_atom_vdw[ rsd2.atom_type_index(j) ] );
				if ( dis2 < bump_dsq ) {
					// E += vdw_scale_factor_ * weights[vdw] * ( ( dis2 - bump_dsq ) **2 ) / bump_dsq
					Real const dE_dr_over_r = vdw_scale_factor_ * weights[ vdw ] * 4.0 * ( dis2 - bump_dsq ) / bump_dsq;
					Vector const f1( i_xyz.cross( j_xyz ) );
					F1 += dE_dr_over_r * f1;
					F2 += dE_dr_over_r * f2;
				}
			}
		} // are rsd1 and rsd2 bonded?

	} // loop over nbrs of rsd1

}



/// @brief VDW_Energy distance cutoff
Distance
VDW_Energy::atomic_interaction_cutoff() const
{
	return 6.0; /// now subtracted off 3.0 from cutoffs in centroid params files
	//return 0.0; /// since all the cutoffs for centroid mode are rolled into the cendist check
}

/// @brief VDW_Energy
void
VDW_Energy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}
core::Size
VDW_Energy::version() const
{
	return 1; // Initial versioning
}



}
}
}
