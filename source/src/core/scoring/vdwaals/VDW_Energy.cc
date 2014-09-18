// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/vdwaals/VDW_Energy.cc
/// @brief  Low resolution (centroid) repulsive energy
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/vdwaals/VDW_Energy.hh>
#include <core/scoring/vdwaals/VDW_EnergyCreator.hh>

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

#include <core/scoring/trie/RotamerTrie.hh>
#include <core/scoring/trie/trie.functions.hh>

#include <core/scoring/vdwaals/VDWTrie.hh>

#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>
#include <core/scoring/etable/etrie/TrieCountPairAll.hh>
#include <core/scoring/etable/etrie/TrieCountPairGeneric.hh>

// Project headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/prof.hh>
#include <core/id/AtomID.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

namespace core {
namespace scoring {
namespace vdwaals {

static thread_local basic::Tracer TR( "core.scoring.vdwaals.VDW_Energy" );

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
VDW_Energy::VDW_Energy( methods::EnergyMethodOptions const & options ):
	parent( new VDW_EnergyCreator ),
	atom_vdw_( ScoringManager::get_instance()->get_AtomVDW( options.atom_vdw_atom_type_set_name() ) ),
	atom_type_set_name_( options.atom_vdw_atom_type_set_name() ),
	vdw_scale_factor_( 0.8 ) // hack from rosetta++
{
}


/// clone
methods::EnergyMethodOP
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

void
VDW_Energy::prepare_rotamers_for_packing(
	pose::Pose const &,
	conformation::RotamerSetBase & set
) const
{
	VDWRotamerTrieOP trie = create_rotamer_trie( set );
	set.store_trie( methods::vdw_method, trie );
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
	//	basic::ProfileThis doit( basic::VDW_ENERGY );
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
// 		using namespace ObjexxFCL::format;
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
	//	basic::ProfileThis doit( basic::VDW_ENERGY );
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

void
VDW_Energy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	assert( set1.resid() != set2.resid() );

	using namespace methods;
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table1( energy_table );
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table2( energy_table );

	temp_table1 = 0; temp_table2 = 0;

	VDWTrieEvaluator vdw_eval( *this, 6, hydrogen_interaction_cutoff2_, weights[ vdw ] );

	VDWRotamerTrieCOP trie1( static_cast< trie::RotamerTrieBase const * > ( set1.get_trie( vdw_method )() ));
	VDWRotamerTrieCOP trie2( static_cast< trie::RotamerTrieBase const * > ( set2.get_trie( vdw_method )() ));

	// figure out which trie countPairFunction needs to be used for this set
	trie::TrieCountPairBaseOP cp = get_count_pair_function_trie( set1, set2, pose );

	trie1->trie_vs_trie( *trie2, *cp, vdw_eval, temp_table1, temp_table2 );
	energy_table += temp_table1;
}

trie::TrieCountPairBaseOP
VDW_Energy::get_count_pair_function_trie(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	core::pose::Pose const & pose
) const
{
	conformation::Residue const & res1( pose.residue( set1.resid() ) );
	conformation::Residue const & res2( pose.residue( set2.resid() ) );
	trie::RotamerTrieBaseCOP trie1( static_cast< trie::RotamerTrieBase const * > ( set1.get_trie( methods::vdw_method )() ));
	trie::RotamerTrieBaseCOP trie2( static_cast< trie::RotamerTrieBase const * > ( set2.get_trie( methods::vdw_method )() ));

	return get_count_pair_function_trie( res1, res2, trie1, trie2 );
}

trie::TrieCountPairBaseOP
VDW_Energy::get_count_pair_function_trie(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	trie::RotamerTrieBaseCOP trie1,
	trie::RotamerTrieBaseCOP trie2
) const
{
	using namespace etable::count_pair;
	using namespace etable::etrie;

	CPResidueConnectionType connection = CountPairFactory::determine_residue_connection( res1, res2 );
	Size conn1 = trie1->get_count_pair_data_for_residue( res2.seqpos() );
	Size conn2 = trie2->get_count_pair_data_for_residue( res1.seqpos() );
	if ( connection == CP_ONE_BOND ) {
		return new VDWTrieCountPair1B( conn1, conn2 );
	} else if ( connection == CP_NO_BONDS ) {
		return new TrieCountPairAll;
	} else {
		return new TrieCountPairGeneric( res1, res2, conn1, conn2 );
	}


	return 0;
}

core::Size
VDW_Energy::version() const
{
	return 1; // Initial versioning
}

void
VDW_Energy::calculate_hydrogen_interaction_cutoff()
{
	TR.Debug << "calculating hydrogen_interaction_cutoff" << std::endl;
	// get the relevant atomset
	chemical::AtomTypeSet const & atom_set
		( *chemical::ChemicalManager::get_instance()->atom_type_set( atom_vdw_.atom_type_set_name() ) );

	hydrogen_interaction_cutoff2_ = 0;
	Size which_ii(0), which_jj(0);
	for ( core::Size ii = 1; ii <= atom_set.n_atomtypes(); ++ii ) {
		if ( atom_set[ ii ].is_hydrogen() ) {
			for ( core::Size jj = ii; jj <= atom_set.n_atomtypes(); ++jj ) {
				if ( atom_set[ jj ].is_hydrogen() ) {
					Real iijj_interaction_dist = atom_vdw_(ii)[jj];
					if ( iijj_interaction_dist > hydrogen_interaction_cutoff2_ ) {
						hydrogen_interaction_cutoff2_ = iijj_interaction_dist * iijj_interaction_dist;
						which_ii = ii;
						which_jj = jj;
					}
				}
			}
		}
	}
	if ( which_ii != 0 ) {
		TR.Debug << "Found hydrogen interaction radius 2 of " << hydrogen_interaction_cutoff2_ << "A^2 (" << std::sqrt( hydrogen_interaction_cutoff2_ );
		TR.Debug << " A) between atoms " << atom_set[ which_ii ].atom_type_name() << " and " << atom_set[ which_jj ].atom_type_name() << std::endl;
	} else {
		TR.Debug << "Did not find any hydrogen atoms in atom type set " << atom_set.name() << std::endl;
	}

}

VDWRotamerTrieOP
VDW_Energy::create_rotamer_trie(
	conformation::RotamerSetBase const & rotset
) const
{
	using namespace trie;
	using namespace etable::etrie;

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamerset( rotset ) );
	VDWAtom at;
	if ( cpdata_map.has_pseudobonds() ||
			cpdata_map.max_connpoints_for_residue() > 1 ||
			cpdata_map.n_entries() > 3 ) {
		CountPairDataGeneric cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		CountPairData_1_1 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 2 ) {
		CountPairData_1_2 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 3 ) {
		CountPairData_1_3 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, atomic_interaction_cutoff() );
	} else {
		/// As of 10/21, all count pair data combinations should be covered. This code should not execute.
		std::cerr << "Unsupported number of residue connections in trie construction." << std::endl;
		utility_exit();
		return 0;
	}
}


VDWTrieEvaluator::VDWTrieEvaluator(
	VDW_Energy const & vdw,
	Real const atomic_interaction_cutoff,
	Real const hydrogen_interaction_cutoff2,
	Real const vdw_weight
) :
	vdw_( vdw ),
	atomic_interaction_cutoff_( atomic_interaction_cutoff ),
	hydrogen_interaction_cutoff2_( hydrogen_interaction_cutoff2 ),
	vdw_weight_( vdw_weight )
{}

Distance
VDWTrieEvaluator::atomic_interaction_cutoff() const
{
	return atomic_interaction_cutoff_;
}

Real
VDWTrieEvaluator::hydrogen_interaction_cutoff2() const
{
	return hydrogen_interaction_cutoff2_;
}

Real
VDWTrieEvaluator::vdw_weight() const
{
	return vdw_weight_;
}

}
}
}
