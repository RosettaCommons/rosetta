// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  molecular mechanics lj energy
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/methods/MMLJEnergyIntra.hh>
#include <core/scoring/methods/MMLJEnergyIntraCreator.hh>
#include <core/scoring/mm/MMLJScore.hh>

// Project headers
#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoringManager.hh>

#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/MinimizationData.hh>
#include <basic/datacache/CacheableData.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/NeighborList.tmpl.hh>

#include <core/scoring/Energies.hh>

#include <core/kinematics/MinimizerMapBase.hh>


// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <iostream>

#include <core/scoring/mm/MMLJEnergyTable.hh>
#include <core/scoring/mm/MMLJLibrary.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the MMTorsionEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MMLJEnergyIntraCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MMLJEnergyIntra );
}

ScoreTypes
MMLJEnergyIntraCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( mm_lj_intra_atr );
	sts.push_back( mm_lj_intra_rep );
	return sts;
}

MMLJEnergyIntra::MMLJEnergyIntra() :
	parent( EnergyMethodCreatorOP( new MMLJEnergyIntraCreator ) ),
	potential_( scoring::ScoringManager::get_instance()->get_MMLJEnergyTable() )
{}

/// clone
EnergyMethodOP
MMLJEnergyIntra::clone() const
{
	return EnergyMethodOP( new MMLJEnergyIntra() );
}

void
MMLJEnergyIntra::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & /*minmap*/,
	basic::datacache::BasicDataCache &,
	ResSingleMinimizationData & res_data_cache
) const {

	//std::cout << "DEBUG: In MMLJEnergyIntra setup_for_minimizing_for_residue()" << std::endl;

	if ( pose.energies().use_nblist_auto_update() ) return;

	// update the existing nblist if it's already present in the res_data_cache object
	core::scoring::ResidueNblistDataOP nbdata( utility::pointer::static_pointer_cast< core::scoring::ResidueNblistData > ( res_data_cache.get_data( mm_lj_intra_nblist ) ) );
	if ( ! nbdata ) nbdata = ResidueNblistDataOP( new ResidueNblistData );

	// get a reference to the MMLJLibrary we are using
	core::scoring::mm::MMLJLibrary const & library = potential_.mm_lj_score().mm_lj_library();

	etable::count_pair::CountPairFunctionCOP count_pair = get_intrares_countpair( rsd, pose, sfxn );

	nbdata->initialize( rsd, count_pair, library.nblist_dis2_cutoff_XX(), library.nblist_dis2_cutoff_XH(), library.nblist_dis2_cutoff_HH() );
	res_data_cache.set_data( mm_lj_intra_nblist, nbdata );
}



void
MMLJEnergyIntra::residue_pair_energy(
	conformation::Residue const & /*rsd1*/,
	conformation::Residue const & /*rsd2*/,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & /*emap*/
) const
{}


void
MMLJEnergyIntra::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & scrfxn,
	EnergyMap & emap
) const
{
	//std::cout << "DEBUG: In MMLJEnergyIntra eval_intrares_energy()" << std::endl;

	using namespace chemical;
	using namespace conformation;
	using namespace etable;
	using namespace count_pair;

	Real total_rep(0.0), total_atr(0.0);

	// get residue info
	ResidueType const & rsdtype( rsd.type() );
	Size const & rsdnatom = rsdtype.natoms();

	// get count pair function
	//CountPairFunctionOP cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn( (*this).get_intrares_countpair( rsd, pose, scrfxn ) );

	// iterate over all pairs of atom in the residue
	for ( int i = 1, i_end = rsdnatom; i <= i_end; ++i ) {
		conformation::Atom const & atom1( rsd.atom(i) );
		for ( int j = i + 1, j_end = rsdnatom; j <= j_end; ++j ) {
			conformation::Atom const & atom2( rsd.atom(j) );

			Real weight(1.0); // unused
			Size path_dist(0);

			// ask count pair if we should score it
			if ( !cpfxn->count( i, j, weight, path_dist ) ) continue;

			// calc dist
			Real dist_squared( atom1.xyz().distance_squared( atom2.xyz() ) );
			// calc energy
			Real rep(0), atr(0);
			potential_.score( rsdtype.atom( i ).mm_atom_type_index(), rsdtype.atom( j ).mm_atom_type_index(), path_dist, dist_squared, rep, atr );
			total_rep += rep;
			total_atr += atr;
		}
	}

	// add energies to emap
	emap[ mm_lj_intra_rep ] += total_rep;
	emap[ mm_lj_intra_atr ] += total_atr;
}


void
MMLJEnergyIntra::eval_intrares_energy_ext(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & res_data_cache,
	pose::Pose const & pose,
	ScoreFunction const & /*sfxn*/,
	EnergyMap & emap
) const {
	//std::cout << "DEBUG: In MMLJEnergyIntra eval_intrares_energy_ext()" << std::endl;

	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( utility::pointer::dynamic_pointer_cast< ResidueNblistData const > ( res_data_cache.get_data( mm_lj_intra_nblist ) ) );
	ResidueNblistData const & nblist( static_cast< ResidueNblistData const & > ( res_data_cache.get_data_ref( mm_lj_intra_nblist ) ) );

	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		conformation::Atom const & atom1( rsd.atom( neighbs[ ii ].atomno1() ) );
		conformation::Atom const & atom2( rsd.atom( neighbs[ ii ].atomno2() ) );

		// get squared distance
		Real dist_squared( atom1.xyz().distance_squared( atom2.xyz() ) );

		// get pat distance b/w atoms
		Size path_dist( neighbs[ ii ].path_dist() );

		// calc score
		Real rep(0), atr(0);
		potential_.score( atom1.mm_type(), atom2.mm_type(), path_dist, dist_squared, rep, atr );

		// add energies to emap
		emap[ mm_lj_intra_rep ] += rep;
		emap[ mm_lj_intra_atr ] += atr;
	}
}


void
MMLJEnergyIntra::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & res_data_cache,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {

	//std::cout << "DEBUG: In MMLJEnergyIntra eval_intrares_derivatives()" << std::endl;

	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( utility::pointer::dynamic_pointer_cast< ResidueNblistData const > ( res_data_cache.get_data( mm_lj_intra_nblist ) ) );
	ResidueNblistData const & nblist( static_cast< ResidueNblistData const & > ( res_data_cache.get_data_ref( mm_lj_intra_nblist ) ) );

	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		conformation::Atom const & atom1( rsd.atom( neighbs[ ii ].atomno1() ) );
		conformation::Atom const & atom2( rsd.atom( neighbs[ ii ].atomno2() ) );

		// initialize f1 and f2
		Vector f1( atom1.xyz().cross( atom2.xyz() ) );
		Vector f2( atom1.xyz() - atom2.xyz() );

		// get squared distance
		Real dist_squared( atom1.xyz().distance_squared( atom2.xyz() ) );

		// get pat distance b/w atoms
		Size path_dist( neighbs[ ii ].path_dist() );

		// calc deriv // make sure we are calculating dE_dR_over_r (SEE NOTE BELOW)
		Real drep(0), datr(0);
		potential_.deriv_score( atom1.mm_type(), atom2.mm_type(), path_dist, dist_squared, drep, datr );
		Real deriv( weights[ mm_lj_intra_rep ] * drep + weights[ mm_lj_intra_atr ] * datr ); // DOUG DOUG DOUG NOTE: THIS MIGHT NEED TO BE DIVIED BY SQRT(R) LOOK AT MMLJSCORE

		// if non-zero add to vectors
		if ( deriv != 0.0 ) {
			atom_derivs[ neighbs[ ii ].atomno1() ].f1() += deriv * f1;
			atom_derivs[ neighbs[ ii ].atomno1() ].f2() += deriv * f2;
			atom_derivs[ neighbs[ ii ].atomno2() ].f1() += -1 * deriv * f1;
			atom_derivs[ neighbs[ ii ].atomno2() ].f2() += -1 * deriv * f2;
		}
	}
}

/// @brief MMLJEnergy does not have an atomic interation threshold
Distance
MMLJEnergyIntra::atomic_interaction_cutoff() const
{
	// this will probably screw up other stuff, but it might not since etable goes to zero at 5.5
	return 7.0;
}

/// @brief MMLJEnergy is context independent; indicates that no context graphs are required
void
MMLJEnergyIntra::indicate_required_context_graphs(utility::vector1< bool > & ) const
{}

/// @brief required for neighbor list and to be more lke the ETable
etable::count_pair::CountPairFunctionCOP
MMLJEnergyIntra::get_count_pair_function(
	Size res1,
	Size res2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) const
{
	return get_count_pair_function( pose.residue( res1 ), pose.residue( res2 ), pose, sfxn );
}

/// @brief required for neighbor list and to be more lke the ETable
etable::count_pair::CountPairFunctionCOP
MMLJEnergyIntra::get_count_pair_function(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*sfxn*/
) const
{
	using namespace etable;
	using namespace count_pair;

	return CountPairFunctionOP( CountPairFactory::create_count_pair_function( res1, res2, CP_CROSSOVER_3 ) );
}

/// @brief required for neighbor list and to be more lke the ETable
etable::count_pair::CountPairFunctionOP
MMLJEnergyIntra::get_intrares_countpair(
	conformation::Residue const & res,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*sfxn*/
) const
{
	using namespace etable;
	using namespace count_pair;

	return CountPairFunctionOP( CountPairFactory::create_intrares_count_pair_function( res, CP_CROSSOVER_3 ) );
}
core::Size
MMLJEnergyIntra::version() const
{
	return 1; // Initial versioning
}


} // namespace methods
} // namespace scoring
} // namespace core
