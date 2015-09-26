// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

/// @details this method heinously const-casts the Pose reference it's given
/// this is because we need to shuck a nblist into the Pose's Energies
/// because the ResSingleMinimizationData isn't working for a NeighborListOP?
void
MMLJEnergyIntra::setup_for_minimizing_for_residue(
	conformation::Residue const & /*rsd*/,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & minmap,
	ResSingleMinimizationData & res_data_cache
) const {
	
	// AMW: oh maybe the issue here is that we keep on calling it
	// and putting the nblist in the pose over and over and this debug_assert
	// that is failing is saying "hey there is already a nblist in here"
	
	if ( !pose.energies().use_nblist() ) return;
	
	//Return right now so we don't set the nblist again?
	if ( pose.energies().use_nblist_of_type( EnergiesCacheableDataType::MM_LJ_INTRA_NBLIST ) ) return;
	
	// stash our nblist inside the pose's energies object
	// get a reference to the MMLJLibrary we are using
	core::scoring::mm::MMLJLibrary const & library = potential_.mm_lj_score().mm_lj_library();
	
	// setup the atom-atom nblist
	NeighborListOP nblist( new NeighborList(
		minmap.domain_map(),
		library.nblist_dis2_cutoff_XX(),
		library.nblist_dis2_cutoff_XH(),
		library.nblist_dis2_cutoff_HH()) );
	
	if ( pose.energies().use_nblist_auto_update() ) {
		// setting this to one angstrom fudge factor
		nblist->set_auto_update( 1 ); // MAGIC NUMBER
	}
	
	// this partially becomes the MMLJEnergy classes's responsibility
	nblist->setup( pose, sfxn, *this );
	
	res_data_cache.set_data( mm_lj_intra_nblist, utility::pointer::dynamic_pointer_cast< basic::datacache::CacheableData/*COP*/ >( NeighborListDataOP( new NeighborListData( nblist ) ) ) );
	//std::cout << "In MMLJEnergyIntra setup_for_minimizing_for_residue()" << std::endl;
}
	
void
MMLJEnergyIntra::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const
{
	// Could this be conflicting somehow?
	return;
	// taken for the most part from BaseETableEnergy.hh

	if ( pose.energies().use_nblist() ) {
		// stash our nblist inside the pose's energies object
		Energies & energies( pose.energies() );

		// get a reference to the MMLJLibrary we are using
		core::scoring::mm::MMLJLibrary const & library = potential_.mm_lj_score().mm_lj_library();

		// setup the atom-atom nblist
		NeighborListOP nblist( new NeighborList(
			min_map.domain_map(),
			library.nblist_dis2_cutoff_XX(),
			library.nblist_dis2_cutoff_XH(),
			library.nblist_dis2_cutoff_HH()) );

		if ( pose.energies().use_nblist_auto_update() ) {
			// setting this to one angstrom fudge factor
			nblist->set_auto_update( 1 ); // MAGIC NUMBER
		}
		// this partially becomes the MMLJEnergy classes's responsibility
		nblist->setup( pose, sfxn, *this );

		energies.set_nblist( EnergiesCacheableDataType::MM_LJ_INTRA_NBLIST, nblist );
		NeighborList nblist_blah( energies.nblist( EnergiesCacheableDataType::MM_LJ_INTRA_NBLIST ) );
		//std::cout << "In MMLJEnergyIntra setup_for_minimizing()" << std::endl;
	}
}


bool
MMLJEnergyIntra::defines_intrares_energy( EnergyMap const & /*weights*/ ) const
{
	return true;
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
			if ( cpfxn->count( i, j, weight, path_dist ) ) {
				// calc dist
				Real dist_squared( atom1.xyz().distance_squared( atom2.xyz() ) );
				// calc energy
				Real rep(0), atr(0);
				potential_.score( rsdtype.atom( i ).mm_atom_type_index(), rsdtype.atom( j ).mm_atom_type_index(), path_dist, dist_squared, rep, atr );
				total_rep += rep;
				total_atr += atr;
				//      std::cout << "INTRA"
				//          << " RSD: "  << std::setw(22) << rsdtype.name()
				//          << " ATM1: " << std::setw(4)  << rsdtype.mm_atom_name(i)
				//          << " ATM2: " << std::setw(4)  << rsdtype.mm_atom_name(j)
				//          << " PDST: " << std::setw(4)  << path_dist
				//          << " DIST: " << std::setw(8)  << atom1.xyz().distance( atom2.xyz() )
				//          << " WGHT: " << std::setw(4)  << weight
				//          << " ENER: " << std::setw(8)  << energy
				//          << std::endl;

			}
		}
	}

	// add energies to emap
	emap[ mm_lj_intra_rep ] += total_rep;
	emap[ mm_lj_intra_atr ] += total_atr;
}

void
MMLJEnergyIntra::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const & /* domain_map*/,
	ScoreFunction const & /*sfxn*/,
	EnergyMap const & /*weights*/,
	Vector & F1,
	Vector & F2
) const
{
	// get atom1 from id
	conformation::Atom const & atom1( pose.residue( id.rsd() ).atom( id.atomno() ) );

	// check to see if we are using the neighbor list
	if ( pose.energies().use_nblist() ) {

		// get atom1's neighbors
		scoring::Energies const & energies(   pose.energies() );
		scoring::NeighborList const & nblist( energies.nblist( EnergiesCacheableDataType::MM_LJ_INTRA_NBLIST ) );
		scoring::AtomNeighbors const & nbrs(  nblist.atom_neighbors( id ) );

		// iterate over all of atom1's neighbors using the neighbor list
		for ( scoring::AtomNeighbors::const_iterator iter=nbrs.begin(), iter_end=nbrs.end(); iter != iter_end; ++iter ) {
			scoring::AtomNeighbor const & nbr( *iter );
			conformation::Atom const & atom2( pose.residue( nbr.rsd() ).atom( nbr.atomno() ) );

			// calculate f1 and f2
			Vector f1( atom1.xyz().cross( atom2.xyz() ) );
			Vector f2( atom1.xyz() - atom2.xyz() );

			// get count pair weight from neighbor list
			// Real const cp_weight( nbr.weight() );

			// get distance
			Real dist_squared( atom1.xyz().distance_squared( atom2.xyz() ) );

			// get path distance
			Size path_dist( nbr.path_dist() );

			// calc deriv
			Real drep(0), datr(0);
			potential_.deriv_score( atom1.mm_type(), atom2.mm_type(), path_dist, dist_squared, drep, datr );
			Real deriv( drep + datr );

			if ( deriv != 0.0 ) {
				F1 += deriv * f1;
				F2 += deriv * f2;
			}
		}
	} else {
		utility_exit_with_message("non-nblist minimize!");
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
	if ( !pose.energies().use_nblist() ) {
		utility_exit_with_message("non-nblist minimize!");
	}
	
	scoring::NeighborList const & nblist( *utility::pointer::dynamic_pointer_cast< const NeighborListData >( res_data_cache.get_data( mm_lj_intra_nblist ) )->nblist() );

	for ( Size atomno = 1; atomno <= rsd.type().natoms(); ++atomno ) {
		// get atom1 from id
		conformation::Atom const & atom1( rsd.atom( atomno ) );
		scoring::AtomNeighbors const & nbrs(  nblist.atom_neighbors( rsd.seqpos(), atomno ) );

		// iterate over all of atom1's neighbors using the neighbor list
		for ( scoring::AtomNeighbors::const_iterator iter=nbrs.begin(), iter_end=nbrs.end(); iter != iter_end; ++iter ) {
			scoring::AtomNeighbor const & nbr( *iter );
			conformation::Atom const & atom2( pose.residue( nbr.rsd() ).atom( nbr.atomno() ) );
				
			// calculate f1 and f2
			Vector f1( atom1.xyz().cross( atom2.xyz() ) );
			Vector f2( atom1.xyz() - atom2.xyz() );
	
			// get count pair weight from neighbor list
			// Real const cp_weight( nbr.weight() );
	
			// get distance
			Real dist_squared( atom1.xyz().distance_squared( atom2.xyz() ) );
	
			// get path distance
			Size path_dist( nbr.path_dist() );
	
			// calc deriv
			Real drep(0), datr(0);
			potential_.deriv_score( atom1.mm_type(), atom2.mm_type(), path_dist, dist_squared, drep, datr );
			Real deriv( weights[ mm_lj_intra_rep ] * drep + weights[ mm_lj_intra_atr ] * datr );

			if ( deriv != 0.0 ) {
				atom_derivs[ atomno ].f1() += deriv * f1;
				atom_derivs[ atomno ].f2() += deriv * f2;
			}
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
