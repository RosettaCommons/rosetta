// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MMLJEnergyInter.hh
/// @brief  molecular mechanics lj energy
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/methods/MMLJEnergyInter.hh>
#include <core/scoring/methods/MMLJEnergyInterCreator.hh>
#include <core/scoring/mm/MMLJScore.hh>

#include <core/scoring/mm/mmtrie/MMEnergyTableAtom.hh>

// Project headers
#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

#include <core/scoring/NeighborList.hh>
#include <core/scoring/NeighborList.tmpl.hh>

#include <core/scoring/Energies.hh>

#include <core/kinematics/MinimizerMapBase.hh>

#include <basic/prof.hh>

#include <core/scoring/ScoringManager.hh>

#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/trie/CPDataCorrespondence.hh>
#include <core/scoring/trie/RotamerDescriptor.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <core/scoring/trie/TrieCollection.hh>
#include <core/scoring/trie/TrieCountPairBase.hh>
#include <core/scoring/trie/trie.functions.hh>

#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/TrieCountPair1BC3.hh>
#include <core/scoring/etable/etrie/TrieCountPairAll.hh>

//#include <core/pack/rotamer_set/RotamerSetFactory.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

#include <utility/numbers.hh>

// C++ headers
#include <iostream>

#include <core/scoring/mm/MMLJLibrary.hh>
#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Project serialization headers
#include <core/scoring/trie/RotamerTrie.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>

#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the MMTorsionEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MMLJEnergyInterCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MMLJEnergyInter );
}

ScoreTypes
MMLJEnergyInterCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( mm_lj_inter_atr );
	sts.push_back( mm_lj_inter_rep );
	return sts;
}

MMLJEnergyInter::MMLJEnergyInter() :
	parent( EnergyMethodCreatorOP( new MMLJEnergyInterCreator ) ),
	potential_( scoring::ScoringManager::get_instance()->get_MMLJEnergyTable() )
{}

/// clone
EnergyMethodOP
MMLJEnergyInter::clone() const
{
	return EnergyMethodOP( new MMLJEnergyInter() );
}

void
MMLJEnergyInter::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const
{
	// taken for the most part from BaseETableEnergy.hh

	if ( !pose.energies().use_nblist() ) return;

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

	energies.set_nblist( EnergiesCacheableDataType::MM_LJ_INTER_NBLIST, nblist );
	//std::cout << "In MMLJEnergyInter setup_for_minimizing()" << std::endl;
}

// The MMLJEnergyInter method stores a vector of rotamer trie objects in the Energies
// object for use in rapid rotamer/background energy calculations.  Overrides default
// do-nothing behavior.
void
MMLJEnergyInter::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	using namespace trie;

	TrieCollectionOP tries( new TrieCollection );
	tries->total_residue( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		// Do not compute energy for virtual residues.
		if ( pose.residue_type(ii).aa() == core::chemical::aa_vrt ) continue;

		RotamerTrieBaseOP one_rotamer_trie = create_rotamer_trie( pose.residue( ii ), pose );
		tries->trie( ii, one_rotamer_trie );
	}
	pose.energies().data().set( EnergiesCacheableDataType::MM_LJ_TRIE_COLLECTION, tries );

}

/// @brief create a rotamer trie for a particular set, deciding upon the kind of count pair data that
/// needs to be contained by the trie.
///
trie::RotamerTrieBaseOP
MMLJEnergyInter::create_rotamer_trie(
	conformation::RotamerSetBase const & rotset,
	pose::Pose const & // will be need to create tries for disulfides
) const
{
	using namespace trie;
	using namespace etable::etrie;
	using namespace core::scoring::mm::mmtrie;

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamerset( rotset ) );
	if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		MMEnergyTableAtom at; CountPairData_1_1 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, hydrogen_interaction_cutoff()  );
	} else if ( cpdata_map.n_entries() == 2 ) {
		MMEnergyTableAtom at; CountPairData_1_2 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, hydrogen_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 3 ) {
		MMEnergyTableAtom at; CountPairData_1_3 cpdat;
		return create_trie( rotset, at, cpdat, cpdata_map, hydrogen_interaction_cutoff() );
	} else {
		std::cerr << "Unsupported number of residue connections in trie construction." << std::endl;
		utility_exit();
		return nullptr;
	}
}

trie::RotamerTrieBaseOP
MMLJEnergyInter::create_rotamer_trie(
	conformation::Residue const & res,
	pose::Pose const & // will be need to create tries for disulfides
) const
{
	using namespace trie;
	using namespace etable::etrie;
	using namespace core::scoring::mm::mmtrie;

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamer( res ) );
	if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		MMEnergyTableAtom at; CountPairData_1_1 cpdat;
		return create_trie( res, at, cpdat, cpdata_map, hydrogen_interaction_cutoff()  );
	} else if ( cpdata_map.n_entries() == 2 ) {
		MMEnergyTableAtom at; CountPairData_1_2 cpdat;
		return create_trie( res, at, cpdat, cpdata_map, hydrogen_interaction_cutoff() );
	} else if ( cpdata_map.n_entries() == 3 ) {
		MMEnergyTableAtom at; CountPairData_1_3 cpdat;
		return create_trie( res, at, cpdat, cpdata_map, hydrogen_interaction_cutoff() );
	} else {
		std::cerr << "Unsupported number of residue connections in trie construction." << std::endl;
		utility_exit();
		return nullptr;
	}
}


// @brief Creates a rotamer trie for the input set of rotamers and stores the trie
// in the rotamer set.
void
MMLJEnergyInter::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & set
) const
{
	trie::RotamerTrieBaseOP rottrie = create_rotamer_trie( set, pose );
	set.store_trie( methods::mm_lj_energy_inter_method, rottrie );
}

// @brief Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
MMLJEnergyInter::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	using namespace trie;

	trie::RotamerTrieBaseOP one_rotamer_trie = create_rotamer_trie( pose.residue( resid ), pose );

	// grab non-const & of the cached tries and replace resid's trie with a new one.
	auto & trie_collection
		( static_cast< TrieCollection & > (pose.energies().data().get( EnergiesCacheableDataType::MM_LJ_TRIE_COLLECTION )));
	trie_collection.trie( resid, one_rotamer_trie );
}


bool
MMLJEnergyInter::defines_intrares_energy( EnergyMap const & ) const
{
	return false;
}

void
MMLJEnergyInter::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scrfxn,
	EnergyMap & emap
) const
{
	using namespace chemical;
	using namespace conformation;
	using namespace etable;
	using namespace count_pair;

	Real total_rep( 0.0 ), total_atr( 0.0 );

	// get residue info
	ResidueType const & rsd1type( rsd1.type() );
	ResidueType const & rsd2type( rsd2.type() );
	Size const & rsd1natom = rsd1type.natoms();
	Size const & rsd2natom = rsd2type.natoms();

	// get count pair function
	//CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionCOP cpfxn( (*this).get_count_pair_function( rsd1, rsd2, pose, scrfxn ) );

	// iterate over all pairs of atom between the residues
	for ( int i = 1, i_end = rsd1natom; i <= i_end; ++i ) {
		conformation::Atom const & atom1( rsd1.atom(i) );
		for ( int j = 1, j_end = rsd2natom; j <= j_end; ++j ) {
			conformation::Atom const & atom2( rsd2.atom(j) );

			Real weight(1.0); // unused
			Size path_dist(0);

			// ask count pair if we should score it
			if ( ! cpfxn->count( i, j, weight, path_dist ) ) continue;

			// calc dist
			Real dist_squared( atom1.xyz().distance_squared( atom2.xyz() ) );
			// calc energy
			Real rep(0), atr(0);
			potential_.score( rsd1type.atom( i ).mm_atom_type_index(), rsd2type.atom( j ).mm_atom_type_index(), path_dist, dist_squared, rep, atr );

			if ( utility::isnan(rep) ) {
				std::cout << "REP NAN REP NAN REP NAN" << std::endl;
				rep = 0;
				potential_.score( rsd1type.atom( i ).mm_atom_type_index(), rsd2type.atom( j ).mm_atom_type_index(), path_dist, dist_squared, rep, atr );
			}
			if ( utility::isnan(atr) ) {
				std::cout << "ATR NAN ATR NAN ATR NAN" << std::endl;
				atr = 0;
				potential_.score( rsd1type.atom( i ).mm_atom_type_index(), rsd2type.atom( j ).mm_atom_type_index(), path_dist, dist_squared, rep, atr );
			}

			total_rep += rep;
			total_atr += atr;

			//     std::cout << "INTER"
			//          << " RSD1: " << std::setw(22) << rsd1type.name()
			//          << " ATM1: " << std::setw(4)  << rsd1type.atom(i).mm_atom_name()
			//          << " RSD2: " << std::setw(22) << rsd2type.name()
			//          << " ATM2: " << std::setw(4)  << rsd2type.atom(j).mm_atom_name()
			//          << " PDST: " << std::setw(3)  << path_dist
			//          << " DIST: " << std::setw(8)  << std::sqrt(dist_squared)
			//          << " WGHT: " << std::setw(4)  << weight
			//          << " REP: " << std::setw(8)  << rep
			//          << " ATR: " << std::setw(8)  << atr
			//          << std::endl;
		}
	}

	// add energies to emap
	emap[ mm_lj_inter_rep ] += total_rep;
	emap[ mm_lj_inter_atr ] += total_atr;
}

void
MMLJEnergyInter::eval_atom_derivative(
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
	if ( !pose.energies().use_nblist() ) {
		utility_exit_with_message("non-nblist minimize!");
	} else {

		// get atom1's neighbors
		scoring::Energies const & energies( pose.energies() );
		scoring::NeighborList const & nblist( energies.nblist( EnergiesCacheableDataType::MM_LJ_INTRA_NBLIST ) );
		scoring::AtomNeighbors const & nbrs( nblist.atom_neighbors( id ) );

		// iterate over all of atom1's neighbors using the neighbor list
		for ( auto const & nbr : nbrs ) {
			conformation::Atom const & atom2( pose.residue( nbr.rsd() ).atom( nbr.atomno() ) );

			// calculate f1 and f2
			Vector f1( atom1.xyz().cross( atom2.xyz() ) );
			Vector f2( atom1.xyz() - atom2.xyz() );

			// get count pair weight from neighbor list
			//Real const cp_weight( nbr.weight() );

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
	}
}

void
MMLJEnergyInter::eval_intrares_energy(
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap &
) const {}

void
MMLJEnergyInter::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & /*weights*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	debug_assert( set1.resid() != set2.resid() );

	using namespace methods;
	using namespace trie;
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table1( energy_table );
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table2( energy_table );

	temp_table1 = 0; temp_table2 = 0;

	//  // save weight information so that its available during tvt execution
	//  elec_weight_ = weights[ fa_elec ];

	RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( mm_lj_energy_inter_method ) ));
	RotamerTrieBaseCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( mm_lj_energy_inter_method ) ));

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp = get_count_pair_function_trie( set1, set2, pose, sfxn );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	trie1->trie_vs_trie( *trie2, *cp, *this, temp_table1, temp_table2 );

	/// add in the energies calculated by the tvt alg.
	energy_table += temp_table1;
	//std::cout << "FINISHED evaluate_rotamer_pair_energies" << std::endl;
}

void
MMLJEnergyInter::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & /*weights*/,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	using namespace methods;
	using namespace trie;

	// allocate space for the trie-vs-trie algorithm
	utility::vector1< core::PackerEnergy > temp_vector1( set.num_rotamers(), 0.0 );
	utility::vector1< core::PackerEnergy > temp_vector2( set.num_rotamers(), 0.0 );

	RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< core::scoring::trie::RotamerTrieBase const > ( set.get_trie( mm_lj_energy_inter_method ) ) );
	RotamerTrieBaseCOP trie2 = ( static_cast< TrieCollection const & >
		( pose.energies().data().get( EnergiesCacheableDataType::MM_LJ_TRIE_COLLECTION )) ).trie( residue.seqpos() );

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp = get_count_pair_function_trie( pose.residue( set.resid() ), residue, trie1, trie2, pose, sfxn );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	trie1->trie_vs_path( *trie2, *cp, *this, temp_vector1, temp_vector2 );

	/// add in the energies calculated by the tvt alg.
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		energy_vector[ ii ] += temp_vector1[ ii ];
	}
	//std::cout << "FINISHED evaluate_rotamer_background_energies" << std::endl;
}


/// @brief MMLJEnergy does not have an atomic interation threshold
Distance
MMLJEnergyInter::atomic_interaction_cutoff() const
{
	// this will probably screw up other stuff, but it might not since etable goes to zero at 5.5
	return 7.0;
}

/// @brief MMLJEnergy is context independent; indicates that no context graphs are required
void
MMLJEnergyInter::indicate_required_context_graphs(utility::vector1< bool > & ) const
{}

/// @brief required for neighbor list and to be more lke the ETable
etable::count_pair::CountPairFunctionCOP
MMLJEnergyInter::get_count_pair_function(
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
MMLJEnergyInter::get_count_pair_function(
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
MMLJEnergyInter::get_intrares_countpair(
	conformation::Residue const & res,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*sfxn*/
) const
{
	using namespace etable;
	using namespace count_pair;

	return CountPairFunctionOP( CountPairFactory::create_intrares_count_pair_function( res, CP_CROSSOVER_3 ) );
}

/// @brief figure out the trie count pair function to use
/// Need to refactor this so that the decision "what kind of count pair behavior should I use" can be decoupled
/// from class instantiation, and therefore shared between the creation of the trie count pair classes and the regular
/// count pair classes
trie::TrieCountPairBaseOP
MMLJEnergyInter::get_count_pair_function_trie(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) const
{
	conformation::Residue const & res1( pose.residue( set1.resid() ) );
	conformation::Residue const & res2( pose.residue( set2.resid() ) );

	trie::RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( methods::mm_lj_energy_inter_method ) ));
	trie::RotamerTrieBaseCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( methods::mm_lj_energy_inter_method ) ));

	return get_count_pair_function_trie( res1, res2, trie1, trie2, pose, sfxn );
}


trie::TrieCountPairBaseOP
MMLJEnergyInter::get_count_pair_function_trie(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	trie::RotamerTrieBaseCOP trie1,
	trie::RotamerTrieBaseCOP trie2,
	pose::Pose const &,
	ScoreFunction const &
) const
{
	using namespace etable::count_pair;
	using namespace trie;
	using namespace etable::etrie;

	TrieCountPairBaseOP tcpfxn;

	CPResidueConnectionType connection( CountPairFactory::determine_residue_connection( res1, res2 ) );
	Size conn1 = trie1->get_count_pair_data_for_residue( res2.seqpos() );
	Size conn2 = trie2->get_count_pair_data_for_residue( res1.seqpos() );

	if ( connection != CP_NO_BONDS ) {
		tcpfxn = TrieCountPairBaseOP( new TrieCountPair1BC3( conn1, conn2 ) );
	} else {
		tcpfxn = TrieCountPairBaseOP( new TrieCountPairAll );
	}
	return tcpfxn;

}

void
MMLJEnergyInter::bump_energy_full(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	EnergyMap tbemap;
	MMLJEnergyInter::residue_pair_energy( rsd1, rsd2, pose, sfxn, tbemap );
	emap[ mm_lj_inter_rep ] += tbemap[ mm_lj_inter_rep ];
	emap[ mm_lj_inter_atr ] += tbemap[ mm_lj_inter_atr ];
}

void
MMLJEnergyInter::bump_energy_backbone(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	// this is not correct, it should iterate over the backbone atoms only
	// but I wanna see if it works
	EnergyMap tbemap;
	MMLJEnergyInter::residue_pair_energy( rsd1, rsd2, pose, sfxn, tbemap );
	emap[ mm_lj_inter_rep ] += tbemap[ mm_lj_inter_rep ];
	emap[ mm_lj_inter_atr ] += tbemap[ mm_lj_inter_atr ];
}
core::Size
MMLJEnergyInter::version() const
{
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
// typedef core::scoring::trie::RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, core::scoring::etable::etrie::CountPairDataGeneric > MMEtableRotTrieGeneric;
//
// SAVE_AND_LOAD_SERIALIZABLE( MMEtableRotTrieGeneric );
// CEREAL_REGISTER_TYPE( MMEtableRotTrieGeneric )

typedef core::scoring::trie::RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, core::scoring::etable::etrie::CountPairData_1_1 > MMEtableRotTrie11; SAVE_AND_LOAD_SERIALIZABLE( MMEtableRotTrie11 );
CEREAL_REGISTER_TYPE( MMEtableRotTrie11 )

typedef core::scoring::trie::RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, core::scoring::etable::etrie::CountPairData_1_2 > MMEtableRotTrie12; SAVE_AND_LOAD_SERIALIZABLE( MMEtableRotTrie12 );
CEREAL_REGISTER_TYPE( MMEtableRotTrie12 )

typedef core::scoring::trie::RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, core::scoring::etable::etrie::CountPairData_1_3 > MMEtableRotTrie13; SAVE_AND_LOAD_SERIALIZABLE( MMEtableRotTrie13 );
CEREAL_REGISTER_TYPE( MMEtableRotTrie13 )


CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_methods_MMLJEnergyInter )
#endif // SERIALIZATION
