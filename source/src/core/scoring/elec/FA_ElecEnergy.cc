// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/FA_ElecEnergy.cc
/// @brief  Electrostatic energy with a distance-dependant dielectric
/// @author Phil Bradley
/// @author Andrew Leaver-Fay
/// @author Modified by James Gleixner and Liz Kellogg
/// @author Modified by Vikram K. Mulligan (vmullig@uw.edu) -- added data caching.

// Unit headers
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/elec/FA_ElecEnergyCreator.hh>

// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/scoring/hbonds/HBondSet.hh>

// Project headers
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/pose/Pose.hh>

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
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>

#include <core/scoring/etable/etrie/TrieCountPair1BC4.hh>
#include <core/scoring/etable/etrie/TrieCountPairAll.hh>
#include <core/scoring/etable/etrie/TrieCountPairNone.hh>
#include <core/scoring/etable/etrie/TrieCountPairGeneric.hh>

#include <core/scoring/hbonds/hbonds.hh>

#include <core/conformation/residue_datacache.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#ifdef SERIALIZATION
// Project serialization headers
#include <core/scoring/trie/RotamerTrie.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>

#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static THREAD_LOCAL basic::Tracer TR( "core.scoring.elec.FA_ElecEnergy" );

/////////////////////////////////////////////////////////////////////////////////////////
///
/// Hacky (hence the name) implementation of distance dependant dielectric electrostatics,
/// with near and far distance cutoffs.


//     alternatives: WARSHEL (from ligand.cc)
//     E = 322.0637*q1*q2/r/e(r)
//     if ( r < 3 ) e(r) = 16.55
//     else         e(r) = 1 + 60*(1-exp(-0.1*r))
//     Warshel, A. Russell, S. T., Q. Rev. Biophys., 1984, 17, 283-422
//
//
//     sigmoidal dielectric: (hingerty 1985)
//
//     e(r) = D - (D-D0)/2 [ (rS)**2 + 2rS + 2]exp(-rS)
//     with eg:
//     D = 78, D0 = 1, S = 0.3565 (rouzina&bloomfield)
//     D = 80, D0 = 4, S = 0.4 (rohs)

namespace core {
namespace scoring {
namespace elec {


/// @details This must return a fresh instance of the FA_ElecEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
FA_ElecEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new FA_ElecEnergy( options ) );
}

ScoreTypes
FA_ElecEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_elec );
	sts.push_back( fa_elec_bb_bb );
	sts.push_back( fa_elec_bb_sc );
	sts.push_back( fa_elec_sc_sc );
	sts.push_back( fa_intra_elec );
	return sts;
}
////////////////////////////////////////////////////////////////////////////
FA_ElecEnergy::FA_ElecEnergy( methods::EnergyMethodOptions const & options ):
	parent( methods::EnergyMethodCreatorOP( new FA_ElecEnergyCreator ) ),
	coulomb_( options ),
	exclude_protein_protein_( options.exclude_protein_protein_fa_elec() ),
	exclude_monomer_( options.exclude_monomer_fa_elec() ),
	exclude_DNA_DNA_( options.exclude_DNA_DNA() ),
	eval_intrares_ST_only_( options.eval_intrares_elec_ST_only() )
{
	initialize();
}


////////////////////////////////////////////////////////////////////////////
FA_ElecEnergy::FA_ElecEnergy( FA_ElecEnergy const & src ):
	parent( src ),
	coulomb_( src.coulomb() ),
	exclude_protein_protein_( src.exclude_protein_protein_ ),
	exclude_monomer_( src.exclude_monomer_ ),
	exclude_DNA_DNA_( src.exclude_DNA_DNA_ ),
	eval_intrares_ST_only_( src.eval_intrares_ST_only_ )
{
	initialize();
}


void
FA_ElecEnergy::initialize() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	coulomb_.initialize();

	// read countpair tables from DB
	use_cp_rep_ = option[ score::elec_representative_cp ]() || option[ score::elec_representative_cp_flip ]();
	flip_cp_rep_ = option[ score::elec_representative_cp_flip ]();
	if ( use_cp_rep_ ) {
		get_cp_tables(); //Lazily loaded once and only once.
	}
}

void
FA_ElecEnergy::get_cp_tables() {
	cp_rep_map_byname_ = core::scoring::ScoringManager::get_instance()->get_cp_rep_map_byname();
}

core::Size
FA_ElecEnergy::get_countpair_representative_atom(
	core::chemical::ResidueType const & restype,
	core::Size atm_i
) const {
	if ( !use_cp_rep_ ) return atm_i;

	// for now ...
	if ( !restype.is_protein() ) return atm_i;

	std::map< chemical::ResidueType const *, std::map<core::Size,core::Size> >::const_iterator iter = cp_rep_map_.find( & restype );
	if ( iter == cp_rep_map_.end() ) {
		// make parameters
		std::map<core::Size,core::Size> rsd_map;

		std::map< std::string, std::map<std::string,std::string> >::const_iterator name_iter = cp_rep_map_byname_->find( restype.name3() );

		if ( name_iter == cp_rep_map_byname_->end() ) {
			TR.Trace << "Warning!  Unable to find countpair representatives for restype " << restype.name3() << std::endl;
		} else {
			std::map<std::string,std::string> const & atms = name_iter->second;
			for ( std::map<std::string,std::string>::const_iterator atom_iter = atms.begin(), atom_iter_end = atms.end(); atom_iter!=atom_iter_end; ++atom_iter ) {
				if ( restype.has(atom_iter->first) && restype.has(atom_iter->second) ) {
					core::Size idx1 = restype.atom_index(atom_iter->first);
					core::Size idx2 = restype.atom_index(atom_iter->second);
					rsd_map.insert( std::make_pair(idx1,idx2) );
				} else {
					TR.Trace << "Warning!  Unable to find atompair " << atom_iter->first << "," << atom_iter->second  << " for " << restype.name3() << " (" << restype.name() << ")" << std::endl;
				}
			}
		}

		cp_rep_map_[ & restype ] = rsd_map;

		iter = cp_rep_map_.find( & restype );
	}

	std::map<core::Size,core::Size> const &mapping_i = iter->second;
	std::map<core::Size,core::Size>::const_iterator iter_map_i = mapping_i.find(atm_i);
	if ( iter_map_i == mapping_i.end() ) {
		return atm_i;
	} else {
		return iter_map_i->second;
	}
}


/// clone
methods::EnergyMethodOP
FA_ElecEnergy::clone() const
{
	return methods::EnergyMethodOP( new FA_ElecEnergy( *this ) );
}

void
FA_ElecEnergy::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_nres_mono(pose);

	if ( ! pose.energies().use_nblist() ) return;

	// stash our nblist inside the pose's energies object
	Energies & energies( pose.energies() );

	// setup the atom-atom nblist
	NeighborListOP nblist;
	Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? option[ run::nblist_autoupdate_narrow ] : 1.5;
	Real const XX = coulomb().max_dis() + 2 * tolerated_motion;
	nblist = NeighborListOP( new NeighborList( min_map.domain_map(), XX*XX, XX*XX, XX*XX) );
	if ( pose.energies().use_nblist_auto_update() ) {
		nblist->set_auto_update( tolerated_motion );
	}
	// this partially becomes the EtableEnergy classes's responsibility
	nblist->setup( pose, sfxn, *this);
	energies.set_nblist( EnergiesCacheableDataType::ELEC_NBLIST, nblist );
}


void
FA_ElecEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & /*sfxn*/ ) const
{
	set_nres_mono(pose); // why?
	pose.update_residue_neighbors();
}


void
FA_ElecEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const
{
	set_nres_mono(pose);
	pose.update_residue_neighbors();
	if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::ELEC_NBLIST ) );
		nblist.prepare_for_scoring( pose, scfxn, *this );
	}
}


// The FA_ElectEnergy method stores a vector of rotamer trie objects in the Energies
// object for use in rapid rotamer/background energy calculations.  Overrides default
// do-nothing behavior.
void
FA_ElecEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	using namespace trie;

	set_nres_mono(pose);

	TrieCollectionOP tries( new TrieCollection );
	tries->total_residue( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		// Do not compute energy for virtual residues.
		if ( pose.residue(ii).aa() == core::chemical::aa_vrt ) continue;

		// ensure that the resmaps are initialized for this
		if ( use_cp_rep_ ) get_countpair_representative_atom(pose.residue( ii ).type(),1);

		RotamerTrieBaseOP one_rotamer_trie = create_rotamer_trie( pose.residue( ii ), pose );
		tries->trie( ii, one_rotamer_trie );
	}
	pose.energies().data().set( EnergiesCacheableDataType::ELEC_TRIE_COLLECTION, tries );
}

// @brief Creates a rotamer trie for the input set of rotamers and stores the trie
// in the rotamer set.
void
FA_ElecEnergy::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & set
) const
{
	// ensure that the resmaps are initialized for everyrestype in the set
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		get_countpair_representative_atom( set.rotamer( ii )->type(),1);
	}

	trie::RotamerTrieBaseOP rottrie = create_rotamer_trie( set, pose );
	set.store_trie( methods::elec_method, rottrie );
}


// @brief Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
FA_ElecEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	conformation::Residue const & rsd( pose.residue( resid ) );

	using namespace trie;
	trie::RotamerTrieBaseOP one_rotamer_trie = create_rotamer_trie( rsd, pose );

	// grab non-const & of the cached tries and replace resid's trie with a new one.
	TrieCollection & trie_collection
		( static_cast< TrieCollection & > (pose.energies().data().get( EnergiesCacheableDataType::ELEC_TRIE_COLLECTION )));
	trie_collection.trie( resid, one_rotamer_trie );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
FA_ElecEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & /*sf*/,
	EnergyMap & emap
) const
{
	if ( pose.energies().use_nblist() ) return;
	using namespace etable::count_pair;

	Real score(0.0);

	Real attached_h_max_dis2 = hydrogen_interaction_cutoff2();

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;

	// NULL if no info
	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		// assuming only a single bond right now -- generalizing to arbitrary topologies
		// also assuming crossover of 4, should be closest (?) to classic rosetta
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		Real d2;
		for ( Size ii = 1; ii <= rsd1.nheavyatoms(); ++ii ) {
			Size ii_rep = get_countpair_representative_atom( rsd1.type(), ii );
			for ( Size jj = 1; jj <= rsd2.nheavyatoms(); ++jj ) {
				Size jj_rep = get_countpair_representative_atom( rsd2.type(), jj );

				Real weight=1.0;
				Size path_dist( 0 );
				if ( cpfxn->count( ii_rep, jj_rep, weight, path_dist ) ) {
					score += score_atom_pair( rsd1, rsd2, ii, jj, emap, weight, d2 );
				} else {
					d2 = rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) );
				}

				if ( d2 > attached_h_max_dis2 ) continue;

				Size ii_hatbegin( rsd1.attached_H_begin()[ ii ] ), ii_hatend( rsd1.attached_H_end()[ ii ] );
				Size jj_hatbegin( rsd2.attached_H_begin()[ jj ] ), jj_hatend( rsd2.attached_H_end()[ jj ] );
				for ( Size kk = ii_hatbegin; kk <= ii_hatend; ++kk ) {
					Size kk_rep = get_countpair_representative_atom( rsd1.type(), kk );

					weight = 1.0;
					path_dist = 0;
					if ( cpfxn->count( kk_rep, jj_rep, weight, path_dist ) ) {
						score += score_atom_pair( rsd1, rsd2, kk, jj, emap, weight, d2 );
					}
				}
				for ( Size kk = jj_hatbegin; kk <= jj_hatend; ++kk ) {
					Size kk_rep = get_countpair_representative_atom( rsd2.type(), kk );

					weight = 1.0;
					path_dist = 0;
					if ( cpfxn->count( ii_rep, kk_rep, weight, path_dist ) ) {
						score += score_atom_pair( rsd1, rsd2, ii, kk, emap, weight, d2 );
					}
				}
				for ( Size kk = ii_hatbegin; kk <= ii_hatend; ++kk ) {
					Size kk_rep = get_countpair_representative_atom( rsd1.type(), kk );
					for ( Size ll = jj_hatbegin; ll <= jj_hatend; ++ll ) {
						Size ll_rep = get_countpair_representative_atom( rsd2.type(), ll );

						weight = 1.0;
						path_dist = 0;
						if ( cpfxn->count( kk_rep, ll_rep, weight, path_dist ) ) {
							score += score_atom_pair( rsd1, rsd2, kk, ll, emap, weight, d2 );
						}
					}
				}
			}
		}
	} else {
		Real d2;
		for ( Size ii = 1; ii <= rsd1.nheavyatoms(); ++ii ) {
			for ( Size jj = 1; jj <= rsd2.nheavyatoms(); ++jj ) {
				Real weight=1.0;
				score += score_atom_pair( rsd1, rsd2, ii, jj, emap, weight, d2 );

				if ( d2 > attached_h_max_dis2 ) continue;

				Size ii_hatbegin( rsd1.attached_H_begin()[ ii ] ), ii_hatend( rsd1.attached_H_end()[ ii ] );
				Size jj_hatbegin( rsd2.attached_H_begin()[ jj ] ), jj_hatend( rsd2.attached_H_end()[ jj ] );
				for ( Size kk = ii_hatbegin; kk <= ii_hatend; ++kk ) {
					score += score_atom_pair( rsd1, rsd2, kk, jj, emap, weight, d2 );
				}
				for ( Size kk = jj_hatbegin; kk <= jj_hatend; ++kk ) {
					score += score_atom_pair( rsd1, rsd2, ii, kk, emap, weight, d2 );
				}
				for ( Size kk = ii_hatbegin; kk <= ii_hatend; ++kk ) {
					for ( Size ll = jj_hatbegin; ll <= jj_hatend; ++ll ) {
						score += score_atom_pair( rsd1, rsd2, kk, ll, emap, weight, d2 );
					}
				}
			}
		}
	}
	emap[ fa_elec ] += score;
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;
}

void
FA_ElecEnergy::eval_intrares_energy(
	conformation::Residue const &rsd,
	pose::Pose const &/*pose*/,
	ScoreFunction const &sf,
	EnergyMap &emap ) const {
	if ( sf.get_weight( fa_intra_elec ) == 0 ) return;
	using namespace etable::count_pair;

	Real score(0.0);

	if ( eval_intrares_ST_only() &&
			!(rsd.aa() == chemical::aa_ser || rsd.aa() == chemical::aa_thr ||
			rsd.aa() == chemical::aa_dse || rsd.aa() == chemical::aa_dth )
			) return;

	Size iOG=0, iHG=0, iN=0;
	utility::vector1< Size > iHs;
	if ( eval_intrares_ST_only() ) {
		iN = rsd.atom_index( "N" );
		if ( rsd.is_lower_terminus() ) {
			debug_assert( rsd.has("1H") && rsd.has("2H") && rsd.has("3H") );
			iHs.resize(3);
			iHs[1] = rsd.atom_index( "1H" );
			iHs[2] = rsd.atom_index( "2H" );
			iHs[3] = rsd.atom_index( "3H" );
		} else {
			iHs.push_back( rsd.atom_index( "H" ) );
		}
		if ( rsd.aa() == chemical::aa_ser || rsd.aa() == chemical::aa_dse ) {
			debug_assert( rsd.has("HG") );
			iOG = rsd.atom_index( "OG" );
			iHG = rsd.atom_index( "HG" );
		} else { // threonine
			debug_assert( rsd.has("HG1") );
			iOG = rsd.atom_index( "OG1" );
			iHG = rsd.atom_index( "HG1" );
		}
	}

	// assuming only a single bond right now -- generalizing to arbitrary topologies
	// also assuming crossover of 4, should be closest (?) to classic rosetta
	CountPairFunctionOP cpfxn =
		CountPairFactory::create_intrares_count_pair_function( rsd, CP_CROSSOVER_4 );
	Real d2;
	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		Size ii_rep = get_countpair_representative_atom( rsd.type(), ii );
		for ( Size jj = ii+1; jj <= rsd.natoms(); ++jj ) {
			Size jj_rep = get_countpair_representative_atom( rsd.type(), jj );

			if ( eval_intrares_ST_only() ) {
				if ( !( ( ii == iN && jj == iOG ) || ( ii == iN && jj == iHG ) ||
						( ii == iOG && iHs.contains(jj) ) || ( iHs.contains(ii) && jj == iHG ) )
						) continue;
			}

			Real weight=1.0;

			Size path_dist( 0 );
			if ( cpfxn->count( ii_rep, jj_rep, weight, path_dist ) ) {
				score += score_atom_pair( rsd, rsd, ii, jj, emap, weight, d2 );
			}
		}
	}

	emap[ fa_intra_elec ] += score;
	//std::cout << rsd.seqpos() << ' ' << score << std::endl;
}


bool
FA_ElecEnergy::minimize_in_whole_structure_context( pose::Pose const & pose ) const
{
	return pose.energies().use_nblist_auto_update();
}


bool
FA_ElecEnergy::defines_score_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool res_moving_wrt_eachother
) const
{
	if ( rsd1.seqpos() == rsd2.seqpos() ) {
		return false;
	} else if ( exclude_protein_protein_ && rsd1.is_protein() && rsd2.is_protein() ) {
		return false;
	} else if ( exclude_monomer_ && monomer_test( rsd1.seqpos(), rsd2.seqpos()) ) {
		return false;
	} else if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) {
		return false;
	}

	return res_moving_wrt_eachother;
}


bool
FA_ElecEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}


void
FA_ElecEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const & /*sf*/,
	EnergyMap & emap
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( elec_pair_nblist ) ));
	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( elec_pair_nblist ) ) );
	Real dsq, score( 0.0 );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		score += score_atom_pair( rsd1, rsd2, neighbs[ ii ].atomno1(), neighbs[ ii ].atomno2(), emap, neighbs[ ii ].weight(), dsq );
	}
	emap[ fa_elec ] += score;
}

void
FA_ElecEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const & /*rsd*/,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*scorefxn*/,
	kinematics::MinimizerMapBase const & /*min_map*/,
	ResSingleMinimizationData & /*resdata*/
) const
{
}

void
FA_ElecEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & /*scorefxn*/,
	kinematics::MinimizerMapBase const & /*minmap*/,
	ResSingleMinimizationData const & /*res1data*/,
	ResSingleMinimizationData const & /*res2data*/,
	ResPairMinimizationData & pair_data
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( pose.energies().use_nblist_auto_update() ) return;

	etable::count_pair::CountPairFunctionCOP count_pair = get_count_pair_function( rsd1, rsd2 );
	debug_assert( rsd1.seqpos() < rsd2.seqpos() );

	// update the existing nblist if it's already present in the min_data object
	ResiduePairNeighborListOP nblist(
		utility::pointer::static_pointer_cast< core::scoring::ResiduePairNeighborList > ( pair_data.get_data( elec_pair_nblist ) ));
	if ( ! nblist ) nblist = ResiduePairNeighborListOP( new ResiduePairNeighborList );

	Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
	Real const XX2 = std::pow( coulomb().max_dis() + 2*tolerated_narrow_nblist_motion, 2 );

	if ( !use_cp_rep_ ) {
		nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair );
	} else {
		// ensure that the resmaps are initialized for this pair
		get_countpair_representative_atom(rsd1.type(),1);
		get_countpair_representative_atom(rsd2.type(),1);

		nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair, cp_rep_map_[&(rsd1.type())], cp_rep_map_[&(rsd2.type())] );
	}

	pair_data.set_data( elec_pair_nblist, nblist );
}

bool
FA_ElecEnergy::requires_a_setup_for_scoring_for_residue_opportunity_during_minimization( pose::Pose const & ) const
{
	return true;
}

void
FA_ElecEnergy::setup_for_scoring_for_residue(
	conformation::Residue const & /*rsd*/,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*sfxn*/,
	ResSingleMinimizationData & /*resdata*/
) const
{
}

bool
FA_ElecEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const &  ) const
{
	return true;
}

void
FA_ElecEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data,
	basic::datacache::BasicDataCache &
) const
{
	setup_for_scoring_for_residue( rsd, pose, sfxn, min_data );
}

void
FA_ElecEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & /*res1data*/,
	ResSingleMinimizationData const & /*res2data*/,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( elec_pair_nblist ) ));

	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( elec_pair_nblist ) ) );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	weight_triple wtrip;
	setup_weight_triple( weights, wtrip );
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		Size at1 = neighbs[ ii ].atomno1();
		Size at2 = neighbs[ ii ].atomno2();
		Vector const & atom1xyz( rsd1.xyz( at1 ) );
		Vector const & atom2xyz( rsd2.xyz( at2 ) );

		Real const at1_charge( rsd1.atomic_charge( at1 ) );
		Real const at2_charge( rsd2.atomic_charge( at2 ) );

		Vector f1( 0.0 ), f2 = ( atom1xyz - atom2xyz );
		Real const dis2( f2.length_squared() );
		Real dE_dr_over_r = neighbs[ ii ].weight() * coulomb().eval_dfa_elecE_dr_over_r( dis2, at1_charge, at2_charge );
		Real sfxn_weight = elec_weight(
			rsd1.atom_is_backbone( at1 ),
			rsd2.atom_is_backbone( at2 ),
			wtrip );

		if ( dE_dr_over_r != 0.0 ) {
			f1 = atom1xyz.cross( atom2xyz );
			Vector f1s = dE_dr_over_r * sfxn_weight * f1;
			Vector f2s = dE_dr_over_r * sfxn_weight * f2;
			r1_atom_derivs[ at1 ].f1() += f1s;
			r1_atom_derivs[ at1 ].f2() += f2s;
			r2_atom_derivs[ at2 ].f1() -= f1s;
			r2_atom_derivs[ at2 ].f2() -= f2s;
		}
	}
}


void
FA_ElecEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & /*min_data*/,
	pose::Pose const & /*pose*/,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	using namespace etable::count_pair;
	if ( weights[ fa_intra_elec ] == 0 ) return;
	Real sfxn_weight = weights[fa_intra_elec];

	if ( eval_intrares_ST_only() &&
			!(rsd.aa() == chemical::aa_ser || rsd.aa() == chemical::aa_thr ||
			rsd.aa() == chemical::aa_dse || rsd.aa() == chemical::aa_dth )
			) return;

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_intrares_count_pair_function( rsd, CP_CROSSOVER_4 );

	Size iN=0, iHG=0, iOG=0;
	utility::vector1< Size > iHs;
	if ( eval_intrares_ST_only() ) {
		iN = rsd.atom_index( "N" );
		if ( rsd.is_lower_terminus() ) {
			debug_assert( rsd.has("1H") && rsd.has("2H") && rsd.has("3H") );
			iHs.resize(3);
			iHs[1] = rsd.atom_index( "1H" );
			iHs[2] = rsd.atom_index( "2H" );
			iHs[3] = rsd.atom_index( "3H" );
		} else {
			iHs.push_back( rsd.atom_index( "H" ) );
		}
		if ( rsd.aa() == chemical::aa_ser || rsd.aa() == chemical::aa_dse ) {
			debug_assert( rsd.has("HG") );
			iOG = rsd.atom_index( "OG" );
			iHG = rsd.atom_index( "HG" );
		} else { // threonine/D-threonine
			debug_assert( rsd.has("HG1") );
			iOG = rsd.atom_index( "OG1" );
			iHG = rsd.atom_index( "HG1" );
		}
	}

	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		Size ii_rep = get_countpair_representative_atom( rsd.type(), ii );
		for ( Size jj = ii+1; jj <= rsd.natoms(); ++jj ) {
			if ( eval_intrares_ST_only() ) {
				if ( !( ( ii == iN && jj == iOG ) || ( ii == iN && jj == iHG ) ||
						( ii == iOG && iHs.contains(jj)) || ( iHs.contains(ii) && jj == iHG ) )
						) continue;
			}

			Size jj_rep = get_countpair_representative_atom( rsd.type(), jj );
			Vector const & atom1xyz( rsd.xyz( ii ) );
			Vector const & atom2xyz( rsd.xyz( jj ) );

			Real const at1_charge( rsd.atomic_charge( ii ) );
			Real const at2_charge( rsd.atomic_charge( jj ) );

			Vector f2 = ( atom1xyz - atom2xyz );
			Real const dis2( f2.length_squared() );
			Real dE_dr_over_r = 0;
			Size path_dist( 0 );
			Real weight=1.0;
			if ( !cpfxn->count( ii_rep, jj_rep, weight, path_dist ) ) continue;
			dE_dr_over_r = weight * coulomb().eval_dfa_elecE_dr_over_r( dis2, at1_charge, at2_charge );

			if ( dE_dr_over_r == 0.0 ) continue;

			Vector f1 = atom1xyz.cross( atom2xyz );
			f1 *= dE_dr_over_r * sfxn_weight;
			f2 *= dE_dr_over_r * sfxn_weight;
			atom_derivs[ ii ].f1() += f1;
			atom_derivs[ ii ].f2() += f2;
			atom_derivs[ jj ].f1() -= f1;
			atom_derivs[ jj ].f2() -= f2;
		}
	}
}



/// @details for use only with the nblist auto-update algorithm
void
FA_ElecEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,// domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using namespace etable::count_pair;
	if ( ! pose.energies().use_nblist_auto_update() ) return;

	// what is my charge?
	Size const i( atom_id.rsd() );
	Size const ii( atom_id.atomno() );
	conformation::Residue const & irsd( pose.residue( i ) );
	Real const ii_charge( irsd.atomic_charge( ii ) );

	if ( ii_charge == 0.0 ) return;

	Vector const & ii_xyz( irsd.xyz(ii) );
	bool const ii_isbb( irsd.atom_is_backbone( ii ) );

	debug_assert( pose.energies().use_nblist() );
	NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::ELEC_NBLIST ) );
	AtomNeighbors const & nbrs( nblist.atom_neighbors(i,ii) );

	weight_triple wtrip;
	setup_weight_triple( weights, wtrip );

	for ( scoring::AtomNeighbors::const_iterator it2=nbrs.begin(),
			it2e=nbrs.end(); it2 != it2e; ++it2 ) {
		scoring::AtomNeighbor const & nbr( *it2 );
		Size const j( nbr.rsd() );
		Size const jj( nbr.atomno() );
		conformation::Residue const & jrsd( pose.residue( j ) );

		Real const jj_charge( jrsd.atomic_charge(jj) );
		if ( jj_charge == 0.0 ) continue; /// should prune out such atoms when constructing the neighborlist!
		Vector const & jj_xyz( jrsd.xyz( jj ) );
		Vector f2 = ( ii_xyz - jj_xyz );
		Real const dis2( f2.length_squared() );
		Real const dE_dr_over_r = nbr.weight() * coulomb().eval_dfa_elecE_dr_over_r( dis2, ii_charge, jj_charge );
		if ( dE_dr_over_r == 0.0 ) continue;

		Real sfxn_weight = elec_weight( ii_isbb, jrsd.atom_is_backbone( jj ), wtrip );
		Vector f1 = ii_xyz.cross( jj_xyz );
		f1 *= dE_dr_over_r * sfxn_weight;
		f2 *= dE_dr_over_r * sfxn_weight;
		F1 += f1;
		F2 += f2;
	}
}

void
FA_ElecEnergy::backbone_backbone_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*scorefxn*/,
	EnergyMap & emap
) const
{
	using namespace etable::count_pair;
	using namespace chemical;

	Real score(0.0);

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;

	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		// assuming only a single bond right now -- generalizing to arbitrary topologies
		// also assuming crossover of 4, should be closest (?) to classic rosetta
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		AtomIndices const & rsd1_bb_atoms( rsd1.type().all_bb_atoms() );
		AtomIndices const & rsd2_bb_atoms( rsd2.type().all_bb_atoms() );

		for ( Size ii=1, ii_end = rsd1_bb_atoms.size(); ii<= ii_end; ++ii ) {
			Size const i = rsd1_bb_atoms[ ii ];
			Size i_rep = get_countpair_representative_atom( rsd1.type(), i );
			Vector const & i_xyz( rsd1.xyz(i) );
			Real const i_charge( rsd1.atomic_charge(i) );
			if ( i_charge == 0.0 ) continue;
			for ( Size jj=1, jj_end = rsd2_bb_atoms.size(); jj<= jj_end; ++jj ) {
				Size const j = rsd2_bb_atoms[ jj ];
				Size j_rep = get_countpair_representative_atom( rsd2.type(), j );
				Real const j_charge( rsd2.atomic_charge(j) );
				if ( j_charge == 0.0 ) continue;
				Real weight(1.0);
				Size path_dist( 0 );
				if ( cpfxn->count( i_rep, j_rep, weight, path_dist ) ) {
					score += weight *
						coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), j_charge);
				}
			}
		}

	} else {
		// no countpair!
		AtomIndices const & rsd1_bb_atoms( rsd1.type().all_bb_atoms() );
		AtomIndices const & rsd2_bb_atoms( rsd2.type().all_bb_atoms() );
		for ( Size ii=1, ii_end = rsd1_bb_atoms.size(); ii<= ii_end; ++ii ) {
			Size const i = rsd1_bb_atoms[ ii ];
			Vector const & i_xyz( rsd1.xyz(i) );
			Real const i_charge( rsd1.atomic_charge(i) );
			if ( i_charge == 0.0 ) continue;
			for ( Size jj=1, jj_end = rsd2_bb_atoms.size(); jj<= jj_end; ++jj ) {
				Size const j = rsd2_bb_atoms[ jj ];
				Real const j_charge( rsd2.atomic_charge(j) );
				if ( j_charge == 0.0 ) continue;
				score += coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), j_charge );
			}
		}
	}
	emap[ fa_elec_bb_bb ] += score;
	emap[ fa_elec ] += score;
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;
}

void
FA_ElecEnergy::backbone_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*scorefxn*/,
	EnergyMap & emap
) const
{
	using namespace etable::count_pair;
	using namespace chemical;

	Real score(0.0);

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;

	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		// assuming only a single bond right now -- generalizing to arbitrary topologies
		// also assuming crossover of 4, should be closest (?) to classic rosetta
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		AtomIndices const & rsd1_bb_atoms( rsd1.type().all_bb_atoms() );
		AtomIndices const & rsd2_sc_atoms( rsd2.type().all_sc_atoms() );

		for ( Size ii=1, ii_end = rsd1_bb_atoms.size(); ii<= ii_end; ++ii ) {
			Size const i = rsd1_bb_atoms[ ii ];
			Size i_rep = get_countpair_representative_atom( rsd1.type(), i );
			Vector const & i_xyz( rsd1.xyz(i) );
			Real const i_charge( rsd1.atomic_charge(i) );
			if ( i_charge == 0.0 ) continue;
			for ( Size jj=1, jj_end = rsd2_sc_atoms.size(); jj<= jj_end; ++jj ) {
				Size const j = rsd2_sc_atoms[ jj ];
				Size j_rep = get_countpair_representative_atom( rsd2.type(), j );
				Real const j_charge( rsd2.atomic_charge(j) );
				if ( j_charge == 0.0 ) continue;
				Real weight(1.0);
				Size path_dist( 0 );
				if ( cpfxn->count( i_rep, j_rep, weight, path_dist ) ) {
					score += weight *
						coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), j_charge);
				}
			}
		}

	} else {
		// no countpair!
		AtomIndices const & rsd1_bb_atoms( rsd1.type().all_bb_atoms() );
		AtomIndices const & rsd2_sc_atoms( rsd2.type().all_sc_atoms() );
		for ( Size ii=1, ii_end = rsd1_bb_atoms.size(); ii<= ii_end; ++ii ) {
			Size const i = rsd1_bb_atoms[ ii ];
			Vector const & i_xyz( rsd1.xyz(i) );
			Real const i_charge( rsd1.atomic_charge(i) );
			if ( i_charge == 0.0 ) continue;
			for ( Size jj=1, jj_end = rsd2_sc_atoms.size(); jj<= jj_end; ++jj ) {
				Size const j = rsd2_sc_atoms[ jj ];
				Real const j_charge( rsd2.atomic_charge(j) );
				if ( j_charge == 0.0 ) continue;
				score += coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), j_charge );
			}
		}
	}
	emap[ fa_elec_bb_sc ] += score;
	emap[ fa_elec ] += score;
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;

}


void
FA_ElecEnergy::sidechain_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*scorefxn*/,
	EnergyMap & emap
) const
{
	using namespace etable::count_pair;
	using namespace chemical;

	Real score(0.0);

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;

	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		// assuming only a single bond right now -- generalizing to arbitrary topologies
		// also assuming crossover of 4, should be closest (?) to classic rosetta
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		AtomIndices const & rsd1_sc_atoms( rsd1.type().all_sc_atoms() );
		AtomIndices const & rsd2_sc_atoms( rsd2.type().all_sc_atoms() );

		for ( Size ii=1, ii_end = rsd1_sc_atoms.size(); ii<= ii_end; ++ii ) {
			Size const i = rsd1_sc_atoms[ ii ];
			Size i_rep = get_countpair_representative_atom( rsd1.type(), i );
			Vector const & i_xyz( rsd1.xyz(i) );
			Real const i_charge( rsd1.atomic_charge(i) );
			if ( i_charge == 0.0 ) continue;
			for ( Size jj=1, jj_end = rsd2_sc_atoms.size(); jj<= jj_end; ++jj ) {
				Size const j = rsd2_sc_atoms[ jj ];
				Size j_rep = get_countpair_representative_atom( rsd2.type(), j );
				Real const j_charge( rsd2.atomic_charge(j) );
				if ( j_charge == 0.0 ) continue;
				Real weight(1.0);
				Size path_dist( 0 );
				if ( cpfxn->count( i_rep, j_rep, weight, path_dist ) ) {
					score += weight *
						coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), j_charge);
				}
			}
		}

	} else {
		// no countpair!
		AtomIndices const & rsd1_sc_atoms( rsd1.type().all_sc_atoms() );
		AtomIndices const & rsd2_sc_atoms( rsd2.type().all_sc_atoms() );
		for ( Size ii=1, ii_end = rsd1_sc_atoms.size(); ii<= ii_end; ++ii ) {
			Size const i = rsd1_sc_atoms[ ii ];
			Vector const & i_xyz( rsd1.xyz(i) );
			Real const i_charge( rsd1.atomic_charge(i) );
			if ( i_charge == 0.0 ) continue;
			for ( Size jj=1, jj_end = rsd2_sc_atoms.size(); jj<= jj_end; ++jj ) {
				Size const j = rsd2_sc_atoms[ jj ];
				Real const j_charge( rsd2.atomic_charge(j) );
				if ( j_charge == 0.0 ) continue;
				score += coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), j_charge );
			}
		}
	}
	emap[ fa_elec_sc_sc ] += score;
	emap[ fa_elec ] += score;
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;

}

void
FA_ElecEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	if ( ! pose.energies().use_nblist() || ! pose.energies().use_nblist_auto_update() ) return;

	EnergyMap tbenergy_map;
	// add in contributions from the nblist atom-pairs
	NeighborList const & nblist
		( pose.energies().nblist( EnergiesCacheableDataType::ELEC_NBLIST ) );

	nblist.check_domain_map( pose.energies().domain_map() );
	utility::vector1< conformation::Residue const * > resvect;
	resvect.reserve( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		resvect.push_back( & pose.residue( ii ) );
	}

	Real bb_sc_scores[ 3 ] = {0.0, 0.0, 0.0};
	Real total_score( 0.0 );
	for ( Size i=1, i_end = pose.size(); i<= i_end; ++i ) {
		conformation::Residue const & ires( *resvect[i] );
		for ( Size ii=1, ii_end=ires.natoms(); ii<= ii_end; ++ii ) {
			AtomNeighbors const & nbrs( nblist.upper_atom_neighbors(i,ii) );
			int ii_isbb = ires.atom_is_backbone( ii );
			for ( AtomNeighbors::const_iterator nbr_iter=nbrs.begin(),
					nbr_end=nbrs.end(); nbr_iter!= nbr_end; ++nbr_iter ) {
				AtomNeighbor const & nbr( *nbr_iter );

				Size const  j( nbr.rsd() );
				Size const jj( nbr.atomno() );

				conformation::Residue const & jres( *resvect[j] );
				int jj_isbb = jres.atom_is_backbone( jj );

				debug_assert( ii_isbb + jj_isbb >= 0 && ii_isbb + jj_isbb < 3 );  //?

				Real score = nbr.weight() *
					coulomb().eval_atom_atom_fa_elecE( ires.xyz(ii), ires.atomic_charge(ii), jres.xyz(jj), jres.atomic_charge(jj) );

				bb_sc_scores[ ii_isbb + jj_isbb ] += score;
				total_score += score;
			}
		}
	}
	totals[ fa_elec ] += total_score;
	totals[ fa_elec_bb_bb ] += bb_sc_scores[ 2 ];
	totals[ fa_elec_bb_sc ] += bb_sc_scores[ 1 ];
	totals[ fa_elec_sc_sc ] += bb_sc_scores[ 0 ];
}

void
FA_ElecEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & , //weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	debug_assert( set1.resid() != set2.resid() );

	// Since a rotamer set may include multiple residue types,
	// we'll make our decision based on what's currently in the Pose.
	if ( ! defines_score_for_residue_pair(pose.residue(set1.resid()), pose.residue(set2.resid()), true) ) return;

	using namespace methods;
	using namespace trie;
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table1( energy_table );
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table2( energy_table );

	temp_table1 = 0; temp_table2 = 0;

	// save weight information so that its available during tvt execution
	core::Real wt_bb_bb = sfxn.weights()[ fa_elec ] + sfxn.weights()[ fa_elec_bb_bb ];
	core::Real wt_bb_sc = sfxn.weights()[ fa_elec ] + sfxn.weights()[ fa_elec_bb_sc ];
	core::Real wt_sc_bb = sfxn.weights()[ fa_elec ] + sfxn.weights()[ fa_elec_bb_sc ];
	core::Real wt_sc_sc = sfxn.weights()[ fa_elec ] + sfxn.weights()[ fa_elec_sc_sc ];

	/// this will later retrieve a stored rotamer trie from inside the set;
	RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( elec_method ) ));
	RotamerTrieBaseCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( elec_method ) ));

	//fpd get rid of mutable data, use evaluator instead
	electrie::ElecTrieEvaluator eleceval(
		wt_bb_bb, wt_bb_sc, wt_sc_bb, wt_sc_sc, *this );

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp = get_count_pair_function_trie( set1, set2, pose, sfxn );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	trie1->trie_vs_trie( *trie2, *cp, eleceval, temp_table1, temp_table2 );

	/// add in the energies calculated by the tvt alg.
	energy_table += temp_table1;
	//std::cout << "FINISHED evaluate_rotamer_pair_energies" << std::endl;

	// There should be a way to turn this on without recompiling...
	// debug
	/*
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table3( energy_table );
	temp_table3 = 0;
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set1.num_rotamers(); ii <= ii_end; ++ii ) {
	for ( Size jj = 1, jj_end = set2.num_rotamers(); jj <= jj_end; ++jj ) {
	emap.zero();
	residue_pair_energy( *set1.rotamer( ii ), *set2.rotamer( jj ), pose, sfxn, emap );
	temp_table3( jj, ii ) += weights.dot( emap );
	if ( std::abs( temp_table1( jj, ii ) - temp_table3( jj, ii )) > 0.001 ) {
	std::cout << "FA_ElecE: Residues " << set1.resid() << " & " << set2.resid() << " rotamers: " << ii << " & " << jj;
	std::cout << " tvt/reg discrepancy: tvt= " <<  temp_table1( jj, ii ) << " reg= " << temp_table3( jj, ii );
	std::cout << " delta: " << temp_table1( jj, ii ) - temp_table3( jj, ii ) << std::endl;
	}
	}
	}
	std::cout << "Finished RPE calcs for residues " << set1.resid() << " & " << set2.resid() << std::endl;
	*/
}

void
FA_ElecEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & , //weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	// Since a rotamer set may include multiple residue types,
	// we'll make our decision based on what's currently in the Pose.
	if ( ! defines_score_for_residue_pair( pose.residue(set.resid()), residue, true) ) return;

	using namespace methods;
	using namespace trie;

	// allocate space for the trie-vs-trie algorithm
	utility::vector1< core::PackerEnergy > temp_vector1( set.num_rotamers(), 0.0 );
	utility::vector1< core::PackerEnergy > temp_vector2( set.num_rotamers(), 0.0 );

	// save weight information so that its available during tvt execution
	core::Real wt_bb_bb = sfxn.weights()[ fa_elec ] + sfxn.weights()[ fa_elec_bb_bb ];
	core::Real wt_bb_sc = sfxn.weights()[ fa_elec ] + sfxn.weights()[ fa_elec_bb_sc ];
	core::Real wt_sc_bb = sfxn.weights()[ fa_elec ] + sfxn.weights()[ fa_elec_bb_sc ];
	core::Real wt_sc_sc = sfxn.weights()[ fa_elec ] + sfxn.weights()[ fa_elec_sc_sc ];

	RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set.get_trie( elec_method ) ));
	RotamerTrieBaseCOP trie2 = ( static_cast< TrieCollection const & >
		( pose.energies().data().get( EnergiesCacheableDataType::ELEC_TRIE_COLLECTION )) ).trie( residue.seqpos() );

	if ( trie2 == NULL ) return;

	//fpd get rid of mutable data, use evaluator instead
	electrie::ElecTrieEvaluator eleceval(
		wt_bb_bb, wt_bb_sc, wt_sc_bb, wt_sc_sc, *this );

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp = get_count_pair_function_trie( pose.residue( set.resid() ), residue, trie1, trie2, pose, sfxn );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	trie1->trie_vs_path( *trie2, *cp, eleceval, temp_vector1, temp_vector2 );

	//std::cout << "FINISHED evaluate_rotamer_background_energies" << std::endl;
	/// add in the energies calculated by the tvt alg.
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		energy_vector[ ii ] += temp_vector1[ ii ];
	}

	//debug
	/*
	utility::vector1< Energy > temp_vector3( energy_vector.size(), 0.0f );
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
	emap.zero();
	residue_pair_energy( *set.rotamer( ii ), residue, pose, sfxn, emap );
	temp_vector3[ ii ] += weights.dot( emap );
	if ( std::abs( temp_vector1[ ii ] - temp_vector3[ ii ]) > 0.001 ) {
	std::cout << "FA_ElecE: Residues " << set.resid() << " & " << residue.seqpos() << " rotamers: " << ii << " & bg";
	std::cout << " tvt/reg discrepancy: tvt= " <<  temp_vector1[ ii ] << " reg= " << temp_vector3[ ii ];
	std::cout << " delta: " << temp_vector1[ ii ] - temp_vector3[ ii ] << std::endl;
	}
	}
	std::cout << "Finished Rotamer BG calcs for residues " << set.resid() << " & " << residue.seqpos() << std::endl;
	*/

}

// hpark: add in explicit enumeration with intra-res portion,
// not using trie wouldn't slow down much?
void
FA_ElecEnergy::evaluate_rotamer_intrares_energies(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< core::PackerEnergy > & energies
) const {
	//TR << "evaluate rotamer intra" << std::endl;
	if ( std::abs( sfxn.get_weight( scoring::fa_intra_elec ) ) < 1.0e-9 ) return;

	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		EnergyMap emap; emap[ fa_intra_elec ] = 0.0;
		conformation::Residue const & ii_rotamer( *set.rotamer( ii ));
		eval_intrares_energy( ii_rotamer, pose, sfxn, emap );
		//if( std::abs( emap[ fa_intra_elec ] ) > 1.0e-9 ){
		// TR << "ii/energy: " << ii << " " << emap[ fa_intra_elec ] << std::endl;
		//}
		energies[ ii ] += emap[ fa_intra_elec ];
	}
}

// could be updated base
bool
FA_ElecEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	if ( std::abs( weights[ fa_intra_elec ] ) > 1.0e-9 ) {
		return true;
	}
	return false;
}


/// @brief FA_ElecEnergy distance cutoff
///
/// Reports the maximum heavy atom/heavy atom distance at which two residues have a non-zero fa_elec interaction energy.
Distance
FA_ElecEnergy::atomic_interaction_cutoff() const
{
	return hydrogen_interaction_cutoff();
}

etable::count_pair::CountPairFunctionCOP
FA_ElecEnergy::get_intrares_countpair(
	conformation::Residue const &res,
	pose::Pose const &,
	ScoreFunction const &
) const
{
	using namespace etable::count_pair;
	return CountPairFactory::create_intrares_count_pair_function( res, CP_CROSSOVER_4 );
}

etable::count_pair::CountPairFunctionCOP
FA_ElecEnergy::get_count_pair_function(
	Size const res1,
	Size const res2,
	pose::Pose const & pose,
	ScoreFunction const &
) const
{
	using namespace etable::count_pair;

	//fd not sure if this is necessary....
	//if ( res1 == res2 ) {
	// return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );
	//}

	conformation::Residue const & rsd1( pose.residue( res1 ) );
	conformation::Residue const & rsd2( pose.residue( res2 ) );
	return get_count_pair_function( rsd1, rsd2 );
}


etable::count_pair::CountPairFunctionCOP
FA_ElecEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	using namespace etable::count_pair;

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );

	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		return CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	}
	return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairAll ) );

}


/// @brief FA_ElecEnergy
void
FA_ElecEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}


void
FA_ElecEnergy::setup_weight_triple(
	EnergyMap const & weights,
	weight_triple & wttrip
) const
{
	wttrip.wbb_bb_ = weights[ fa_elec ] + weights[ fa_elec_bb_bb ];
	wttrip.wbb_sc_ = weights[ fa_elec ] + weights[ fa_elec_bb_sc ];
	wttrip.wsc_bb_ = weights[ fa_elec ] + weights[ fa_elec_bb_sc ];
	wttrip.wsc_sc_ = weights[ fa_elec ] + weights[ fa_elec_sc_sc ];
}

// create an elec trie rotamer descriptor from a single rotamer
template <class CPDAT>
void
create_rotamer_descriptor(
	conformation::Residue const & res,
	trie::CPDataCorrespondence const &cpdata_map,
	std::map< chemical::ResidueType const *, std::map<core::Size,core::Size> > const &cp_reps,
	trie::RotamerDescriptor< electrie::ElecAtom, CPDAT > & rotamer_descriptor
)
{
	using namespace trie;
	using namespace electrie;

	rotamer_descriptor.natoms( res.natoms() );

	Size count_added_atoms = 0;
	for ( Size jj = 1; jj <= res.nheavyatoms(); ++jj ) {
		ElecAtom newatom;
		CPDAT cpdata;
		std::map< chemical::ResidueType const *, std::map<core::Size,core::Size> >::const_iterator it=cp_reps.find( &(res.type()));
		core::Size jj_rep = jj;
		if ( it != cp_reps.end() ) jj_rep=lookup_cp_map( it->second, jj );
		initialize_cpdata_for_atom( cpdata, jj_rep, res, cpdata_map );

		newatom.atom( res.atom(jj) );
		newatom.is_hydrogen( false );
		newatom.isbb( res.atom_is_backbone( jj ) );
		newatom.charge( res.atomic_charge(jj) );

		RotamerDescriptorAtom< ElecAtom, CPDAT > rdatom( newatom, cpdata );
		rotamer_descriptor.atom( ++count_added_atoms, rdatom );

		for ( Size kk = res.attached_H_begin( jj ),
				kk_end = res.attached_H_end( jj );
				kk <= kk_end; ++kk ) {
			ElecAtom newhatom;
			newhatom.atom( res.atom(kk) );
			newhatom.is_hydrogen( true );
			newhatom.isbb( res.atom_is_backbone( kk ) );
			newhatom.charge( res.atomic_charge( kk ) );

			CPDAT hcpdata;
			std::map< chemical::ResidueType const *, std::map<core::Size,core::Size> >::const_iterator it=cp_reps.find( &(res.type()));
			core::Size kk_rep = kk;
			if ( it != cp_reps.end() ) kk_rep=lookup_cp_map( it->second, kk );
			initialize_cpdata_for_atom( hcpdata, kk_rep, res, cpdata_map );

			RotamerDescriptorAtom< ElecAtom, CPDAT> hrdatom( newhatom, hcpdata );
			rotamer_descriptor.atom( ++count_added_atoms, hrdatom );
		}
	}
}


/// @brief create a rotamer trie for a particular set, deciding upon the kind of count pair data that
/// needs to be contained by the trie.
///
trie::RotamerTrieBaseOP
FA_ElecEnergy::create_rotamer_trie(
	conformation::RotamerSetBase const & rotset,
	pose::Pose const & // will be need to create tries for disulfides
) const
{
	using namespace trie;
	using namespace etable::etrie;
	using namespace electrie;

	trie::RotamerTrieBaseOP retval;

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamerset( rotset ) );
	if ( cpdata_map.has_pseudobonds() ||
			cpdata_map.max_connpoints_for_residue() > 1 ||
			cpdata_map.n_entries() > 3 ) {
		utility::vector1< RotamerDescriptor< electrie::ElecAtom, CountPairDataGeneric > > rotamer_descriptors( rotset.num_rotamers() );
		for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
			create_rotamer_descriptor( *rotset.rotamer( ii ), cpdata_map, cp_rep_map_, rotamer_descriptors[ ii ] );
			rotamer_descriptors[ ii ].rotamer_id( ii );
		}
		sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );
		retval = electrie::ElecRotamerTrieOP( new RotamerTrie< electrie::ElecAtom, CountPairDataGeneric >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		utility::vector1< RotamerDescriptor< ElecAtom, CountPairData_1_1 > > rotamer_descriptors( rotset.num_rotamers() );
		for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
			create_rotamer_descriptor( *rotset.rotamer( ii ), cpdata_map, cp_rep_map_, rotamer_descriptors[ ii ] );
			rotamer_descriptors[ ii ].rotamer_id( ii );
		}
		sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );
		retval = electrie::ElecRotamerTrieOP( new RotamerTrie< electrie::ElecAtom, CountPairData_1_1 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 2 ) {
		utility::vector1< RotamerDescriptor< ElecAtom, CountPairData_1_2 > > rotamer_descriptors( rotset.num_rotamers() );
		for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
			create_rotamer_descriptor( *rotset.rotamer( ii ), cpdata_map, cp_rep_map_, rotamer_descriptors[ ii ] );
			rotamer_descriptors[ ii ].rotamer_id( ii );
		}
		sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );
		retval = electrie::ElecRotamerTrieOP( new RotamerTrie< electrie::ElecAtom, CountPairData_1_2 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 3 ) {
		utility::vector1< RotamerDescriptor< ElecAtom, CountPairData_1_3 > > rotamer_descriptors( rotset.num_rotamers() );
		for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
			create_rotamer_descriptor( *rotset.rotamer( ii ), cpdata_map, cp_rep_map_, rotamer_descriptors[ ii ] );
			rotamer_descriptors[ ii ].rotamer_id( ii );
		}
		sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );
		retval = electrie::ElecRotamerTrieOP( new RotamerTrie< electrie::ElecAtom, CountPairData_1_3 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else {
		utility_exit_with_message( "Unknown residue connection in FA_ElecEnergy::create_rotamer_trie");
	}

	for ( Size ii = 1; ii <= cpdata_map.n_entries(); ++ii ) {
		retval->set_resid_2_connection_entry( cpdata_map.resid_for_entry( ii ), ii );
	}
	return retval;
}

/// @details Create a one-residue rotamer trie
trie::RotamerTrieBaseOP
FA_ElecEnergy::create_rotamer_trie(
	conformation::Residue const & res,
	pose::Pose const & /*pose*/
) const
{
	using namespace trie;
	using namespace etable::etrie;
	using namespace electrie;

	trie::RotamerTrieBaseOP retval;

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamer( res ) );
	if ( cpdata_map.has_pseudobonds() ||
			cpdata_map.max_connpoints_for_residue() > 1 ||
			cpdata_map.n_entries() > 3 ) {
		utility::vector1< RotamerDescriptor< ElecAtom, CountPairDataGeneric > > rotamer_descriptors( 1 );
		create_rotamer_descriptor( res, cpdata_map, cp_rep_map_, rotamer_descriptors[ 1 ] );
		rotamer_descriptors[ 1 ].rotamer_id( 1 );
		retval = electrie::ElecRotamerTrieOP( new RotamerTrie< ElecAtom, CountPairDataGeneric >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		utility::vector1< RotamerDescriptor< ElecAtom, CountPairData_1_1 > > rotamer_descriptors( 1 );
		create_rotamer_descriptor( res, cpdata_map, cp_rep_map_, rotamer_descriptors[ 1 ] );
		rotamer_descriptors[ 1 ].rotamer_id( 1 );
		retval = electrie::ElecRotamerTrieOP( new RotamerTrie< ElecAtom, CountPairData_1_1 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 2 ) {
		utility::vector1< RotamerDescriptor< ElecAtom, CountPairData_1_2 > > rotamer_descriptors( 1 );
		create_rotamer_descriptor( res, cpdata_map, cp_rep_map_, rotamer_descriptors[ 1 ] );
		rotamer_descriptors[ 1 ].rotamer_id( 1 );
		retval = electrie::ElecRotamerTrieOP( new RotamerTrie< ElecAtom, CountPairData_1_2 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 3 ) {
		utility::vector1< RotamerDescriptor< ElecAtom, CountPairData_1_3 > > rotamer_descriptors( 1 );
		create_rotamer_descriptor( res, cpdata_map, cp_rep_map_, rotamer_descriptors[ 1 ] );
		rotamer_descriptors[ 1 ].rotamer_id( 1 );
		retval = electrie::ElecRotamerTrieOP( new RotamerTrie< ElecAtom, CountPairData_1_3 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else {
		utility_exit_with_message( "Unknown residue connection in FA_ElecEnergy::create_rotamer_trie");
	}

	for ( Size ii = 1; ii <= cpdata_map.n_entries(); ++ii ) {
		retval->set_resid_2_connection_entry( cpdata_map.resid_for_entry( ii ), ii );
	}
	return retval;
}

/// @brief figure out the trie count pair function to use
/// Need to refactor this so that the decision "what kind of count pair behavior should I use" can be decoupled
/// from class instantiation, and therefore shared between the creation of the trie count pair classes and the regular
/// count pair classes
trie::TrieCountPairBaseOP
FA_ElecEnergy::get_count_pair_function_trie(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) const
{
	conformation::Residue const & res1( pose.residue( set1.resid() ) );
	conformation::Residue const & res2( pose.residue( set2.resid() ) );

	trie::RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( methods::elec_method ) ));
	trie::RotamerTrieBaseCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( methods::elec_method ) ));

	return get_count_pair_function_trie( res1, res2, trie1, trie2, pose, sfxn );
}


trie::TrieCountPairBaseOP
FA_ElecEnergy::get_count_pair_function_trie(
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
	if ( ! defines_score_for_residue_pair(res1, res2, true) ) return trie::TrieCountPairBaseOP( new TrieCountPairNone() );

	/// code needs to be added here to deal with multiple bonds (and psuedubonds!) between residues,
	/// but ultimately, this code is incompatible with designing both disulfides and non-disulfies
	/// at the same residue...
	CPResidueConnectionType connection = CountPairFactory::determine_residue_connection( res1, res2 );
	Size conn1 = trie1->get_count_pair_data_for_residue( res2.seqpos() );
	Size conn2 = trie2->get_count_pair_data_for_residue( res1.seqpos() );

	if ( connection == CP_ONE_BOND ) {
		tcpfxn = TrieCountPairBaseOP( new TrieCountPair1BC4( conn1, conn2 ) );
	} else if ( connection == CP_NO_BONDS ) {
		tcpfxn = TrieCountPairBaseOP( new TrieCountPairAll );
	} else {
		tcpfxn = TrieCountPairBaseOP( new TrieCountPairGeneric( res1, res2, conn1, conn2 ) );
	}
	return tcpfxn;

}


inline
Real
FA_ElecEnergy::score_atom_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Size const at1,
	Size const at2,
	EnergyMap & emap,
	Real const cpweight,
	Real & d2
) const
{

	Real energy;
	energy = cpweight * coulomb().eval_atom_atom_fa_elecE(
		rsd1.xyz(at1), rsd1.atomic_charge(at1),
		rsd2.xyz(at2), rsd2.atomic_charge(at2), d2);

	if ( rsd1.atom_is_backbone(at1) && rsd2.atom_is_backbone(at2) ) {
		emap[ fa_elec_bb_bb ] += energy;
	} else if ( !rsd1.atom_is_backbone(at1) && ! rsd2.atom_is_backbone(at2) ) {
		emap[ fa_elec_sc_sc ] += energy;
	} else {
		emap[ fa_elec_bb_sc ] += energy;
	}

	return energy;
}

core::Size
FA_ElecEnergy::version() const
{
	return 1; // Initial versioning
}


void
FA_ElecEnergy::set_nres_mono(
	core::pose::Pose const & pose
) const {
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue(i).is_upper_terminus() ) {
			nres_monomer_ = i;
			//std::cerr << "nres_monomer_ " << i << std::endl;
			return;
		}
	}
}


bool
FA_ElecEnergy::monomer_test(
	Size irsd,
	Size jrsd
) const {
	return (irsd <= nres_monomer_ && jrsd <= nres_monomer_ ) ||
		(irsd >  nres_monomer_ && jrsd >  nres_monomer_ );
}



} // namespace elec
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

typedef core::scoring::trie::RotamerTrie< core::scoring::elec::electrie::ElecAtom, core::scoring::etable::etrie::CountPairDataGeneric > EtrieRotTrieGeneric;

SAVE_AND_LOAD_SERIALIZABLE( EtrieRotTrieGeneric );
CEREAL_REGISTER_TYPE( EtrieRotTrieGeneric )

typedef core::scoring::trie::RotamerTrie< core::scoring::elec::electrie::ElecAtom, core::scoring::etable::etrie::CountPairData_1_1 > EtrieRotTrie11;
SAVE_AND_LOAD_SERIALIZABLE( EtrieRotTrie11 );
CEREAL_REGISTER_TYPE( EtrieRotTrie11 )

typedef core::scoring::trie::RotamerTrie< core::scoring::elec::electrie::ElecAtom, core::scoring::etable::etrie::CountPairData_1_2 > EtrieRotTrie12;
SAVE_AND_LOAD_SERIALIZABLE( EtrieRotTrie12 );
CEREAL_REGISTER_TYPE( EtrieRotTrie12 )

typedef core::scoring::trie::RotamerTrie< core::scoring::elec::electrie::ElecAtom, core::scoring::etable::etrie::CountPairData_1_3 > EtrieRotTrie13;
SAVE_AND_LOAD_SERIALIZABLE( EtrieRotTrie13 );
CEREAL_REGISTER_TYPE( EtrieRotTrie13 )


CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_elec_FA_ElecEnergy )
#endif
