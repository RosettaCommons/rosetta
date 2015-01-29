// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/HBondEnergy.fwd.hh
/// @brief  Hydrogen bond energy method forward declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay



// Unit Headers
#include <core/scoring/hbonds/HBondEnergy.hh>
#include <core/scoring/hbonds/HBondEnergyCreator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/scoring/hbonds/hbtrie/HBAtom.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.hh>
#include <core/scoring/hbonds/hbtrie/HBondTrie.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBCountPairFunction.hh>

#include <core/scoring/methods/LK_BallEnergy.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/trie/TrieCollection.hh>
#include <core/scoring/trie/RotamerTrieBase.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <core/scoring/trie/RotamerDescriptor.hh>

#include <core/scoring/dssp/Dssp.hh>  // for helix-len-dep hbond strength

#include <core/scoring/methods/EnergyMethodOptions.hh>

// membrane specific initialization
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneInfo.hh>

#include <basic/datacache/BasicDataCache.hh>

#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableData.hh>

// Numeric Headers
#include <numeric/numeric.functions.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/id/types.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <utility/vector1.hh>
#include <boost/unordered_map.hpp>

namespace core {
namespace scoring {
namespace hbonds {

class HBondResidueMinData;
class HBondResPairMinData;

typedef utility::pointer::shared_ptr< HBondResidueMinData >       HBondResidueMinDataOP;
typedef utility::pointer::shared_ptr< HBondResidueMinData const > HBondResidueMinDataCOP;

typedef utility::pointer::shared_ptr< HBondResPairMinData >       HBondResPairMinDataOP;
typedef utility::pointer::shared_ptr< HBondResPairMinData const > HBondResPairMinDataCOP;

/// @brief A class to hold data for the HBondEnergy class used in
/// score and derivative evaluation.
class HBondResidueMinData : public basic::datacache::CacheableData
{
public:
	HBondResidueMinData() : bb_don_avail_( true ), bb_acc_avail_( true ) {}
	virtual ~HBondResidueMinData() {}

	virtual basic::datacache::CacheableDataOP clone() const
	{ return basic::datacache::CacheableDataOP( new HBondResidueMinData( *this ) ); }

	void set_bb_don_avail( bool setting ) { bb_don_avail_ = setting; }
	void set_bb_acc_avail( bool setting ) { bb_acc_avail_ = setting; }

	bool bb_don_avail() const { return bb_don_avail_; }
	bool bb_acc_avail() const { return bb_acc_avail_; }

	void set_natoms( Size setting ) { natoms_ = setting; }
	Size natoms() const { return natoms_; }

	void set_nneighbors( Size setting ) { nneighbors_ = setting; }
	Size nneighbors() const { return nneighbors_; }

private:
	Size natoms_;
	Size nneighbors_;

	bool bb_don_avail_;
	bool bb_acc_avail_;
};

class HBondResPairMinData : public basic::datacache::CacheableData
{
public:
	HBondResPairMinData();
	virtual ~HBondResPairMinData() {}
	virtual basic::datacache::CacheableDataOP clone() const
	{ return basic::datacache::CacheableDataOP( new HBondResPairMinData( *this ) ); }

	void set_res1_data( HBondResidueMinDataCOP );
	void set_res2_data( HBondResidueMinDataCOP );

	HBondResidueMinData const & res1_data() const { return *res1_dat_; }
	HBondResidueMinData const & res2_data() const { return *res2_dat_; }

	//void update_natoms();

	void clear_hbonds();
	void add_hbond( HBond const & );

private:

	HBondResidueMinDataCOP res1_dat_;
	HBondResidueMinDataCOP res2_dat_;

};

HBondResPairMinData::HBondResPairMinData() {}
void HBondResPairMinData::set_res1_data( HBondResidueMinDataCOP dat ) { res1_dat_ = dat; }
void HBondResPairMinData::set_res2_data( HBondResidueMinDataCOP dat ) { res2_dat_ = dat; }

static thread_local basic::Tracer tr( "core.scoring.hbonds.HbondEnergy" );

/// @details This must return a fresh instance of the HBondEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
HBondEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new HBondEnergy( options.hbond_options() ) );
}

ScoreTypes
HBondEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( hbond_lr_bb );
	sts.push_back( hbond_sr_bb );
	sts.push_back( hbond_bb_sc );
	sts.push_back( hbond_sr_bb_sc );
	sts.push_back( hbond_lr_bb_sc );
	sts.push_back( hbond_sc );
	sts.push_back( hbond_intra ); //Currently affects only RNA.
	sts.push_back( hbond );
	return sts;
}


/// ctor
HBondEnergy::HBondEnergy( HBondOptions const & opts ):
	parent( methods::EnergyMethodCreatorOP( methods::EnergyMethodCreatorOP( new HBondEnergyCreator ) ) ),
	options_( HBondOptionsCOP( HBondOptionsOP( new HBondOptions( opts ) ) )),
	database_( HBondDatabase::get_database(opts.params_database_tag()) )
{}

/// copy ctor
HBondEnergy::HBondEnergy( HBondEnergy const & src ):
	parent( src ),
	options_( HBondOptionsCOP( HBondOptionsOP( new HBondOptions( *src.options_) ) )) ,
	database_( src.database_)
{}

HBondEnergy::~HBondEnergy() {}

/// clone
methods::EnergyMethodOP
HBondEnergy::clone() const
{
	return methods::EnergyMethodOP( new HBondEnergy( *this ) );
}

///
void
HBondEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	using EnergiesCacheableDataType::HBOND_SET;
	using EnergiesCacheableDataType::HBOND_TRIE_COLLECTION;

	pose.update_residue_neighbors();
	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( *options_ ) );

	// membrane object initialization
	if (options_->Mbhbond()) {
		Membrane_FAPotential const & memb_potential_ = ScoringManager::get_instance()->get_Membrane_FAPotential();
		memb_potential_.compute_fa_projection( pose );
		normal_ = MembraneEmbed_from_pose( pose ).normal();
		center_ = MembraneEmbed_from_pose( pose ).center();
		thickness_ = Membrane_FAEmbed_from_pose( pose ).thickness();
		steepness_ = Membrane_FAEmbed_from_pose( pose ).steepness();
	}

	// membrane framework object initialization
	if ( options_->mphbond() ) {

		// Initialize membrane specific parameters
		normal_ = pose.conformation().membrane_info()->membrane_normal();
		center_ = pose.conformation().membrane_info()->membrane_center();
		thickness_ = pose.conformation().membrane_info()->membrane_thickness();
		steepness_ = pose.conformation().membrane_info()->membrane_steepness();

	}

	hbond_set->setup_for_residue_pair_energies( pose );
	pose.energies().data().set( HBOND_SET, hbond_set );

	using namespace trie;
	using namespace hbtrie;

	TrieCollectionOP tries( new TrieCollection );
	tries->total_residue( pose.total_residue() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		// Do not compute energy for virtual residues.
		if ( pose.residue(ii).aa() == core::chemical::aa_vrt) continue;

		HBondRotamerTrieOP one_rotamer_trie = create_rotamer_trie( pose.residue( ii ), pose );
		tries->trie( ii, one_rotamer_trie );
	}
	pose.energies().data().set( HBOND_TRIE_COLLECTION, tries );
}

void
HBondEnergy::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & set
) const
{
	using namespace hbtrie;

	HBondRotamerTrieOP rottrie = create_rotamer_trie( set, pose );
	//std::cout << "--------------------------------------------------" << std::endl << " HBROTAMER TRIE: " << set.resid() << std::endl;
	//rottrie->print();
	set.store_trie( methods::hbond_method, rottrie );
}

// Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
HBondEnergy::update_residue_for_packing( pose::Pose & pose, Size resid ) const
{
	using namespace trie;
	using namespace hbtrie;
	using EnergiesCacheableDataType::HBOND_TRIE_COLLECTION;

	HBondRotamerTrieOP one_rotamer_trie = create_rotamer_trie( pose.residue( resid ), pose );

	// grab non-const & of the cached tries and replace resid's trie with a new one.
	TrieCollection & trie_collection( static_cast< TrieCollection & > (pose.energies().data().get( HBOND_TRIE_COLLECTION )));
	trie_collection.trie( resid, one_rotamer_trie );
}



///
void
HBondEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	pose.update_residue_neighbors();
	HBondSetOP hbond_set( new hbonds::HBondSet( *options_, pose.total_residue() ) );

	// membrane object initialization
	if (options_->Mbhbond()) {
		Membrane_FAPotential const & memb_potential_ = ScoringManager::get_instance()->get_Membrane_FAPotential();
		memb_potential_.compute_fa_projection( pose );
		normal_ = MembraneEmbed_from_pose( pose ).normal();
		center_ = MembraneEmbed_from_pose( pose ).center();
		thickness_ = Membrane_FAEmbed_from_pose( pose ).thickness();
		steepness_ = Membrane_FAEmbed_from_pose( pose ).steepness();
	}

	// membrane framework supported membrane object initialization
	if ( options_->mphbond() ) {

		// Initialize membrane parameters
		normal_ = pose.conformation().membrane_info()->membrane_normal();
		center_ = pose.conformation().membrane_info()->membrane_center();
		thickness_ = pose.conformation().membrane_info()->membrane_thickness();
		steepness_ = pose.conformation().membrane_info()->membrane_steepness();

	}

	//fpd  we need secstruct info in some cases ... don't change while minimizing
	//fpd      MUST be called before setup_for_residue_pair_energies
	if (options_->length_dependent_srbb() && ! pose.energies().use_nblist()) {
		dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );
	}

	hbond_set->setup_for_residue_pair_energies( pose );

	/// During minimization, keep the set of bb/bb hbonds "fixed" by using the old boolean values.
	if ( pose.energies().use_nblist() && pose.energies().data().has( HBOND_SET ) ) {
		HBondSet const & existing_set = static_cast< HBondSet const & > (pose.energies().data().get( HBOND_SET ));
		hbond_set->copy_bb_donor_acceptor_arrays( existing_set );
	}
	pose.energies().data().set( HBOND_SET, hbond_set );

}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @details Note that this only evaluates sc-sc and sc-bb energies unless options_->decompose_bb_hb_into_pair_energies
/// is set to true, in which case, this function also evaluates bb-bb energies.
/// Note also that this function enforces the bb/sc hbond exclusion rule.
void
HBondEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	if ( rsd1.seqpos() == rsd2.seqpos() ) return;
	if ( options_->exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( HBOND_SET )));

	// this only works because we have already called
	// hbond_set->setup_for_residue_pair_energies( pose )

	// Non-pairwise additive exclusion rules:
	// exclude backbone-backbone hbond if set in options_ (if, say, they were pre-computed)
	// exclude backbone-sidechain hbond if backbone-backbone hbond already in hbond_set*
	// exclude sidechain-backbone hbond if backbone-backbone hbond already in hbond_set*
	// * these two rules are only enforced as long as bb_donor_acceptor_check is "true"

	// mjo Historically, if this exclusion rule is not
	// enforced--accoring to Brian Kuhlman--"Serines are put up and down
	// helices".  According to John Karanicolas, amide acceptors have a
	// local energy minima where one lone pair moves in line with the
	// base-acceptor, and the other lone pair is delocalized in between
	// the base and acceptor atoms.  In this configuration it is
	// energetically disfavorable to make multiple hbonds with the
	// acceptor.

	// NOTE: "bsc" -> acc=bb don=sc
	//       "scb" -> don=sc don=bb

	bool exclude_bsc = false, exclude_scb = false;
	if (rsd1.is_protein()) exclude_scb = options_->bb_donor_acceptor_check() && hbond_set.don_bbg_in_bb_bb_hbond(rsd1.seqpos());
	if (rsd2.is_protein()) exclude_bsc = options_->bb_donor_acceptor_check() && hbond_set.acc_bbg_in_bb_bb_hbond(rsd2.seqpos());

	//pba membrane dependent hbond potential hack
	if ( options_->Mbhbond() && options_->mphbond() ) {

		identify_hbonds_1way_membrane(
			*database_,
			rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd2.seqpos()),
			false /*calculate_derivative*/,
			!options_->decompose_bb_hb_into_pair_energies(), exclude_bsc, exclude_scb, false,
			*options_,
			emap,
			pose);

		exclude_bsc = exclude_scb = false;
		if (rsd2.is_protein()) exclude_scb = options_->bb_donor_acceptor_check() && hbond_set.don_bbg_in_bb_bb_hbond(rsd2.seqpos());
		if (rsd1.is_protein()) exclude_bsc = options_->bb_donor_acceptor_check() && hbond_set.acc_bbg_in_bb_bb_hbond(rsd1.seqpos());

		identify_hbonds_1way_membrane(
			*database_,
			rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			!options_->decompose_bb_hb_into_pair_energies(), exclude_bsc, exclude_scb, false,
			*options_,
			emap,
			pose);

	  } else {

		//ss-dep weights
		SSWeightParameters ssdep;
		ssdep.ssdep_ = options_->length_dependent_srbb();
		ssdep.l_ = options_->length_dependent_srbb_lowscale();
		ssdep.h_ = options_->length_dependent_srbb_highscale();
		ssdep.len_l_ = options_->length_dependent_srbb_minlength();
		ssdep.len_h_ = options_->length_dependent_srbb_maxlength();
		Real ssdep_weight_factor = get_ssdep_weight(rsd1, rsd2, pose, ssdep);

		identify_hbonds_1way(
			*database_,
			rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd2.seqpos()),
			false /*calculate_derivative*/,
			!options_->decompose_bb_hb_into_pair_energies(), exclude_bsc, exclude_scb, false,
			*options_,
			emap, num_hbonds_, ssdep_weight_factor);

		exclude_bsc = exclude_scb = false;
		if (rsd2.is_protein()) exclude_scb = options_->bb_donor_acceptor_check() && hbond_set.don_bbg_in_bb_bb_hbond(rsd2.seqpos());
		if (rsd1.is_protein()) exclude_bsc = options_->bb_donor_acceptor_check() && hbond_set.acc_bbg_in_bb_bb_hbond(rsd1.seqpos());

		identify_hbonds_1way(
			*database_,
			rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			!options_->decompose_bb_hb_into_pair_energies(), exclude_bsc, exclude_scb, false,
			*options_,
			emap, num_hbonds_, ssdep_weight_factor);
		}

	//std::cout << std::endl << num_hbonds_.size() << std::endl;
}

bool
HBondEnergy::defines_score_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	bool res_moving_wrt_eachother
) const
{
	return res_moving_wrt_eachother;
}

bool
HBondEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

bool
HBondEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

/// @details computes the residue pair energy during minimization; this includes bb/bb energies
/// as opposed to the standard residue_pair_energy interface, which does not include bb/bb energies.
/// On the other hand, this interface presumes that no new bb/bb hydrogen bonds are formed during
/// the course of minimization, and no existing bb/bb hydrogen bonds are lost.
/// Note: this function does not directly enforce the bb/sc exclusion rule logic, but rather,
/// takes the boolean "bb_don_avail" and "bb_acc_avail" data stored in the pairdata object
void
HBondEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & pairdata,
	pose::Pose const & pose, //pba
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( rsd1.xyz( rsd1.nbr_atom() ).distance_squared( rsd2.xyz( rsd2.nbr_atom() ) )
		> std::pow( rsd1.nbr_radius() + rsd2.nbr_radius() + atomic_interaction_cutoff(), 2 ) ) return;

debug_assert( utility::pointer::dynamic_pointer_cast< HBondResPairMinData const > ( pairdata.get_data( hbond_respair_data ) ));
	HBondResPairMinData const & hb_pair_dat( static_cast< HBondResPairMinData const & > ( pairdata.get_data_ref( hbond_respair_data ) ));

	// membrane dependent hbond potential hack
	if ( options_->Mbhbond() && options_->mphbond() ) {

		{ // scope        /// 1st == evaluate hbonds with donor atoms on rsd1
		/// case A: sc is acceptor, bb is donor && res2 is the acceptor residue -> look at the donor availability of residue 1
		bool exclude_scb( ! hb_pair_dat.res1_data().bb_don_avail() );
		/// case B: bb is acceptor, sc is donor && res2 is the acceptor residue -> look at the acceptor availability of residue 2
		bool exclude_bsc( ! hb_pair_dat.res2_data().bb_acc_avail() );

		identify_hbonds_1way_membrane(
			*database_,
			rsd1, rsd2,
			hb_pair_dat.res1_data().nneighbors(), hb_pair_dat.res2_data().nneighbors(),
			false /*calculate_derivative*/,
			false, exclude_bsc, exclude_scb, false,
			*options_,
			emap, pose);
		}

		{ // scope
		/// 2nd == evaluate hbonds with donor atoms on rsd2
		/// case A: sc is acceptor, bb is donor && res1 is the acceptor residue -> look at the donor availability of residue 2
		bool exclude_scb( ! hb_pair_dat.res2_data().bb_don_avail() );
		/// case B: bb is acceptor, sc is donor && res1 is the acceptor residue -> look at the acceptor availability of residue 1
		bool exclude_bsc( ! hb_pair_dat.res1_data().bb_acc_avail() );

		identify_hbonds_1way_membrane(
			*database_,
			rsd2, rsd1,
			hb_pair_dat.res2_data().nneighbors(), hb_pair_dat.res1_data().nneighbors(),
			false /*calculate_derivative*/,
			false, exclude_bsc, exclude_scb, false,
			*options_,
			emap, pose);
		}

	} else {

		//ss-dep weights
		SSWeightParameters ssdep;
		ssdep.ssdep_ = options_->length_dependent_srbb();
		ssdep.l_ = options_->length_dependent_srbb_lowscale();
		ssdep.h_ = options_->length_dependent_srbb_highscale();
		ssdep.len_l_ = options_->length_dependent_srbb_minlength();
		ssdep.len_h_ = options_->length_dependent_srbb_maxlength();
		Real ssdep_weight_factor = get_ssdep_weight(rsd1, rsd2, pose, ssdep);

		{ // scope
		/// 1st == evaluate hbonds with donor atoms on rsd1
		/// case A: sc is acceptor, bb is donor && res2 is the acceptor residue -> look at the donor availability of residue 1
		bool exclude_scb( ! hb_pair_dat.res1_data().bb_don_avail() );
		/// case B: bb is acceptor, sc is donor && res2 is the acceptor residue -> look at the acceptor availability of residue 2
		bool exclude_bsc( ! hb_pair_dat.res2_data().bb_acc_avail() );

		identify_hbonds_1way(
			*database_,
			rsd1, rsd2,
			hb_pair_dat.res1_data().nneighbors(), hb_pair_dat.res2_data().nneighbors(),
			false /*calculate_derivative*/,
			false, exclude_bsc, exclude_scb, false,
			*options_,
			emap, ssdep_weight_factor);
		}

		{ // scope
		/// 2nd == evaluate hbonds with donor atoms on rsd2
		/// case A: sc is acceptor, bb is donor && res1 is the acceptor residue -> look at the donor availability of residue 2
		bool exclude_scb( ! hb_pair_dat.res2_data().bb_don_avail() );
		/// case B: bb is acceptor, sc is donor && res1 is the acceptor residue -> look at the acceptor availability of residue 1
		bool exclude_bsc( ! hb_pair_dat.res1_data().bb_acc_avail() );

		identify_hbonds_1way(
			*database_,
			rsd2, rsd1,
			hb_pair_dat.res2_data().nneighbors(), hb_pair_dat.res1_data().nneighbors(),
			false /*calculate_derivative*/,
			false, exclude_bsc, exclude_scb, false,
			*options_,
			emap, ssdep_weight_factor);
		}

	}

}

/// @details Note that this function helps enforce the bb/sc exclusion rule by setting the donor and acceptor availability
/// for backbone donors and acceptors.  If the backbone-sidechain-exclusion rule is not being enforced, then this function
/// marks all donors and acceptors as being available.  If it is being enforced, then is uses the hbondset functions
/// don_bbg_in_bb_bb_hbond and acc_bbg_in_bb_bb_hbond.  The decisions made in this function impact the evaluation
/// of energies in the above residue_pair_energy_ext method.
void
HBondEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData & res_data_cache
) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	HBondSet const & hbondset = static_cast< HBondSet const & > (pose.energies().data().get( HBOND_SET ));
	HBondResidueMinDataOP hbresdata( 0 );
	if ( res_data_cache.get_data( hbond_res_data ) ) {
		hbresdata = utility::pointer::static_pointer_cast< core::scoring::hbonds::HBondResidueMinData > ( res_data_cache.get_data( hbond_res_data ) );
		// assume that bb-don-avail and bb-acc-avail are already initialized
	} else {
		hbresdata = HBondResidueMinDataOP( new HBondResidueMinData );
		hbresdata->set_nneighbors( hbondset.nbrs( rsd.seqpos() ) );
		if ( rsd.is_protein() ) {
			hbresdata->set_bb_don_avail( options_->bb_donor_acceptor_check() ? ! hbondset.don_bbg_in_bb_bb_hbond( rsd.seqpos() ) : true );
			hbresdata->set_bb_acc_avail( options_->bb_donor_acceptor_check() ? ! hbondset.acc_bbg_in_bb_bb_hbond( rsd.seqpos() ) : true );
		}
		res_data_cache.set_data( hbond_res_data, hbresdata );
	}
	hbresdata->set_natoms( rsd.natoms() );
}


void
HBondEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const & res1_data_cache,
	ResSingleMinimizationData const & res2_data_cache,
	ResPairMinimizationData & data_cache
) const
{
	HBondResPairMinDataOP hbpairdat;
	if ( data_cache.get_data( hbond_respair_data ) ) {
		// assume that hbpairdat has already been pointed at its two residues, and that it may need to update
		// its size based on a change to the number of atoms from a previous initialization.
	debug_assert( utility::pointer::dynamic_pointer_cast< HBondResPairMinData > ( data_cache.get_data( hbond_respair_data ) ));
		hbpairdat = utility::pointer::static_pointer_cast< core::scoring::hbonds::HBondResPairMinData > ( data_cache.get_data( hbond_respair_data ) );
	} else {
	debug_assert( utility::pointer::dynamic_pointer_cast< HBondResidueMinData const > ( res1_data_cache.get_data( hbond_res_data ) ));
	debug_assert( utility::pointer::dynamic_pointer_cast< HBondResidueMinData const > ( res2_data_cache.get_data( hbond_res_data ) ));

		hbpairdat = HBondResPairMinDataOP( new HBondResPairMinData );
		hbpairdat->set_res1_data( utility::pointer::static_pointer_cast< core::scoring::hbonds::HBondResidueMinData const > ( res1_data_cache.get_data( hbond_res_data ) ));
		hbpairdat->set_res2_data( utility::pointer::static_pointer_cast< core::scoring::hbonds::HBondResidueMinData const > ( res2_data_cache.get_data( hbond_res_data ) ));
		data_cache.set_data( hbond_respair_data, hbpairdat );
	}
	//hbpairdat->update_natoms();

}

bool
HBondEnergy::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return false;
}


/// @details Triplication of the loops that iterate across hbond donors and acceptors
/// Find all hbonds for a pair of residues and add those found hbonds to the
/// hb_pair_dat object; these hbonds will be used for derivative evaluation, so evaluate
/// the F1/F2 derivative vectors now.  This function respects the exclude bsc and scb variables
/// to avoid hbonds.  Called by setup_for_derivatives_for_residue_pair
/// pba modified
void
HBondEnergy::hbond_derivs_1way(
	EnergyMap const & weights,
	HBondSet const & hbond_set,
	HBondDatabaseCOP database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	Real const ssdep_weight_factor,
	// output
	utility::vector1< DerivVectorPair > & don_atom_derivs,
	utility::vector1< DerivVectorPair > & acc_atom_derivs
) const
{
	EnergyMap emap;

	bool is_intra_res = (don_rsd.seqpos() == acc_rsd.seqpos());
	if ( is_intra_res && !calculate_intra_res_hbonds( don_rsd, hbond_set.hbond_options() ) ) return;

	// <f1,f2> -- derivative vectors
	HBondDerivs deriv;

	for ( chemical::AtomIndices::const_iterator
			hnum  = don_rsd.Hpos_polar().begin(),	hnume = don_rsd.Hpos_polar().end();
			hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		Vector const & hatm_xyz( don_rsd.atom(hatm).xyz() );
		Vector const & datm_xyz( don_rsd.atom(datm).xyz() );

		for ( chemical::AtomIndices::const_iterator
				anum  = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end();
				anum != anume; ++anum ) {

			Size const aatm( *anum );

			if ( acc_rsd.atom_is_backbone(aatm)){
				if ( ! datm_is_bb && exclude_bsc ) continue; // if the donor is sc, the acceptor bb, and exclude_b(a)sc(d)
			} else {
				if (datm_is_bb && exclude_scb) continue; // if the donor is bb, the acceptor sc, and exclude_sc(a)b(d)
			}

			// rough filter for existance of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalTuple hbe_type( datm, don_rsd, aatm, acc_rsd);

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
		debug_assert( base2 > 0 && base != base2 );

			hb_energy_deriv( *database, *options_, hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, true /*eval deriv*/, deriv);

			if (unweighted_energy >= MAX_HB_ENERGY) continue;

			//pba buggy?
			Real weighted_energy = // evn weight * weight-set[ hbtype ] weight
				(! hbond_set.hbond_options().use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, hbond_set.hbond_options() )) *
				hb_eval_type_weight( hbe_type.eval_type(), weights, is_intra_res, hbond_set.hbond_options().put_intra_into_total() );
			weighted_energy *= ssdep_weight_factor;

			// membrane specific correction
			if ( options_->Mbhbond() && options_->mphbond()  ) {
				weighted_energy = get_membrane_depth_dependent_weight(normal_, center_, thickness_,
				steepness_, don_nb, acc_nb, hatm_xyz, acc_rsd.atom(aatm ).xyz()) *
					hb_eval_type_weight( hbe_type.eval_type(), weights, is_intra_res, hbond_set.hbond_options().put_intra_into_total() );
			}

			/// Only accumulate h and acc derivs for now; soon, don, h, acc, a-base and abase2 derivs
			don_atom_derivs[ datm ].f1() += weighted_energy * deriv.don_deriv.f1();
			don_atom_derivs[ datm ].f2() += weighted_energy * deriv.don_deriv.f2();
			don_atom_derivs[ hatm ].f1() += weighted_energy * deriv.h_deriv.f1();
			don_atom_derivs[ hatm ].f2() += weighted_energy * deriv.h_deriv.f2();
			acc_atom_derivs[ aatm ].f1() += weighted_energy * deriv.acc_deriv.f1();
			acc_atom_derivs[ aatm ].f2() += weighted_energy * deriv.acc_deriv.f2();

			// ring-acceptor derivative assignment logic is tricky
			assign_abase_derivs( *options_, acc_rsd, aatm, hbe_type, deriv.abase_deriv, weighted_energy, acc_atom_derivs );


			acc_atom_derivs[ base2 ].f1() += weighted_energy * deriv.abase2_deriv.f1();
			acc_atom_derivs[ base2 ].f2() += weighted_energy * deriv.abase2_deriv.f2();

		} // loop over acceptors
	} // loop over donors

}

void
HBondEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{

	if ( calculate_intra_res_hbonds( rsd, *options_ ) ) {
		using EnergiesCacheableDataType::HBOND_SET;
		HBondSet const & hbond_set = static_cast< HBondSet const & > (pose.energies().data().get( HBOND_SET ));
		bool const exclude_scb( false );
		hbond_derivs_1way( weights, hbond_set, database_, rsd, rsd, 1, 1, exclude_scb, exclude_scb, 1, atom_derivs, atom_derivs );
	}
}

void
HBondEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // stores the hbond database in the cached hbond set
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	if ( rsd1.xyz( rsd1.nbr_atom() ).distance_squared( rsd2.xyz( rsd2.nbr_atom() ) )
		> std::pow( rsd1.nbr_radius() + rsd2.nbr_radius() + atomic_interaction_cutoff(), 2 ) ) return;

	/// Iterate across all acceptor and donor atom pairs for these two residues, and write down the hydrogen bonds
	/// that are formed.

	using EnergiesCacheableDataType::HBOND_SET;

	HBondSet const & hbondset = static_cast< HBondSet const & > (pose.energies().data().get( HBOND_SET ));

debug_assert( utility::pointer::dynamic_pointer_cast< HBondResPairMinData const > ( min_data.get_data( hbond_respair_data ) ));
	HBondResPairMinData const & hb_pair_dat = static_cast< HBondResPairMinData const & > ( min_data.get_data_ref( hbond_respair_data ) );

	Size const rsd1nneighbs( hb_pair_dat.res1_data().nneighbors() );
	Size const rsd2nneighbs( hb_pair_dat.res2_data().nneighbors() );

	//ss-dep weights
	SSWeightParameters ssdep;
	ssdep.ssdep_ = options_->length_dependent_srbb();
	ssdep.l_ = options_->length_dependent_srbb_lowscale();
	ssdep.h_ = options_->length_dependent_srbb_highscale();
	ssdep.len_l_ = options_->length_dependent_srbb_minlength();
	ssdep.len_h_ = options_->length_dependent_srbb_maxlength();
	Real ssdep_weight_factor = get_ssdep_weight(rsd1, rsd2, pose, ssdep);

	{ // scope
	/// 1st == find hbonds with donor atoms on rsd1
	/// case A: sc is acceptor, bb is donor && res2 is the acceptor residue -> look at the donor availability of residue 1
	bool exclude_scb( ! hb_pair_dat.res1_data().bb_don_avail() );
	/// case B: bb is acceptor, sc is donor && res2 is the acceptor residue -> look at the acceptor availability of residue 2
	bool exclude_bsc( ! hb_pair_dat.res2_data().bb_acc_avail() );

	hbond_derivs_1way( weights, hbondset, database_, rsd1, rsd2, rsd1nneighbs, rsd2nneighbs, exclude_bsc, exclude_scb, ssdep_weight_factor, r1_atom_derivs, r2_atom_derivs );
	}

	{ // scope
	/// 2nd == evaluate hbonds with donor atoms on rsd2
	/// case A: sc is acceptor, bb is donor && res1 is the acceptor residue -> look at the donor availability of residue 2
	bool exclude_scb( ! hb_pair_dat.res2_data().bb_don_avail() );
	/// case B: bb is acceptor, sc is donor && res1 is the acceptor residue -> look at the acceptor availability of residue 1
	bool exclude_bsc( ! hb_pair_dat.res1_data().bb_acc_avail() );


	hbond_derivs_1way( weights, hbondset, database_, rsd2, rsd1, rsd2nneighbs, rsd1nneighbs, exclude_bsc, exclude_scb, ssdep_weight_factor, r2_atom_derivs, r1_atom_derivs );
	}

}

void
HBondEnergy::backbone_backbone_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap & emap
) const
{
	/// if we're including bb/bb energies in the energy maps, then they need to be calculated
	/// in the backbone_backbone_energy so that:
	/// residue_pair_energy = backbone_backbone_energy(r1,r2) + backbone_sidechain_energy(r1,r2) +
	///                       backbone_sidechain_energy(r2,r1) + sidechain_sidechain_energy(r1,r2)

	if ( ! options_->decompose_bb_hb_into_pair_energies() ) return;

	using EnergiesCacheableDataType::HBOND_SET;

	if ( rsd1.seqpos() == rsd2.seqpos() ) return;
	if ( options_->exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >( pose.energies().data().get( HBOND_SET )));

	// membrane dependent hbond potential hack
	if ( options_->Mbhbond() && options_->mphbond() ) {

		identify_hbonds_1way_membrane(
			*database_,
			rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			false, true, true, true, /* calc bb_bb, don't calc bb_sc, sc_bb, or sc_sc */
			*options_,
			emap, pose);

		identify_hbonds_1way_membrane(
			*database_,
			rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			false, true, true, true, /* calc bb_bb, don't calc bb_sc, sc_bb, or sc_sc */
			*options_,
			emap, pose);

	} else {

		identify_hbonds_1way(
			*database_,
			rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			false, true, true, true, /* calc bb_bb, don't calc bb_sc, sc_bb, or sc_sc */
			*options_,
			emap);

		identify_hbonds_1way(
			*database_,
			rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			false, true, true, true, /* calc bb_bb, don't calc bb_sc, sc_bb, or sc_sc */
			*options_,
			emap);

	}
}

/// @details Note: this function enforces the bb/sc hbond exclusion rule
void
HBondEnergy::backbone_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap & emap
) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	if ( rsd1.seqpos() == rsd2.seqpos() ) return;
	if ( options_->exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & > ( pose.energies().data().get( HBOND_SET )));

	// this only works because we have already called
	// hbond_set->setup_for_residue_pair_energies( pose )

	// membrane dependent hbond potential hack
	if ( options_->Mbhbond() && options_->mphbond() ) {

		/// If we're enforcing the bb/sc exclusion rule, and residue1 is a protein residue, and if residue 1's backbone-donor group
		/// is already participating in a bb/bb hbond, then do not evaluate the identify_hbonds_1way_membrane function
		if ( !options_->bb_donor_acceptor_check() || ! rsd1.is_protein() || ! hbond_set.don_bbg_in_bb_bb_hbond(rsd1.seqpos())) {
			identify_hbonds_1way_membrane(
				*database_,
				rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd2.seqpos()),
				false /*calculate_derivative*/,
				true, true, false, true,
				*options_,
				emap, pose);
		}

		/// If we're enforcing the bb/sc exclusion rule, and residue1 is a protein residue, and if residue 1's backbone-acceptor group
		/// is already participating in a bb/bb hbond, then do not evaluate the identify_hbonds_1way_membrane function
		if ( !options_->bb_donor_acceptor_check() || ! rsd1.is_protein() || ! hbond_set.acc_bbg_in_bb_bb_hbond(rsd1.seqpos())) {
			identify_hbonds_1way_membrane(
				*database_,
				rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
				false /*calculate_derivative*/,
				true, false, true, true,
				*options_,
				emap, pose);
		}

	} else {
		/// If we're enforcing the bb/sc exclusion rule, and residue1 is a protein residue, and if residue 1's backbone-donor group
		/// is already participating in a bb/bb hbond, then do not evaluate the identify_hbonds_1way function
		if ( !options_->bb_donor_acceptor_check() || !rsd1.is_protein() || !hbond_set.don_bbg_in_bb_bb_hbond(rsd1.seqpos())) {
			identify_hbonds_1way(
				*database_,
				rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd2.seqpos()),
				false /*calculate_derivative*/,
				true, true, false, true,
				*options_,
				emap);
		}

		/// If we're enforcing the bb/sc exclusion rule, and residue1 is a protein residue, and if residue 1's backbone-acceptor group
		/// is already participating in a bb/bb hbond, then do not evaluate the identify_hbonds_1way function
		if ( !options_->bb_donor_acceptor_check() || !rsd1.is_protein() || !hbond_set.acc_bbg_in_bb_bb_hbond(rsd1.seqpos())) {
			 identify_hbonds_1way(
				*database_,
				rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
				false /*calculate_derivative*/,
				true, false, true, true,
				*options_,
				emap);
		}
	}
}

void
HBondEnergy::sidechain_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap & emap
) const
{
	Size nbrs1 = pose.energies().tenA_neighbor_graph().get_node( rsd1.seqpos() )->num_neighbors_counting_self_static();
	Size nbrs2 = pose.energies().tenA_neighbor_graph().get_node( rsd2.seqpos() )->num_neighbors_counting_self_static();

	// membrane dependent hbond potential hack
	if ( options_->Mbhbond() && options_->mphbond() ) {

		identify_hbonds_1way_membrane(
			*database_,
			rsd1, rsd2, nbrs1, nbrs2,
			false /*calculate_derivative*/,
			true, true, true, false, *options_, emap,
			pose);
		identify_hbonds_1way_membrane(
			*database_,
			rsd2, rsd1, nbrs2, nbrs1,
			false /*calculate_derivative*/,
			true, true, true, false, *options_, emap,
			pose);

	} else {

		identify_hbonds_1way(
			*database_,
			rsd1, rsd2, nbrs1, nbrs2,
			false /*calculate_derivative*/,
			true, true, true, false, *options_, emap);
		identify_hbonds_1way(
			*database_,
			rsd2, rsd1, nbrs2, nbrs1,
			false /*calculate_derivative*/,
			true, true, true, false, *options_, emap);
	}

}


void
HBondEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
debug_assert( set1.resid() != set2.resid() );

	if ( options_->exclude_DNA_DNA() && pose.residue( set1.resid() ).is_DNA() && pose.residue( set2.resid() ).is_DNA() ) return;

	using namespace methods;
	using namespace hbtrie;
	using namespace trie;
	using EnergiesCacheableDataType::HBOND_SET;

	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table1( energy_table );
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table2( energy_table );

	temp_table1 = 0; temp_table2 = 0;

	// save weight information so that its available during tvt execution
	// and also the neighbor counts for the two residues.
	weights_ = weights;
	res1_ = set1.resid();
	res2_ = set2.resid();

	rotamer_seq_sep_ = pose.residue( set2.resid() ).polymeric_oriented_sequence_distance( pose.residue( set1.resid() ) );

	if ( true ) { // super_hacky
		hbonds::HBondSet const & hbond_set
			( static_cast< hbonds::HBondSet const & >
				( pose.energies().data().get( HBOND_SET )));
		res1_nb_ = hbond_set.nbrs( set1.resid() );
		res2_nb_ = hbond_set.nbrs( set2.resid() );
	} else {
		res1_nb_ = pose.energies().tenA_neighbor_graph().get_node( set1.resid() )->num_neighbors_counting_self();
		res2_nb_ = pose.energies().tenA_neighbor_graph().get_node( set2.resid() )->num_neighbors_counting_self();
	}
	HBondRotamerTrieCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( hbond_method ) ));
	HBondRotamerTrieCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( hbond_method ) ));

	//prepare_for_residue_pair( set1.resid(), set2.resid(), pose );
//debug_assert( rep_scoretype() == fa_rep || rep_scoretype() == coarse_fa_rep );

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp( new HBCountPairFunction );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	//std::cout << "BEGIN TVT " << set1.resid() << " and " << set2.resid() << std::endl;
	trie1->trie_vs_trie( *trie2, *cp, *this, temp_table1, temp_table2 );
	//std::cout << "END TVT " << set1.resid() << " and " << set2.resid() << std::endl;

	/// add in the energies calculated by the tvt alg.
	energy_table += temp_table1;
	//std::cout << "FINISHED evaluate_rotamer_pair_energies" << std::endl;

	/*
	// debug
	//using namespace pack;

	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table3( energy_table );
	temp_table3 = 0;
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set1.num_rotamers(); ii <= ii_end; ++ii ) {
		for ( Size jj = 1, jj_end = set2.num_rotamers(); jj <= jj_end; ++jj ) {
			emap.zero();
			residue_pair_energy( *set1.rotamer( ii ), *set2.rotamer( jj ), pose, sfxn, emap );
			temp_table3( jj, ii ) += weights.dot( emap );
			if ( std::abs( temp_table1( jj, ii ) - temp_table3( jj, ii )) > 0.001 ) {
				std::cout << "Residues " << set1.resid() << " & " << set2.resid() << " rotamers: " << ii << " & " << jj;
				std::cout << " tvt/reg discrepancy: tvt= " <<  temp_table1( jj, ii ) << " reg= " << temp_table3( jj, ii );
				std::cout << " delta: " << temp_table1( jj, ii ) - temp_table3( jj, ii ) << std::endl;
			}
		}
	}
	std::cout << "Finished RPE calcs for residues " << set1.resid() << " & " << set2.resid() << std::endl;
	*/
}



//@brief overrides default rotamer/background energy calculation and uses
// the trie-vs-trie algorithm instead
void
HBondEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & ,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	using namespace methods;
	using namespace hbtrie;
	using namespace trie;
	using EnergiesCacheableDataType::HBOND_SET;
	using EnergiesCacheableDataType::HBOND_TRIE_COLLECTION;

	if ( options_->exclude_DNA_DNA() && residue.is_DNA() && pose.residue( set.resid() ).is_DNA() ) return;

	// allocate space for the trie-vs-trie algorithm
	utility::vector1< core::PackerEnergy > temp_vector1( set.num_rotamers(), 0.0 );
	utility::vector1< core::PackerEnergy > temp_vector2( set.num_rotamers(), 0.0 );

	// save weight information so that its available during tvt execution
	weights_ = weights;
	res1_ = set.resid();      //Added by Parin Sripakdeevong (sripakpa@stanford.edu)
	res2_ = residue.seqpos(); //Added by Parin Sripakdeevong (sripakpa@stanford.edu)
	rotamer_seq_sep_ = pose.residue( residue.seqpos() ).polymeric_oriented_sequence_distance( pose.residue( set.resid() ) );

	if ( true ) { // super_hacky
		hbonds::HBondSet const & hbond_set
			( static_cast< hbonds::HBondSet const & >
				( pose.energies().data().get( HBOND_SET )));
		res1_nb_ = hbond_set.nbrs( set.resid() );
		res2_nb_ = hbond_set.nbrs( residue.seqpos() );
	} else {
		res1_nb_ = pose.energies().tenA_neighbor_graph().get_node( set.resid() )->num_neighbors_counting_self();
		res2_nb_ = pose.energies().tenA_neighbor_graph().get_node( residue.seqpos() )->num_neighbors_counting_self();
	}

	HBondRotamerTrieCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set.get_trie( hbond_method ) ) );
	HBondRotamerTrieCOP trie2 = ( static_cast< TrieCollection const & >
																( pose.energies().data().get( HBOND_TRIE_COLLECTION )) ).trie( residue.seqpos() );

	//fpd
	if (trie2 == NULL) return;

	TrieCountPairBaseOP cp( new HBCountPairFunction );

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

	/*
	//debug
	utility::vector1< Energy > temp_vector3( energy_vector.size(), 0.0f );
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		emap.zero();
		residue_pair_energy( *set.rotamer( ii ), residue, pose, sfxn, emap );
		temp_vector3[ ii ] += weights.dot( emap );
		if ( std::abs( temp_vector1[ ii ] - temp_vector3[ ii ]) > 0.001 ) {
			std::cout << "Residues " << set.resid() << " & " << residue.seqpos() << " rotamers: " << ii << " & bg";
			std::cout << " tvt/reg discrepancy: tvt= " << temp_vector1[ ii ] << " reg= " << temp_vector3[ ii ];
			std::cout << " delta: " << temp_vector1[ ii ] - temp_vector3[ ii ] << std::endl;
		}
	}
	std::cout << "Finished Rotamer BG calcs for residues " << set.resid() << " & " << residue.seqpos() << std::endl;
	*/

}


///
void
HBondEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	/// Don't add in bb/bb hbond energies during minimization
	if ( pose.energies().use_nblist() ) return;

	if (options_->decompose_bb_hb_into_pair_energies()) return;

	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >( pose.energies().data().get( HBOND_SET )));

	EnergyMap hbond_emap(totals);

	// the current logic is that we fill the hbond set with backbone
	// hbonds only at the beginning of scoring. this is done to setup
	// the bb-bb hbond exclusion logic. so the hbondset should only
	// include bb-bb hbonds.
	// but see get_hb_don_chem_type in hbonds_geom.cc -- that only
	// classifies protein backbone donors as backbone, and the energy
	// accumulation by type is influenced by that via HBeval_lookup
	//
	// the important thing is that there's no double counting, which
	// is I think true since both fill_hbond_set and get_rsd-rsd-energy
	// use atom_is_backbone to check...
//debug_assert( std::abs( bb_scE ) < 1e-3 && std::abs( scE ) < 1e-3 );


	// this is to replicate buggy behavior regarding protein-backbone -- dna-backbone hbonds
	Real original_bb_sc    = totals[ hbond_bb_sc ];
	Real original_sr_bb_sc = totals[ hbond_sr_bb_sc ];
	Real original_lr_bb_sc = totals[ hbond_lr_bb_sc ];
	Real original_sc       = totals[ hbond_sc ];
	Real original_intra    = totals[ hbond_intra ];
	// end replicate

	get_hbond_energies( hbond_set, totals );

	// begin replicate
	totals[ hbond_bb_sc ]    = original_bb_sc;
	totals[ hbond_sr_bb_sc ] = original_sr_bb_sc;
	totals[ hbond_lr_bb_sc ] = original_lr_bb_sc;
	totals[ hbond_sc ]       = original_sc;
	totals[ hbond_intra ]    = original_intra;
	// end replicate

}

///@brief HACK!  MAX_R defines the maximum donorH to acceptor distance.
// The atomic_interaction_cutoff method is meant to return the maximum distance
// between two *heavy atoms* for them to have a zero interaction energy.
// I am currently assuming a 1.35 A maximum distance between a hydrogen and the
// heavy atom it is bound to, stealing this number from the CYS.params file since
// the HG in CYS is much further from it's SG than aliphatic hydrogens are from their carbons.
// This is a bad idea.  Someone come up with a way to fix this!
//
// At 4.35 A interaction cutoff, the hbond energy function is incredibly short ranged!
Distance
HBondEnergy::atomic_interaction_cutoff() const
{
	return MAX_R + 1.35; // MAGIC NUMBER
}

/// @brief the atomic interaction cutoff and the hydrogen interaction cutoff are the same.
Real
HBondEnergy::hydrogen_interaction_cutoff2() const
{
	return (MAX_R + 1.35) * ( MAX_R + 1.35 );
}


///@brief HBondEnergy is context sensitive
void
HBondEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & context_graphs_required
) const
{
	context_graphs_required[ ten_A_neighbor_graph ] = true;
}

bool
HBondEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	return ( weights[hbond_intra] > 0.0 || weights[hbond] > 0.0 );
}


void
HBondEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & ,
	EnergyMap & emap
) const
{
	if ( calculate_intra_res_hbonds( rsd, *options_ ) ){
		TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
		Size const rsd_nb = tenA_neighbor_graph.get_node( rsd.seqpos() )->num_neighbors_counting_self_static();
		identify_intra_res_hbonds( *database_, rsd, rsd_nb, *options_, emap);
	}
}

void
create_rotamer_descriptor(
	conformation::Residue const & res,
	hbonds::HBondOptions const & options,
	hbonds::HBondSet const & hbond_set,
	trie::RotamerDescriptor< hbtrie::HBAtom, hbtrie::HBCPData > & rotamer_descriptor
)
{
	using namespace trie;
	using namespace hbtrie;

	Size const resid( res.seqpos() );

	Size n_to_add(0);
	utility::vector1< int > add_to_trie( res.natoms(), 0 );
	for ( Size jj = 1; jj <= res.natoms(); ++jj ) {
		if ( res.atom_type_set()[ res.atom( jj ).type() ].is_acceptor() ) {
			//std::cout << "acc=" << jj << " ";
			add_to_trie[ jj ] = 1;
		} //else if ( res.atom_type_set()[ res.atom( jj ).type() ].is_hydrogen() &&
			else if ( res.atom_is_hydrogen( jj ) &&
				res.atom_type_set()[ res.atom( res.type().atom_base( jj ) ).type() ].is_donor()) {
			add_to_trie[ jj ] = 1;
			add_to_trie[ res.type().atom_base( jj ) ] = 1;
			//std::cout << "don=" << jj << " donb= " << res.type().atom_base( jj ) << " ";
		}
	}
	//std::cout << "rotamer trie for residue: " << res.type().name() << std::endl;
	for ( Size jj = 1; jj <= res.natoms(); ++jj ) {
		if ( add_to_trie[ jj ] == 1 ){
			++n_to_add;
			//std::cout << jj << " ";
		}
	}
	//std::cout << std::endl;

	if ( n_to_add == 0 ) {
		// What happens if there are NO hydrogen bonding atoms?  It would be nice if we didn't have to create
		// a trie at all, but I'm pretty sure the indexing logic requires that we have a place-holder rotamer.
		add_to_trie[ 1 ] = 1; ++n_to_add;
	}
	rotamer_descriptor.natoms( n_to_add );

	Size count_added_atoms = 0;
	for ( Size jj = 1; jj <= res.nheavyatoms(); ++jj ) {
		if ( add_to_trie[ jj ] == 0 ) continue;

		HBAtom newatom;
		HBCPData cpdata;

		newatom.xyz( res.atom(jj).xyz() );
		newatom.base_xyz( res.xyz( res.atom_base( jj )) );
		newatom.is_hydrogen( false );
		newatom.is_backbone( res.atom_is_backbone( jj ) ) ;

		//Following preserves hbond_sc, hbond_bb_sc, as preferred by other developers for protein/DNA.
		newatom.is_protein( res.is_protein() );
		newatom.is_dna( res.is_DNA() );

		if ( res.atom_type_set()[ res.atom( jj ).type() ].is_acceptor() ) {

			newatom.hb_chem_type( get_hb_acc_chem_type( jj, res ));
			//newatom.orientation_vector( create_acc_orientation_vector( res, jj ));
			//newatom.base_xyz( res.xyz( res.atom_base( jj )) );
			newatom.base2_xyz( res.xyz( res.abase2( jj )) );

			cpdata.is_sc( ! res.type().atom_is_backbone( jj ) );

			/// Count-pair data is responsible for enforcing the sc/bb hbond exclusion rule.
			/// If we're not using the rule, set "avoid_sc_hbonds" to false.
			cpdata.avoid_sc_hbonds( options.bb_donor_acceptor_check() &&
				! cpdata.is_sc() &&
				res.type().is_protein() &&
				hbond_set.acc_bbg_in_bb_bb_hbond( resid ) );
		}


		RotamerDescriptorAtom< HBAtom, HBCPData > rdatom( newatom, cpdata );
		rotamer_descriptor.atom( ++count_added_atoms, rdatom );

		for ( Size kk = res.attached_H_begin( jj ),
				kk_end = res.attached_H_end( jj );
				kk <= kk_end; ++kk ) {
			if ( add_to_trie[ kk ] == 0 ) continue;

			HBAtom newhatom;
			newhatom.xyz( res.atom(kk).xyz() );
			newhatom.base_xyz( res.xyz( res.atom_base( kk )) );
			newhatom.base2_xyz( Vector( 0.0, 0.0, 0.0 ) );
			//newhatom.orientation_vector( create_don_orientation_vector( res, kk ));
			newhatom.hb_chem_type( get_hb_don_chem_type( res.atom_base(kk), res ));
			newhatom.is_hydrogen( true );
			newhatom.is_backbone( res.atom_is_backbone( kk ) ) ;

			//Following preserves hbond_sc, hbond_bb_sc, as preferred by other developers for protein/DNA.
			newhatom.is_protein( res.is_protein() );
			newhatom.is_dna( res.is_DNA() );

			HBCPData hcpdata;
			hcpdata.is_sc( ! res.type().atom_is_backbone( kk ) );

			/// Count-pair data is responsible for enforcing the sc/bb hbond exclusion rule.
			/// If we're not using the rule, set "avoid_sc_hbonds" to false.
			hcpdata.avoid_sc_hbonds( options.bb_donor_acceptor_check() &&
				! hcpdata.is_sc() &&
				res.type().is_protein() &&
				hbond_set.don_bbg_in_bb_bb_hbond( resid ) );

			RotamerDescriptorAtom< HBAtom, HBCPData > hrdatom( newhatom, hcpdata );
			rotamer_descriptor.atom( ++count_added_atoms, hrdatom );

		}
	}
}

hbtrie::HBondRotamerTrieOP
HBondEnergy::create_rotamer_trie(
	conformation::RotamerSetBase const & rotset,
	pose::Pose const & pose
) const
{
	using namespace trie;
	using namespace hbtrie;
	using EnergiesCacheableDataType::HBOND_SET;

	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
		( pose.energies().data().get( HBOND_SET ) ) );

	//Size const resid( rotset.resid() );

	utility::vector1< RotamerDescriptor< HBAtom, HBCPData > > rotamer_descriptors( rotset.num_rotamers() );

	for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
		conformation::ResidueCOP ii_rotamer( rotset.rotamer( ii ) );
		create_rotamer_descriptor( *ii_rotamer, *options_, hbond_set, rotamer_descriptors[ ii ] );
		rotamer_descriptors[ ii ].rotamer_id( ii );
	}

	sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );

	return hbtrie::HBondRotamerTrieOP( new RotamerTrie< HBAtom, HBCPData >( rotamer_descriptors, atomic_interaction_cutoff()) );

}

hbtrie::HBondRotamerTrieOP
HBondEnergy::create_rotamer_trie(
	conformation::Residue const & res,
	pose::Pose const & pose
) const
{
	using namespace trie;
	using namespace hbtrie;
	using EnergiesCacheableDataType::HBOND_SET;

	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
		( pose.energies().data().get( HBOND_SET ) ) );

	utility::vector1< RotamerDescriptor< HBAtom, HBCPData > > rotamer_descriptors( 1 );

	create_rotamer_descriptor( res, *options_, hbond_set, rotamer_descriptors[ 1 ] );
	rotamer_descriptors[ 1 ].rotamer_id( 1 );

	return hbtrie::HBondRotamerTrieOP( new RotamerTrie< HBAtom, HBCPData >( rotamer_descriptors, atomic_interaction_cutoff()) );

}


///@brief code to evaluate a hydrogen bond energy for the trie that
/// didn't belong in the header itself -- it certainly does enough work
/// such that inlining it would not likely produce a speedup.
Energy
HBondEnergy::drawn_out_heavyatom_hydrogenatom_energy(
	hbtrie::HBAtom const & at1, // atom 1 is the heavy atom, the acceptor
	hbtrie::HBAtom const & at2, // atom 2 is the hydrogen atom, the donor
	bool flipped // is at1 from residue 1?
) const
{

	// When acc and don are both polymers and on the same chain:
	// ss = acc.seqpos - don.seqpos
	int ss = (flipped ? -rotamer_seq_sep_ : rotamer_seq_sep_);
	HBEvalTuple hbe_type = hbond_evaluation_type( at2, 0,   // donor atom
		at1, ss); // acceptor atom
	Energy hbenergy;
	//hb_energy_deriv_u( database_, hbe_type, at2.xyz(), at2.xyz() /*apl -- donor atom coordinate goes here, but is only used for derivatives */,
	//	at2.orientation_vector(),
	//	at1.xyz(), at1.xyz() /* apl -- acceptor-base coordinate goes here, but is only used for derivatives */,
	//	at1.orientation_vector(),
	//	Vector(-1.0,-1.0,-1.0), // abase2 xyz -- this is now wrong
	//	hbenergy,
	//	false /*evaluate_derivative*/, DUMMY_DERIVS );

	hb_energy_deriv( *database_, *options_, hbe_type,
		at2.base_xyz(), at2.xyz(), // donor heavy atom, donor hydrogen,
		at1.xyz(), at1.base_xyz(), at1.base2_xyz(), // acceptor, acceptor base, acceptor base2
		hbenergy, false, DUMMY_DERIVS );
	//std::cout << "drawn_out_heavyatom_hydrogenatom_energy " << hbenergy << " " << at1 << " " << at2 << std::endl;

	if ( hbenergy >= MAX_HB_ENERGY ) return 0.0; // no hbond

	//std::cout << "drawn_out_heavyatom_hydrogenatom_energy " << hbenergy << " " << at1 << " " << at2 << std::endl;

	Real envweight( 1.0 );
	if ( options_->use_hb_env_dep() ){
		envweight = ( flipped ? get_environment_dependent_weight( hbe_type, res2_nb_, res1_nb_, *options_ ) :
			get_environment_dependent_weight( hbe_type, res1_nb_, res2_nb_, *options_ ));
	}
	//pba membrane specific correction
	if ( options_->Mbhbond() && options_->mphbond() ) {

		Real membrane_depth_dependent_weight( 1.0 );
		membrane_depth_dependent_weight = ( flipped ? get_membrane_depth_dependent_weight(normal_, center_, thickness_,
			steepness_, res2_nb_, res1_nb_, at2.xyz(), at1.xyz()) : get_membrane_depth_dependent_weight(normal_,
			center_, thickness_, steepness_, res2_nb_, res1_nb_, at2.xyz(), at1.xyz()) );

		envweight = membrane_depth_dependent_weight;
	}

	Real weighted_energy(hb_eval_type_weight(hbe_type.eval_type(), weights_, res1_==res2_) * hbenergy * envweight);

	//std::cout << "weighted energy: " << weighted_energy << " " << hbenergy << " " << envweight << std::endl;
	return weighted_energy;
}


/// @details
/// Version 2: 2011-06-27 Fixing chi2 SER/THR and chi3 TYR derivatives when they act as acceptors.
core::Size
HBondEnergy::version() const
{
	//return 1; // Initial versioning
	return 2;
}

} // hbonds
} // scoring
} // core
