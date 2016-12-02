// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LK_BallEnergy.cc
/// @brief  Orientation dependent variant of the LK Solvation using
/// @author Phil Bradley

#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>


// Unit headers
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <core/scoring/lkball/LK_BallEnergyCreator.hh>
#include <core/scoring/lkball/LK_BallInfo.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/residue_datacache.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/constraints/AngleConstraint.hh>

#include <basic/options/util.hh> // HACK
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/prof.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

#include <utility/vector1.functions.hh>

//trie
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/trie/CPDataCorrespondence.hh>
#include <core/scoring/trie/RotamerDescriptor.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <core/scoring/trie/TrieCollection.hh>
#include <core/scoring/trie/trie.functions.hh>

#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>

#include <core/scoring/etable/etrie/TrieCountPair1BC4.hh>
#include <core/scoring/etable/etrie/TrieCountPair1BC3.hh>
#include <core/scoring/etable/etrie/TrieCountPairAll.hh>
#include <core/scoring/etable/etrie/TrieCountPairNone.hh>
#include <core/scoring/etable/etrie/TrieCountPairGeneric.hh>


#ifdef WIN32
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#endif

namespace core {
namespace scoring {
namespace lkball {


static THREAD_LOCAL basic::Tracer TR("core.scoring.methods.LK_BallEnergy");


/// @details This must return a fresh instance of the LK_ball class,
/// never an instance already in use
methods::EnergyMethodOP
LK_BallEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new LK_BallEnergy( options ) );
}

ScoreTypes
LK_BallEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( lk_ball_iso );
	sts.push_back( lk_ball );
	sts.push_back( lk_ball_wtd );
	sts.push_back( lk_ball_bridge );
	sts.push_back( lk_ball_bridge_uncpl );
	return sts;
}

/// inline retrieval functions here:
inline
LKB_RotamerSetInfo const &
retrieve_lkb_rotamer_set_info( conformation::RotamerSetBase const & set ) {
	return static_cast< LKB_RotamerSetInfo const & >
		( set.data().get( conformation::RotamerSetCacheableDataType::LK_BALL_ROTAMER_SET_INFO ) );
}

inline
LKB_ResidueInfo const &
retrieve_lkb_residue_info( pose::Pose const & pose, Size const seqpos ) {
	debug_assert( seqpos && seqpos <= ( static_cast< LKB_PoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::LK_BALL_POSE_INFO ) )).size() );
	return ( static_cast< LKB_PoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::LK_BALL_POSE_INFO ) )[ seqpos ] );
}

inline
LKB_ResidueInfo &
retrieve_nonconst_lkb_residue_info( pose::Pose & pose, Size const seqpos ) {
	debug_assert( seqpos && seqpos <= ( static_cast< LKB_PoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::LK_BALL_POSE_INFO ) )).size() );
	return ( static_cast< LKB_PoseInfo & >
		( pose.data().get( pose::datacache::CacheableDataType::LK_BALL_POSE_INFO ) )[ seqpos ] );
}


class LKB_ResPairMinData : public basic::datacache::CacheableData {
public:
	LKB_ResPairMinData();
	virtual ~LKB_ResPairMinData() {}
	virtual basic::datacache::CacheableDataOP clone() const
	{ return basic::datacache::CacheableDataOP( new LKB_ResPairMinData( *this ) ); }

	void
	initialize(
		LKB_ResidueInfoCOP res1_data,
		LKB_ResidueInfoCOP res2_data
	);

	LKB_ResidueInfo const & res1_data() const { return *res1_data_; }
	LKB_ResidueInfo const & res2_data() const { return *res2_data_; }

	bool
	initialized() const { return initialized_; }


private:

	LKB_ResidueInfoCOP res1_data_;
	LKB_ResidueInfoCOP res2_data_;

	bool initialized_;
};

typedef utility::pointer::shared_ptr< LKB_ResPairMinData >       LKB_ResPairMinDataOP;
typedef utility::pointer::shared_ptr< LKB_ResPairMinData const > LKB_ResPairMinDataCOP;

LKB_ResPairMinData::LKB_ResPairMinData():
	initialized_( false )
{}

void
LKB_ResPairMinData::initialize(
	LKB_ResidueInfoCOP res1_data,
	LKB_ResidueInfoCOP res2_data
)
{
	initialized_ = true;
	res1_data_ = res1_data;
	res2_data_ = res2_data;
}


/////////////////////////////////////// mindata retrieval functions
inline
LKB_ResPairMinData &
retrieve_nonconst_lkb_pairdata(
	ResPairMinimizationData & pairdata
)
{
	LKB_ResPairMinDataOP lkb_pairdata(0);
	if ( pairdata.get_data( lkb_respair_data ) ) {
		debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResPairMinData > ( pairdata.get_data( lkb_respair_data )));
		lkb_pairdata = utility::pointer::static_pointer_cast< LKB_ResPairMinData > ( pairdata.get_data( lkb_respair_data ));
	} else {
		lkb_pairdata = LKB_ResPairMinDataOP( new LKB_ResPairMinData );
		pairdata.set_data( lkb_respair_data, lkb_pairdata );
	}
	return *lkb_pairdata;
}

/////////////////////////////////////// mindata retrieval functions
inline
LKB_ResPairMinData const &
retrieve_lkb_pairdata(
	ResPairMinimizationData const & pairdata
)
{
	debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResPairMinData const > ( pairdata.get_data( lkb_respair_data )));
	return ( static_cast< LKB_ResPairMinData const & > ( pairdata.get_data_ref( lkb_respair_data ) ) );
}
/////////////////////////////////////// mindata retrieval functions
inline
LKB_ResidueInfo &
retrieve_nonconst_lkb_resdata(
	ResSingleMinimizationData & resdata
)
{
	LKB_ResidueInfoOP lkb_resdata( 0 );
	if ( resdata.get_data( lkb_res_data ) ) {
		debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResidueInfo > ( resdata.get_data( lkb_res_data )));
		lkb_resdata = utility::pointer::static_pointer_cast< LKB_ResidueInfo > ( resdata.get_data( lkb_res_data ) );
	} else {
		lkb_resdata = LKB_ResidueInfoOP( new LKB_ResidueInfo );
		resdata.set_data( lkb_res_data, lkb_resdata );
	}
	return *lkb_resdata;
}
/////////////////////////////////////// mindata retrieval functions
inline
LKB_ResidueInfo const &
retrieve_lkb_resdata(
	ResSingleMinimizationData const & resdata
)
{
	debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResidueInfo const > ( resdata.get_data( lkb_res_data )));
	return ( static_cast< LKB_ResidueInfo const & > ( resdata.get_data_ref( lkb_res_data ) ) );
}
/////////////////////////////////////// mindata retrieval functions
inline
LKB_ResidueInfoCOP
retrieve_lkb_resdata_ptr(
	ResSingleMinimizationData const & resdata
)
{
	debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResidueInfo const > ( resdata.get_data( lkb_res_data )));
	return ( utility::pointer::static_pointer_cast< LKB_ResidueInfo const > ( resdata.get_data( lkb_res_data ) ) );
}


LK_BallEnergy::LK_BallEnergy( methods::EnergyMethodOptions const & options ):
	parent             ( methods::EnergyMethodCreatorOP( new LK_BallEnergyCreator ) ),
	etable_            ( ScoringManager::get_instance()->etable( options ).lock() ),
	solv1_             ( ScoringManager::get_instance()->etable( options ).lock()->solv1()),
	solv2_             ( ScoringManager::get_instance()->etable( options ).lock()->solv2()),
	dsolv1_            ( ScoringManager::get_instance()->etable( options ).lock()->dsolv1()),
	etable_bins_per_A2_( ScoringManager::get_instance()->etable( options ).lock()->get_bins_per_A2() ),
	slim_etable_       ( basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ] ),
	use_intra_dna_cp_crossover_4_( true ),
	ramp_width_A2_     ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_ramp_width_A2 ]() ),
	overlap_width_A2_  ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_overlap_width_A2 ]() ),
	multi_water_fade_  ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_water_fade ]() ),
	lkbridge_angle_widthscale_ ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_bridge_angle_widthscale ]() ),
	save_bridging_waters_ (false)
{
	// overlap target length
	core::Real cos_overlap_base_angle = -1.0/3.0;
	core::Real heavyatom_water_len2 = 2.65; // MUST MATCH LK_BallInfo.cc!!
	if ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_donor ].user() ) {
		heavyatom_water_len2 = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_donor ]();
	}
	overlap_target_len_ = std::sqrt( heavyatom_water_len2 + heavyatom_water_len2 - 2 * heavyatom_water_len2 * cos_overlap_base_angle );

	setup_d2_bounds();
	if ( !slim_etable_ ) { runtime_assert( solv1_.size()>0 ); }

	core::Real eps = 0.1;
	lkb_max_dis_ = std::max( etable_->max_dis(), 2*heavyatom_water_len2 + std::sqrt( overlap_width_A2_ ) + eps );
	lkb_max_dis2_ = lkb_max_dis_*lkb_max_dis_;
	fasol_max_dis2_ = etable_->get_safe_max_dis2();
}

Distance
LK_BallEnergy::atomic_interaction_cutoff() const
{
	return lkb_max_dis_;
}


/// clone
methods::EnergyMethodOP
LK_BallEnergy::clone() const
{
	return methods::EnergyMethodOP( new LK_BallEnergy( *this ) );
}


LK_BallEnergy::LK_BallEnergy( LK_BallEnergy const & src ):
	ContextIndependentTwoBodyEnergy( src ),
	etable_(src.etable_),
	solv1_( src.solv1_ ),
	solv2_( src.solv2_ ),
	dsolv1_( src.dsolv1_ ),
	etable_bins_per_A2_( src.etable_bins_per_A2_ ),
	lkb_max_dis_( src.lkb_max_dis_ ),
	lkb_max_dis2_( src.lkb_max_dis2_ ),
	fasol_max_dis2_( src.fasol_max_dis2_ ),
	slim_etable_( src.slim_etable_ ),
	use_intra_dna_cp_crossover_4_( src.use_intra_dna_cp_crossover_4_ ),
	ramp_width_A2_     ( src.ramp_width_A2_ ),
	overlap_width_A2_     ( src.overlap_width_A2_ ),
	multi_water_fade_  ( src.multi_water_fade_ ),
	lkbridge_angle_widthscale_ ( src.lkbridge_angle_widthscale_ ),
	overlap_target_len_  ( src.overlap_target_len_ ),
	save_bridging_waters_ (src.save_bridging_waters_)
{
	setup_d2_bounds();
}


void
compute_and_store_pose_waters(
	pose::Pose & pose
)
{
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
	// using namespace core::pack::rotamer_set; // WaterPackingInfo
	LKB_PoseInfoOP info( new LKB_PoseInfo() );
	for ( Size i=1; i<= pose.size(); ++i ) {
		info->append( LKB_ResidueInfoOP( new LKB_ResidueInfo( pose.residue(i) ) ) );
	}
	pose.data().set( pose::datacache::CacheableDataType::LK_BALL_POSE_INFO, info );
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// New minimization interface code ////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
LK_BallEnergy::minimize_in_whole_structure_context( pose::Pose const &pose ) const
{
	//return false;
	return pose.energies().use_nblist_auto_update();
}


void
LK_BallEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &, // scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	ResSingleMinimizationData & resdata
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	LKB_ResidueInfo & info( retrieve_nonconst_lkb_resdata( resdata ) );
	info.initialize( rsd.type() );
	info.build_waters( rsd );
}

void
LK_BallEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const &, // rsd1,
	conformation::Residue const &, // rsd2,
	pose::Pose const & pose,
	ScoreFunction const &, //scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	ResSingleMinimizationData const & res1data,
	ResSingleMinimizationData const & res2data,
	ResPairMinimizationData & pairdata
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	LKB_ResPairMinData & lkb_pairdata( retrieve_nonconst_lkb_pairdata( pairdata ) );
	lkb_pairdata.initialize( retrieve_lkb_resdata_ptr( res1data ),
		retrieve_lkb_resdata_ptr( res2data ) );
}


bool
LK_BallEnergy::requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & ) const
{
	return true;
}

void
LK_BallEnergy::setup_for_scoring_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &, // sfxn,
	ResSingleMinimizationData & resdata
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	LKB_ResidueInfo & info( retrieve_nonconst_lkb_resdata( resdata ) );
	if ( !info.matches_residue_type( rsd.type() ) ) {
		std::cout << "LK_BallEnergy::setup_for_scoring_for_residue:: lkb-info mismatch: " << info.residue_type().name() << ' ' <<
			rsd.type().name() << std::endl;
		info.initialize( rsd.type() );
	}
	info.build_waters( rsd ); // already initialized in setup for minimizing for rsd
}

bool
LK_BallEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const &  ) const
{
	return true;
}

void
LK_BallEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const
{
	/// just compute water locations
	setup_for_scoring_for_residue( rsd, pose, sfxn, min_data );
}


void
LK_BallEnergy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
	compute_and_store_pose_waters( pose ); // could check task and do only some

	//fpd trie
	using namespace trie;
	using namespace lkbtrie;
	TrieCollectionOP tries( new TrieCollection );
	tries->total_residue( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		// Do not compute energy for virtual residues.
		if ( pose.residue(ii).aa() == core::chemical::aa_vrt ) continue;

		LKBRotamerTrieOP one_rotamer_trie = create_rotamer_trie( pose.residue( ii ), pose );
		tries->trie( ii, one_rotamer_trie );
	}
	pose.energies().data().set( EnergiesCacheableDataType::LKB_TRIE_COLLECTION, tries );
}


void
LK_BallEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	/// update waters for residue that has changed during packing, eg in rotamer trials
	/// need to double-check the current logic on this...
	// retrieve_nonconst_lkb_residue_info( pose, resid ).build_waters( pose.residue( resid ) );
	LKB_ResidueInfo & info( retrieve_nonconst_lkb_residue_info( pose, resid ) );
	conformation::Residue const & rsd( pose.residue( resid ) );
	if ( !info.matches_residue_type( rsd.type() ) ) {
		std::cout << "LK_BallEnergy::update_residue_for_packing:: lkb-info mismatch: " << info.residue_type().name() << ' ' <<
			rsd.type().name() << std::endl;
		info.initialize( rsd.type() );
	}
	info.build_waters( rsd );

	//fpd trie
	using namespace trie;
	using namespace lkbtrie;
	LKBRotamerTrieOP one_rotamer_trie = create_rotamer_trie( pose.residue( resid ), pose );

	// grab non-const & of the cached tries and replace resid's trie with a new one.
	TrieCollection & trie_collection
		( static_cast< TrieCollection & > (pose.energies().data().get( EnergiesCacheableDataType::LKB_TRIE_COLLECTION )));
	trie_collection.trie( resid, one_rotamer_trie );
}


void
LK_BallEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const & //scfxn
) const
{
	pose.update_residue_neighbors();
	compute_and_store_pose_waters( pose );

	if ( save_bridging_waters_ ) {
		bridging_waters_.clear(  );
	}
}

void
LK_BallEnergy::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & rotamer_set
) const
{
	//TR.Trace << "prepare_rotamers_for_packing: " << rotamer_set.num_rotamers() << std::endl;

	/// create a rotamer set info object
	LKB_RotamerSetInfoOP info( new LKB_RotamerSetInfo );

	for ( Size n=1; n<= rotamer_set.num_rotamers(); ++n ) {
		conformation::ResidueOP rot( rotamer_set.nonconst_rotamer(n) );
		LKB_ResidueInfoOP rotinfo( new LKB_ResidueInfo( *rot ) );
		rot->nonconst_data_ptr()->set( conformation::residue_datacache::LK_BALL_INFO, rotinfo->clone() ); // DataCache::set() does not clone
		info->append( rotinfo );
	}

	rotamer_set.data().set( conformation::RotamerSetCacheableDataType::LK_BALL_ROTAMER_SET_INFO, info );

	//fpd trie
	using namespace trie;
	using namespace lkbtrie;
	LKBRotamerTrieOP rottrie = create_rotamer_trie( rotamer_set, pose );
	//std::cout << "--------------------------------------------------" << std::endl << " HBROTAMER TRIE: " << set.resid() << std::endl;
	//rottrie->print();
	rotamer_set.store_trie( methods::lkball_method, rottrie );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
LK_BallEnergy::setup_d2_bounds()
{
	Real const h2o_radius( 1.4 );
	chemical::AtomTypeSet const & atom_set( *(etable_->atom_set().lock() ) );
	d2_low_.resize( atom_set.n_atomtypes() );
	for ( Size i=1; i<= atom_set.n_atomtypes(); ++i ) {
		chemical::AtomType const & atype( atom_set[ i ] );
		Real const d2_high( numeric::square( h2o_radius + atype.lj_radius() ) ); // was 3.0 * 3.0
		d2_low_[ i ] = std::max( 0.0, d2_high - ramp_width_A2_ );                // was 2.1 * 2.1
		TR.Trace << "d2_low_high: " << atype.name() << ' ' <<
			std::sqrt( d2_low_[i] ) << ' ' << std::sqrt( d2_high ) << std::endl;
	}

	/// which atomtypes are "charged"? Narg, Nlys, OOC (for splitting into a separate score)

	atom_type_is_charged_.resize( atom_set.n_atomtypes() );
	for ( Size i=1; i<= atom_set.n_atomtypes(); ++i ) {
		std::string const & name( atom_set[ i ].name() );
		debug_assert( Size( atom_set.atom_type_index( name ) ) == i ); // sanity
		if ( name == "Narg" || name == "Nlys" || name == "OOC" ) atom_type_is_charged_[ i ] = true;
		else atom_type_is_charged_[ i ] = false;
		TR.Trace << "atom_type_is_charged_: " << name << ' ' << atom_type_is_charged_[i] << std::endl;
	}

	lk_ball_prefactor_.clear();
	lk_ball_prefactor_.resize( atom_set.n_atomtypes(), 1.0 );
}


/// @details  Stolen from LK_SigmoidalFunc in lk_hack
/// d2_delta = d2 - d2_low
///
Real
LK_BallEnergy::eval_lk_fraction( Real const d2_delta, Real const width ) const
{
	debug_assert( d2_delta >= -0.001 && d2_delta <= width + 0.001 );
	static Real const inv_range( 1.0 / width );
	Real const xprime( inv_range * d2_delta );
	return ( 1 - xprime*xprime ) * ( 1 - xprime*xprime );
}

Real
LK_BallEnergy::eval_d_lk_fraction_dr_over_r( Real const d2_delta, Real const width ) const
{
	debug_assert( d2_delta >= -0.001 && d2_delta <= width + 0.001 );
	static Real const inv_range( 1.0 / width );
	Real const xprime( inv_range * d2_delta );
	return -8.0 * inv_range * ( 1 - xprime * xprime ) * xprime;  //?? why 8 (and not 4)
}


//// helper
Real
LK_BallEnergy::get_lk_fractional_contribution_for_single_water(
	Vector const & atom2_xyz,
	Size const atom2_type,
	Vector const & atom1_water
) const
{
	//Real max_frac( 0.0 ), frac, d2;
	Real const d2_low( d2_low_[ atom2_type ] );

	// find the closest water:
	Real const d2_delta = atom2_xyz.distance_squared( atom1_water ) - d2_low;
	if ( d2_delta > ramp_width_A2_ ) return 0.0;
	else if ( d2_delta < 0.0 ) return 1.0;
	else return eval_lk_fraction( d2_delta, ramp_width_A2_ );
}

// get water-water bridging contribution
Real
LK_BallEnergy::get_lkbr_fractional_contribution(
	Vector const & atom1_base,
	Vector const & atom2_base,
	Vectors const & atom1_waters,
	Vectors const & atom2_waters,
	Real const & lk_desolvation_sum,
	Real const & lkbr_wt,
	Real const & lkbr_uncpl_wt
) const {
	utility::vector1< numeric::xyzVector<core::Real> > dweighted_water_ddi;
	Real weighted_water_d2_delta(0.0), angleterm_lkbr(0.0), pointterm_lkbr(0.0), d_angleterm_lkbr_dr(0.0);
	return get_lkbr_fractional_contribution(
		atom1_base, atom2_base, atom1_waters, atom2_waters,
		dweighted_water_ddi, weighted_water_d2_delta, pointterm_lkbr, angleterm_lkbr, d_angleterm_lkbr_dr, lk_desolvation_sum, lkbr_wt, lkbr_uncpl_wt, false);
}


// get water-water bridging contribution
// may want to split out deriv-computing and non-deriv computing version since it is a bit more expensive than in non-br lk case
Real
LK_BallEnergy::get_lkbr_fractional_contribution(
	Vector const & atom1_base,
	Vector const & atom2_base,
	Vectors const & atom1_waters,
	Vectors const & atom2_waters,
	utility::vector1< numeric::xyzVector<core::Real> > & d_weighted_d2_d_di,  // derivative of weighted distance w.r.t. water 1 movement
	Real & weighted_d2_water_delta,
	Real & pointterm_lkbr,
	Real & angleterm_lkbr,
	Real & d_angleterm_lkbr_dr,
	Real const & lk_desolvation_sum,
	Real const & lkbr_wt,
	Real const & lkbr_uncpl_wt,
	bool compute_derivs /*=true*/
) const {
	pointterm_lkbr = 0;
	angleterm_lkbr = 0;
	d_angleterm_lkbr_dr = 0;
	weighted_d2_water_delta = 0.0;
	if ( compute_derivs ) {
		d_weighted_d2_d_di.clear();
		d_weighted_d2_d_di.resize( atom1_waters.size(), numeric::xyzVector<core::Real>(0.0,0.0,0.0) );
	}

	if ( atom1_waters.empty() || atom2_waters.empty() ) return 0.0;

	// softmax of all water combinations
	Real d2_delta, d2_delta_wt;
	for ( Size idx1 = 1; idx1 <= atom1_waters.size(); ++idx1 ) {
		for ( Size idx2 = 1; idx2 <= atom2_waters.size(); ++idx2 ) {
			numeric::xyzVector< core::Real > delta_ij = atom1_waters[idx1] - atom2_waters[idx2];
			d2_delta = delta_ij.length_squared();
			d2_delta_wt = exp( -d2_delta/multi_water_fade_ );
			weighted_d2_water_delta += d2_delta_wt;

			if ( compute_derivs ) {
				d_weighted_d2_d_di[idx1] += d2_delta_wt * delta_ij;
			}

			if ( save_bridging_waters_ ) {
				Real score = 0.0, lkbridge_frac = 0.0;
				if ( d2_delta < overlap_width_A2_ ) lkbridge_frac = d2_delta < 0.0 ? Real( 1.0 ) : eval_lk_fraction( d2_delta, overlap_width_A2_ );
				core::Real lkbr_uncpl_score = lkbr_uncpl_wt*lkbridge_frac;
				core::Real lkbr_score = lkbr_wt*lk_desolvation_sum*lkbridge_frac;
				score = ( lkbr_uncpl_score + lkbr_score );
				if ( score < -1e-4 ) {
					bridging_waters_.push_back( ScoredBridgingWater ( lkbr_uncpl_score, lkbr_score, 0.5*(atom1_waters[idx1] + atom2_waters[idx2]) ) );
				}
			}
		}
	}

	if ( compute_derivs ) {
		for ( Size idx = 1; idx <= atom1_waters.size(); ++idx ) {
			d_weighted_d2_d_di[idx] /= weighted_d2_water_delta;
		}
	}

	weighted_d2_water_delta = -multi_water_fade_ * log( weighted_d2_water_delta );

	Real frac( 0.0 );
	if ( weighted_d2_water_delta < overlap_width_A2_ ) {
		frac = ( weighted_d2_water_delta < 0.0 ? Real( 1.0 ) : eval_lk_fraction( weighted_d2_water_delta, overlap_width_A2_ ) );
	}

	pointterm_lkbr = frac;

	// now multiply by base angle (assumption! .. all water dists are equal!)
	if ( lkbridge_angle_widthscale_!=0 ) {
		core::Real act_len = (atom1_base-atom2_base).length();
		core::Real base_atom_delta = std::abs( overlap_target_len_-act_len );
		if ( base_atom_delta > overlap_width_A2_ ) {
			frac = 0;
			angleterm_lkbr = d_angleterm_lkbr_dr = 0;
		} else {
			core::Real angle_constraint_width = lkbridge_angle_widthscale_ * overlap_width_A2_;
			angleterm_lkbr = eval_lk_fraction( base_atom_delta, angle_constraint_width);

			frac *= angleterm_lkbr;

			if ( compute_derivs ) {
				core::Real base_atom_delta_sign = (overlap_target_len_ > act_len )?1.0:-1.0;
				d_angleterm_lkbr_dr = base_atom_delta_sign*eval_d_lk_fraction_dr_over_r(base_atom_delta, angle_constraint_width);
			}
		}
	} else {
		angleterm_lkbr = 1.0;
		d_angleterm_lkbr_dr = 0;
	}
	return frac;
}

/// @note  closest_water may be set to 0 upon return if none of the waters are within ramp_width_A2_
Real
LK_BallEnergy::get_lk_fractional_contribution(
	Vector const & atom2_xyz,
	Size const atom2_type,
	Vectors const & atom1_waters,
	utility::vector1< Real > & d_weighted_d2_d_di,  // per water contribution
	Real & weighted_d2_water_delta
) const
{
	Real const d2_low( d2_low_[ atom2_type ] );

	// softmax of closest water
	d_weighted_d2_d_di.resize( atom1_waters.size() );

	weighted_d2_water_delta = 0.0;
	Real d2_delta;
	for ( Size idx = 1; idx <= atom1_waters.size(); ++idx ) {
		d2_delta = atom2_xyz.distance_squared( atom1_waters[idx] ) - d2_low;
		d_weighted_d2_d_di[idx] = exp( -d2_delta/multi_water_fade_ );
		weighted_d2_water_delta += d_weighted_d2_d_di[idx];
		//TR << "d2_delta = " << d2_delta << std::endl;
	}

	for ( Size idx = 1; idx <= atom1_waters.size(); ++idx ) {
		d_weighted_d2_d_di[idx] /= weighted_d2_water_delta;
	}

	weighted_d2_water_delta = -multi_water_fade_ * log( weighted_d2_water_delta );
	//TR << "weighted_d2_water_delta = " << weighted_d2_water_delta << std::endl;

	Real frac( 0.0 );
	if ( weighted_d2_water_delta < ramp_width_A2_ ) {
		frac = ( weighted_d2_water_delta < 0.0 ? Real( 1.0 ) : eval_lk_fraction( weighted_d2_water_delta, ramp_width_A2_ ) );
	}
	//TR << "frac = " << frac << std::endl;

	return frac;
}


//fpd  this function largely duplicates the code above
//fpd  this version is called during scoring, while the above is called during minimizaiton
//fpd  the overhead of creating/deleting the per-water derivative vector is significant
//fpd     during packing, so we avoid that here
Real
LK_BallEnergy::get_lk_fractional_contribution(
	Vector const & atom2_xyz,
	Size const atom2_type,
	Vectors const & atom1_waters
) const
{
	Real const d2_low( d2_low_[ atom2_type ] );

	// softmax of closest water
	Real weighted_d2_water_delta = 0.0;
	Real d2_delta;
	for ( Size idx = 1; idx <= atom1_waters.size(); ++idx ) {
		d2_delta = atom2_xyz.distance_squared( atom1_waters[idx] ) - d2_low;
		weighted_d2_water_delta += exp( -d2_delta/multi_water_fade_ );
	}

	weighted_d2_water_delta = -multi_water_fade_ * log( weighted_d2_water_delta );

	Real frac( 0.0 );
	if ( weighted_d2_water_delta < ramp_width_A2_ ) {
		frac = ( weighted_d2_water_delta < 0.0 ? Real( 1.0 ) : eval_lk_fraction( weighted_d2_water_delta, ramp_width_A2_ ) );
	}
	return frac;
}

/// This guy is used during scoring if we are not minimizing, or inside linmem_ig packing and things like that...
void
LK_BallEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &sf,
	EnergyMap & emap
) const
{
	/// there might be data stashed in these residues if we came through certain packing routes
	using conformation::residue_datacache::LK_BALL_INFO;
	bool rsd1_cached = rsd1.data_ptr() != 0 && rsd1.data_ptr()->has( LK_BALL_INFO );
	bool rsd2_cached = rsd2.data_ptr() != 0 && rsd2.data_ptr()->has( LK_BALL_INFO );

	if ( rsd1_cached && rsd2_cached ) {
		residue_pair_energy( rsd1,
			*( dynamic_cast< LKB_ResidueInfo const * >( rsd1.data_ptr()->get_raw_const_ptr( LK_BALL_INFO ))),
			rsd2,
			*( dynamic_cast< LKB_ResidueInfo const * >( rsd2.data_ptr()->get_raw_const_ptr( LK_BALL_INFO ))),
			sf, emap );
	} else if ( rsd1_cached ) {
		LKB_ResidueInfo const & info2( retrieve_lkb_residue_info( pose, rsd2.seqpos() ) );
		residue_pair_energy( rsd1,
			*( dynamic_cast< LKB_ResidueInfo const * >( rsd1.data_ptr()->get_raw_const_ptr( LK_BALL_INFO ))),
			rsd2, info2, sf, emap );
	} else if ( rsd2_cached ) {
		LKB_ResidueInfo const & info1( retrieve_lkb_residue_info( pose, rsd1.seqpos() ) );
		residue_pair_energy( rsd1, info1, rsd2,
			*( dynamic_cast< LKB_ResidueInfo const * >( rsd2.data_ptr()->get_raw_const_ptr( LK_BALL_INFO ))),
			sf, emap );
	} else {
		LKB_ResidueInfo const & info1( retrieve_lkb_residue_info( pose, rsd1.seqpos() ) );
		LKB_ResidueInfo const & info2( retrieve_lkb_residue_info( pose, rsd2.seqpos() ) );
		residue_pair_energy( rsd1, info1, rsd2, info2, sf, emap );
	}
}

void
LK_BallEnergy::sidechain_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &, // pose,
	ScoreFunction const &sf,
	EnergyMap & emap
) const
{
	using conformation::residue_datacache::LK_BALL_INFO;

	/// if we got here we should have come through packing...
	debug_assert( rsd1.data_ptr() != 0 &&
		rsd1.data_ptr()->get_const_ptr( LK_BALL_INFO ) != 0 &&
		dynamic_cast< LKB_ResidueInfo const * > ( rsd1.data_ptr()->get_raw_const_ptr( LK_BALL_INFO )));
	debug_assert( rsd2.data_ptr() != 0 &&
		rsd2.data_ptr()->get_const_ptr( LK_BALL_INFO ) != 0 &&
		dynamic_cast< LKB_ResidueInfo const * > ( rsd2.data_ptr()->get_raw_const_ptr( LK_BALL_INFO )));

	residue_pair_energy( rsd1,
		*( dynamic_cast< LKB_ResidueInfo const * >( rsd1.data_ptr()->get_raw_const_ptr( LK_BALL_INFO ))),
		rsd2,
		*( dynamic_cast< LKB_ResidueInfo const * >( rsd2.data_ptr()->get_raw_const_ptr( LK_BALL_INFO ))),
		sf, emap );
}

// this guy is used during minimization
bool
LK_BallEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}
void
LK_BallEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & pairdata,
	pose::Pose const & pose,
	ScoreFunction const &sf,
	EnergyMap & emap
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	LKB_ResPairMinData const & lkb_pairdata( retrieve_lkb_pairdata( pairdata ) );

	residue_pair_energy( rsd1, lkb_pairdata.res1_data(), rsd2, lkb_pairdata.res2_data(), sf, emap );
}


/// helper function for outsiders
Real
LK_BallEnergy::calculate_lk_desolvation_of_single_atom_by_residue(
	Size const atom1,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
)
{
	using namespace etable::count_pair;
	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	// setup residue information
	Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
	Size const atom1_type_index( rsd1.atom( atom1 ).type() );

	Real total_desolvation( 0.0 );
	for ( Size atom2=1; atom2<= rsd2.nheavyatoms(); ++atom2 ) {
		Vector const & atom2_xyz( rsd2.xyz( atom2 ) );

		Real cp_weight = 1.0;
		Size pathdist;
		if ( ! cpfxn->count( atom1, atom2, cp_weight, pathdist ) ) continue;

		Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

		if ( ( d2 >= fasol_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

		Real lk_desolvation_of_atom1_by_atom2;
		if ( slim_etable_ ) {
			Real lk_desolvation_of_atom2_by_atom1;
			// not sure the order is correct here:
			etable_->analytic_lk_energy( rsd1.atom( atom1 ), rsd2.atom( atom2 ), lk_desolvation_of_atom1_by_atom2,
				lk_desolvation_of_atom2_by_atom1 );
			lk_desolvation_of_atom1_by_atom2 *= cp_weight;
		} else {
			// setup for solvation Etable lookups
			Size const atom2_type_index( rsd2.atom( atom2 ).type() );
			Real const d2_bin = d2 * etable_bins_per_A2_;
			int disbin = static_cast< int >( d2_bin ) + 1;
			Real frac = d2_bin - ( disbin - 1 );
			int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

			lk_desolvation_of_atom1_by_atom2 = cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] );
		}
		total_desolvation += lk_desolvation_of_atom1_by_atom2;
	}
	return total_desolvation;
}

/// get the lk-ball desolvation of atom1 by atom2, and the unoriented lk desolvation of atom1 by atom2

void
LK_BallEnergy::calculate_lk_ball_atom_energies(
	Size const atom1,
	conformation::Residue const & rsd1,
	Vectors const & atom1_waters,
	Size const atom2,
	conformation::Residue const & rsd2,
	Real & lk_desolvation_of_atom1_by_atom2,
	Real & lk_ball_desolvation_of_atom1_by_atom2 // includes lk-fraction
) const
{
	using namespace etable::count_pair;

	// initialize values
	lk_desolvation_of_atom1_by_atom2 = 0.0;
	lk_ball_desolvation_of_atom1_by_atom2 = 0.0;

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	Real cp_weight = 1.0; Size pathdist;
	if ( !cpfxn->count( atom1, atom2, cp_weight, pathdist ) ) return; /// NO COUNTPAIR

	Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
	Vector const & atom2_xyz( rsd2.xyz( atom2 ) );

	Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

	if ( ( d2 >= fasol_max_dis2_) || ( d2 == Real(0.0) ) ) return; // TOO FARAWAY (OR SAME ATOM?)

	Size const atom2_type_index( rsd2.atom( atom2 ).type() );
	if ( slim_etable_ ) {
		Real lk_desolvation_of_atom2_by_atom1;
		// note sure the order is correct here:
		etable_->analytic_lk_energy( rsd1.atom( atom1 ), rsd2.atom( atom2 ), lk_desolvation_of_atom1_by_atom2,
			lk_desolvation_of_atom2_by_atom1 );
		lk_desolvation_of_atom1_by_atom2 *= cp_weight;
	} else {
		// setup for solvation Etable lookups
		Size const atom1_type_index( rsd1.atom( atom1 ).type() );

		Real const d2_bin = d2 * etable_bins_per_A2_;
		int disbin = static_cast< int >( d2_bin ) + 1;
		Real frac = d2_bin - ( disbin - 1 );
		int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

		lk_desolvation_of_atom1_by_atom2 = ( cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] ) );
	}
	lk_ball_desolvation_of_atom1_by_atom2 = lk_desolvation_of_atom1_by_atom2 *
		get_lk_fractional_contribution( atom2_xyz, atom2_type_index, atom1_waters );
}


void
LK_BallEnergy::calculate_lk_ball_atom_energies_cp(
	Size const atom1,
	conformation::Residue const & rsd1,
	Vectors const & atom1_waters,
	Size const atom2,
	conformation::Residue const & rsd2,
	etable::count_pair::CPCrossoverBehavior const & cp_crossover,
	Real & lk_desolvation_of_atom1_by_atom2,
	Real & lk_ball_desolvation_of_atom1_by_atom2 // includes lk-fraction
) const
{
	using namespace etable::count_pair;

	// initialize values
	lk_desolvation_of_atom1_by_atom2 = 0.0;
	lk_ball_desolvation_of_atom1_by_atom2 = 0.0;

	CountPairFunctionOP cpfxn
		( rsd1.seqpos() == rsd2.seqpos() ?
		CountPairFactory::create_intrares_count_pair_function( rsd1, cp_crossover ) :
		CountPairFactory::create_count_pair_function( rsd1, rsd2, cp_crossover ) );

	if ( rsd1.seqpos() == rsd2.seqpos() && atom1 == atom2 ) return;

	Real cp_weight = 1.0; Size pathdist;
	if ( !cpfxn->count( atom1, atom2, cp_weight, pathdist ) ) return; /// NO COUNTPAIR

	Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
	Vector const & atom2_xyz( rsd2.xyz( atom2 ) );

	Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

	if ( ( d2 >= fasol_max_dis2_) || ( d2 == Real(0.0) ) ) return; // TOO FARAWAY (OR SAME ATOM?)

	Size const atom2_type_index( rsd2.atom( atom2 ).type() );
	if ( slim_etable_ ) {
		Real lk_desolvation_of_atom2_by_atom1;
		// note sure the order is correct here:
		etable_->analytic_lk_energy( rsd1.atom( atom1 ), rsd2.atom( atom2 ), lk_desolvation_of_atom1_by_atom2,
			lk_desolvation_of_atom2_by_atom1 );
		lk_desolvation_of_atom1_by_atom2 *= cp_weight;
	} else {
		// setup for solvation Etable lookups
		Size const atom1_type_index( rsd1.atom( atom1 ).type() );

		Real const d2_bin = d2 * etable_bins_per_A2_;
		int disbin = static_cast< int >( d2_bin ) + 1;
		Real frac = d2_bin - ( disbin - 1 );
		int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

		lk_desolvation_of_atom1_by_atom2 = ( cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] ) );
	}

	lk_ball_desolvation_of_atom1_by_atom2 = lk_desolvation_of_atom1_by_atom2 *
		get_lk_fractional_contribution( atom2_xyz, atom2_type_index, atom1_waters );
}


Real
LK_BallEnergy::calculate_lk_desolvation_of_single_atom_by_residue_no_count_pair(
	Size const atom1,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
)
{
	//using namespace etable::count_pair;

	// setup residue information
	Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
	Size const atom1_type_index( rsd1.atom( atom1 ).type() );

	Real total_desolvation( 0.0 );
	for ( Size atom2=1; atom2<= rsd2.nheavyatoms(); ++atom2 ) {
		Vector const & atom2_xyz( rsd2.xyz( atom2 ) );

		Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

		if ( ( d2 >= fasol_max_dis2_) || ( d2 < 1e-3 ) ) continue; // exclude self...

		Real lk_desolvation_of_atom1_by_atom2;
		if ( slim_etable_ ) {
			Real lk_desolvation_of_atom2_by_atom1;
			// note sure the order is correct here:
			etable_->analytic_lk_energy( rsd1.atom( atom1 ), rsd2.atom( atom2 ), lk_desolvation_of_atom1_by_atom2,
				lk_desolvation_of_atom2_by_atom1 );
		} else {
			// setup for solvation Etable lookups
			Size const atom2_type_index( rsd2.atom( atom2 ).type() );
			Real const d2_bin = d2 * etable_bins_per_A2_;
			int disbin = static_cast< int >( d2_bin ) + 1;
			Real frac = d2_bin - ( disbin - 1 );
			int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

			lk_desolvation_of_atom1_by_atom2 = ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ];
		}
		total_desolvation += lk_desolvation_of_atom1_by_atom2;
	}
	return total_desolvation;
}

void
LK_BallEnergy::eval_desolvation_derivs_no_count_pair(
	Real const d2,
	Size const atom1,
	conformation::Residue const & rsd1,
	Size const atom2,
	conformation::Residue const & rsd2,
	Real & atom1_lk_desolvation_by_atom2_deriv,
	Real & atom2_lk_desolvation_by_atom1_deriv
)
{
	if ( ( d2 >= fasol_max_dis2_) || ( d2 == Real(0.0) ) ) {
		atom1_lk_desolvation_by_atom2_deriv = 0.0;
		atom2_lk_desolvation_by_atom1_deriv = 0.0;
		return;
	}

	if ( slim_etable_ ) {
		Real inv_dis;
		etable_->analytic_lk_derivatives( rsd1.atom( atom1 ), rsd2.atom( atom2 ),
			atom1_lk_desolvation_by_atom2_deriv, atom2_lk_desolvation_by_atom1_deriv, inv_dis );
	} else {
		Real const d2_bin( d2 * etable_bins_per_A2_ );
		Size const disbin( static_cast< int >( d2_bin ) + 1 );
		Real const frac( d2_bin - ( disbin - 1 ) );

		/// this index into solv1 or dsolv1 should represent desolvation of atom1 by atom2
		int const linear_index12( dsolv1_.index( disbin, rsd2.atom( atom2 ).type(), rsd1.atom( atom1 ).type() ) );
		/// this index into solv1 or dsolv1 should represent desolvation of atom2 by atom1
		int const linear_index21( dsolv1_.index( disbin, rsd1.atom( atom1 ).type(), rsd2.atom( atom2 ).type() ) );

		atom1_lk_desolvation_by_atom2_deriv = dsolv1_[ linear_index12 ] * ( 1-frac ) + dsolv1_[ linear_index12 + 1 ] * frac;
		atom2_lk_desolvation_by_atom1_deriv = dsolv1_[ linear_index21 ] * ( 1-frac ) + dsolv1_[ linear_index21 + 1 ] * frac;
	}
}

/// HACKY HELPER FXN TAKEN DIRECTLY FROM BaseEtableEnergy.tmpl.hh
/// IF YOU CHANGE THE BEHAVIOR THERE YOU SHOULD CHANGE HERE AS WELL
///
scoring::etable::count_pair::CPCrossoverBehavior
determine_crossover_behavior(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	bool const use_intra_dna_cp_crossover_4
)
{
	using namespace scoring::etable::count_pair;
	// Ask "are these two residues polymers and do they share a polymeric bond?"
	if ( res1.is_polymer_bonded(res2) && res2.is_polymer_bonded(res1) ) {
		if ( ( !( res1.is_protein() && res2.is_protein() ) ) &&
				( !( use_intra_dna_cp_crossover_4 && res1.is_DNA() && res2.is_DNA() ) ) &&
				( !( res1.is_RNA() && res2.is_RNA() ) ) ) {
			return CP_CROSSOVER_3;
		} else {
			return CP_CROSSOVER_4; // peptide bond w/ or w/o rama, but definately w/o mm_twist
		}
	} else if ( res1.seqpos() == res2.seqpos() ) {
		// logic for controlling intra-residue count pair behavior goes here; for now, default to crossover 3
		return CP_CROSSOVER_3;
	} else {
		return CP_CROSSOVER_3; // e.g. disulfides where seqsep != 1
	}
}

void
LK_BallEnergy::accumulate_single_atom_contributions(
	Size const,
	Size const,
	Vectors const & atom1_waters,
	utility::vector1< Real > const & atom1_wts,
	conformation::Residue const &,
	Size const atom2_type_index,
	Vector const & atom2_xyz,
	Real const lk_desolvation_of_atom1_by_atom2,
	EnergyMap & emap
) const
{
	/// for first checkin to trunk, take old approach of only counting things for atoms with waters attached
	/// more logic here in the blab branch, hence existence as separate fxn
	if ( atom1_waters.empty() ) return;

	emap[ lk_ball_iso ] += lk_desolvation_of_atom1_by_atom2;
	Real const lk_desolvation_of_atom1_by_atom2_lkb
		( lk_desolvation_of_atom1_by_atom2 *
		get_lk_fractional_contribution( atom2_xyz, atom2_type_index, atom1_waters ) );
	emap[ lk_ball ] += lk_desolvation_of_atom1_by_atom2_lkb;
	emap[ lk_ball_wtd ] += ( atom1_wts[1] * lk_desolvation_of_atom1_by_atom2 +
		atom1_wts[2] * lk_desolvation_of_atom1_by_atom2_lkb );
}

void
LK_BallEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	LKB_ResidueInfo const & rsd1_info,
	conformation::Residue const & rsd2,
	LKB_ResidueInfo const & rsd2_info,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	//PROF_START( basic::LK_BALL_RESIDUE_PAIR_ENERGY );

	using namespace etable::count_pair;
	//bool const verbose( false );
	utility::vector1< Vectors > const & rsd1_waters( rsd1_info.waters() );
	utility::vector1< Vectors > const & rsd2_waters( rsd2_info.waters() );

	utility::vector1< utility::vector1< Real > > const & rsd1_atom_wts( rsd1_info.atom_weights() );
	utility::vector1< utility::vector1< Real > > const & rsd2_atom_wts( rsd2_info.atom_weights() );

	CPCrossoverBehavior crossover = determine_crossover_behavior( rsd1, rsd2, use_intra_dna_cp_crossover_4_ );
	CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, crossover );

	bool use_lkbr = (sfxn.get_weight( core::scoring::lk_ball_bridge )!=0);
	bool use_lkbr_uncpl = (sfxn.get_weight( core::scoring::lk_ball_bridge_uncpl)!=0);

	// setup residue information
	for ( Size atom1=1; atom1<= rsd1.nheavyatoms(); ++atom1 ) {
		Vectors const & atom1_waters( rsd1_waters[ atom1 ] );
		Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
		Size const atom1_type_index( rsd1.atom( atom1 ).type() );
		utility::vector1< Real > const & atom1_weights( rsd1_atom_wts[atom1] );

		for ( Size atom2=1; atom2<= rsd2.nheavyatoms(); ++atom2 ) {
			Vectors const & atom2_waters( rsd2_waters[ atom2 ] );
			Vector const & atom2_xyz( rsd2.xyz( atom2 ) );

			if ( atom1_waters.empty() && atom2_waters.empty() ) continue;

			utility::vector1< Real > const & atom2_weights( rsd2_atom_wts[atom2] );

			Real cp_weight = 1.0;
			Size pathdist;
			if ( ! cpfxn->count( atom1, atom2, cp_weight, pathdist ) ) continue;

			Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

			if ( d2 == Real(0.0) ) continue; // sanity check
			if ( !use_lkbr_uncpl && d2 >= fasol_max_dis2_ ) continue;
			if ( d2 >= lkb_max_dis2_ ) continue;

			Real lk_desolvation_of_atom1_by_atom2 = 0.0, lk_desolvation_of_atom2_by_atom1 = 0.0;

			// only lk_ball_bridge_uncpl goes beyond fasol_max_dis2_
			if ( d2 < fasol_max_dis2_ ) {
				Size const atom2_type_index( rsd2.atom( atom2 ).type() );
				if ( slim_etable_ ) {
					etable_->analytic_lk_energy( rsd1.atom( atom1 ), rsd2.atom( atom2 ), lk_desolvation_of_atom1_by_atom2,
						lk_desolvation_of_atom2_by_atom1 );
					lk_desolvation_of_atom1_by_atom2 *= cp_weight;
					lk_desolvation_of_atom2_by_atom1 *= cp_weight;
				} else {
					// setup for solvation Etable lookups
					Real const d2_bin = d2 * etable_bins_per_A2_;
					int disbin = static_cast< int >( d2_bin ) + 1;
					Real frac = d2_bin - ( disbin - 1 );
					int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

					lk_desolvation_of_atom1_by_atom2 = cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] );
					lk_desolvation_of_atom2_by_atom1 = cp_weight * ( ( 1. - frac ) * solv2_[ l1 ] + frac * solv2_[ l1+1 ] );
				}

				accumulate_single_atom_contributions( atom1, atom1_type_index, atom1_waters, atom1_weights,
					rsd1, atom2_type_index, atom2_xyz,
					lk_desolvation_of_atom1_by_atom2, emap );

				accumulate_single_atom_contributions( atom2, atom2_type_index, atom2_waters, atom2_weights,
					rsd2, atom1_type_index, atom1_xyz,
					lk_desolvation_of_atom2_by_atom1, emap );
			}

			// fpd - get lk_ball_bridge
			if ( atom1_waters.empty() || atom2_waters.empty() ) continue;
			if ( !use_lkbr_uncpl && !use_lkbr ) continue;

			Real lk_desolvation_sum = lk_desolvation_of_atom1_by_atom2+lk_desolvation_of_atom2_by_atom1;
			core::Real lkbr_wt = sfxn.get_weight( core::scoring::lk_ball_bridge );
			core::Real lkbr_uncpl_wt = sfxn.get_weight( core::scoring::lk_ball_bridge_uncpl );
			core::Real lkbridge_frac = get_lkbr_fractional_contribution(
				atom1_xyz, atom2_xyz,
				atom1_waters, atom2_waters,
				lk_desolvation_sum,
				lkbr_wt, lkbr_uncpl_wt );

			emap[ lk_ball_bridge ] += (lk_desolvation_of_atom1_by_atom2+lk_desolvation_of_atom2_by_atom1) * lkbridge_frac;
			emap[ lk_ball_bridge_uncpl ] += lkbridge_frac;
		} // atom2
	} // atom1

	//PROF_STOP( basic::LK_BALL_RESIDUE_PAIR_ENERGY );
}



/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////


void
LK_BallEnergy::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const &
) const
{
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
	pose.update_residue_neighbors();
	compute_and_store_pose_waters( pose );
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
}


/// @note  atom2 is desolvating atom1. atom1_waters is non-empty
/// @note  Pretend that atom1 is the atom whose derivs are being calculated. weight_factor may include -1 term
/// to switch the order...
///
void
LK_BallEnergy::sum_deriv_contributions_for_heavyatom_pair_one_way(
	Size const heavyatom1,
	conformation::Residue const & rsd1,
	LKB_ResidueInfo const & rsd1_info,
	Size const heavyatom2,
	conformation::Residue const & rsd2,
	LKB_ResidueInfo const & rsd2_info,
	EnergyMap const & weights,
	Real const weight_factor,
	Real const d2,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const
{
	Vectors const & atom1_waters = rsd1_info.waters()[heavyatom1];
	utility::vector1< Real > const & atom1_wts = rsd1_info.atom_weights()[heavyatom1];
	Vectors const & atom2_waters = rsd2_info.waters()[heavyatom2];

	if ( atom1_waters.empty() ) return;

	if ( atom1_wts.size() != 2 ) {
		std::cerr << "LK_BallEnergy: bad atom1_wts " << rsd1.name() << ' ' << rsd1.atom_name( heavyatom1 ) << ' ' <<
			atom1_wts.size() << std::endl;
		utility_exit_with_message("bad atom1_wts: "+rsd1.name()+" "+rsd1.atom_name(heavyatom1));
	}

	Real lk_deriv, lk_score, other_lk_score;  // other_lk score needed for bridging water

	if ( slim_etable_ ) {
		Real other_lk_deriv, inv_dis;
		etable_->analytic_lk_energy( rsd1.atom( heavyatom1 ), rsd2.atom( heavyatom2 ), lk_score, other_lk_score );
		etable_->analytic_lk_derivatives( rsd1.atom( heavyatom1 ), rsd2.atom( heavyatom2 ), lk_deriv, other_lk_deriv, inv_dis );
	} else {
		Real const d2_bin( d2 * etable_bins_per_A2_ );
		Size const disbin( static_cast< int >( d2_bin ) + 1 );
		Real const frac( d2_bin - ( disbin - 1 ) );

		/// this index into solv1 or dsolv1 should represent desolvation of heavyatom1 by atom2
		Size const heavyatom1_type_index( rsd1.atom( heavyatom1 ).type() );
		int const linear_index( dsolv1_.index( disbin, rsd2.atom( heavyatom2 ).type(), heavyatom1_type_index ) );

		lk_deriv = dsolv1_[ linear_index ] * ( 1-frac ) + dsolv1_[ linear_index+1 ] * frac;
		lk_score =  solv1_[ linear_index ] * ( 1-frac ) +  solv1_[ linear_index+1 ] * frac;
		other_lk_score =  solv2_[ linear_index ] * ( 1-frac ) +  solv2_[ linear_index+1 ] * frac;
	}

	// get weights
	Real const // include contributions from 3 score types:  lk_ball_iso, lk_ball, and lk_ball_wtd
		lk_ball_iso_weight  ( weights[ lk_ball_iso ] + atom1_wts[1] * weights[ lk_ball_wtd ] ), // HACKING
		lk_ball_aniso_weight( weights[ lk_ball     ] + atom1_wts[2] * weights[ lk_ball_wtd ] ),
		lk_ball_bridge_weight( weights[ lk_ball_bridge ] ),
		lk_ball_bridge_uncpl_weight( weights[ lk_ball_bridge_uncpl ] );

	Vector const & heavyatom1_xyz( rsd1.xyz( heavyatom1 ) );
	Vector const & atom2_xyz( rsd2.xyz( heavyatom2 ) );
	Vector f1( heavyatom1_xyz.cross( atom2_xyz ) ), f2( heavyatom1_xyz - atom2_xyz );
	Real const inv_dis( 1.0 / std::sqrt( d2 ) );

	bool skip_heavyatom = ( d2 > fasol_max_dis2_ );

	// now the lk_ball deriv
	//
	// lk_ball = fraction * lk
	// so d lk_ball = fraction * ( d lk ) + ( d fraction ) * lk
	// fpd with softmax:
	//              = fraction * (d lk / d dsum * d dsum/d di ) + ( d fraction / d dsum * d dsum/d di ) * lk
	if ( !skip_heavyatom ) {
		Real weighted_water_d2(0);
		utility::vector1< core::Real > dwwd2_ddi;
		Real const lk_fraction(
			get_lk_fractional_contribution( atom2_xyz, rsd2.atom( heavyatom2 ).type(), atom1_waters, dwwd2_ddi, weighted_water_d2 ) );

		// (1) the derivs for the parts that don't involve waters:
		{
			Real dE_dr_over_r = weight_factor * lk_ball_iso_weight * lk_deriv * inv_dis;
			dE_dr_over_r += weight_factor * lk_ball_aniso_weight * lk_fraction * lk_deriv * inv_dis;
			r1_at_derivs[heavyatom1].f1() += dE_dr_over_r * f1;
			r1_at_derivs[heavyatom1].f2() += dE_dr_over_r * f2;
			r2_at_derivs[heavyatom2].f1() -= dE_dr_over_r * f1;
			r2_at_derivs[heavyatom2].f2() -= dE_dr_over_r * f2;
		}

		// (2) waters-atom2 interaction
		// note that there's no derivative unless we're in the ramping zone:
		// we only do this if this term hasn't already been captured by one of our dependent hydrogens
		// we assume that for heavyatoms with dependent polar hydrogens, every water belongs to one of the hydrogens
		if ( weighted_water_d2 < ramp_width_A2_ && weighted_water_d2 > 0.0 ) {
			for ( Size i=1; i<=atom1_waters.size(); ++i ) {
				Vector const & atom1_water_xyz( atom1_waters[ i ] );

				// update f1 and f2 to reflect water-atom2 as the interaction
				Vector f1w = atom1_water_xyz.cross( atom2_xyz );
				Vector f2w = atom1_water_xyz - atom2_xyz;

				// what is the derivative of the lk_fraction term wrt r?
				Real const dE_dr_over_r
					( weight_factor * lk_ball_aniso_weight * lk_score * dwwd2_ddi[i] * eval_d_lk_fraction_dr_over_r( weighted_water_d2, ramp_width_A2_ ) );

				// Derivatives for the desolvating atom
				r2_at_derivs[heavyatom2].f1() -= dE_dr_over_r * f1w;
				r2_at_derivs[heavyatom2].f2() -= dE_dr_over_r * f2w;

				// derivatives for the desolvated atoms
				// dR/datom1 = dR/dwater * dwater/datom1
				//core::Real dwater = f2.length();
				//Real denom = dwater*dwater*dwater;
				{
					Size atom1 = rsd1_info.get_water_builder( rsd1 , heavyatom1 )[i].atom1();
					Vector const & r1_atom1_xyz( rsd1.xyz( atom1 ) );

					numeric::xyzMatrix< Real >const & dwater_datom1 = rsd1_info.atom1_derivs()[heavyatom1][i];
					Vector dwater_datom1x ( dwater_datom1(1,1), dwater_datom1(2,1), dwater_datom1(3,1) );
					Vector dwater_datom1y ( dwater_datom1(1,2), dwater_datom1(2,2), dwater_datom1(3,2) );
					Vector dwater_datom1z ( dwater_datom1(1,3), dwater_datom1(2,3), dwater_datom1(3,3) );
					Vector dRdatom( f2w.dot( dwater_datom1x ), f2w.dot( dwater_datom1y ), f2w.dot( dwater_datom1z ) );

					Vector f2t = dRdatom;
					Vector f1t = r1_atom1_xyz.cross( r1_atom1_xyz - dRdatom );

					r1_at_derivs[atom1].f1() += dE_dr_over_r * f1t;
					r1_at_derivs[atom1].f2() += dE_dr_over_r * f2t;
				}
				{
					Size atom2 = rsd1_info.get_water_builder( rsd1 , heavyatom1 )[i].atom2();
					Vector const & r1_atom2_xyz( rsd1.xyz( atom2 ) );

					numeric::xyzMatrix< Real >const & dwater_datom2 = rsd1_info.atom2_derivs()[heavyatom1][i];
					Vector dwater_datom2x ( dwater_datom2(1,1), dwater_datom2(2,1), dwater_datom2(3,1) );
					Vector dwater_datom2y ( dwater_datom2(1,2), dwater_datom2(2,2), dwater_datom2(3,2) );
					Vector dwater_datom2z ( dwater_datom2(1,3), dwater_datom2(2,3), dwater_datom2(3,3) );
					Vector dRdatom( f2w.dot( dwater_datom2x ), f2w.dot( dwater_datom2y ), f2w.dot( dwater_datom2z ) );

					Vector f2t = dRdatom;
					Vector f1t = r1_atom2_xyz.cross( r1_atom2_xyz - dRdatom );

					r1_at_derivs[atom2].f1() += dE_dr_over_r * f1t;
					r1_at_derivs[atom2].f2() += dE_dr_over_r * f2t;
				}
				{
					Size atom3 = rsd1_info.get_water_builder( rsd1 , heavyatom1 )[i].atom3();
					Vector const & r1_atom3_xyz( rsd1.xyz( atom3 ) );

					numeric::xyzMatrix< Real >const & dwater_datom3 = rsd1_info.atom3_derivs()[heavyatom1][i];
					Vector dwater_datom3x ( dwater_datom3(1,1), dwater_datom3(2,1), dwater_datom3(3,1) );
					Vector dwater_datom3y ( dwater_datom3(1,2), dwater_datom3(2,2), dwater_datom3(3,2) );
					Vector dwater_datom3z ( dwater_datom3(1,3), dwater_datom3(2,3), dwater_datom3(3,3) );
					Vector dRdatom( f2w.dot( dwater_datom3x ), f2w.dot( dwater_datom3y ), f2w.dot( dwater_datom3z ) );

					Vector f2t = dRdatom;
					Vector f1t = r1_atom3_xyz.cross( r1_atom3_xyz - dRdatom );

					r1_at_derivs[atom3].f1() += dE_dr_over_r * f1t;
					r1_at_derivs[atom3].f2() += dE_dr_over_r * f2t;
				}
			}
		}
	}

	// (3) waters-waters interaction
	// only if lk_ball_bridge is turned on
	// unlike the other terms, this is _symmetric_ -- instead of A desolv B and B desolv A, we just get A/B water overlap computed once
	//   thus, we only compute the derivatives for the A component here
	if ( atom2_waters.empty() ) return;
	if ( lk_ball_bridge_weight==0 && lk_ball_bridge_uncpl_weight==0 ) return;

	Real bridging_lk_factor = lk_ball_bridge_uncpl_weight + lk_ball_bridge_weight * (lk_score+other_lk_score);
	Real weighted_d2_water_delta(0), angleterm_lkbr(0), pointterm_lkbr(0), d_angleterm_lkbr_dr(0);
	utility::vector1< numeric::xyzVector<core::Real> > d_weighted_d2_d_di;

	Real lk_desolvation_sum = 0.0, lkbr_wt = 0.0, lkbr_uncpl_wt = 0.0;
	Real const lkbr_fraction(
		get_lkbr_fractional_contribution(
		heavyatom1_xyz, atom2_xyz, atom1_waters, atom2_waters, d_weighted_d2_d_di, weighted_d2_water_delta, pointterm_lkbr, angleterm_lkbr, d_angleterm_lkbr_dr, lk_desolvation_sum, lkbr_wt, lkbr_uncpl_wt ) );

	// A: change in LK (used as a scalefactor) on _bridge but not _bridge_uncpl
	if ( lk_ball_bridge_weight != 0 ) {
		Real const dE_dr_over_r( weight_factor * lk_ball_bridge_weight * lkbr_fraction * lk_deriv * inv_dis );
		r1_at_derivs[heavyatom1].f1() += dE_dr_over_r * f1;
		r1_at_derivs[heavyatom1].f2() += dE_dr_over_r * f2;
		r2_at_derivs[heavyatom2].f1() -= dE_dr_over_r * f1;
		r2_at_derivs[heavyatom2].f2() -= dE_dr_over_r * f2;
	}

	// A': water-angle potential (also used as a scalefactor)
	if ( lkbridge_angle_widthscale_!=0 ) {
		Real const dE_dr_over_r( -0.25*weight_factor * bridging_lk_factor * pointterm_lkbr * d_angleterm_lkbr_dr * inv_dis );
		r1_at_derivs[heavyatom1].f1() += dE_dr_over_r * f1;
		r1_at_derivs[heavyatom1].f2() += dE_dr_over_r * f2;
		r2_at_derivs[heavyatom2].f1() -= dE_dr_over_r * f1;
		r2_at_derivs[heavyatom2].f2() -= dE_dr_over_r * f2;
	}

	// B: change in water positions
	if ( weighted_d2_water_delta < overlap_width_A2_ && weighted_d2_water_delta > 0.0 ) {
		for ( Size i=1; i<=atom1_waters.size(); ++i ) {
			// what is the derivative of the lkbr_fraction term wrt r?
			numeric::xyzVector<core::Real> const dE_dwi
				( weight_factor * bridging_lk_factor * d_weighted_d2_d_di[i] * angleterm_lkbr * eval_d_lk_fraction_dr_over_r( weighted_d2_water_delta, overlap_width_A2_ ) );

			// derivatives for the desolvated atoms
			// dR/datom1 = dR/dwater * dwater/datom1
			{
				Size atom1 = rsd1_info.get_water_builder( rsd1 , heavyatom1 )[i].atom1();
				Vector const & r1_atom1_xyz( rsd1.xyz( atom1 ) );

				numeric::xyzMatrix< Real >const & dwater_datom1 = rsd1_info.atom1_derivs()[heavyatom1][i];
				Vector dwater_datom1x ( dwater_datom1(1,1), dwater_datom1(2,1), dwater_datom1(3,1) );
				Vector dwater_datom1y ( dwater_datom1(1,2), dwater_datom1(2,2), dwater_datom1(3,2) );
				Vector dwater_datom1z ( dwater_datom1(1,3), dwater_datom1(2,3), dwater_datom1(3,3) );
				Vector dRdatom( dE_dwi.dot( dwater_datom1x ), dE_dwi.dot( dwater_datom1y ), dE_dwi.dot( dwater_datom1z ) );

				Vector f2t = dRdatom;
				Vector f1t = r1_atom1_xyz.cross( r1_atom1_xyz - dRdatom );

				r1_at_derivs[atom1].f1() += f1t;
				r1_at_derivs[atom1].f2() += f2t;
			}
			{
				Size atom2 = rsd1_info.get_water_builder( rsd1 , heavyatom1 )[i].atom2();
				Vector const & r1_atom2_xyz( rsd1.xyz( atom2 ) );

				numeric::xyzMatrix< Real >const & dwater_datom2 = rsd1_info.atom2_derivs()[heavyatom1][i];
				Vector dwater_datom2x ( dwater_datom2(1,1), dwater_datom2(2,1), dwater_datom2(3,1) );
				Vector dwater_datom2y ( dwater_datom2(1,2), dwater_datom2(2,2), dwater_datom2(3,2) );
				Vector dwater_datom2z ( dwater_datom2(1,3), dwater_datom2(2,3), dwater_datom2(3,3) );
				Vector dRdatom( dE_dwi.dot( dwater_datom2x ), dE_dwi.dot( dwater_datom2y ), dE_dwi.dot( dwater_datom2z ) );

				Vector f2t = dRdatom;
				Vector f1t = r1_atom2_xyz.cross( r1_atom2_xyz - dRdatom );

				r1_at_derivs[atom2].f1() += f1t;
				r1_at_derivs[atom2].f2() += f2t;
			}
			{
				Size atom3 = rsd1_info.get_water_builder( rsd1 , heavyatom1 )[i].atom3();
				Vector const & r1_atom3_xyz( rsd1.xyz( atom3 ) );

				numeric::xyzMatrix< Real >const & dwater_datom3 = rsd1_info.atom3_derivs()[heavyatom1][i];
				Vector dwater_datom3x ( dwater_datom3(1,1), dwater_datom3(2,1), dwater_datom3(3,1) );
				Vector dwater_datom3y ( dwater_datom3(1,2), dwater_datom3(2,2), dwater_datom3(3,2) );
				Vector dwater_datom3z ( dwater_datom3(1,3), dwater_datom3(2,3), dwater_datom3(3,3) );
				Vector dRdatom( dE_dwi.dot( dwater_datom3x ), dE_dwi.dot( dwater_datom3y ), dE_dwi.dot( dwater_datom3z ) );

				Vector f2t = dRdatom;
				Vector f1t = r1_atom3_xyz.cross( r1_atom3_xyz - dRdatom );

				r1_at_derivs[atom3].f1() += f1t;
				r1_at_derivs[atom3].f2() += f2t;
			}
		}
	}
}

/// @note  Assumes that atom1 is the "moving" atom, ie the atom for which eval_atom_derivative was called
/// @note  Calculates the water positions for atom2 if d2 < safe_max_dis2
void
LK_BallEnergy::sum_deriv_contributions_for_heavyatom_pair(
	Real const d2,
	Size const heavyatom1,
	conformation::Residue const & rsd1,
	LKB_ResidueInfo const & rsd1_info,
	Size const heavyatom2,
	conformation::Residue const & rsd2,
	LKB_ResidueInfo const & rsd2_info,
	pose::Pose const &,
	EnergyMap const & weights,
	Real const cp_weight,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const
{
	sum_deriv_contributions_for_heavyatom_pair_one_way(
		heavyatom1, rsd1, rsd1_info, heavyatom2, rsd2, rsd2_info, weights,
		cp_weight, d2, r1_at_derivs, r2_at_derivs );

	sum_deriv_contributions_for_heavyatom_pair_one_way(
		heavyatom2, rsd2, rsd2_info, heavyatom1, rsd1, rsd1_info, weights,
		cp_weight, d2, r2_at_derivs, r1_at_derivs );
}

void
LK_BallEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /* context_graphs_required */ ) const
{}


void
LK_BallEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & , //weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	debug_assert( set1.resid() != set2.resid() );

	using namespace methods;
	using namespace trie;
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table1( energy_table );
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table2( energy_table );

	temp_table1 = 0; temp_table2 = 0;

	RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie(  methods::lkball_method ) ));
	RotamerTrieBaseCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie(  methods::lkball_method ) ));

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp = get_count_pair_function_trie( set1, set2, pose, sfxn );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	lkbtrie::LKBTrieEvaluator lkbeval(
		sfxn.get_weight(lk_ball), sfxn.get_weight(lk_ball_iso), sfxn.get_weight(lk_ball_wtd),
		sfxn.get_weight(lk_ball_bridge), sfxn.get_weight(lk_ball_bridge_uncpl),
		*this, etable_
	);
	trie1->trie_vs_trie( *trie2, *cp, lkbeval, temp_table1, temp_table2 );

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
	std::cout << "lkballE: Residues " << set1.resid() << " & " << set2.resid() << " rotamers: " << ii << " & " << jj;
	std::cout << " tvt/reg discrepancy: tvt= " <<  temp_table1( jj, ii ) << " reg= " << temp_table3( jj, ii );
	std::cout << " delta: " << temp_table1( jj, ii ) - temp_table3( jj, ii ) << std::endl;
	}
	}
	}
	std::cout << "Finished RPE calcs for residues " << set1.resid() << " & " << set2.resid() << std::endl;
	*/
}

void
LK_BallEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & , //weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	using namespace methods;
	using namespace trie;
	using namespace lkbtrie;

	// allocate space for the trie-vs-trie algorithm
	utility::vector1< core::PackerEnergy > temp_vector1( set.num_rotamers(), 0.0 );
	utility::vector1< core::PackerEnergy > temp_vector2( set.num_rotamers(), 0.0 );

	RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set.get_trie( methods::lkball_method ) ));
	RotamerTrieBaseCOP trie2 = ( static_cast< TrieCollection const & >
		( pose.energies().data().get( EnergiesCacheableDataType::LKB_TRIE_COLLECTION )) ).trie( rsd.seqpos() );

	if ( trie2 == NULL ) return;

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp = get_count_pair_function_trie( pose.residue( set.resid() ), rsd, trie1, trie2, pose, sfxn );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	lkbtrie::LKBTrieEvaluator lkbeval(
		sfxn.get_weight(lk_ball), sfxn.get_weight(lk_ball_iso), sfxn.get_weight(lk_ball_wtd),
		sfxn.get_weight(lk_ball_bridge), sfxn.get_weight(lk_ball_bridge_uncpl),
		*this, etable_
	);
	trie1->trie_vs_path( *trie2, *cp, lkbeval, temp_vector1, temp_vector2 );

	/// add in the energies calculated by the tvt alg.
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		energy_vector[ ii ] += temp_vector1[ ii ];
	}
}


void
LK_BallEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & res1data,
	ResSingleMinimizationData const & res2data,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
	debug_assert( r1_at_derivs.size() >= rsd1.natoms() );
	debug_assert( r2_at_derivs.size() >= rsd2.natoms() );

	// retrieve some info
	LKB_ResidueInfo const & rsd1_info( retrieve_lkb_resdata( res1data ) );
	LKB_ResidueInfo const & rsd2_info( retrieve_lkb_resdata( res2data ) );

	ResiduePairNeighborList const & nblist =
		static_cast< ResiduePairNeighborList const & > (min_data.get_data_ref( etab_pair_nblist ));

	bool use_lkbr_uncpl = (weights[core::scoring::lk_ball_bridge_uncpl]!=0);

	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		Size const heavyatom1( neighbs[ ii ].atomno1() ), heavyatom2( neighbs[ ii ].atomno2() );
		if ( rsd1.atom_is_hydrogen( heavyatom1 ) || rsd2.atom_is_hydrogen( heavyatom2 ) ) continue;
		Real const cp_weight( neighbs[ ii ].weight() );

		Real const d2( rsd1.xyz( heavyatom1 ).distance_squared( rsd2.xyz( heavyatom2 ) ) );

		if ( d2 == Real(0.0) ) continue; // sanity check
		if ( !use_lkbr_uncpl && d2 >= fasol_max_dis2_ ) continue;
		if ( d2 >= lkb_max_dis2_ ) continue;

		// fpd ... new version works at the heavyatom level
		sum_deriv_contributions_for_heavyatom_pair(
			d2, heavyatom1, rsd1, rsd1_info, heavyatom2, rsd2, rsd2_info, pose, weights, cp_weight, r1_at_derivs, r2_at_derivs );
	}
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
}


//
//fpd this function is only called during nblist auto update
void
LK_BallEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &sf,
	EnergyMap & totals
) const
{
	if ( !pose.energies().use_nblist() || !pose.energies().use_nblist_auto_update() ) return;

	EnergyMap tbenergy_map;
	// add in contributions from the nblist atom-pairs
	NeighborList const & nblist ( pose.energies().nblist( EnergiesCacheableDataType::ETABLE_NBLIST ) );
	nblist.check_domain_map( pose.energies().domain_map() );

	/// Trick to avoid calls to Conformation::residue()
	utility::vector1< conformation::Residue const * > resvect;
	resvect.reserve( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		resvect.push_back( & pose.residue( ii ) );
	}

	for ( Size i=1, i_end = pose.size(); i<= i_end; ++i ) {
		conformation::Residue const & rsd1( *resvect[i] );
		LKB_ResidueInfo const & rsd1_info( retrieve_lkb_residue_info( pose, rsd1.seqpos() ) );
		utility::vector1< Vectors > const & rsd1_waters( rsd1_info.waters() );
		utility::vector1< utility::vector1< Real > > const & rsd1_atom_wts( rsd1_info.atom_weights() );

		for ( Size ii=1, ii_end=rsd1.nheavyatoms(); ii<= ii_end; ++ii ) {
			Vectors const & atom1_waters( rsd1_waters[ ii ] );
			Vector const & atom1_xyz( rsd1.xyz( ii ) );
			Size const atom1_type_index( rsd1.atom( ii ).type() );
			utility::vector1< Real > const & atom1_weights( rsd1_atom_wts[ii] );

			AtomNeighbors const & nbrs( nblist.upper_atom_neighbors(i,ii) );
			for ( AtomNeighbors::const_iterator nbr=nbrs.begin(),
					nbr_end=nbrs.end(); nbr!= nbr_end; ++nbr ) {
				Size const  j( nbr->rsd() );
				Size const jj( nbr->atomno() );

				conformation::Residue const & rsd2( *resvect[j] );

				if ( rsd2.atom_is_hydrogen(jj) ) continue;
				if ( i == j ) continue;

				LKB_ResidueInfo const & rsd2_info( retrieve_lkb_residue_info( pose, rsd2.seqpos() ) );
				utility::vector1< Vectors > const & rsd2_waters( rsd2_info.waters() );
				utility::vector1< utility::vector1< Real > > const & rsd2_atom_wts( rsd2_info.atom_weights() );
				utility::vector1< Real > const & atom2_weights( rsd2_atom_wts[jj] );

				Vectors const & atom2_waters( rsd2_waters[ jj ] );
				Vector const & atom2_xyz( rsd2.xyz( jj ) );

				if ( atom1_waters.empty() && atom2_waters.empty() ) continue;

				Real const cp_weight( nbr->weight_func()*nbr->weight() );  // fpd is this correct (same CP as etable)

				Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

				Real lk_desolvation_of_atom1_by_atom2, lk_desolvation_of_atom2_by_atom1;
				Size const atom2_type_index( rsd2.atom( jj ).type() );
				if ( slim_etable_ ) {
					etable_->analytic_lk_energy( rsd1.atom( ii ), rsd2.atom( jj ), lk_desolvation_of_atom1_by_atom2,
						lk_desolvation_of_atom2_by_atom1 );
					lk_desolvation_of_atom1_by_atom2 *= cp_weight;
					lk_desolvation_of_atom2_by_atom1 *= cp_weight;

				} else {
					// setup for solvation Etable lookups
					Real const d2_bin = d2 * etable_bins_per_A2_;
					int disbin = static_cast< int >( d2_bin ) + 1;
					Real frac = d2_bin - ( disbin - 1 );
					int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

					lk_desolvation_of_atom1_by_atom2 = cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] );
					lk_desolvation_of_atom2_by_atom1 = cp_weight * ( ( 1. - frac ) * solv2_[ l1 ] + frac * solv2_[ l1+1 ] );
				}

				if ( atom1_waters.size() > 0 ) {
					accumulate_single_atom_contributions( ii, atom1_type_index, atom1_waters, atom1_weights,
						rsd1, atom2_type_index, atom2_xyz,
						lk_desolvation_of_atom1_by_atom2, tbenergy_map );
				}

				if ( atom2_waters.size() > 0 ) {
					accumulate_single_atom_contributions( jj, atom2_type_index, atom2_waters, atom2_weights,
						rsd2, atom1_type_index, atom1_xyz,
						lk_desolvation_of_atom2_by_atom1, tbenergy_map );
				}

				if ( sf.get_weight( core::scoring::lk_ball_bridge ) != 0 || sf.get_weight( core::scoring::lk_ball_bridge_uncpl) != 0 ) {
					Real lk_desolvation_sum = lk_desolvation_of_atom1_by_atom2+lk_desolvation_of_atom2_by_atom1;
					Real lkbr_wt = sf.get_weight( core::scoring::lk_ball_bridge );
					Real lkbr_uncpl_wt = sf.get_weight( core::scoring::lk_ball_bridge_uncpl );
					Real lkbridge_frac = get_lkbr_fractional_contribution(
						atom1_xyz, atom2_xyz,
						atom1_waters, atom2_waters,
						lk_desolvation_sum,
						lkbr_wt, lkbr_uncpl_wt );
					tbenergy_map[ lk_ball_bridge ] += (lk_desolvation_of_atom1_by_atom2+lk_desolvation_of_atom2_by_atom1) * lkbridge_frac;
					tbenergy_map[ lk_ball_bridge_uncpl ] += lkbridge_frac;
				}
			}
		}
	}
	totals += tbenergy_map;
}


//
//fpd  this function is only called during during minimization when nblist auto update is on
//fpd   there is some redundancy in this function as some terms are computed multiple times
void
LK_BallEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & , //sfxn, // needed for non-nblist minimization
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	Size const idresid = id.rsd();

	conformation::Residue const & rsd1( pose.residue(idresid) );
	Size const atom1( id.atomno() );

	NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::ETABLE_NBLIST ) );

	bool use_lkbr_uncpl = (weights[core::scoring::lk_ball_bridge_uncpl]!=0);

	LKB_ResidueInfo const & rsd1_info( retrieve_lkb_residue_info( pose, rsd1.seqpos() ) );
	if ( !rsd1_info.has_waters() ) return;
	scoring::AtomNeighbors const & nbrs( nblist.atom_neighbors(idresid,atom1) );

	if ( !rsd1.atom_is_hydrogen( atom1 ) ) {
		for ( scoring::AtomNeighbors::const_iterator it2=nbrs.begin(), it2e=nbrs.end(); it2 != it2e; ++it2 ) {
			// part 1: all heavyatoms desolvated by this atom
			scoring::AtomNeighbor const & nbr( *it2 );
			Size const nbrresid = nbr.rsd();
			Real const cp_weight( nbr.weight() );  // do not use nbr->weight_func() here
			if ( nbrresid == idresid ) continue; // no intra

			conformation::Residue const & rsd2( pose.residue(nbrresid) );
			Size const heavyatom2( nbr.atomno() );
			if ( rsd2.atom_is_hydrogen( heavyatom2 ) ) continue;

			LKB_ResidueInfo const & rsd2_info( retrieve_lkb_residue_info( pose, nbrresid ) );
			if ( rsd2_info.waters()[heavyatom2].size() == 0 ) continue;

			Real const d2( rsd1.xyz( atom1 ).distance_squared( rsd2.xyz( heavyatom2 ) ) );
			if ( d2 == Real(0.0) ) continue; // sanity check
			if ( !use_lkbr_uncpl && d2 >= fasol_max_dis2_ ) continue;
			if ( d2 >= lkb_max_dis2_ ) continue;

			utility::vector1< DerivVectorPair > r1_at_derivs(rsd1.natoms()), r2_at_derivs(rsd2.natoms());
			sum_deriv_contributions_for_heavyatom_pair_one_way(
				heavyatom2, rsd2, rsd2_info, atom1, rsd1, rsd1_info, weights, cp_weight, d2, r2_at_derivs, r1_at_derivs );

			F1 += r1_at_derivs[atom1].f1();
			F2 += r1_at_derivs[atom1].f2();
		}
	}

	// part 2: _all_ atoms desolvating this atom
	for ( Size i=1; i<=rsd1.nheavyatoms(); ++i ) {
		Size nwaters = rsd1_info.waters()[i].size();
		WaterBuilders const & builder = rsd1_info.get_water_builder( rsd1 , i );

		bool need_to_calculate = false;
		for ( Size j=1; j<=nwaters && !need_to_calculate; ++j ) {
			if ( builder[j].atom1() == atom1 || builder[j].atom2() == atom1 || builder[j].atom3() == atom1 ) {
				need_to_calculate = true;
			}
		}

		if ( !need_to_calculate ) continue;

		scoring::AtomNeighbors const & nbrs2( nblist.atom_neighbors(idresid,i) );

		for ( scoring::AtomNeighbors::const_iterator it2=nbrs2.begin(), it2e=nbrs2.end(); it2 != it2e; ++it2 ) {
			scoring::AtomNeighbor const & nbr( *it2 );
			Size const nbrresid = nbr.rsd();
			Real const cp_weight( nbr.weight() );  // do not use nbr->weight_func() here

			if ( nbrresid == idresid ) continue; // no intra

			conformation::Residue const & rsd2( pose.residue(nbrresid) );
			LKB_ResidueInfo const & rsd2_info( retrieve_lkb_residue_info( pose, nbrresid ) );

			Size const heavyatom2( nbr.atomno() );
			if ( rsd2.atom_is_hydrogen( heavyatom2 ) ) continue;
			Real const d2( rsd1.xyz( i ).distance_squared( rsd2.xyz( heavyatom2 ) ) );

			if ( ( d2 >= lkb_max_dis2_ ) || ( d2 == Real(0.0) ) ) continue; // no contribution
			if ( rsd1_info.waters()[i].size() == 0 ) continue;

			utility::vector1< DerivVectorPair > r1_at_derivs(rsd1.natoms()), r2_at_derivs(rsd2.natoms());
			sum_deriv_contributions_for_heavyatom_pair_one_way(
				i, rsd1, rsd1_info, heavyatom2, rsd2, rsd2_info, weights, cp_weight, d2, r1_at_derivs, r2_at_derivs );

			F1 += r1_at_derivs[atom1].f1();
			F2 += r1_at_derivs[atom1].f2();
		}
	}
}

// create a lkb trie rotamer descriptor from a single rotamer
template <class CPDAT>
void
create_rotamer_descriptor(
	conformation::Residue const & res,
	LKB_ResidueInfo const &lkb_resinfo,
	trie::CPDataCorrespondence const & cpdata_map,
	trie::RotamerDescriptor< lkbtrie::LKBAtom, CPDAT > & rotamer_descriptor
)
{
	using namespace trie;
	using namespace lkbtrie;

	rotamer_descriptor.natoms( res.nheavyatoms() );

	Size count_added_atoms = 0;
	for ( Size jj = 1; jj <= res.nheavyatoms(); ++jj ) {
		LKBAtom newatom;
		CPDAT cpdata;
		initialize_cpdata_for_atom( cpdata, jj, res, cpdata_map );

		newatom.atom( res.atom(jj) );
		newatom.waters( lkb_resinfo.waters()[jj] );
		newatom.atom_weights( lkb_resinfo.atom_weights()[jj] );

		RotamerDescriptorAtom< LKBAtom, CPDAT > rdatom( newatom, cpdata );
		rotamer_descriptor.atom( ++count_added_atoms, rdatom );
	}
}


lkbtrie::LKBRotamerTrieOP
LK_BallEnergy::create_rotamer_trie(
	conformation::RotamerSetBase const & rotset,
	pose::Pose const & /*pose*/
) const
{
	using namespace trie;
	using namespace lkbtrie;
	using namespace etable::etrie;

	LKB_RotamerSetInfo const & info( retrieve_lkb_rotamer_set_info( rotset ) );

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamerset( rotset ) );

	lkbtrie::LKBRotamerTrieOP retval;

	if ( cpdata_map.has_pseudobonds() ||
			cpdata_map.max_connpoints_for_residue() > 1 ||
			cpdata_map.n_entries() > 3 ) {
		utility::vector1< RotamerDescriptor< LKBAtom, CountPairDataGeneric > > rotamer_descriptors( rotset.num_rotamers() );
		for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
			create_rotamer_descriptor( *rotset.rotamer( ii ), info[ii], cpdata_map, rotamer_descriptors[ ii ] );
			rotamer_descriptors[ ii ].rotamer_id( ii );
		}
		sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );
		retval = lkbtrie::LKBRotamerTrieOP( new RotamerTrie< LKBAtom, CountPairDataGeneric >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		utility::vector1< RotamerDescriptor< LKBAtom, CountPairData_1_1 > > rotamer_descriptors( rotset.num_rotamers() );
		for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
			create_rotamer_descriptor( *rotset.rotamer( ii ), info[ii], cpdata_map, rotamer_descriptors[ ii ] );
			rotamer_descriptors[ ii ].rotamer_id( ii );
		}
		sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );
		retval = lkbtrie::LKBRotamerTrieOP( new RotamerTrie< LKBAtom, CountPairData_1_1 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 2 ) {
		utility::vector1< RotamerDescriptor< LKBAtom, CountPairData_1_2 > > rotamer_descriptors( rotset.num_rotamers() );
		for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
			create_rotamer_descriptor( *rotset.rotamer( ii ), info[ii], cpdata_map, rotamer_descriptors[ ii ] );
			rotamer_descriptors[ ii ].rotamer_id( ii );
		}
		sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );
		retval = lkbtrie::LKBRotamerTrieOP( new RotamerTrie< LKBAtom, CountPairData_1_2 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 3 ) {
		utility::vector1< RotamerDescriptor< LKBAtom, CountPairData_1_3 > > rotamer_descriptors( rotset.num_rotamers() );
		for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
			create_rotamer_descriptor( *rotset.rotamer( ii ), info[ii], cpdata_map, rotamer_descriptors[ ii ] );
			rotamer_descriptors[ ii ].rotamer_id( ii );
		}
		sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );
		retval = lkbtrie::LKBRotamerTrieOP( new RotamerTrie< LKBAtom, CountPairData_1_3 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else {
		utility_exit_with_message( "Unknown residue connection in LK_BallEnergy::create_rotamer_trie");
	}

	for ( Size ii = 1; ii <= cpdata_map.n_entries(); ++ii ) {
		retval->set_resid_2_connection_entry( cpdata_map.resid_for_entry( ii ), ii );
	}

	return retval;
}

lkbtrie::LKBRotamerTrieOP
LK_BallEnergy::create_rotamer_trie(
	conformation::Residue const & res,
	pose::Pose const & pose
) const
{
	using namespace trie;
	using namespace lkbtrie;
	using namespace etable::etrie;

	LKB_ResidueInfo const & info( retrieve_lkb_residue_info( pose, res.seqpos() ) );

	CPDataCorrespondence cpdata_map( create_cpdata_correspondence_for_rotamer( res ) );

	lkbtrie::LKBRotamerTrieOP retval;

	if ( cpdata_map.has_pseudobonds() ||
			cpdata_map.max_connpoints_for_residue() > 1 ||
			cpdata_map.n_entries() > 3 ) {
		utility::vector1< RotamerDescriptor< LKBAtom, CountPairDataGeneric > > rotamer_descriptors( 1 );
		create_rotamer_descriptor( res, info, cpdata_map, rotamer_descriptors[ 1 ] );
		rotamer_descriptors[ 1 ].rotamer_id( 1 );
		retval = lkbtrie::LKBRotamerTrieOP( new RotamerTrie< LKBAtom, CountPairDataGeneric >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 1 || cpdata_map.n_entries() == 0 /* HACK! */ ) {
		utility::vector1< RotamerDescriptor< LKBAtom, CountPairData_1_1 > > rotamer_descriptors( 1 );
		create_rotamer_descriptor( res, info, cpdata_map, rotamer_descriptors[ 1 ] );
		rotamer_descriptors[ 1 ].rotamer_id( 1 );
		retval = lkbtrie::LKBRotamerTrieOP( new RotamerTrie< LKBAtom, CountPairData_1_1 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 2 ) {
		utility::vector1< RotamerDescriptor< LKBAtom, CountPairData_1_2 > > rotamer_descriptors( 1 );
		create_rotamer_descriptor( res, info, cpdata_map, rotamer_descriptors[ 1 ] );
		rotamer_descriptors[ 1 ].rotamer_id( 1 );
		return lkbtrie::LKBRotamerTrieOP( new RotamerTrie< LKBAtom, CountPairData_1_2 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else if ( cpdata_map.n_entries() == 3 ) {
		utility::vector1< RotamerDescriptor< LKBAtom, CountPairData_1_3 > > rotamer_descriptors( 1 );
		create_rotamer_descriptor( res, info, cpdata_map, rotamer_descriptors[ 1 ] );
		rotamer_descriptors[ 1 ].rotamer_id( 1 );
		retval = lkbtrie::LKBRotamerTrieOP( new RotamerTrie< LKBAtom, CountPairData_1_3 >( rotamer_descriptors, atomic_interaction_cutoff()) );
	} else {
		utility_exit_with_message( "Unknown residue connection in LK_BallEnergy::create_rotamer_trie");
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
LK_BallEnergy::get_count_pair_function_trie(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) const
{
	conformation::Residue const & res1( pose.residue( set1.resid() ) );
	conformation::Residue const & res2( pose.residue( set2.resid() ) );

	trie::RotamerTrieBaseCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( methods::lkball_method ) ));
	trie::RotamerTrieBaseCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( methods::lkball_method ) ));

	return get_count_pair_function_trie( res1, res2, trie1, trie2, pose, sfxn );
}


trie::TrieCountPairBaseOP
LK_BallEnergy::get_count_pair_function_trie(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	trie::RotamerTrieBaseCOP trie1,
	trie::RotamerTrieBaseCOP trie2,
	pose::Pose const & ,
	ScoreFunction const &
) const
{
	using namespace etable::count_pair;
	using namespace trie;
	using namespace etable::etrie;

	TrieCountPairBaseOP tcpfxn;

	CPResidueConnectionType connection = CountPairFactory::determine_residue_connection( res1, res2 );
	Size conn1 = trie1->get_count_pair_data_for_residue( res2.seqpos() );
	Size conn2 = trie2->get_count_pair_data_for_residue( res1.seqpos() );

	if ( connection == CP_ONE_BOND ) {
		CPCrossoverBehavior crossover = determine_crossover_behavior( res1, res2, use_intra_dna_cp_crossover_4_ );
		switch ( crossover ) {
		case CP_CROSSOVER_3 :
			tcpfxn = TrieCountPairBaseOP( new TrieCountPair1BC3( conn1, conn2 ) );
			break;
		case CP_CROSSOVER_4 :
			tcpfxn = TrieCountPairBaseOP( new TrieCountPair1BC4( conn1, conn2 ) );
			break;
		default :
			utility_exit();
			break;
		}
	} else if ( connection == CP_MULTIPLE_BONDS_OR_PSEUDOBONDS ) {
		CPCrossoverBehavior crossover = determine_crossover_behavior( res1, res2, use_intra_dna_cp_crossover_4_ );

		TrieCountPairGenericOP cpgen( new TrieCountPairGeneric( res1, res2, conn1, conn2 ) );
		if ( crossover == CP_CROSSOVER_3 ) {
			cpgen->crossover( 3 );
		} else if ( crossover == CP_CROSSOVER_4 ) {
			cpgen->crossover( 4 );
		} else {
			utility_exit();
		}
		cpgen->hard_crossover( false );
		tcpfxn = cpgen;
	} else {
		tcpfxn = TrieCountPairBaseOP( new TrieCountPairAll );
	}
	return tcpfxn;
}

core::Size
LK_BallEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}

