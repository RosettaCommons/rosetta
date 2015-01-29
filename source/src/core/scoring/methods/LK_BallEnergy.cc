// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/LK_BallEnergy.hh
/// @brief  LK Solvation using hemisphere culling class declaration
/// @author David Baker
/// @author Andrew Leaver-Fay
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/GenBornPotential.fwd.hh>


// Unit headers
#include <core/scoring/methods/LK_BallEnergy.hh>
#include <core/scoring/methods/LK_BallEnergyCreator.hh>
#include <core/scoring/methods/LK_BallInfo.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>
//#include <core/pack/rotamer_set/RotamerSet.hh>
//#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
// #include <core/io/pdb/pose_io.hh> // HACK
// #include <fstream> // HACK

#include <core/scoring/constraints/AngleConstraint.hh>

#include <basic/options/util.hh> // HACK
#include <basic/options/keys/dna.OptionKeys.gen.hh>

#include <basic/prof.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

#include <utility/vector1.functions.hh> // HACK

#ifdef WIN32
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#endif

/// LAZY using
//using core::pack::rotamer_set::RotamerSet;

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the LK_hack class,
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
	return sts;
}



static thread_local basic::Tracer TR( "core.scoring.methods.LK_BallEnergy" );

Real const LK_BallEnergy::ramp_width_A2_( 5.0 );


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



class LKB_ResPairMinData : public basic::datacache::CacheableData
{
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
// 	void set_res1_data( LKB_ResidueInfoCOP );
// 	void set_res2_data( LKB_ResidueInfoCOP );

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
	debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResPairMinData > ( pairdata.get_data( lkb_respair_data ) ));
		lkb_pairdata = utility::pointer::static_pointer_cast< LKB_ResPairMinData > ( pairdata.get_data( lkb_respair_data ) );
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
debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResPairMinData const > ( pairdata.get_data( lkb_respair_data ) ) );
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
	debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResidueInfo > ( resdata.get_data( lkb_res_data ) ) );
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
debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResidueInfo const > ( resdata.get_data( lkb_res_data ) ) );
	return ( static_cast< LKB_ResidueInfo const & > ( resdata.get_data_ref( lkb_res_data ) ) );
}
/////////////////////////////////////// mindata retrieval functions
inline
LKB_ResidueInfoCOP
retrieve_lkb_resdata_ptr(
												ResSingleMinimizationData const & resdata
												)
{
debug_assert( utility::pointer::dynamic_pointer_cast< LKB_ResidueInfo const > ( resdata.get_data( lkb_res_data ) ) );
	return ( utility::pointer::static_pointer_cast< core::scoring::methods::LKB_ResidueInfo const > ( resdata.get_data( lkb_res_data ) ) );
}





/// HACKING //////////////////////////
// void
// LK_BallEnergy::setup_hack()
// {
// 	using namespace options;
// 	using namespace OptionKeys::dna::specificity::lk_ball_hack;
// 	if ( option[ lk_ball_positions ].user() ) {
// 		positions_ = option[ lk_ball_positions ]();
// 	}
// 	include_all_dna_ = false; // wait on this until we can do some dna relaxes and estimate reference energies
// }

// bool
// LK_BallEnergy::include_residue( conformation::Residue const & rsd ) const {
// 	return ( ( include_all_dna_ && rsd.is_DNA() ) ||
// 					 ( utility::has_element( positions_, rsd.seqpos() ) ) );
// }

// void
// LK_BallEnergy::add_my_score_types()
// {
// 	add_score_type( lk_ball_iso );
// 	add_score_type( lk_ball );
// 	add_score_type( lk_polar );
// 	add_score_type( lk_polar_nw );
// 	add_score_type( lk_nonpolar );
// 	add_score_type( lk_charged );

// 	add_score_type( lk_ball_xd );
// 	add_score_type( lk_polar_xd );
// 	add_score_type( lk_polar_nw_xd );
// 	add_score_type( lk_nonpolar_xd );


// 	add_score_type( lk_ball_dd );
// 	add_score_type( lk_polar_dd );
// 	add_score_type( lk_polar_nw_dd );
// 	add_score_type( lk_nonpolar_dd );
// }

LK_BallEnergy::LK_BallEnergy( EnergyMethodOptions const & options ):
	parent             ( methods::EnergyMethodCreatorOP( new LK_BallEnergyCreator ) ),
	etable_            (  ScoringManager::get_instance()->etable( options.etable_type() ) ),
	// TODO: locking once of casting to etable::EtableCOP would be better
	solv1_             (  ScoringManager::get_instance()->etable( options.etable_type() ).lock()->solv1()),
	solv2_             (  ScoringManager::get_instance()->etable( options.etable_type() ).lock()->solv2()),
	dsolv1_            (  ScoringManager::get_instance()->etable( options.etable_type() ).lock()->dsolv1()),
	safe_max_dis2_     (  ScoringManager::get_instance()->etable( options.etable_type() ).lock()->get_safe_max_dis2() ),
	etable_bins_per_A2_(  ScoringManager::get_instance()->etable( options.etable_type() ).lock()->get_bins_per_A2() ),
	use_intra_dna_cp_crossover_4_( true )
{
	setup_d2_bounds();
}

// LK_BallEnergy::LK_BallEnergy( etable::Etable const & etable_in):
// 	etable_(etable_in),
// 	solv1_(etable_in.solv1()),
// 	solv2_(etable_in.solv2()),
// 	dsolv1_( etable_in.dsolv1() ),
// 	safe_max_dis2_( etable_in.get_safe_max_dis2() ),
// 	etable_bins_per_A2_( etable_in.get_bins_per_A2())
// // 	d2_low_ ( 2.0 * 2.0 ),
// // 	d2_high_( 3.1 * 3.1 )
// {
// 	//setup_hack();
// 	add_my_score_types();

// 	setup_d2_bounds();
// // 	add_score_type( lk_ball_iso );
// // 	add_score_type( lk_ball_iso_nw );
// }


Distance
LK_BallEnergy::atomic_interaction_cutoff() const
{
	etable::EtableCOP etable_op( etable_ );
	return etable_op->max_dis();
}


/// clone
EnergyMethodOP
LK_BallEnergy::clone() const
{
	return EnergyMethodOP( new LK_BallEnergy( *this ) );
}


LK_BallEnergy::LK_BallEnergy( LK_BallEnergy const & src ):
	ContextIndependentTwoBodyEnergy( src ),
	etable_(src.etable_),
	solv1_( src.solv1_ ),
	solv2_( src.solv2_ ),
	dsolv1_( src.dsolv1_ ),
	safe_max_dis2_( src.safe_max_dis2_ ),
	etable_bins_per_A2_( src.etable_bins_per_A2_ ),
	use_intra_dna_cp_crossover_4_( src.use_intra_dna_cp_crossover_4_ )
{
	setup_d2_bounds();
}


///
void
compute_and_store_pose_waters(
															pose::Pose & pose
															)
{
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
//	using namespace core::pack::rotamer_set; // WaterPackingInfo
	LKB_PoseInfoOP info( new LKB_PoseInfo() );
	for ( Size i=1; i<= pose.total_residue(); ++i ) info->append( LKB_ResidueInfoOP( new LKB_ResidueInfo( pose, pose.residue(i) ) ) );
	pose.data().set( pose::datacache::CacheableDataType::LK_BALL_POSE_INFO, info );
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// New minimization interface code ////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
LK_BallEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const
{
	return false;
}



void
LK_BallEnergy::setup_for_minimizing_for_residue(
																								conformation::Residue const & rsd,
																								pose::Pose const &, // pose,
																								ScoreFunction const &, // scorefxn,
																								kinematics::MinimizerMapBase const &, // min_map,
																								ResSingleMinimizationData & resdata
																								) const
{
	LKB_ResidueInfo & info( retrieve_nonconst_lkb_resdata( resdata ) );
	info.initialize( rsd.type() );
	info.build_waters( rsd );
}

void
LK_BallEnergy::setup_for_minimizing_for_residue_pair(
																										 conformation::Residue const &, // rsd1,
																										 conformation::Residue const &, // rsd2,
																										 pose::Pose const &,
																										 ScoreFunction const &, //scorefxn,
																										 kinematics::MinimizerMapBase const &, // min_map,
																										 ResSingleMinimizationData const & res1data,
																										 ResSingleMinimizationData const & res2data,
																										 ResPairMinimizationData & pairdata
																										 ) const
{

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
																						 pose::Pose const &,// pose,
																						 ScoreFunction const &, // sfxn,
																						 ResSingleMinimizationData & resdata
																						 ) const
{
	LKB_ResidueInfo & info( retrieve_nonconst_lkb_resdata( resdata ) );
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


///
void
LK_BallEnergy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
	pose.update_residue_neighbors();
	compute_and_store_pose_waters( pose ); // could check task and do only some
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
}

void
LK_BallEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	/// update waters for residue that has changed during packing, eg in rotamer trials
	/// need to double-check the current logic on this...
	retrieve_nonconst_lkb_residue_info( pose, resid ).build_waters( pose.residue( resid ) );
}

///
void
LK_BallEnergy::setup_for_scoring(
																 pose::Pose & pose,
																 ScoreFunction const & //scfxn
) const
{
	// silly hack:
	// if ( std::abs( scfxn.get_weight( lk_ball_wtd ) ) > 1e-3 ) {
	// 	runtime_assert( basic::options::option[ basic::options::OptionKeys::dna::specificity::no_SP3_acceptor_waters ] );
	// }

	pose.update_residue_neighbors();
	compute_and_store_pose_waters( pose );
}

void
LK_BallEnergy::prepare_rotamers_for_packing(
																						pose::Pose const & pose,
																						conformation::RotamerSetBase & rotamer_set
																						) const
{
//	using namespace pack::rotamer_set;
	//TR.Trace << "prepare_rotamers_for_packing: " << rotamer_set.num_rotamers() << std::endl;

	/// create a rotamer set info object
	LKB_RotamerSetInfoOP info( new LKB_RotamerSetInfo );

	for ( Size n=1; n<= rotamer_set.num_rotamers(); ++n ) {
		info->append( LKB_ResidueInfoOP( new LKB_ResidueInfo( pose, *( rotamer_set.rotamer( n ) ) ) ) );
	}

	rotamer_set.data().set( conformation::RotamerSetCacheableDataType::LK_BALL_ROTAMER_SET_INFO, info );

}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
LK_BallEnergy::setup_d2_bounds()
{
	Real const h2o_radius( 1.4 );
	etable::EtableCOP etable_op( etable_ );
	chemical::AtomTypeSetCOP atom_set_op( etable_op->atom_set() );
	chemical::AtomTypeSet const & atom_set( *atom_set_op );
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

// 	if ( options::option[ options::OptionKeys::dna::specificity::lk_ball_prefactor ].user() ) {
// 		utility::vector1< std::string > const lkbp_string
// 			( options::option[ options::OptionKeys::dna::specificity::lk_ball_prefactor ] );
// 		if ( lkbp_string.size()%2 != 0 ) utility_exit_with_message("lk_ball_prefactor should have even length");
// 		for ( Size ii=0; ii< lkbp_string.size()/2; ++ii ) {
// 			std::string const atom_type_name( lkbp_string[ 2*ii+1 ] );
// 			Real const prefactor( float_of( lkbp_string[ 2*ii+2 ] ) );
// 			lk_ball_prefactor_[ atom_set.atom_type_index( atom_type_name ) ] = prefactor;
// 			TR.Trace << "lk_ball_prefactor: " << atom_type_name << ' ' << prefactor << std::endl;
// 		}
// 	}
}


/// @details  Stolen from LK_SigmoidalFunc in lk_hack
/// d2_delta = d2 - d2_low
///
Real
LK_BallEnergy::eval_lk_fraction( Real const d2_delta ) const
{
debug_assert( d2_delta >= -0.001 && d2_delta <= ramp_width_A2_ + 0.001 );
	static Real const inv_range( 1.0 / ramp_width_A2_ );
	Real const xprime( inv_range * d2_delta );
	return ( 1 - xprime*xprime ) * ( 1 - xprime*xprime );
}

Real
LK_BallEnergy::eval_d_lk_fraction_dr_over_r( Real const d2_delta ) const
{
debug_assert( d2_delta >= -0.001 && d2_delta <= ramp_width_A2_ + 0.001 );
	static Real const inv_range( 1.0 / ramp_width_A2_ );
	Real const xprime( inv_range * d2_delta );
	return -8.0 * inv_range * ( 1 - xprime * xprime ) * xprime;

	//Real const xprime_squared( numeric::square( inv_range * ( d2 - d2_low_ ) ) );
	//return -8.0 * ( 1 - xprime_squared ) * xprime_squared;

}

/// @note  closest_water may be set to 0 upon return if none of the waters are within ramp_width_A2_
Real
LK_BallEnergy::get_lk_fractional_contribution(
																							Vector const & atom2_xyz,
																							Size const atom2_type,
																							Vectors const & atom1_waters,
																							Size & closest_water,
																							Real & closest_water_d2_delta
																							) const
{
	//Real max_frac( 0.0 ), frac, d2;
	Real const d2_low( d2_low_[ atom2_type ] );

	// find the closest water:
	closest_water = 0;
	closest_water_d2_delta = 100.0;
	Real d2_delta;
	for ( Vectors::const_iterator water= atom1_waters.begin(), water_end= atom1_waters.end();
				water != water_end; ++water ) {
		d2_delta = atom2_xyz.distance_squared( *water ) - d2_low;
		if ( d2_delta > ramp_width_A2_ ) continue; // frac is 0.0
		else if ( d2_delta < closest_water_d2_delta ) {
			closest_water_d2_delta = d2_delta;
			closest_water = water - atom1_waters.begin() + 1; // hack...
			if ( d2_delta < 0.0 ) break; // no point -- frac is already 1.0 at this point
		}
	}

	/// now compute the fraction
	Real frac( 0.0 );
	if ( closest_water ) {
		frac = ( closest_water_d2_delta < 0.0 ? Real( 1.0 ) : eval_lk_fraction( closest_water_d2_delta ) );
	}
	return frac;
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
	else return eval_lk_fraction( d2_delta );
}


Real
LK_BallEnergy::get_lk_fractional_contribution(
																							Vector const & atom2_xyz,
																							Size const atom2_type,
																							Vectors const & atom1_waters
																							) const
{
	Size closest_water(0);
	Real closest_water_d2_delta(0.0);
	return get_lk_fractional_contribution( atom2_xyz, atom2_type, atom1_waters, closest_water, closest_water_d2_delta );
}



/// This guy is used during scoring if we are not minimizing
void
LK_BallEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	LKB_ResidueInfo const & info1( retrieve_lkb_residue_info( pose, rsd1.seqpos() ) );
	LKB_ResidueInfo const & info2( retrieve_lkb_residue_info( pose, rsd2.seqpos() ) );
	residue_pair_energy( rsd1, info1, rsd2, info2, emap );
	//residue_pair_energy( rsd1, info1.waters(), rsd2, info2.waters(), emap );
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
																				pose::Pose const &,// pose,
																				ScoreFunction const &,
																				EnergyMap & emap
																				) const
{
	LKB_ResPairMinData const & lkb_pairdata( retrieve_lkb_pairdata( pairdata ) );

	residue_pair_energy( rsd1, lkb_pairdata.res1_data(), rsd2, lkb_pairdata.res2_data(), emap );
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
		if ( cpfxn->count( atom1, atom2, cp_weight, pathdist ) ) {
			Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			// setup for solvation Etable lookups
			Size const atom2_type_index( rsd2.atom( atom2 ).type() );
			Real const d2_bin = d2 * etable_bins_per_A2_;
			int	disbin = static_cast< int >( d2_bin ) + 1;
			Real	frac = d2_bin - ( disbin - 1 );
			int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

			Real const lk_desolvation_of_atom1_by_atom2
				( cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] ) );

			total_desolvation += lk_desolvation_of_atom1_by_atom2;
		}
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

	if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) return; // TOO FARAWAY (OR SAME ATOM?)

	// setup for solvation Etable lookups
	Size const atom1_type_index( rsd1.atom( atom1 ).type() );
	Size const atom2_type_index( rsd2.atom( atom2 ).type() );

	Real const d2_bin = d2 * etable_bins_per_A2_;
	int	disbin = static_cast< int >( d2_bin ) + 1;
	Real	frac = d2_bin - ( disbin - 1 );
	int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

	lk_desolvation_of_atom1_by_atom2 = ( cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] ) );

	Real dummy_real;
	Size dummy_size;
	lk_ball_desolvation_of_atom1_by_atom2 = lk_desolvation_of_atom1_by_atom2 *
		get_lk_fractional_contribution( atom2_xyz, atom2_type_index, atom1_waters, dummy_size, dummy_real );
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

	if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) return; // TOO FARAWAY (OR SAME ATOM?)

	// setup for solvation Etable lookups
	Size const atom1_type_index( rsd1.atom( atom1 ).type() );
	Size const atom2_type_index( rsd2.atom( atom2 ).type() );

	Real const d2_bin = d2 * etable_bins_per_A2_;
	int	disbin = static_cast< int >( d2_bin ) + 1;
	Real	frac = d2_bin - ( disbin - 1 );
	int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

	lk_desolvation_of_atom1_by_atom2 = ( cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] ) );

	Real dummy_real;
	Size dummy_size;
	lk_ball_desolvation_of_atom1_by_atom2 = lk_desolvation_of_atom1_by_atom2 *
		get_lk_fractional_contribution( atom2_xyz, atom2_type_index, atom1_waters, dummy_size, dummy_real );
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

		if ( ( d2 >= safe_max_dis2_) || ( d2 < 1e-3 ) ) continue; // exclude self...

		// setup for solvation Etable lookups
		Size const atom2_type_index( rsd2.atom( atom2 ).type() );
		Real const d2_bin = d2 * etable_bins_per_A2_;
		int	disbin = static_cast< int >( d2_bin ) + 1;
		Real	frac = d2_bin - ( disbin - 1 );
		int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

		Real const lk_desolvation_of_atom1_by_atom2
			( ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] ) );

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
	if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) {
		atom1_lk_desolvation_by_atom2_deriv = 0.0;
		atom2_lk_desolvation_by_atom1_deriv = 0.0;
		return;
	}

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

/// HACKY HELPER FXN TAKEN DIRECTLY FROM BaseEtableEnergy.tmpl.hh
/// IF YOU CHANGE THE BEHAVIOR THERE YOU SHOULD CHANGE HERE AS WELL
///
scoring::etable::count_pair::CPCrossoverBehavior
determine_crossover_behavior(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	bool const use_intra_dna_cp_crossover_4
	//pose::Pose const &
	//ScoreFunction const & sfxn
)
{
	using namespace scoring::etable::count_pair;
	// maybe should ask "are these two residues polymers and do they share a polymeric bond"
//debug_assert( !sfxn.has_zero_weight( mm_twist ) );
	if ( res1.polymeric_sequence_distance(res2) == 1 ) {
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
	}else {
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
	if ( !atom1_waters.empty() ) {
		emap[ lk_ball_iso ] += lk_desolvation_of_atom1_by_atom2;
		Size dummy_size;
		Real dummy_real;
		Real const lk_desolvation_of_atom1_by_atom2_lkb
			( lk_desolvation_of_atom1_by_atom2 *
				get_lk_fractional_contribution( atom2_xyz, atom2_type_index, atom1_waters, dummy_size, dummy_real ) );
		emap[ lk_ball ] += lk_desolvation_of_atom1_by_atom2_lkb;
		emap[ lk_ball_wtd ] += ( atom1_wts[1] * lk_desolvation_of_atom1_by_atom2 +
														 atom1_wts[2] * lk_desolvation_of_atom1_by_atom2_lkb );
	}
}



void
LK_BallEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	LKB_ResidueInfo const & rsd1_info,
	conformation::Residue const & rsd2,
	LKB_ResidueInfo const & rsd2_info,
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

	// the old way:
	//
	//CountPairFunctionOP cpfxn =
	//	CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	CPCrossoverBehavior crossover = determine_crossover_behavior( rsd1, rsd2, use_intra_dna_cp_crossover_4_ ); //, pose ); //, sfxn );
	CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, crossover );

	// setup residue information
	for ( Size atom1=1; atom1<= rsd1.nheavyatoms(); ++atom1 ) {
		Vectors const & atom1_waters( rsd1_waters[ atom1 ] );
		Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
		Size const atom1_type_index( rsd1.atom( atom1 ).type() );
		//chemical::AtomType const & atom1_type( rsd1.atom_type( atom1 ) );
		utility::vector1< Real > const & atom1_weights( rsd1_atom_wts[atom1] );

		for ( Size atom2=1; atom2<= rsd2.nheavyatoms(); ++atom2 ) {
			Vectors const & atom2_waters( rsd2_waters[ atom2 ] );
			Vector const & atom2_xyz( rsd2.xyz( atom2 ) );
			//chemical::AtomType const & atom2_type( rsd2.atom_type( atom2 ) );

			if ( atom1_waters.empty() && atom2_waters.empty() ) continue;

			utility::vector1< Real > const & atom2_weights( rsd2_atom_wts[atom2] );

			Real cp_weight = 1.0; Size pathdist;
			if ( cpfxn->count( atom1, atom2, cp_weight, pathdist ) ) {
				Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

				if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

				// setup for solvation Etable lookups
				Size const atom2_type_index( rsd2.atom( atom2 ).type() );
				Real const d2_bin = d2 * etable_bins_per_A2_;
				int	disbin = static_cast< int >( d2_bin ) + 1;
				Real	frac = d2_bin - ( disbin - 1 );
				int const l1 = solv1_.index( disbin, atom2_type_index, atom1_type_index );

				Real lk_desolvation_of_atom1_by_atom2
					( cp_weight * ( ( 1. - frac ) * solv1_[ l1 ] + frac * solv1_[ l1+1 ] ) );

				Real lk_desolvation_of_atom2_by_atom1
					( cp_weight * ( ( 1. - frac ) * solv2_[ l1 ] + frac * solv2_[ l1+1 ] ) );

				accumulate_single_atom_contributions( atom1, atom1_type_index, atom1_waters, atom1_weights,
																							rsd1, atom2_type_index, atom2_xyz,
																							lk_desolvation_of_atom1_by_atom2, emap );

				accumulate_single_atom_contributions( atom2, atom2_type_index, atom2_waters, atom2_weights,
																							rsd2, atom1_type_index, atom1_xyz,
																							lk_desolvation_of_atom2_by_atom1, emap );

			} // count pair
		} // atom2
	} // atom1


	//PROF_STOP( basic::LK_BALL_RESIDUE_PAIR_ENERGY );
}


///I AM HAVING A SLEEPOVER





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


/**

	 Note that we calculate the lk_ball_iso derivative as well as the lk_ball derivative...

	 Derivatives are only included for heavyatoms.

	 For a non-polar heavyatom, the derivs are for all polar atoms that it desolvates

	 For a polar heavyatom, derivs are for all polar atoms that it desolvates as well as all atoms it's being
	 desolvated by.

	 Given an atom desolvating a polar atom:
	 * the lk_ball_iso deriv is the standard lk deriv, but make sure we use the correct array! (see LK_hack code)

	 * the lk_ball score = wt * lk_polar, so the derivs have two components. One looks like wt * lk_polar deriv
	   contribution. The other looks like the lk_polar term * the derivative of the wt. The derivative of the wt
		 is found by getting the closest water, taking the derivative of the wt term wrt distance and using f1/f2
		 contributions for the desolvating atom xyz and the water xyz.

 **/

/// @details  Compute the f1 and f2 vectors for the derivative of the lk_fraction term for atom2 desolvating atom1
///
/// @note f1 and f2 are zeroed and filled within this function, rather than being accumulated into
/// @note We pretend that atom1 is the "moving" atom
/// Should eliminate code duplication with next routine... currently this is for outside use.

Real
LK_BallEnergy::eval_lk_ball_fraction_deriv(
																					 Vector const & atom2_xyz,
																					 Size const atom2_type,
																					 Vectors const & atom1_waters,
																					 bool const evaluate_deriv,
																					 Vector & f1,
																					 Vector & f2
																					 ) const
{
	///
	Size closest_water(0);
	Real closest_water_d2_delta( 100.0);
	Real const lk_fraction( get_lk_fractional_contribution( atom2_xyz, atom2_type, atom1_waters,
																													closest_water, closest_water_d2_delta ) );

	if ( evaluate_deriv ) {
		f1.clear(); f2.clear();

		// water1-atom2 interaction
		// note that there's no derivative unless we're in the ramping zone:
		if ( closest_water_d2_delta < ramp_width_A2_ && closest_water_d2_delta > 0.0 ) {
			Vector const & atom1_water_xyz( atom1_waters[ closest_water ] );
			// what is the derivative of the lk_fraction term wrt r?
			Real const dE_dr_over_r( eval_d_lk_fraction_dr_over_r( closest_water_d2_delta ) );
			f1 = dE_dr_over_r * ( atom1_water_xyz.cross( atom2_xyz ) );
			f2 = dE_dr_over_r * ( atom1_water_xyz - atom2_xyz );
		}
	}
	return lk_fraction;
}




/// @note  atom2 is desolvating atom1. atom1_waters is non-empty
/// @note  Pretend that atom1 is the atom whose derivs are being calculated. weight_factor may include -1 term
/// to switch the order...
///
void
LK_BallEnergy::sum_deriv_contributions_for_atom_pair_one_way(
																														 Size const atom1,
																														 conformation::Residue const & rsd1,
																														 Vectors const & atom1_waters,
																														 utility::vector1< Real > const & atom1_wts,
																														 Size const atom2,
																														 conformation::Residue const & rsd2,
																														 EnergyMap const & weights,
																														 Real const weight_factor,
																														 Real const d2,
																														 Vector & F1,
																														 Vector & F2
																														 ) const
{
	if ( atom1_waters.empty() ) return;
	if ( atom1_wts.size() != 2 ) {
		std::cerr << "LK_BallEnergy: bad atom1_wts " << rsd1.name() << ' ' << rsd1.atom_name( atom1 ) << ' ' <<
			atom1_wts.size() << std::endl;
		utility_exit_with_message("bad atom1_wts: "+rsd1.name()+" "+rsd1.atom_name(atom1));
	}


	Real const d2_bin( d2 * etable_bins_per_A2_ );
	Size const disbin( static_cast< int >( d2_bin ) + 1 );
	Real const frac( d2_bin - ( disbin - 1 ) );

	/// this index into solv1 or dsolv1 should represent desolvation of atom1 by atom2
	Size const atom1_type_index( rsd1.atom( atom1 ).type() );
	int const linear_index( dsolv1_.index( disbin, rsd2.atom( atom2 ).type(), atom1_type_index ) );

	Real const lk_deriv( dsolv1_[ linear_index ] * ( 1-frac ) + dsolv1_[ linear_index+1 ] * frac );
	Real const lk_score(  solv1_[ linear_index ] * ( 1-frac ) +  solv1_[ linear_index+1 ] * frac );

	//// first for the direct lk interaction between atom1 and atom2 /////////////////////////////////////////////
	Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
	Vector const & atom2_xyz( rsd2.xyz( atom2 ) );
	Vector f1( atom1_xyz.cross( atom2_xyz ) ), f2( atom1_xyz - atom2_xyz );

	Real const inv_dis( 1.0 / std::sqrt( d2 ) );

	Real const // include contributions from 3 score types:  lk_ball_iso, lk_ball, and lk_ball_wtd
		lk_ball_iso_weight  ( weights[ lk_ball_iso ] + atom1_wts[1] * weights[ lk_ball_wtd ] ), // HACKING
		lk_ball_aniso_weight( weights[ lk_ball     ] + atom1_wts[2] * weights[ lk_ball_wtd ] );

	{ // the derivs for the parts that don't involve waters:
		Real const dE_dr_over_r( weight_factor * lk_ball_iso_weight * lk_deriv * inv_dis );
		F1 += dE_dr_over_r * f1;
		F2 += dE_dr_over_r * f2;
	}


	// now the lk_ball deriv
	//
	// lk_ball = fraction * lk
	// so d lk_ball = fraction * ( d lk ) + ( d fraction ) * lk
	//
	Size closest_water(0);
	Real closest_water_d2_delta( 100.0);
	Real const lk_fraction( get_lk_fractional_contribution( atom2_xyz, rsd2.atom( atom2 ).type(), atom1_waters,
																													closest_water, closest_water_d2_delta ) );

	// 1st term -- atom1-atom2 interaction
	{
		Real const dE_dr_over_r( weight_factor * lk_ball_aniso_weight * lk_fraction * lk_deriv * inv_dis );
		F1 += dE_dr_over_r * f1;
		F2 += dE_dr_over_r * f2;
	}

	// 2nd term -- water1-atom2 interaction
	// note that there's no derivative unless we're in the ramping zone:
	if ( closest_water_d2_delta < ramp_width_A2_ && closest_water_d2_delta > 0.0 ) {
		Vector const & atom1_water_xyz( atom1_waters[ closest_water ] );
		// update f1 and f2 to reflect water-atom2 as the interaction
		f1 = atom1_water_xyz.cross( atom2_xyz );
		f2 = atom1_water_xyz - atom2_xyz;
		// what is the derivative of the lk_fraction term wrt r?
		Real const dE_dr_over_r
			( weight_factor * lk_ball_aniso_weight * lk_score * eval_d_lk_fraction_dr_over_r( closest_water_d2_delta ) );
		F1 += dE_dr_over_r * f1;
		F2 += dE_dr_over_r * f2;
	}
}

/// @note  Assumes that atom1 is the "moving" atom, ie the atom for which eval_atom_derivative was called
/// @note  Calculates the water positions for atom2 if d2 < safe_max_dis2
void
LK_BallEnergy::sum_deriv_contributions_for_atom_pair(
																										 Size const atom1,
																										 conformation::Residue const & rsd1,
																										 LKB_ResidueInfo const & rsd1_info,
																										 Size const atom2,
																										 conformation::Residue const & rsd2,
																										 LKB_ResidueInfo const & rsd2_info,
																										 pose::Pose const &,
																										 EnergyMap const & weights,
																										 Real const cp_weight,
																										 Vector & F1,
																										 Vector & F2
																										 ) const
{
	Real const d2( rsd1.xyz( atom1 ).distance_squared( rsd2.xyz( atom2 ) ) );

	if ( ( d2 >= safe_max_dis2_ ) || ( d2 == Real(0.0) ) ) return; // no contribution

	Vectors const & atom1_waters( rsd1_info.waters()[ atom1 ] );
	Vectors const & atom2_waters( rsd2_info.waters()[ atom2 ] );
	utility::vector1< Real > const & atom1_wts( rsd1_info.atom_weights()[ atom1 ] );
	utility::vector1< Real > const & atom2_wts( rsd2_info.atom_weights()[ atom2 ] );

	sum_deriv_contributions_for_atom_pair_one_way( atom1, rsd1, atom1_waters, atom1_wts, atom2, rsd2, weights,
																								 cp_weight, d2, F1, F2 );

	sum_deriv_contributions_for_atom_pair_one_way( atom2, rsd2, atom2_waters, atom2_wts, atom1, rsd1, weights,
																								 -1.0 * cp_weight, d2, F1, F2 );

}


	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;

// void
// LK_BallEnergy::eval_atom_derivative(
// 	id::AtomID const & id,
// 	pose::Pose const & pose,
// 	kinematics::DomainMap const &, // domain_map,
// 	ScoreFunction const & /*sfxn*/, // needed for non-nblist minimization
// 	EnergyMap const & weights,
// 	Vector & F1,
// 	Vector & F2
// ) const
// {
// 	//PROF_START( basic::LK_BALL_EVAL_ATOM_DERIVATIVE );

// 	Size const idresid = id.rsd();
// 	Size const atom1( id.atomno() );
// 	conformation::Residue const & rsd1( pose.residue( idresid ) );
// 	if ( atom1 > rsd1.nheavyatoms() ) {
// 		//PROF_STOP( basic::LK_BALL_EVAL_ATOM_DERIVATIVE );
// 		return; // NO HYDROGEN INTERACTIONS
// 	}

//  	LKB_ResidueInfo const & rsd1_info( retrieve_lkb_residue_info( pose, idresid ) );

// debug_assert( pose.energies().use_nblist() );

// 	scoring::AtomNeighbors const & nbrs
// 		( pose.energies().nblist( basic::ETABLE_NBLIST ).atom_neighbors( id ) );

// 	Vector f1,f2;
// 	for ( scoring::AtomNeighbors::const_iterator it2=nbrs.begin(),
// 					it2e=nbrs.end(); it2 != it2e; ++it2 ) {
// 		scoring::AtomNeighbor const & nbr( *it2 );
// 		Real const cp_weight( nbr.weight() );
// 		if ( Size( nbr.atomno() ) > pose.residue( nbr.rsd() ).nheavyatoms() ) continue; // NO HYDROGEN INTERACTIONS
// 		LKB_ResidueInfo const & rsd2_info( retrieve_lkb_residue_info( pose, nbr.rsd() ) );
// 		sum_contributions_for_atom_pair( atom1, rsd1, rsd1_info, nbr.atomno(), pose.residue( nbr.rsd() ), rsd2_info, pose,
// 																		 weights, cp_weight, F1, F2 );
// 	}
// 	//PROF_STOP( basic::LK_BALL_EVAL_ATOM_DERIVATIVE );
// }



void
LK_BallEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /* context_graphs_required */ ) const
{}



void
LK_BallEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const &, // pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	using namespace conformation;
	using namespace numeric;

	//TR.Trace << "rotamer_pair" << std::endl;

	//PROF_START( basic::LK_BALL_ROTAMER_PAIR_ENERGIES );

	LKB_RotamerSetInfo const & info1( retrieve_lkb_rotamer_set_info( set1 ) );
	LKB_RotamerSetInfo const & info2( retrieve_lkb_rotamer_set_info( set2 ) );


	for ( Size ii = 1; ii <= set1.get_n_residue_types(); ++ii ) {
		Size const ii_offset = set1.get_residue_type_begin( ii );
		for ( Size jj = 1; jj <= set2.get_n_residue_types(); ++jj ) {
			Size const jj_offset = set2.get_residue_type_begin( jj );
			if ( !info1[ ii_offset ].has_waters() && !info2[ jj_offset ].has_waters() ) continue;

			for ( Size kk = 1, kke = set1.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
				Size const kk_rot_id = ii_offset + kk - 1;
				for ( Size ll = 1, lle = set2.get_n_rotamers_for_residue_type( jj ); ll <= lle; ++ll ) {
					Size const ll_rot_id = jj_offset + ll - 1;
					EnergyMap emap;
					residue_pair_energy( *set1.rotamer( kk_rot_id ), info1[ kk_rot_id ],
						*set2.rotamer( ll_rot_id ), info2[ ll_rot_id ], emap );

					energy_table( ll_rot_id, kk_rot_id ) += static_cast< core::PackerEnergy >( weights.dot( emap ) );
				}
			}
		}
	}
	//PROF_STOP( basic::LK_BALL_ROTAMER_PAIR_ENERGIES );
}

void
LK_BallEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{

	//PROF_START( basic::LK_BALL_ROTAMER_BACKGROUND_ENERGIES );
	//TR.Trace << "rotamer_background!" << std::endl;

	using conformation::Residue;
	LKB_ResidueInfo const & rsd_info( retrieve_lkb_residue_info( pose, rsd.seqpos() ) );
	LKB_RotamerSetInfo const & set_info( retrieve_lkb_rotamer_set_info( set ) );

	for ( Size ii = 1; ii <= set.get_n_residue_types(); ++ii ) {
		Size const ii_offset = set.get_residue_type_begin( ii );

		if ( !set_info[ ii_offset ].has_waters() && !rsd_info.has_waters() ) continue;

		for ( Size kk = 1, kke = set.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
			Size const kk_rot_id = ii_offset + kk - 1;

			EnergyMap emap;
			residue_pair_energy( *set.rotamer( kk_rot_id ), set_info[ kk_rot_id ], rsd, rsd_info, emap );

			/**{ // check that this agrees with the old way!
				// explicitly build the residue_info's now
				LKB_ResidueInfo const info1_redo( *set.rotamer( kk_rot_id ) );
				LKB_ResidueInfo const info2_redo( rsd );
				Real lk_ball_energy_redo( 0.0 ), lk_ball_iso_energy_redo( 0.0 );
				residue_pair_energy( *set.rotamer( kk_rot_id ), info1_redo.waters(),
														 rsd, info2_redo.waters(),
														 lk_ball_energy_redo, lk_ball_iso_energy_redo );
			debug_assert( std::abs( lk_ball_energy - lk_ball_energy_redo ) +
								std::abs( lk_ball_iso_energy - lk_ball_iso_energy_redo ) < 1e-3 );

								}**/

			energy_vector[ kk_rot_id ] += static_cast< core::PackerEnergy >( weights.dot( emap ) );
		} // kk - rotamers for residue types
	} // ii - residue types for rotamer set

	//PROF_STOP( basic::LK_BALL_ROTAMER_BACKGROUND_ENERGIES );
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
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
debug_assert( r1_at_derivs.size() >= rsd1.natoms() );
debug_assert( r2_at_derivs.size() >= rsd2.natoms() );

	// retrieve some info
 	// LKB_ResidueInfo const & rsd1_info( retrieve_lkb_residue_info( pose, rsd1.seqpos() ) );
 	// LKB_ResidueInfo const & rsd2_info( retrieve_lkb_residue_info( pose, rsd2.seqpos() ) );
	LKB_ResidueInfo const & rsd1_info( retrieve_lkb_resdata( res1data ) );
	LKB_ResidueInfo const & rsd2_info( retrieve_lkb_resdata( res2data ) );

debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( etab_pair_nblist ) ));
	ResiduePairNeighborList const & nblist
		( static_cast< ResiduePairNeighborList const & > (min_data.get_data_ref( etab_pair_nblist )) );

	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		Size const atom1( neighbs[ ii ].atomno1() ), atom2( neighbs[ ii ].atomno2() );
		if ( rsd1.atom_is_hydrogen( atom1 ) || rsd2.atom_is_hydrogen( atom2 ) ) continue; // NO HYDROGEN INTXNS
		Real const cp_weight( neighbs[ ii ].weight() );

		Vector f1(0,0,0),f2(0,0,0);
		sum_deriv_contributions_for_atom_pair( atom1, rsd1, rsd1_info, atom2, rsd2, rsd2_info, pose, weights,
																					 cp_weight, f1, f2 );

		r1_at_derivs[ neighbs[ ii ].atomno1() ].f1() += f1;
		r1_at_derivs[ neighbs[ ii ].atomno1() ].f2() += f2;
		r2_at_derivs[ neighbs[ ii ].atomno2() ].f1() += -1*f1;
		r2_at_derivs[ neighbs[ ii ].atomno2() ].f2() += -1*f2;
	}
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
}

Size
get_parallel_h_for_arg(
	chemical::ResidueType const & rsd_type,
	Size const hatm
)
{
	using namespace chemical;
debug_assert( rsd_type.aa() == aa_arg );
	typedef std::map< ResidueTypeCOP, utility::vector1< Size > > H_Map;
	static H_Map hmap;
	ResidueTypeCOP rsd_type_ptr( rsd_type.get_self_ptr() );
	H_Map::const_iterator iter = hmap.find( rsd_type_ptr );
	if ( iter == hmap.end() ) {
		// first time seeing this arginine residue type
		utility::vector1< Size > parallel_h( rsd_type.natoms(), 0 );
		parallel_h[ rsd_type.atom_index(  "HE"  ) ] = rsd_type.atom_index( "1HH2" );
		parallel_h[ rsd_type.atom_index( "1HH2" ) ] = rsd_type.atom_index(  "HE"  );
		parallel_h[ rsd_type.atom_index( "2HH1" ) ] = rsd_type.atom_index( "2HH2" );
		parallel_h[ rsd_type.atom_index( "2HH2" ) ] = rsd_type.atom_index( "2HH1" );
		hmap.insert( std::make_pair( rsd_type_ptr, parallel_h ) );
		iter = hmap.find( rsd_type_ptr );
	}
	return iter->second[ hatm ];
}


/// hack to improve geometries of hbonds:
/// apply an additional weighting factor which is used in LK_BallEnergy
/// Do this for sidechain SP2 or Ring acceptors
/// Do this for all sidechain donors
///
/// NOTE: also excluding multiple parallel Arg hbonds
///
/// @note  For the deriv calculation, we assume that the deriv passed in has already been evaluated by hb_energy_deriv
/// and thus represents the deriv of the unweighted energy considering the donor as the moving atom.
void
apply_lk_ball_fraction_weight_for_hbonds(
																				 Size const hatm,
																				 conformation::Residue const & don_rsd,
																				 Size const aatm,
																				 conformation::Residue const & acc_rsd,
																				 Vector const & hatm_xyz,
																				 Vector const & datm_xyz,
																				 Real & unweighted_energy,
																				 bool const evaluate_derivative,
																				 hbonds::HBondDerivs & hbderivs,
																				 Real & don_fraction,
																				 Real & acc_fraction
																				 )
{
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
	// params
	bool const no_lk_ball_for_SP2( false ); //option[ OKDS::no_lk_ball_for_SP2 ] );
	bool const exclude_arg_double_parallel_hbonds( true );

// not considering these next two for the time being:
// 	bool const use_soft_lk_fraction_for_hbond_acceptors( false );
// 	bool const use_soft_lk_fraction_for_hbond_acceptors_from_water( false );


	Size const datm( don_rsd.atom_base( hatm ) );

	// static for efficiency
	static Real const optimal_acceptor_water_distance( 2.85 ); /// these should agree with LK_BallInfo.cc, right ??
	static Real const optimal_donor_water_distance   ( 2.85 ); /// the values there are 2.65 !!
	static Real const optimal_acceptor_angle_SP2( numeric::conversions::radians( 120.0 ) );
	static Real const n1_coeff_SP2( -1.0 * optimal_acceptor_water_distance * std::cos( optimal_acceptor_angle_SP2 ) );
	static Real const n2_coeff_SP2(        optimal_acceptor_water_distance * std::sin( optimal_acceptor_angle_SP2 ) );
	static methods::EnergyMethodOptions energy_method_options;
	static methods::LK_BallEnergy lk_ball_energy( energy_method_options );

	Vector f1_don( 0.0 ), f2_don( 0.0 ), f1_acc( 0.0 ), f2_acc( 0.0 );
	don_fraction = acc_fraction = 1.0;

	if ( !don_rsd.atom_is_backbone( datm ) ) {
		/// get ideal water location
		utility::vector1< Vector > waters;
		waters.push_back( datm_xyz + optimal_donor_water_distance * ( hatm_xyz - datm_xyz ).normalized_any() );
		if ( don_rsd.aa() == chemical::aa_arg && exclude_arg_double_parallel_hbonds ) {
			Size const other_hatm( get_parallel_h_for_arg( don_rsd.type(), hatm ) );
			if ( other_hatm ) {
				Vector const & other_datm_xyz( don_rsd.xyz( don_rsd.atom_base( other_hatm ) ) ),
					& other_hatm_xyz( don_rsd.xyz( other_hatm ) ), & aatm_xyz( acc_rsd.xyz( aatm ) );
				Vector const other_water( other_datm_xyz + optimal_donor_water_distance *
																	( other_hatm_xyz - other_datm_xyz ).normalized_any() );
				if ( aatm_xyz.distance_squared( other_water ) < aatm_xyz.distance_squared( waters.front() ) ) {
					// closer to the other hydrogen's water!
//  					std::cout << "excluding arg hbond: " << acc_rsd.name1() << ' ' << acc_rsd.seqpos() << ' ' <<
//  						acc_rsd.atom_name( aatm ) << " is closer to ARG " << don_rsd.seqpos() << ' ' <<
//  						don_rsd.atom_name( other_hatm ) << " than to " << don_rsd.atom_name( hatm ) << std::endl;
					don_fraction = acc_fraction = 0.0;
					waters.clear(); // prevent re-evaluation of don_fraction
					if ( evaluate_derivative ) {
						f1_don.zero();
						f2_don.zero();
					}
				}
			}
		}
		if ( !waters.empty() ) {
			// otherwise donor is an arginine with wrong water closer
			don_fraction = lk_ball_energy.eval_lk_ball_fraction_deriv( acc_rsd.xyz( aatm ),acc_rsd.atom( aatm ).type(),
																																 waters, evaluate_derivative, f1_don, f2_don );
		} else {
		debug_assert( don_rsd.aa() == chemical::aa_arg );
		}
	}


	if ( !acc_rsd.atom_is_backbone( aatm ) ) {
		using namespace chemical;
		Hybridization const & hybrid( acc_rsd.atom_type( aatm ).hybridization() );
		if ( hybrid == SP2_HYBRID && no_lk_ball_for_SP2 && acc_rsd.is_protein() ) {
			//TR.Trace << "Not computing acc_fraction " << acc_rsd.atom_name( aatm ) << std::endl;
		} else if ( hybrid == SP2_HYBRID || hybrid == RING_HYBRID ) {
			/// only cases for which we do this:
			Vector const &   aatm_xyz( acc_rsd.xyz( aatm ) );
			Vector const &  abase_xyz( acc_rsd.xyz( acc_rsd.atom_base( aatm ) ) );
			Vector const & abase2_xyz( acc_rsd.xyz( acc_rsd.abase2( aatm ) ) );
			utility::vector1< Vector > waters;

			if ( hybrid == SP2_HYBRID ) {
				Vector const n1( ( aatm_xyz - abase_xyz ).normalized_any() );
				Vector n2( ( abase2_xyz - abase_xyz ) );
				n2 = ( n2 - n1.dot( n2 ) * n1 ).normalized_any();
			debug_assert( std::abs( n1.dot( n2 ) ) < 1e-3 );
				waters.push_back( aatm_xyz + n1_coeff_SP2 * n1 + n2_coeff_SP2 * n2 );
				waters.push_back( aatm_xyz + n1_coeff_SP2 * n1 - n2_coeff_SP2 * n2 );
				//dump_waters_pdb( acc_rsd, waters, "acc_sp2" );
				//dumped[2] = true;
			} else if ( hybrid == RING_HYBRID ) {
				Vector const n1( ( aatm_xyz - 0.5 * abase_xyz - 0.5 * abase2_xyz ).normalized_any() );
				waters.push_back( aatm_xyz + optimal_acceptor_water_distance * n1 );
				//dump_waters_pdb( acc_rsd, waters, "acc_ring" );
				//dumped[3] = true;
			} else {
				utility_exit_with_message("dont do this hybridization state!");
			}

			/// note that this function thinks that the acceptor is the "moving" atom, hence -1.0 in deriv calcs below
			acc_fraction = lk_ball_energy.eval_lk_ball_fraction_deriv( datm_xyz, don_rsd.atom( datm ).type(), waters,
																																 evaluate_derivative, f1_acc, f2_acc );
		}
	}

	if ( acc_fraction < 0.99 || don_fraction < 0.99 ) {
		// debugging
		using namespace ObjexxFCL::format;
// 		std::cout << "apply_lk_ball_fraction_weight_for_hbonds: " <<
// 			F(9,3,don_fraction) << I(4,don_rsd.seqpos()) << ' ' << don_rsd.name1() << ' '<< don_rsd.atom_name( hatm ) <<
// 			F(9,3,acc_fraction) << I(4,acc_rsd.seqpos()) << ' ' << acc_rsd.name1() << ' '<< acc_rsd.atom_name( aatm ) <<
// 			std::endl;
	}


	{ // we used to have an option to soften the dependence, but that was taken out in porting to trunk to get derivs
		// correctly first

		/// do this BEFORE we multiply unweighted_energy by don_fraction and acc_fraction
		if ( evaluate_derivative ) {
			/// note that don_derivs should probably be attached to h_deriv instead, need to think more about this...
			/// also note that f12_don think the donor is the moving atom, while f12_acc think the acceptor is moving
			/// this is the chain rule, baby:
			hbderivs.don_deriv.f1() =  ( hbderivs.don_deriv.f1() * don_fraction * acc_fraction +
																	 unweighted_energy       * f1_don       * acc_fraction +
																	 unweighted_energy       * don_fraction * -1 * f1_acc   );

			hbderivs.don_deriv.f2() =  ( hbderivs.don_deriv.f2() * don_fraction * acc_fraction +
																	 unweighted_energy       * f2_don       * acc_fraction +
																	 unweighted_energy       * don_fraction * -1 * f2_acc   );

			hbderivs.acc_deriv.f1() =  ( hbderivs.acc_deriv.f1() * don_fraction * acc_fraction +
																	 unweighted_energy       * -1 * f1_don  * acc_fraction +
																	 unweighted_energy       * don_fraction * f1_acc   );

			hbderivs.acc_deriv.f2() =  ( hbderivs.acc_deriv.f2() * don_fraction * acc_fraction +
																	 unweighted_energy       * -1 * f2_don  * acc_fraction +
																	 unweighted_energy       * don_fraction * f2_acc   );

			hbderivs.h_deriv.f1() =  ( hbderivs.h_deriv.f1() * don_fraction * acc_fraction );
			hbderivs.h_deriv.f2() =  ( hbderivs.h_deriv.f2() * don_fraction * acc_fraction );

			hbderivs.abase_deriv.f1() =  ( hbderivs.abase_deriv.f1() * don_fraction * acc_fraction );
			hbderivs.abase_deriv.f2() =  ( hbderivs.abase_deriv.f2() * don_fraction * acc_fraction );

			hbderivs.abase2_deriv.f1() =  ( hbderivs.abase2_deriv.f1() * don_fraction * acc_fraction );
			hbderivs.abase2_deriv.f2() =  ( hbderivs.abase2_deriv.f2() * don_fraction * acc_fraction );

			// the old way: (this was just wrt the donor; we had already negated f1_acc,f2_acc)
			//
			// 			deriv.first  =  ( deriv.first       * don_fraction * acc_fraction +
			// 												unweighted_energy * f1_don       * acc_fraction +
			// 												unweighted_energy * don_fraction * f1_acc       );

			// 			deriv.second =  ( deriv.second      * don_fraction * acc_fraction +
			// 												unweighted_energy * f2_don       * acc_fraction +
			// 												unweighted_energy * don_fraction * f2_acc       );

		}

		/// now apply the new weighting factor
		unweighted_energy = unweighted_energy * don_fraction * acc_fraction;

	} // use a softer damping scheme
	//std::cout << "LK_BallEnergy.cc: " << __LINE__ << std::endl;
}





core::Size
LK_BallEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
// 					if ( verbose ) { // HACKING
// 						Real closest_water_dis2(0.0);
// 						Size closest_water(0);
// 						Real const frac( get_lk_fractional_contribution( atom2_xyz, atom2_type_index, atom1_waters, closest_water,
// 																														 closest_water_dis2));
// 						//std::cout << "LK_BALL nwaters= " << atom1_waters.size() << " distance= " << F(9,3,std::sqrt(d2) ) <<
// 							" lk_ball_iso= " << F(9,3,lk_desolvation_of_atom1_by_atom2)<<
// 							" lk_fraction= " << F(9,3,frac) <<
// 							" closest_water_dis= " << F(9,3,std::sqrt( d2_low_[ atom2_type_index ] + closest_water_dis2 ) ) <<
// 							" desolvated_atom: " << rsd1.name() << ' ' << rsd1.seqpos() << ' ' << rsd1.atom_name( atom1 ) <<
// 							" desolvating_atom: " << rsd2.name() << ' ' << rsd2.seqpos() << ' ' << rsd2.atom_name( atom2 ) <<
// 							std::endl;
// 					}



// 					if ( verbose ) { // HACKING
// 						Real closest_water_dis2(0.0);
// 						Size closest_water(0);
// 						Real const frac( get_lk_fractional_contribution( atom1_xyz, atom1_type_index, atom2_waters, closest_water,
// 																														 closest_water_dis2));
// 						//std::cout << "LK_BALL nwaters= " << atom2_waters.size() << " distance= " << F(9,3,std::sqrt(d2) ) <<
// 							" lk_ball_iso= " << F(9,3,lk_desolvation_of_atom2_by_atom1)<<
// 							" lk_fraction= " << F(9,3,frac) <<
// 							" closest_water_dis= " << F(9,3,std::sqrt( d2_low_[ atom1_type_index ] + closest_water_dis2 ) ) <<
// 							" desolvated_atom: " << rsd2.name() << ' ' << rsd2.seqpos() << ' ' << rsd2.atom_name( atom2 ) <<
// 							" desolvating_atom: " << rsd1.name() << ' ' << rsd1.seqpos() << ' ' << rsd1.atom_name( atom1 ) <<
// 							std::endl;
// 					}




// 				if ( (  atom1_is_polar && lk_desolvation_of_atom1_by_atom2 < -1e-3 ) || // debugging
// 						 ( !atom1_is_polar && lk_desolvation_of_atom1_by_atom2 > 1e-3 ) ) {
// 					//std::cout << "lk_funny_sign: " << lk_desolvation_of_atom1_by_atom2 << ' ' <<
// 						rsd1.atom_name( atom1 ) << ' ' << rsd1.name() << std::endl;
// 				}

// 				if ( (  atom2_is_polar && lk_desolvation_of_atom2_by_atom1 < -1e-3 ) || // debugging
// 						 ( !atom2_is_polar && lk_desolvation_of_atom2_by_atom1 > 1e-3 ) ) {
// 					//std::cout << "lk_funny_sign: " << lk_desolvation_of_atom2_by_atom1 << ' ' <<
// 						rsd2.atom_name( atom2 ) << ' ' << rsd2.name() << std::endl;
// 				}
