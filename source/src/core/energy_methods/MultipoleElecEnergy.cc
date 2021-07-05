// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/MultipoleElecEnergy.cc
/// @brief  Fixed multipole electrostatics using Jay Ponder's approach in Tinker/Amoeba
/// @author Jim Havranek


// Unit headers
#include <core/energy_methods/MultipoleElecEnergy.hh>
#include <core/energy_methods/MultipoleElecEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/MultipoleElecPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/DenseEnergyContainer.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/conformation/Residue.hh>


#include <utility/vector1.hh>

static basic::Tracer TR( "core.energy_methods.MultipoleElecEnergy" );

// Utility headers


// C++
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


namespace core {
namespace energy_methods {



inline
core::scoring::MultipoleElecResidueInfo const &
retrieve_mp_residue_info( pose::Pose const & pose, Size const seqpos ) {
	debug_assert( seqpos && seqpos <= ( static_cast< core::scoring::MultipoleElecPoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) )).size() );
	return ( static_cast< core::scoring::MultipoleElecPoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ).residue_info( seqpos ) );
}

inline
core::scoring::MultipoleElecResidueInfo &
retrieve_nonconst_mp_residue_info( pose::Pose & pose, Size const seqpos ) {
	debug_assert( seqpos && seqpos <= ( static_cast< core::scoring::MultipoleElecPoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) )).size() );
	return ( static_cast< core::scoring::MultipoleElecPoseInfo & >
		( pose.data().get( pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ).residue_info( seqpos ) );
}

class MultipoleElecResPairMinData : public basic::datacache::CacheableData
{
public:
	MultipoleElecResPairMinData();
	~MultipoleElecResPairMinData() override = default;
	basic::datacache::CacheableDataOP clone() const override
	{ return utility::pointer::make_shared< MultipoleElecResPairMinData >( *this ); }

	void
	initialize(
		core::scoring::MultipoleElecResidueInfoCOP res1_data,
		core::scoring::MultipoleElecResidueInfoCOP res2_data
	);

	core::scoring::MultipoleElecResidueInfo const & res1_data() const { return *res1_data_; }
	core::scoring::MultipoleElecResidueInfo const & res2_data() const { return *res2_data_; }

	bool
	initialized() const { return initialized_; }


private:

	core::scoring::MultipoleElecResidueInfoCOP res1_data_;
	core::scoring::MultipoleElecResidueInfoCOP res2_data_;

	bool initialized_{ false };
};

using MultipoleElecResPairMinDataOP = utility::pointer::shared_ptr<MultipoleElecResPairMinData>;
using MultipoleElecResPairMinDataCOP = utility::pointer::shared_ptr<const MultipoleElecResPairMinData>;

MultipoleElecResPairMinData::MultipoleElecResPairMinData()= default;


void
MultipoleElecResPairMinData::initialize(
	core::scoring::MultipoleElecResidueInfoCOP res1_data,
	core::scoring::MultipoleElecResidueInfoCOP res2_data
)
{
	initialized_ = true;
	res1_data_ = res1_data;
	res2_data_ = res2_data;
}

/////////////////////////////////////// mindata retrieval functions

inline
MultipoleElecResPairMinData &
retrieve_nonconst_mp_pairdata(
	core::scoring::ResPairMinimizationData & pairdata
)
{
	MultipoleElecResPairMinDataOP mp_pairdata(nullptr);
	if ( pairdata.get_data( core::scoring::mp_respair_data ) ) {
		debug_assert( utility::pointer::dynamic_pointer_cast< MultipoleElecResPairMinData > ( pairdata.get_data( core::scoring::mp_respair_data ) ));
		mp_pairdata = utility::pointer::static_pointer_cast< MultipoleElecResPairMinData > ( pairdata.get_data( core::scoring::mp_respair_data ) );
	} else {
		mp_pairdata = utility::pointer::make_shared< MultipoleElecResPairMinData >();
		pairdata.set_data( core::scoring::mp_respair_data, mp_pairdata );
	}
	return *mp_pairdata;
}

inline
MultipoleElecResPairMinData const &
retrieve_mp_pairdata(
	core::scoring::ResPairMinimizationData const & pairdata
)
{
	debug_assert( utility::pointer::dynamic_pointer_cast< MultipoleElecResPairMinData const > ( pairdata.get_data( core::scoring::mp_respair_data ) ) );
	return ( static_cast< MultipoleElecResPairMinData const & > ( pairdata.get_data_ref( core::scoring::mp_respair_data ) ) );
}

inline
core::scoring::MultipoleElecResidueInfo &
retrieve_nonconst_mp_resdata(
	core::scoring::ResSingleMinimizationData & resdata
)
{
	using namespace core::scoring;
	MultipoleElecResidueInfoOP mp_resdata( nullptr );
	if ( resdata.get_data( mp_res_data ) ) {
		debug_assert( utility::pointer::dynamic_pointer_cast< MultipoleElecResidueInfo > ( resdata.get_data( mp_res_data ) ) );
		mp_resdata = utility::pointer::static_pointer_cast< MultipoleElecResidueInfo > ( resdata.get_data( mp_res_data ) );
	} else {
		mp_resdata = utility::pointer::make_shared< MultipoleElecResidueInfo >();
		resdata.set_data( mp_res_data, mp_resdata );
	}
	return *mp_resdata;
}

inline
core::scoring::MultipoleElecResidueInfo const &
retrieve_mp_resdata(
	core::scoring::ResSingleMinimizationData const & resdata
)
{
	return ( static_cast< core::scoring::MultipoleElecResidueInfo const & > ( resdata.get_data_ref( core::scoring::mp_res_data ) ) );
}

inline
core::scoring::MultipoleElecResidueInfoCOP
retrieve_mp_resdata_ptr(
	core::scoring::ResSingleMinimizationData const & resdata
)
{
	debug_assert( utility::pointer::dynamic_pointer_cast< core::scoring::MultipoleElecResidueInfo const > ( resdata.get_data( core::scoring::mp_res_data ) ) );
	return ( utility::pointer::static_pointer_cast< core::scoring::MultipoleElecResidueInfo const > ( resdata.get_data( core::scoring::mp_res_data ) ) );
}




/// @details This must return a fresh instance of the MultipoleElecEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
MultipoleElecEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< MultipoleElecEnergy >( options );
}

core::scoring::ScoreTypes
MultipoleElecEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( multipole_elec );
	return sts;
}


MultipoleElecEnergy::MultipoleElecEnergy( MultipoleElecEnergy const & /*src*/ ) = default;


MultipoleElecEnergy::MultipoleElecEnergy( core::scoring::methods::EnergyMethodOptions const & options ):
	parent( utility::pointer::make_shared< MultipoleElecEnergyCreator >() ),
	potential_( core::scoring::ScoringManager::get_instance()->get_MultipoleElecPotential( options ) ),
	exclude_DNA_DNA_( options.exclude_DNA_DNA() )
{}


/// clone
core::scoring::methods::EnergyMethodOP
MultipoleElecEnergy::clone() const
{
	return utility::pointer::make_shared< MultipoleElecEnergy >( *this );
}

bool
MultipoleElecEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size,
	Size
) const
{
	return true;
}


///
void
MultipoleElecEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const & , // residues_repacking,
	utility::vector1< bool > const &
) const
{

	potential_.setup_for_scoring( pose );
	pose.update_residue_neighbors();

}

void
MultipoleElecEnergy::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & rotamer_set
) const
{
	// Need to assign types, etc.
	potential_.get_rotamers_multipole_info( pose, rotamer_set );

	// This will calculate effective Kirkwood radii and cache
	potential_.get_rotamers_effective_radii( pose, rotamer_set );

	// At this point, if we were doing polarization, it would go here.
}


//  void
//  update_residue_for_packing(
//   pose::Pose &,
//   Size /*resid*/ ) const
//  {}
void
MultipoleElecEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	/// update amoeba type information for residue that has changed during packing, eg in rotamer trials
	/// need to double-check the current logic on this...
	potential_.update_residue_for_packing( pose, resid );
}

///
void
MultipoleElecEnergy::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	core::scoring::methods::LongRangeEnergyType const & lr_type( long_range_type() );

	potential_.setup_for_scoring( pose );

	// create a container
	core::scoring::Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;
	} else {
		core::scoring::LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		core::scoring::DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		core::scoring::LREnergyContainerOP new_dec = utility::pointer::make_shared< core::scoring::DenseEnergyContainer >( pose.size(), core::scoring::multipole_elec );
		energies.set_long_range_container( lr_type, new_dec );
	}
}


///
void
MultipoleElecEnergy::setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	core::scoring::methods::LongRangeEnergyType const & lr_type( long_range_type() );

	potential_.setup_for_scoring( pose );

	// create a container
	core::scoring::Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;
	} else {
		core::scoring::LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		core::scoring::DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		core::scoring::LREnergyContainerOP new_dec = utility::pointer::make_shared< core::scoring::DenseEnergyContainer >( pose.size(), core::scoring::multipole_elec );
		energies.set_long_range_container( lr_type, new_dec );
	}

	potential_.calculate_and_store_all_derivs( pose );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
MultipoleElecEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core::scoring;
	//using core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO;
	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	MultipoleElecResidueInfo const & mp1( retrieve_mp_residue_info( pose, rsd1.seqpos() ) );
	MultipoleElecResidueInfo const & mp2( retrieve_mp_residue_info( pose, rsd2.seqpos() ) );

	//TR << "Calculating residue pair energy" << std::endl;

	emap[ core::scoring::multipole_elec ] += potential_.get_res_res_elecE( rsd1, mp1, rsd2, mp2 );
}


/////////////////////////////////
// Minimization specific stuff
/////////////////////////////////

void
MultipoleElecEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose, // pose,
	core::scoring::ScoreFunction const &, // scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	basic::datacache::BasicDataCache &,
	core::scoring::ResSingleMinimizationData & resdata
) const
{
	using namespace core::scoring;
	auto const & mp_info
		( static_cast< MultipoleElecPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) );

	MultipoleElecResidueInfo const & mp_pose_data( mp_info.residue_info( rsd.seqpos() ) );

	MultipoleElecResidueInfoOP mp_resdata( mp_pose_data.copy_clone() );

	resdata.set_data( mp_res_data, mp_resdata );
	// Assign and rotate moment information
	//TR << "Calling assign_residue from setup_from minimizing" << std::endl;
	potential_.assign_residue_amoeba_type( rsd, *mp_resdata );
	potential_.align_residue_multipole_axes( pose, rsd, *mp_resdata );
	// Copy over the induced dipoles, which must be obtained from
	// whole structure relaxation.
}


void
MultipoleElecEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const &, // rsd1,
	conformation::Residue const &, // rsd2,
	pose::Pose const &,
	core::scoring::ScoreFunction const &, //scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	core::scoring::ResSingleMinimizationData const & res1data,
	core::scoring::ResSingleMinimizationData const & res2data,
	core::scoring::ResPairMinimizationData & pairdata
) const
{
	MultipoleElecResPairMinData & mp_pairdata( retrieve_nonconst_mp_pairdata( pairdata ) );
	mp_pairdata.initialize( retrieve_mp_resdata_ptr( res1data ),
		retrieve_mp_resdata_ptr( res2data ) );
}


bool
MultipoleElecEnergy::requires_a_setup_for_scoring_for_residue_opportunity_during_minimization( pose::Pose const & ) const
{
	return true;
}

void
MultipoleElecEnergy::setup_for_scoring_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &, // sfxn,
	core::scoring::ResSingleMinimizationData & resdata
) const
{
	using namespace core::scoring;
	MultipoleElecResidueInfo & info( retrieve_nonconst_mp_resdata( resdata ) );
	// Assign and rotate moment information
	// TR << "Calling assign_residue from setup_for_scoring_for_residue" << std::endl;
	potential_.assign_residue_amoeba_type( rsd, info );
	potential_.align_residue_multipole_axes( pose, rsd, info );
}

bool
MultipoleElecEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const &  ) const
{
	return true;
}


void
MultipoleElecEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::scoring::ResSingleMinimizationData & min_data,
	basic::datacache::BasicDataCache &
) const
{
	/// just compute water locations
	setup_for_scoring_for_residue( rsd, pose, sfxn, min_data );
}


bool
MultipoleElecEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

void
MultipoleElecEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::scoring::ResPairMinimizationData const & , // pairdata,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core::scoring;
	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	//TR << "Calculating residue pair energy ext" << std::endl;

	MultipoleElecResidueInfo const & mp1( retrieve_mp_residue_info( pose, rsd1.seqpos() ) );
	MultipoleElecResidueInfo const & mp2( retrieve_mp_residue_info( pose, rsd2.seqpos() ) );

	emap[ core::scoring::multipole_elec ] += potential_.get_res_res_elecE( rsd1, mp1, rsd2, mp2 );
}

bool
MultipoleElecEnergy::use_extended_intrares_energy_interface() const
{
	return true;
}

void
MultipoleElecEnergy::eval_intrares_energy_ext(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const & data_cache,
	pose::Pose const & ,
	core::scoring::ScoreFunction const & ,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core::scoring;
	if ( exclude_DNA_DNA_ && rsd.is_DNA() ) return;

	//TR << "Calculating intraresidue energy ext" << std::endl;

	MultipoleElecResidueInfo const & mp_info( retrieve_mp_resdata( data_cache ) );
	emap[ core::scoring::multipole_elec ] += potential_.get_res_res_elecE( rsd, mp_info, rsd, mp_info );
}







/////////////////////////////////
// End minimization specific stuff
/////////////////////////////////

void
MultipoleElecEnergy::evaluate_rotamer_intrares_energies(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	utility::vector1< core::PackerEnergy > & energies
) const
{
	using namespace core::scoring;
	using core::conformation::RotamerSetCacheableDataType::MULTIPOLE_ELEC_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set.resid() ).is_DNA() ) return;

	auto const & mp_info
		( set.data().get< MultipoleElecRotamerSetInfo >( MULTIPOLE_ELEC_ROTAMER_SET_INFO ) );

	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {

		//TR << "Calculating rotamer intraresidue energy ext" << std::endl;

		Real const elecE
			( potential_.get_res_res_elecE( *set.rotamer( ii ), mp_info.residue_info( ii ),
			*set.rotamer( ii ), mp_info.residue_info( ii ) ) );

		energies[ ii ] += static_cast< core::PackerEnergy > ( sfxn[ multipole_elec ] * elecE );
	}
}

void
MultipoleElecEnergy::evaluate_rotamer_intrares_energy_maps(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & , // sxfn,
	utility::vector1< core::scoring::EnergyMap > & emaps
) const
{
	using namespace core::scoring;
	using core::conformation::RotamerSetCacheableDataType::MULTIPOLE_ELEC_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set.resid() ).is_DNA() ) return;

	// std::cout << "MultipoleElec rotamer_intrares_energies: " << set.resid() << std::endl;
	auto const & mp_info
		( set.data().get< MultipoleElecRotamerSetInfo >( MULTIPOLE_ELEC_ROTAMER_SET_INFO ) );

	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		//TR << "Calculating rotamer intraresidue energy maps" << std::endl;
		Real const elecE
			( potential_.get_res_res_elecE( *set.rotamer( ii ), mp_info.residue_info( ii ),
			*set.rotamer( ii ), mp_info.residue_info( ii ) ) );

		(emaps[ ii ])[ multipole_elec ] += elecE ;
	}
}

void
MultipoleElecEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & , // sfxn,
	core::scoring::EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	using namespace conformation;
	using namespace numeric;
	using namespace core::scoring;
	using core::conformation::RotamerSetCacheableDataType::MULTIPOLE_ELEC_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set1.resid() ).is_DNA() && pose.residue( set2.resid() ).is_DNA() ) return;

	auto const & mp_info1
		( set1.data().get< MultipoleElecRotamerSetInfo >( MULTIPOLE_ELEC_ROTAMER_SET_INFO ) );

	auto const & mp_info2
		( set2.data().get< MultipoleElecRotamerSetInfo >( MULTIPOLE_ELEC_ROTAMER_SET_INFO ) );

	for ( Size ii = 1; ii <= set1.get_n_residue_types(); ++ii ) {
		Size const ii_offset = set1.get_residue_type_begin( ii );
		Residue const & ii_example_rotamer( *set1.rotamer( ii_offset ));

		Vector const & ii_coord( ii_example_rotamer.nbr_atom_xyz() );
		Real const ii_radius( ii_example_rotamer.nbr_radius() );

		for ( Size jj = 1; jj <= set2.get_n_residue_types(); ++jj ) {
			Size const jj_offset = set2.get_residue_type_begin( jj );
			Residue const & jj_example_rotamer( *set2.rotamer( jj_offset ));

			if ( exclude_DNA_DNA_ && ii_example_rotamer.is_DNA() && jj_example_rotamer.is_DNA() ) continue;

			Vector const & jj_coord( jj_example_rotamer.nbr_atom_xyz() );
			Real const jj_radius( jj_example_rotamer.nbr_radius() );

			if ( ii_coord.distance_squared( jj_coord ) >= std::pow(ii_radius+jj_radius+packing_interaction_cutoff(), 2 ) ) continue;

			for ( Size kk = 1, kke = set1.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
				Size const kk_rot_id = ii_offset + kk - 1;
				for ( Size ll = 1, lle = set2.get_n_rotamers_for_residue_type( jj ); ll <= lle; ++ll ) {
					Size const ll_rot_id = jj_offset + ll - 1;

					//TR << "Calculating rotamer-rotamer pair energy" << std::endl;

					Real const elecE(
						potential_.get_res_res_elecE( *set1.rotamer( kk_rot_id ), mp_info1.residue_info( kk_rot_id ),
						*set2.rotamer( ll_rot_id ), mp_info2.residue_info( ll_rot_id ) ) );

					energy_table( ll_rot_id, kk_rot_id ) +=
						static_cast< core::PackerEnergy >( weights[ core::scoring::multipole_elec ] *  elecE );
				}
			}
		}
	}
}

void
MultipoleElecEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & , // sfxn,
	core::scoring::EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{

	using conformation::Residue;
	using namespace core::scoring;
	using core::conformation::RotamerSetCacheableDataType::MULTIPOLE_ELEC_ROTAMER_SET_INFO;

	auto const & mp_set_info
		( set.data().get< MultipoleElecRotamerSetInfo >( MULTIPOLE_ELEC_ROTAMER_SET_INFO ) );

	MultipoleElecResidueInfo const & mp_rsd_info( retrieve_mp_residue_info( pose, rsd.seqpos() ) );

	for ( Size ii = 1; ii <= set.get_n_residue_types(); ++ii ) {
		Size const ii_offset = set.get_residue_type_begin( ii );
		Residue const & ii_example_rotamer( *set.rotamer( ii_offset ));

		Vector const & ii_coord( ii_example_rotamer.nbr_atom_xyz() );
		Real const ii_radius( ii_example_rotamer.nbr_radius() );

		if ( exclude_DNA_DNA_ && ii_example_rotamer.is_DNA() && rsd.is_DNA() ) continue;

		Vector const & jj_coord( rsd.nbr_atom_xyz() );
		Real const jj_radius( rsd.nbr_radius() );

		if ( ii_coord.distance_squared( jj_coord ) >= std::pow(ii_radius+jj_radius+packing_interaction_cutoff(), 2 ) ) continue;

		for ( Size kk = 1, kke = set.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
			Size const kk_rot_id = ii_offset + kk - 1;

			//TR << "Calculating rotamer-background pair energy" << std::endl;

			Real const elecE(
				potential_.get_res_res_elecE( *set.rotamer( kk_rot_id ), mp_set_info.residue_info( kk_rot_id ),
				rsd, mp_rsd_info ) );
			energy_vector[ kk_rot_id ] += static_cast< core::PackerEnergy > (weights[ core::scoring::multipole_elec ] *  elecE );
		} // kk - rotamers for residue types
	} // ii - residue types for rotamer set
}

void
MultipoleElecEnergy::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & , // sfxn,
	core::scoring::EnergyMap const & ,
	utility::vector1< core::scoring::EnergyMap > & emaps
) const
{

	using conformation::Residue;
	using namespace core::scoring;
	using core::conformation::RotamerSetCacheableDataType::MULTIPOLE_ELEC_ROTAMER_SET_INFO;

	auto const & mp_set_info
		( set.data().get< MultipoleElecRotamerSetInfo >( MULTIPOLE_ELEC_ROTAMER_SET_INFO ) );

	MultipoleElecResidueInfo const & mp_rsd_info( retrieve_mp_residue_info( pose, rsd.seqpos() ) );

	for ( Size ii = 1; ii <= set.get_n_residue_types(); ++ii ) {
		Size const ii_offset = set.get_residue_type_begin( ii );
		Residue const & ii_example_rotamer( *set.rotamer( ii_offset ));

		Vector const & ii_coord( ii_example_rotamer.nbr_atom_xyz() );
		Real const ii_radius( ii_example_rotamer.nbr_radius() );

		if ( exclude_DNA_DNA_ && ii_example_rotamer.is_DNA() && rsd.is_DNA() ) continue;

		Vector const & jj_coord( rsd.nbr_atom_xyz() );
		Real const jj_radius( rsd.nbr_radius() );

		if ( ii_coord.distance_squared( jj_coord ) >= std::pow(ii_radius+jj_radius+packing_interaction_cutoff(), 2 ) ) continue;
		for ( Size kk = 1, kke = set.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
			Size const kk_rot_id = ii_offset + kk - 1;

			//TR << "Calculating rotamer-background pair energy maps" << std::endl;

			Real const elecE
				( potential_.get_res_res_elecE( *set.rotamer( kk_rot_id ), mp_set_info.residue_info( kk_rot_id ),
				rsd, mp_rsd_info ) );
			(emaps[ kk_rot_id ])[ multipole_elec ] += elecE;
		} // kk - rotamers for residue types
	} // ii - residue types for rotamer set
}

/////////////////////////////////////////////////////////////////////////////
///
#ifdef NOTDEF
void
MultipoleElecEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	potential_.eval_atom_derivative( atom_id, weights[ core::scoring::multipole_elec ], pose, domain_map, exclude_DNA_DNA_, F1, F2 );
}
#endif

/// @brief MultipoleElecEnergy requires no context graphs
void
MultipoleElecEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{
}

/// @brief MultipoleElecEnergy does define intraresidue interactions
bool
MultipoleElecEnergy::defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const
{
	return true;
}

void
MultipoleElecEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core::scoring;
	if ( exclude_DNA_DNA_ && rsd.is_DNA() ) return;

	auto const & mp_info
		( static_cast< MultipoleElecPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) ); // SHOULD BE FAST!

	//TR << "Calculating intraresidue energy" << std::endl;

	emap[ core::scoring::multipole_elec ] += potential_.get_res_res_elecE( rsd, mp_info.residue_info( rsd.seqpos() ),
		rsd, mp_info.residue_info( rsd.seqpos() ) );
}

void
MultipoleElecEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const & min_data,
	pose::Pose const & pose,
	core::scoring::EnergyMap const & weights,
	utility::vector1< core::scoring::DerivVectorPair > & atom_derivs
) const
{
	using namespace core::scoring;

	MultipoleElecResidueInfo const & mp( retrieve_mp_resdata( min_data ) );

	Real const factor( weights[ core::scoring::multipole_elec] );

	potential_.eval_residue_pair_derivatives( rsd, rsd, mp, mp, pose, factor,
		atom_derivs, atom_derivs );
}


void
MultipoleElecEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & ,
	conformation::Residue const & ,
	core::scoring::ResSingleMinimizationData const & ,
	core::scoring::ResSingleMinimizationData const & ,
	core::scoring::ResPairMinimizationData const & ,
	pose::Pose const & , // provides context
	core::scoring::EnergyMap const & ,
	utility::vector1< core::scoring::DerivVectorPair > & ,
	utility::vector1< core::scoring::DerivVectorPair > &
) const
{
}


core::Size
MultipoleElecEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
