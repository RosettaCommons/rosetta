// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/VdWTinkerEnergy.cc
/// @brief  VdW treatment using buffered 14-7 approach in Tinker/Amoeba
/// @author Jim Havranek


// Unit headers
#include <core/scoring/methods/VdWTinkerEnergy.hh>
#include <core/scoring/methods/VdWTinkerEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/VdWTinkerPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/DenseEnergyContainer.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>

#include <numeric/xyz.io.hh>

#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.VdWTinkerEnergy" );

// Utility headers


// C++
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


namespace core {
namespace scoring {
namespace methods {

inline
VdWTinkerResidueInfo const &
retrieve_vdw_residue_info( pose::Pose const & pose, Size const seqpos ) {
	assert( seqpos && seqpos <= ( static_cast< VdWTinkerPoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) )).size() );
	return ( static_cast< VdWTinkerPoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) ).residue_info( seqpos ) );
}

inline
VdWTinkerResidueInfo &
retrieve_nonconst_vdw_residue_info( pose::Pose & pose, Size const seqpos ) {
	assert( seqpos && seqpos <= ( static_cast< VdWTinkerPoseInfo const & >
		( pose.data().get( pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) )).size() );
	return ( static_cast< VdWTinkerPoseInfo & >
		( pose.data().get( pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) ).residue_info( seqpos ) );
}

class VdWTinkerResPairMinData : public basic::datacache::CacheableData
{
public:
	VdWTinkerResPairMinData();
	virtual ~VdWTinkerResPairMinData() {}
	virtual basic::datacache::CacheableDataOP clone() const
	{ return basic::datacache::CacheableDataOP( new VdWTinkerResPairMinData( *this ) ); }

	void
	initialize(
		VdWTinkerResidueInfoCOP res1_data,
		VdWTinkerResidueInfoCOP res2_data
	);

	VdWTinkerResidueInfo const & res1_data() const { return *res1_data_; }
	VdWTinkerResidueInfo const & res2_data() const { return *res2_data_; }

	bool
	initialized() const { return initialized_; }


private:

	VdWTinkerResidueInfoCOP res1_data_;
	VdWTinkerResidueInfoCOP res2_data_;

	bool initialized_;
};

typedef utility::pointer::shared_ptr< VdWTinkerResPairMinData >       VdWTinkerResPairMinDataOP;
typedef utility::pointer::shared_ptr< VdWTinkerResPairMinData const > VdWTinkerResPairMinDataCOP;

VdWTinkerResPairMinData::VdWTinkerResPairMinData():
	initialized_( false )
{}


void
VdWTinkerResPairMinData::initialize(
	VdWTinkerResidueInfoCOP res1_data,
	VdWTinkerResidueInfoCOP res2_data
) {
	initialized_ = true;
	res1_data_ = res1_data;
	res2_data_ = res2_data;
}

/////////////////////////////////////// mindata retrieval functions

inline
VdWTinkerResPairMinData &
retrieve_nonconst_vdw_pairdata(
	ResPairMinimizationData & pairdata
) {
	VdWTinkerResPairMinDataOP vdw_pairdata(0);
	if ( pairdata.get_data( vdw_respair_data ) ) {
		assert( utility::pointer::dynamic_pointer_cast< VdWTinkerResPairMinData > ( pairdata.get_data( vdw_respair_data ) ));
		vdw_pairdata = utility::pointer::static_pointer_cast< VdWTinkerResPairMinData > ( pairdata.get_data( vdw_respair_data ) );
	} else {
		vdw_pairdata = VdWTinkerResPairMinDataOP( new VdWTinkerResPairMinData );
		pairdata.set_data( vdw_respair_data, vdw_pairdata );
	}
	return *vdw_pairdata;
}

inline
VdWTinkerResPairMinData const &
retrieve_vdw_pairdata(
	ResPairMinimizationData const & pairdata
) {
	assert( utility::pointer::dynamic_pointer_cast< VdWTinkerResPairMinData const > ( pairdata.get_data( vdw_respair_data ) ) );
	return ( static_cast< VdWTinkerResPairMinData const & > ( pairdata.get_data_ref( vdw_respair_data ) ) );

}

inline
VdWTinkerResidueInfo &
retrieve_nonconst_vdw_resdata(
	ResSingleMinimizationData & resdata
) {
	VdWTinkerResidueInfoOP vdw_resdata( 0 );
	if ( resdata.get_data( vdw_res_data ) ) {
		assert( utility::pointer::dynamic_pointer_cast< VdWTinkerResidueInfo > ( resdata.get_data( vdw_res_data ) ) );
		vdw_resdata = utility::pointer::static_pointer_cast< VdWTinkerResidueInfo > ( resdata.get_data( vdw_res_data ) );
	} else {
		vdw_resdata = VdWTinkerResidueInfoOP( new VdWTinkerResidueInfo );
		resdata.set_data( vdw_res_data, vdw_resdata );
	}
	return *vdw_resdata;
}

inline
VdWTinkerResidueInfo const &
retrieve_vdw_resdata(
	ResSingleMinimizationData const & resdata
) {
	return ( static_cast< VdWTinkerResidueInfo const & > ( resdata.get_data_ref( vdw_res_data ) ) );
}

inline
VdWTinkerResidueInfoCOP
retrieve_vdw_resdata_ptr(
	ResSingleMinimizationData const & resdata
) {
	assert( utility::pointer::dynamic_pointer_cast< VdWTinkerResidueInfo const > ( resdata.get_data( vdw_res_data ) ) );
	return ( utility::pointer::static_pointer_cast< VdWTinkerResidueInfo const > ( resdata.get_data( vdw_res_data ) ) );
}




/// @details This must return a fresh instance of the VdWTinkerEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
VdWTinkerEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new VdWTinkerEnergy( options ) );
}

ScoreTypes
VdWTinkerEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_vdw_tinker );
	return sts;
}


VdWTinkerEnergy::VdWTinkerEnergy( VdWTinkerEnergy const & src ):
	parent( src ),
	potential_( src.potential_ ),
	exclude_DNA_DNA_( src.exclude_DNA_DNA_ )
{}


VdWTinkerEnergy::VdWTinkerEnergy( EnergyMethodOptions const & options ):
	parent( methods::EnergyMethodCreatorOP( new VdWTinkerEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_VdWTinkerPotential() ),
	exclude_DNA_DNA_( options.exclude_DNA_DNA() )
{}


/// clone
EnergyMethodOP
VdWTinkerEnergy::clone() const
{
	return EnergyMethodOP( new VdWTinkerEnergy( *this ) );
}

bool
VdWTinkerEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size,
	Size
) const
{
	return true;
}


///
void
VdWTinkerEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const & , // residues_repacking,
	utility::vector1< bool > const &
) const {
	potential_.setup_for_scoring( pose );
	pose.update_residue_neighbors();

}

void
VdWTinkerEnergy::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & rotamer_set
) const {
	// Need to assign types, etc.
	potential_.get_rotamers_vdw_info( pose, rotamer_set );
}

void
VdWTinkerEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const {
	/// update amoeba type information for residue that has changed during packing, eg in rotamer trials
	/// need to double-check the current logic on this...
	potential_.update_residue_for_packing( pose, resid );
}

///
void
VdWTinkerEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	LongRangeEnergyType const & lr_type( long_range_type() );

	potential_.setup_for_scoring( pose );

	// create a container
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.total_residue() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		LREnergyContainerOP new_dec = LREnergyContainerOP( new DenseEnergyContainer( pose.total_residue(), fa_vdw_tinker ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}


///
void
VdWTinkerEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	LongRangeEnergyType const & lr_type( long_range_type() );

	potential_.setup_for_scoring( pose );

	// create a container
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.total_residue() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		LREnergyContainerOP new_dec = LREnergyContainerOP( new DenseEnergyContainer( pose.total_residue(), fa_vdw_tinker ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
VdWTinkerEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	//using core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO;
	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	VdWTinkerResidueInfo const & mp1( retrieve_vdw_residue_info( pose, rsd1.seqpos() ) );
	VdWTinkerResidueInfo const & mp2( retrieve_vdw_residue_info( pose, rsd2.seqpos() ) );

	emap[ fa_vdw_tinker ] += potential_.get_res_res_vdw( rsd1, mp1, rsd2, mp2 );
}


/////////////////////////////////
// Minimization specific stuff
/////////////////////////////////

void
VdWTinkerEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose, // pose,
	ScoreFunction const &, // scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	ResSingleMinimizationData & resdata
) const {
	VdWTinkerPoseInfo const & vdw_info
		( static_cast< VdWTinkerPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) ) );

	VdWTinkerResidueInfo const & vdw_pose_data( vdw_info.residue_info( rsd.seqpos() ) );

	VdWTinkerResidueInfoOP vdw_resdata( vdw_pose_data.copy_clone() );

	resdata.set_data( vdw_res_data, vdw_resdata );
	// Assign and rotate moment information
	//TR << "Calling assign_residue from setup_from minimizing" << std::endl;
	potential_.assign_residue_amoeba_type( rsd, *vdw_resdata );
	// Copy over the induced dipoles, which must be obtained from
	// whole structure relaxation.
}


void
VdWTinkerEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const &, // rsd1,
	conformation::Residue const &, // rsd2,
	pose::Pose const &,
	ScoreFunction const &, //scorefxn,
	kinematics::MinimizerMapBase const &, // min_map,
	ResSingleMinimizationData const & res1data,
	ResSingleMinimizationData const & res2data,
	ResPairMinimizationData & pairdata
) const {
	VdWTinkerResPairMinData & vdw_pairdata( retrieve_nonconst_vdw_pairdata( pairdata ) );
	vdw_pairdata.initialize( retrieve_vdw_resdata_ptr( res1data ),
		retrieve_vdw_resdata_ptr( res2data ) );
}


bool
VdWTinkerEnergy::requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & ) const
{
	return true;
}

void
VdWTinkerEnergy::setup_for_scoring_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	ScoreFunction const &, // sfxn,
	ResSingleMinimizationData & resdata
) const
{
	VdWTinkerResidueInfo & info( retrieve_nonconst_vdw_resdata( resdata ) );
	// Assign and rotate moment information
	// TR << "Calling assign_residue from setup_for_scoring_for_residue" << std::endl;
	potential_.assign_residue_amoeba_type( rsd, info );
}

bool
VdWTinkerEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const &  ) const
{
	return true;
}


void
VdWTinkerEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const
{
	/// just compute water locations
	setup_for_scoring_for_residue( rsd, pose, sfxn, min_data );
}


bool
VdWTinkerEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

void
VdWTinkerEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & , // pairdata,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	//TR << "Calculating residue pair energy ext" << std::endl;

	VdWTinkerResidueInfo const & mp1( retrieve_vdw_residue_info( pose, rsd1.seqpos() ) );
	VdWTinkerResidueInfo const & mp2( retrieve_vdw_residue_info( pose, rsd2.seqpos() ) );

	emap[ fa_vdw_tinker ] += potential_.get_res_res_vdw( rsd1, mp1, rsd2, mp2 );

}

bool
VdWTinkerEnergy::use_extended_intrares_energy_interface() const
{
	return true;
}

void
VdWTinkerEnergy::eval_intrares_energy_ext(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & data_cache,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap
) const {
	if ( exclude_DNA_DNA_ && rsd.is_DNA() ) return;

	//TR << "Calculating intraresidue energy ext" << std::endl;

	VdWTinkerResidueInfo const & vdw_info( retrieve_vdw_resdata( data_cache ) );
	emap[ fa_vdw_tinker ] += potential_.get_res_res_vdw( rsd, vdw_info, rsd, vdw_info );
}


/////////////////////////////////
// End minimization specific stuff
/////////////////////////////////

void
VdWTinkerEnergy::evaluate_rotamer_intrares_energies(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< core::PackerEnergy > & energies
) const {
	using core::conformation::RotamerSetCacheableDataType::VDWTINKER_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set.resid() ).is_DNA() ) return;

	VdWTinkerRotamerSetInfo const & vdw_info
		( set.data().get< VdWTinkerRotamerSetInfo >( VDWTINKER_ROTAMER_SET_INFO ) );

	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		//TR << "Calculating rotamer intraresidue energy ext" << std::endl;
		Real const vdwE
			( potential_.get_res_res_vdw( *set.rotamer( ii ), vdw_info.residue_info( ii ),
			*set.rotamer( ii ), vdw_info.residue_info( ii ) ) );

		energies[ ii ] += static_cast< core::PackerEnergy > ( sfxn[ fa_vdw_tinker ] * vdwE );
	}
}

void
VdWTinkerEnergy::evaluate_rotamer_intrares_energy_maps(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & , // sxfn,
	utility::vector1< EnergyMap > & emaps
) const
{
	using core::conformation::RotamerSetCacheableDataType::VDWTINKER_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set.resid() ).is_DNA() ) return;

	// std::cout << "VdWTinker rotamer_intrares_energies: " << set.resid() << std::endl;

	VdWTinkerRotamerSetInfo const & vdw_info
		( set.data().get< VdWTinkerRotamerSetInfo >( VDWTINKER_ROTAMER_SET_INFO ) );

	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		//TR << "Calculating rotamer intraresidue energy maps" << std::endl;
		Real const vdwE
			( potential_.get_res_res_vdw( *set.rotamer( ii ), vdw_info.residue_info( ii ),
			*set.rotamer( ii ), vdw_info.residue_info( ii ) ) );

		(emaps[ ii ])[ fa_vdw_tinker ] += vdwE ;
	}
}

void
VdWTinkerEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & , // sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	using namespace conformation;
	using namespace numeric;
	using core::conformation::RotamerSetCacheableDataType::VDWTINKER_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set1.resid() ).is_DNA() && pose.residue( set2.resid() ).is_DNA() ) return;

	VdWTinkerRotamerSetInfo const & vdw_info1
		( set1.data().get< VdWTinkerRotamerSetInfo >( VDWTINKER_ROTAMER_SET_INFO ) );

	VdWTinkerRotamerSetInfo const & vdw_info2
		( set2.data().get< VdWTinkerRotamerSetInfo >( VDWTINKER_ROTAMER_SET_INFO ) );

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
					
					Real const vdwE(
						potential_.get_res_res_vdw( *set1.rotamer( kk_rot_id ), vdw_info1.residue_info( kk_rot_id ),
							*set2.rotamer( ll_rot_id ), vdw_info2.residue_info( ll_rot_id ) ) );
					
					energy_table( ll_rot_id, kk_rot_id ) +=
					static_cast< core::PackerEnergy >( weights[ fa_vdw_tinker ] *  vdwE );
				}
			}
		}
	}
}

void
VdWTinkerEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & , // sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const {
	using conformation::Residue;
	using core::conformation::RotamerSetCacheableDataType::VDWTINKER_ROTAMER_SET_INFO;

	VdWTinkerRotamerSetInfo const & vdw_set_info
		( set.data().get< VdWTinkerRotamerSetInfo >( VDWTINKER_ROTAMER_SET_INFO ) );

	VdWTinkerResidueInfo const & vdw_rsd_info( retrieve_vdw_residue_info( pose, rsd.seqpos() ) );

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
			
			Real const vdwE(
				potential_.get_res_res_vdw( *set.rotamer( kk_rot_id ), vdw_set_info.residue_info( kk_rot_id ),
					rsd, vdw_rsd_info ) );
			energy_vector[ kk_rot_id ] += static_cast< core::PackerEnergy > (weights[ fa_vdw_tinker ] *  vdwE );
		} // kk - rotamers for residue types
	} // ii - residue types for rotamer set
}

void
VdWTinkerEnergy::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & , // sfxn,
	EnergyMap const & ,
	utility::vector1< EnergyMap > & emaps
) const {

	using conformation::Residue;
	using core::conformation::RotamerSetCacheableDataType::VDWTINKER_ROTAMER_SET_INFO;

	VdWTinkerRotamerSetInfo const & vdw_set_info
		( set.data().get< VdWTinkerRotamerSetInfo >( VDWTINKER_ROTAMER_SET_INFO ) );

	VdWTinkerResidueInfo const & vdw_rsd_info( retrieve_vdw_residue_info( pose, rsd.seqpos() ) );

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
			
			Real const vdwE
			( potential_.get_res_res_vdw( *set.rotamer( kk_rot_id ), vdw_set_info.residue_info( kk_rot_id ),
				rsd, vdw_rsd_info ) );
			(emaps[ kk_rot_id ])[ fa_vdw_tinker ] += vdwE;
		} // kk - rotamers for residue types
	} // ii - residue types for rotamer set
}

/////////////////////////////////////////////////////////////////////////////
///
#ifdef NOTDEF
void
VdWTinkerEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	potential_.eval_atom_derivative( atom_id, weights[ fa_vdw_tinker ], pose, domain_map, exclude_DNA_DNA_, F1, F2 );
}
#endif

/// @brief VdWTinkerEnergy requires no context graphs
void
VdWTinkerEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{
}

/// @brief VdWTinkerEnergy does define intraresidue interactions
bool
VdWTinkerEnergy::defines_intrares_energy( EnergyMap const & /*weights*/ ) const
{
	return true;
}

void
VdWTinkerEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if ( exclude_DNA_DNA_ && rsd.is_DNA() ) return;

	VdWTinkerPoseInfo const & vdw_info
		( static_cast< VdWTinkerPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) ) ); // SHOULD BE FAST!

	emap[ fa_vdw_tinker ] += potential_.get_res_res_vdw( rsd, vdw_info.residue_info( rsd.seqpos() ),
		rsd, vdw_info.residue_info( rsd.seqpos() ) );
}

void
VdWTinkerEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {
	VdWTinkerResidueInfo const & mp( retrieve_vdw_resdata( min_data ) );

	Real const factor( weights[ fa_vdw_tinker] );

	potential_.eval_residue_pair_derivatives( rsd, rsd, mp, mp, pose, factor,
		atom_derivs, atom_derivs );
}


void
VdWTinkerEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & data1,
	ResSingleMinimizationData const & data2,
	ResPairMinimizationData const & ,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	VdWTinkerResidueInfo const & mp1( retrieve_vdw_resdata( data1 ) );
	VdWTinkerResidueInfo const & mp2( retrieve_vdw_resdata( data2 ) );

	potential_.eval_residue_pair_derivatives( rsd1, rsd2, mp1, mp2, pose, weights[ fa_vdw_tinker ],
		r1_atom_derivs, r2_atom_derivs );
}


core::Size
VdWTinkerEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
