// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file   devel/facts/FACTSEnergy.cc
// @brief  
// @author Massih Khorvash
// @author: Hahnbeom Park

// Unit headers
#include <core/scoring/facts/FACTSEnergy.hh>
#include <core/scoring/facts/FACTSEnergyCreator.hh>
#include <core/scoring/facts/FACTSPotential.hh>

// Package headers
#include <core/scoring/GenBornPotential.hh>

// Project headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/DenseEnergyContainer.hh>

#include <core/pose/Pose.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>

#include <basic/datacache/CacheableData.hh>
#include <basic/prof.hh>

// Utility headers
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the GenBornEnergy class,
/// never an instance already in use
methods::EnergyMethodOP FACTSEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new FACTSEnergy( options );
}

ScoreTypes FACTSEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( facts_elec );
	sts.push_back( facts_sasa );
	return sts;
}

FACTSEnergy::FACTSEnergy( FACTSEnergy const & src ):
	parent( src ),
	potential_( src.potential_ ),
	exclude_DNA_DNA_( src.exclude_DNA_DNA_ )
{}

FACTSEnergy::FACTSEnergy( EnergyMethodOptions const & options ):
	parent( new FACTSEnergyCreator ),
	potential_( ScoringManager::get_instance()->get_FACTSPotential() ),
	exclude_DNA_DNA_( options.exclude_DNA_DNA() )
{}

/// clone
EnergyMethodOP FACTSEnergy::clone() const
{
	return new FACTSEnergy( *this );
}

bool FACTSEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size,
	Size
) const
{
	return true;
}

///
void
FACTSEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const & residues_repacking,
	utility::vector1< bool > const &
) const
{
	pose.update_residue_neighbors();

	/// build placeholders, compute template born radii, stash in Pose
	potential_.setup_for_packing( pose, residues_repacking );

}

void FACTSEnergy::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	//pack::rotamer_set::RotamerSet & rotamer_set
	conformation::RotamerSetBase & rotamer_set
) const
{
	/// this will stash rotamer born radii in rotamer set
	potential_.get_rotamers_born_radii( pose, rotamer_set );
}

void FACTSEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	/// update born radii for residue that has changed during packing, eg in rotamer trials
	/// need to double-check the current logic on this...
	potential_.update_residue_for_packing( pose, resid );
}

//
void FACTSEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	methods::LongRangeEnergyType const & lr_type( long_range_type() );

	potential_.setup_for_scoring( pose, false );

	// create a container
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		DenseEnergyContainerOP dec( static_cast< DenseEnergyContainer * > ( lrc.get() ) );
		if ( dec->size() != pose.total_residue() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		LREnergyContainerOP new_dec_elec = new DenseEnergyContainer( pose.total_residue(), facts_elec );
		energies.set_long_range_container( lr_type, new_dec_elec );
		LREnergyContainerOP new_dec_sasa = new DenseEnergyContainer( pose.total_residue(), facts_sasa );
		energies.set_long_range_container( lr_type, new_dec_sasa );
	}
}

void FACTSEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	potential_.setup_for_derivatives( pose );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void FACTSEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	FACTSPoseInfo const & facts_info
		( static_cast< FACTSPoseInfo const & >( pose.data().get( pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) ); // SHOULD BE FAST!

	emap[ facts_elec ] += potential_.evaluate_polar_energy( rsd1, facts_info.residue_info( rsd1.seqpos() ),
																													rsd2 );
	if ( rsd1.seqpos() == rsd2.seqpos() ){
		emap[ facts_sasa ] += potential_.evaluate_nonpolar_energy( rsd1, facts_info.residue_info( rsd1.seqpos() ),
																															 rsd2 );
	}
}

void FACTSEnergy::evaluate_rotamer_intrares_energies(
  conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< core::PackerEnergy > & energies
) const
{

	using namespace conformation;
	using namespace numeric;
	using core::conformation::RotamerSetCacheableDataType::FACTS_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set.resid() ).is_DNA() ) return;

	FACTSRotamerSetInfo const & facts_info
		( set.data().get< FACTSRotamerSetInfo >( FACTS_ROTAMER_SET_INFO ) );

	FACTSPoseInfo const & facts_pose_info
		( static_cast< FACTSPoseInfo const & >( pose.data().get( pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) );

	utility::vector1< Real > dBRi1( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dBRi2( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dSAi1( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dSAi2( pose.residue(set.resid()).natoms(), 0.0 );

	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {

		Real const polarE_noncorrect = 
			potential_.evaluate_polar_otf_energy( *set.rotamer( ii ), facts_info.residue_info( ii ),
																						*set.rotamer( ii ), facts_info.residue_info( ii ),
																						dBRi1, dBRi2, false );

		Real const polarE_correct = 0.0;
		//Real const polarE_correct = 
		//		potential_.polar_energy_pack_corrector( pose.residue( set.resid()),
		//																					*set.rotamer( ii ),
		//																					facts_pose_info.residue_info( set.resid() ) );

		Real const polarE = polarE_noncorrect + polarE_correct;
		
		Real const nonpolarE
			( potential_.evaluate_nonpolar_energy( *set.rotamer( ii ), facts_info.residue_info( ii ),
																						 *set.rotamer( ii ) ));

		energies[ ii ] += static_cast< core::PackerEnergy > ( sfxn[ facts_elec ] * polarE 
																													+ sfxn[ facts_sasa ] * nonpolarE);
		
	}
}

void FACTSEnergy::evaluate_rotamer_intrares_energy_maps(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & ,
	utility::vector1< EnergyMap > & emaps
) const
{
	using namespace conformation;
	using namespace numeric;
	using core::conformation::RotamerSetCacheableDataType::FACTS_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set.resid() ).is_DNA() ) return;

	FACTSRotamerSetInfo const & facts_info
		( set.data().get< FACTSRotamerSetInfo >( FACTS_ROTAMER_SET_INFO ) );

	utility::vector1< Real > dBRi1( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dBRi2( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dSAi1( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dSAi2( pose.residue(set.resid()).natoms(), 0.0 );

	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {

		Real const polarE
			( potential_.evaluate_polar_otf_energy( *set.rotamer( ii ), facts_info.residue_info( ii ),
																							*set.rotamer( ii ), facts_info.residue_info( ii ),
																							dBRi1, dBRi2,
																							false ) );
		Real const nonpolarE
			( potential_.evaluate_nonpolar_energy( *set.rotamer( ii ), facts_info.residue_info( ii ),
																						 *set.rotamer( ii ) ));

			(emaps[ ii ])[ facts_elec ] += polarE ;
			(emaps[ ii ])[ facts_sasa ] += nonpolarE ;
	}
}

void FACTSEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	using namespace conformation;
	using namespace numeric;
	using core::conformation::RotamerSetCacheableDataType::FACTS_ROTAMER_SET_INFO;

	//cout << "eval_pairres" << endl;
	if ( exclude_DNA_DNA_ && pose.residue( set1.resid() ).is_DNA() && pose.residue( set2.resid() ).is_DNA() ) return;

	PROF_START( basic::FACTS_ROTAMER_PAIR_ENERGIES );

	FACTSRotamerSetInfo const & facts_info1
		( set1.data().get< FACTSRotamerSetInfo >( FACTS_ROTAMER_SET_INFO ) );

	FACTSRotamerSetInfo const & facts_info2
		( set2.data().get< FACTSRotamerSetInfo >( FACTS_ROTAMER_SET_INFO ) );

	utility::vector1< Real > dBRi1( pose.residue(set1.resid()).natoms(), 0.0 );
	utility::vector1< Real > dBRi2( pose.residue(set2.resid()).natoms(), 0.0 );
	utility::vector1< Real > dSAi1( pose.residue(set1.resid()).natoms(), 0.0 );
	utility::vector1< Real > dSAi2( pose.residue(set2.resid()).natoms(), 0.0 );

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

			if ( ii_coord.distance_squared( jj_coord ) < std::pow(ii_radius+jj_radius+packing_interaction_cutoff(), 2 )) {
				for ( Size kk = 1, kke = set1.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
					Size const kk_rot_id = ii_offset + kk - 1;
					for ( Size ll = 1, lle = set2.get_n_rotamers_for_residue_type( jj ); ll <= lle; ++ll ) {
						Size const ll_rot_id = jj_offset + ll - 1;

						/*
						potential_.evaluate_context_change_for_packing(
 																  pose.residue( set1.resid() ),
																	*set1.rotamer( kk_rot_id ), facts_info1.residue_info( kk_rot_id ),
																	pose.residue( set2.resid() ),
																	*set2.rotamer( ll_rot_id ), facts_info2.residue_info( ll_rot_id ),
																	dBRi1, dBRi2, dSAi1, dSAi2 );
						*/

						Real const polarE
							( potential_.evaluate_polar_otf_energy( *set1.rotamer( kk_rot_id ), facts_info1.residue_info( kk_rot_id ),
																											*set2.rotamer( ll_rot_id ), facts_info2.residue_info( ll_rot_id ),
																											dBRi1, dBRi2,
																											false ));
						Real const nonpolarE
							( potential_.evaluate_nonpolar_energy( *set1.rotamer( kk_rot_id ), facts_info1.residue_info( kk_rot_id ),
																										 *set2.rotamer( ll_rot_id ) ) );

						energy_table( ll_rot_id, kk_rot_id ) += static_cast< core::PackerEnergy >( weights[ facts_elec ] * polarE 
																																											 + weights[ facts_sasa ] * nonpolarE );
						//std::cout << "Pair: " << set1.resid() << " " << set2.resid() << " " << kk_rot_id << " ";
						//std::cout << ll_rot_id << " " << polarE << " " << nonpolarE << std::endl;
					}
				}
			}
		}
	}
	PROF_START( basic::FACTS_ROTAMER_PAIR_ENERGIES );
}

/*
// Below is a lazy version for checking pair term score
// Note that this is also "not exact", since during packing context should be changed, which couldn't be updated on the fly
void FACTSEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	using namespace conformation;
	using namespace numeric;
	using core::conformation::RotamerSetCacheableDataType::FACTS_ROTAMER_SET_INFO;

	//cout << "eval_pairres" << endl;
	if ( exclude_DNA_DNA_ && pose.residue( set1.resid() ).is_DNA() && pose.residue( set2.resid() ).is_DNA() ) return;

	PROF_START( basic::FACTS_ROTAMER_PAIR_ENERGIES );

	FACTSRotamerSetInfo const & facts_info1
		( set1.data().get< FACTSRotamerSetInfo >( FACTS_ROTAMER_SET_INFO ) );

	FACTSRotamerSetInfo const & facts_info2
		( set2.data().get< FACTSRotamerSetInfo >( FACTS_ROTAMER_SET_INFO ) );

	FACTSPoseInfo const & facts_info_pose( static_cast< FACTSPoseInfo const & >
																	 (pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )));

	core::pose::Pose pose_tmp( pose );

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

			if ( ii_coord.distance_squared( jj_coord ) < std::pow(ii_radius+jj_radius+packing_interaction_cutoff(), 2 )) {
				for ( Size kk = 1, kke = set1.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
					Size const kk_rot_id = ii_offset + kk - 1;
					
					pose_tmp.replace_residue( set1.resid(), *set1.rotamer( kk_rot_id), false );

					for ( Size ll = 1, lle = set2.get_n_rotamers_for_residue_type( jj ); ll <= lle; ++ll ) {
						Size const ll_rot_id = jj_offset + ll - 1;

						pose_tmp.replace_residue( set2.resid(), *set2.rotamer( ll_rot_id), false );

						potential_.setup_for_scoring( pose_tmp );
						FACTSResidueInfo const & facts1
							( pose_tmp.data().get< FACTSPoseInfo >( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ).residue_info(set1.resid()) );
						FACTSResidueInfo const & facts2
							( pose_tmp.data().get< FACTSPoseInfo >( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ).residue_info(set2.resid()) );

						//FACTSResidueInfo facts1 = facts_info1.residue_info( kk_rot_id );
						//facts1.initialize( pose_tmp.residue(set1.resid()) );
						//potential_.get_single_rotamer_born_radii( pose_tmp.residue(set1.resid()), pose_tmp, facts_info_pose, facts1 );

						//FACTSResidueInfo facts2 = facts_info2.residue_info( ll_rot_id );
						//facts2.initialize( pose_tmp.residue(set2.resid()) );
						//potential_.get_single_rotamer_born_radii( pose_tmp.residue(set2.resid()), pose_tmp, facts_info_pose, facts2 );

						Real const polarE
							( potential_.evaluate_polar_otf_energy( pose.residue( set1.resid() ), *set1.rotamer( kk_rot_id ), facts1,
																											pose.residue( set2.resid() ), *set2.rotamer( ll_rot_id ), facts2, false ));
						Real const nonpolarE
							( potential_.evaluate_nonpolar_energy( *set1.rotamer( kk_rot_id ), facts1, 
																										 *set2.rotamer( ll_rot_id ) ) );

						energy_table( ll_rot_id, kk_rot_id ) += static_cast< core::PackerEnergy >( weights[ facts_elec ] * polarE 
																																											 + weights[ facts_sasa ] * nonpolarE );
						std::cout << "Pair: " << set1.resid() << " " << set2.resid() << " " << kk_rot_id << " ";
						std::cout << ll_rot_id << " " << polarE << " " << nonpolarE << std::endl;
					}
				}
			}
		}
	}
	PROF_START( basic::FACTS_ROTAMER_PAIR_ENERGIES );
}
*/

void FACTSEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	//cout << "eval_background" << endl;

	PROF_START( basic::FACTS_ROTAMER_BACKGROUND_ENERGIES );

	using conformation::Residue;
	using core::conformation::RotamerSetCacheableDataType::FACTS_ROTAMER_SET_INFO;

	FACTSResidueInfo const & facts_rsd_info
		( pose.data().get< FACTSPoseInfo >( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ).residue_info(rsd.seqpos()));

	FACTSRotamerSetInfo const & facts_set_info
		( set.data().get< FACTSRotamerSetInfo >( FACTS_ROTAMER_SET_INFO ) );

	utility::vector1< Real > dBRi1( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dBRi2( rsd.natoms(), 0.0 );
	utility::vector1< Real > dSAi1( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dSAi2( rsd.natoms(), 0.0 );

	for ( Size ii = 1; ii <= set.get_n_residue_types(); ++ii ) {
		Size const ii_offset = set.get_residue_type_begin( ii );
		Residue const & ii_example_rotamer( *set.rotamer( ii_offset ));

		Vector const & ii_coord( ii_example_rotamer.nbr_atom_xyz() );
		Real const ii_radius( ii_example_rotamer.nbr_radius() );

		if ( exclude_DNA_DNA_ && ii_example_rotamer.is_DNA() && rsd.is_DNA() ) continue;

		Vector const & jj_coord( rsd.nbr_atom_xyz() );
		Real const jj_radius( rsd.nbr_radius() );

		if ( ii_coord.distance_squared( jj_coord ) < std::pow(ii_radius+jj_radius+packing_interaction_cutoff(), 2 )) {
			for ( Size kk = 1, kke = set.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
				Size const kk_rot_id = ii_offset + kk - 1;

				/*
				potential_.evaluate_context_change_for_packing(
 																  pose.residue( set.resid() ),
																	*set.rotamer( kk_rot_id ), facts_set_info.residue_info( kk_rot_id ),
																	rsd,
																	rsd, facts_rsd_info,
																	dBRi1, dBRi2, dSAi1, dSAi2 );
				*/

				Real const polarE
					( potential_.evaluate_polar_otf_energy( *set.rotamer( kk_rot_id ), facts_set_info.residue_info( kk_rot_id ),
																									rsd, facts_rsd_info, 
																									dBRi1, dBRi2, false ) );
				Real const nonpolarE
					( potential_.evaluate_nonpolar_energy( *set.rotamer( kk_rot_id ), facts_set_info.residue_info( kk_rot_id ),
																								 rsd) );
				energy_vector[ kk_rot_id ] += static_cast< core::PackerEnergy > (weights[ facts_elec ] * polarE
																																				 + weights[ facts_sasa ] * nonpolarE );
				//std::cout << "Background: " << ii << " " << kk << " " << kk_rot_id << " " << polarE << " " << nonpolarE << std::endl;
			} // kk - rotamers for residue types
		} // nbrs
	} // ii - residue types for rotamer set
	PROF_STOP( basic::FACTS_ROTAMER_BACKGROUND_ENERGIES );
}

void FACTSEnergy::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & ,
	utility::vector1< EnergyMap > & emaps
) const
{
	using conformation::Residue;
	using core::conformation::RotamerSetCacheableDataType::FACTS_ROTAMER_SET_INFO;

	FACTSResidueInfo const & facts_rsd_info
		( pose.data().get< FACTSPoseInfo >( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ).residue_info(rsd.seqpos()));

	FACTSRotamerSetInfo const & facts_set_info
		( set.data().get< FACTSRotamerSetInfo >( FACTS_ROTAMER_SET_INFO ) );

	utility::vector1< Real > dBRi1( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dBRi2( rsd.natoms(), 0.0 );
	utility::vector1< Real > dSAi1( pose.residue(set.resid()).natoms(), 0.0 );
	utility::vector1< Real > dSAi2( rsd.natoms(), 0.0 );

	for ( Size ii = 1; ii <= set.get_n_residue_types(); ++ii ) {
		Size const ii_offset = set.get_residue_type_begin( ii );
		Residue const & ii_example_rotamer( *set.rotamer( ii_offset ));

		Vector const & ii_coord( ii_example_rotamer.nbr_atom_xyz() );
		Real const ii_radius( ii_example_rotamer.nbr_radius() );

		if ( exclude_DNA_DNA_ && ii_example_rotamer.is_DNA() && rsd.is_DNA() ) continue;

		Vector const & jj_coord( rsd.nbr_atom_xyz() );
		Real const jj_radius( rsd.nbr_radius() );

		if ( ii_coord.distance_squared( jj_coord ) < std::pow(ii_radius+jj_radius+packing_interaction_cutoff(), 2 )) {
			for ( Size kk = 1, kke = set.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
				Size const kk_rot_id = ii_offset + kk - 1;

				Real const polarE
					( potential_.evaluate_polar_otf_energy( *set.rotamer( kk_rot_id ), facts_set_info.residue_info( kk_rot_id ),
																									rsd, facts_rsd_info, 
																									dBRi1, dBRi2,	false) );
				Real const nonpolarE
					( potential_.evaluate_nonpolar_energy( *set.rotamer( kk_rot_id ), facts_set_info.residue_info( kk_rot_id ),
																								 rsd ) );
				(emaps[ kk_rot_id ])[ facts_elec ] += polarE;
				(emaps[ kk_rot_id ])[ facts_sasa ] += nonpolarE;
			} // kk - rotamers for residue types
		} // nbrs
	} // ii - residue types for rotamer set
}


/////////////////////////////////////////////////////////////////////////////
///
void FACTSEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	potential_.eval_atom_polar_derivative( atom_id, weights[ facts_elec ], pose, domain_map, exclude_DNA_DNA_, F1, F2 );
	potential_.eval_atom_nonpolar_derivative( atom_id, weights[ facts_sasa ], pose, domain_map, exclude_DNA_DNA_, F1, F2 );
}

/// @brief FACTSEnergy requires no context graphs
void FACTSEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{
}

/// @brief FACTSEnergy does define intraresidue interactions
bool FACTSEnergy::defines_intrares_energy( EnergyMap const & /*weights*/ ) const
{
	return true;
}

void FACTSEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( exclude_DNA_DNA_ && rsd.is_DNA() ) return;

	FACTSResidueInfo const & facts
		( pose.data().get< FACTSPoseInfo >( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ).residue_info( rsd.seqpos()));
	
	emap[ facts_elec ] += potential_.evaluate_polar_energy( rsd, facts, rsd );
	emap[ facts_sasa ] += potential_.evaluate_nonpolar_energy( rsd, facts, rsd );

}

Size FACTSEnergy::version() const
{
	return 1;
}

}
}
}
