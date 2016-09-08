// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/GenBornEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/GenBornEnergy.hh>
#include <core/scoring/methods/GenBornEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/GenBornPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ContextGraphTypes.hh>

#include <core/scoring/DenseEnergyContainer.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>
#include <basic/prof.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>

#include <utility/vector1.hh>


// Utility headers


// C++
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
/**

This is a reimplementation of Jim Havranek's original rosetta++ Gen Born code.
source files: rosetta++/gb_elec*

**/
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the GenBornEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
GenBornEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new GenBornEnergy( options ) );
}

ScoreTypes
GenBornEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( gb_elec );
	return sts;
}


GenBornEnergy::GenBornEnergy( GenBornEnergy const & src ):
	parent( src ),
	potential_( src.potential_ ),
	exclude_DNA_DNA_( src.exclude_DNA_DNA_ )
{}


GenBornEnergy::GenBornEnergy( EnergyMethodOptions const & options ):
	parent( methods::EnergyMethodCreatorOP( new GenBornEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_GenBornPotential() ),
	exclude_DNA_DNA_( options.exclude_DNA_DNA() )
{}


/// clone
EnergyMethodOP
GenBornEnergy::clone() const
{
	return EnergyMethodOP( new GenBornEnergy( *this ) );
}

bool
GenBornEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size,
	Size
) const
{
	return true;
}


void
GenBornEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const & residues_repacking,
	utility::vector1< bool > const &
) const
{
	pose.update_residue_neighbors();

	/// build placeholders, compute template born radii, stash in Pose
	potential_.setup_for_packing( pose, residues_repacking );

}

void
GenBornEnergy::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & rotamer_set
) const
{
	/// this will stash rotamer born radii in rotamer set
	potential_.get_rotamers_born_radii( pose, rotamer_set );

}


//  void
//  update_residue_for_packing(
//   pose::Pose &,
//   Size /*resid*/ ) const
//  {}
void
GenBornEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	/// update born radii for residue that has changed during packing, eg in rotamer trials
	/// need to double-check the current logic on this...
	potential_.update_residue_for_packing( pose, resid );
}


void
GenBornEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	LongRangeEnergyType const & lr_type( long_range_type() );

	potential_.get_all_born_radii( pose );

	// create a container
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;

	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		LREnergyContainerOP new_dec( new DenseEnergyContainer( pose.size(), gb_elec ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}


void
GenBornEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	potential_.get_all_born_radii( pose );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
GenBornEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	GenBornPoseInfo const & gb_info
		( static_cast< GenBornPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ) ) ); // SHOULD BE FAST!

	emap[ gb_elec ] += potential_.get_res_res_elecE( rsd1, gb_info.residue_info( rsd1.seqpos() ),
		rsd2, gb_info.residue_info( rsd2.seqpos() ) );
}

void
GenBornEnergy::evaluate_rotamer_intrares_energies(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	utility::vector1< core::PackerEnergy > & energies
) const
{
	using namespace conformation;
	using namespace numeric;
	using core::conformation::RotamerSetCacheableDataType::GEN_BORN_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set.resid() ).is_DNA() ) return;


	GenBornRotamerSetInfo const & gb_info
		( set.data().get< GenBornRotamerSetInfo >( GEN_BORN_ROTAMER_SET_INFO ) );

	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {

		Real const elecE
			( potential_.get_res_res_elecE( *set.rotamer( ii ), gb_info.residue_info( ii ),
			*set.rotamer( ii ), gb_info.residue_info( ii ) ) );

		energies[ ii ] += static_cast< core::PackerEnergy > ( sfxn[ gb_elec ] * elecE );
	}
}

void
GenBornEnergy::evaluate_rotamer_intrares_energy_maps(
	conformation::RotamerSetBase const & set,
	pose::Pose const & pose,
	ScoreFunction const & ,
	utility::vector1< EnergyMap > & emaps
) const
{
	using namespace conformation;
	using namespace numeric;
	using core::conformation::RotamerSetCacheableDataType::GEN_BORN_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set.resid() ).is_DNA() ) return;

	// std::cout << "GenBorn rotamer_intrares_energies: " << set.resid() << std::endl;

	GenBornRotamerSetInfo const & gb_info
		( set.data().get< GenBornRotamerSetInfo >( GEN_BORN_ROTAMER_SET_INFO ) );

	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {

		Real const elecE
			( potential_.get_res_res_elecE( *set.rotamer( ii ), gb_info.residue_info( ii ),
			*set.rotamer( ii ), gb_info.residue_info( ii ) ) );

		(emaps[ ii ])[ gb_elec ] += elecE ;
	}
}

void
GenBornEnergy::evaluate_rotamer_pair_energies(
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
	using core::conformation::RotamerSetCacheableDataType::GEN_BORN_ROTAMER_SET_INFO;

	if ( exclude_DNA_DNA_ && pose.residue( set1.resid() ).is_DNA() && pose.residue( set2.resid() ).is_DNA() ) return;

	PROF_START( basic::GEN_BORN_ROTAMER_PAIR_ENERGIES );

	GenBornRotamerSetInfo const & gb_info1
		( set1.data().get< GenBornRotamerSetInfo >( GEN_BORN_ROTAMER_SET_INFO ) );

	GenBornRotamerSetInfo const & gb_info2
		( set2.data().get< GenBornRotamerSetInfo >( GEN_BORN_ROTAMER_SET_INFO ) );


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

					Real const elecE(
						potential_.get_res_res_elecE( *set1.rotamer( kk_rot_id ), gb_info1.residue_info( kk_rot_id ),
						*set2.rotamer( ll_rot_id ), gb_info2.residue_info( ll_rot_id ) ) );

					energy_table( ll_rot_id, kk_rot_id ) += static_cast< core::PackerEnergy >( weights[ gb_elec ] *  elecE );
				}
			}
		}
	}
	PROF_START( basic::GEN_BORN_ROTAMER_PAIR_ENERGIES );
}

void
GenBornEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	PROF_START( basic::GEN_BORN_ROTAMER_BACKGROUND_ENERGIES );

	using conformation::Residue;
	using core::conformation::RotamerSetCacheableDataType::GEN_BORN_ROTAMER_SET_INFO;
	//using core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO;

	GenBornResidueInfo const & gb_rsd_info
		( pose.data().get< GenBornPoseInfo >( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ).residue_info(rsd.seqpos()));
	//GenBornResidueInfo const & gb_rsd_info( retrieve_gen_born_pose_info( pose ).residue_info( rsd.seqpos ) );

	GenBornRotamerSetInfo const & gb_set_info
		( set.data().get< GenBornRotamerSetInfo >( GEN_BORN_ROTAMER_SET_INFO ) );
	//GenBornRotamersInfo const & gb_set_info( retrieve_gen_born_rotamers_info( set ) );

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

			Real const elecE(
				potential_.get_res_res_elecE( *set.rotamer( kk_rot_id ), gb_set_info.residue_info( kk_rot_id ),
				rsd, gb_rsd_info ) );
			energy_vector[ kk_rot_id ] += static_cast< core::PackerEnergy > (weights[ gb_elec ] *  elecE );
		} // kk - rotamers for residue types
	} // ii - residue types for rotamer set
	PROF_STOP( basic::GEN_BORN_ROTAMER_BACKGROUND_ENERGIES );
}

void
GenBornEnergy::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &, // sfxn,
	EnergyMap const & ,
	utility::vector1< EnergyMap > & emaps
) const
{
	// std::cout << "GenBornEnergy: rotamer background: " << set.resid() << ' ' << rsd.seqpos() << std::endl;

	using conformation::Residue;
	using core::conformation::RotamerSetCacheableDataType::GEN_BORN_ROTAMER_SET_INFO;

	GenBornResidueInfo const & gb_rsd_info
		( pose.data().get< GenBornPoseInfo >( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ).residue_info(rsd.seqpos()));
	GenBornRotamerSetInfo const & gb_set_info
		( set.data().get< GenBornRotamerSetInfo >( GEN_BORN_ROTAMER_SET_INFO ) );

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

			Real const elecE
				( potential_.get_res_res_elecE( *set.rotamer( kk_rot_id ), gb_set_info.residue_info( kk_rot_id ),
				rsd, gb_rsd_info ) );
			(emaps[ kk_rot_id ])[ gb_elec ] += elecE;
		} // kk - rotamers for residue types
	} // ii - residue types for rotamer set
}

/////////////////////////////////////////////////////////////////////////////
///
void
GenBornEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	potential_.eval_atom_derivative( atom_id, weights[ gb_elec ], pose, domain_map, exclude_DNA_DNA_, F1, F2 );
}

/// @brief GenBornEnergy requires no context graphs
void
GenBornEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{
}

/// @brief GenBornEnergy does define intraresidue interactions
bool
GenBornEnergy::defines_intrares_energy( EnergyMap const & /*weights*/ ) const
{
	return true;
}

void
GenBornEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( exclude_DNA_DNA_ && rsd.is_DNA() ) return;

	GenBornResidueInfo const & gb
		( pose.data().get< GenBornPoseInfo >( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ).residue_info(rsd.seqpos()));

	emap[ gb_elec ] += potential_.get_res_res_elecE( rsd, gb, rsd, gb );
}
core::Size
GenBornEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
