// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/electron_density_atomwise/ElecDensAtomwiseEnergy.cc
/// @brief  elec_dens_atomwise scoring method implementation
/// @author Fang-Chieh Chou

// Unit headers
#include <core/scoring/electron_density_atomwise/ElectronDensityAtomwise.hh>
#include <core/scoring/electron_density_atomwise/ElecDensAtomwiseEnergy.hh>
#include <core/scoring/electron_density_atomwise/ElecDensAtomwiseEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>


#include <core/conformation/Atom.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/OneToAllEnergyContainer.hh>
#include <core/id/AtomID.hh>
#include <basic/Tracer.hh>

#include <utility/pointer/memory.hh>

using basic::Error;
using basic::Warning;
namespace core {
namespace scoring {
namespace electron_density_atomwise {

using namespace core;
static basic::Tracer TR( "core.scoring.electron_density_atomwise.ElecDensAtomwiseEnergy" );

/// @details This must return a fresh instance of the ElecDensAtomwiseEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ElecDensAtomwiseEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new ElecDensAtomwiseEnergy );
}

ScoreTypes
ElecDensAtomwiseEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( elec_dens_atomwise );
	return sts;
}

methods::LongRangeEnergyType
ElecDensAtomwiseEnergy::long_range_type() const {
	return methods::elec_dens_atomwise_energy;
}

ElecDensAtomwiseEnergy::ElecDensAtomwiseEnergy() :
	parent( methods::EnergyMethodCreatorOP( new ElecDensAtomwiseEnergyCreator ) ) {
	//Load map
	get_density_map();
}

ElecDensAtomwiseEnergy::~ElecDensAtomwiseEnergy() = default;

/// clone
methods::EnergyMethodOP ElecDensAtomwiseEnergy::clone() const {
	return methods::EnergyMethodOP( new ElecDensAtomwiseEnergy( *this ) );
}


void
ElecDensAtomwiseEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	using namespace methods;

	// Do we have a map?
	if ( !get_density_map().isMapLoaded() ) {
		utility_exit_with_message( "Density scoring function called but no map loaded." );
	}

	// make sure the root of the FoldTree is a virtual atom and is followed by a jump
	// if not, emit warning
	kinematics::Edge const &root_edge( *pose.fold_tree().begin() );
	int virt_res_idx = root_edge.start();
	conformation::Residue const &root_res( pose.residue( virt_res_idx ) );
	pose_is_proper = true;

	if ( root_res.aa() != core::chemical::aa_vrt || root_edge.label() < 0 ) {
		pose_is_proper = false;  // we may be able to recover from this some time but for now just exit
		utility_exit_with_message( "Fold tree is not set properly for density scoring!" );
	}

	// create LR energy container
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		OneToAllEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::OneToAllEnergyContainer >( lrc ) );

		// make sure size or root did not change
		if ( dec->size() != pose.size() || dec->fixed() != virt_res_idx ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		TR.Debug << "Creating new one-to-all energy container (" << pose.size() << ")" << std::endl;
		LREnergyContainerOP new_dec( new OneToAllEnergyContainer (
			virt_res_idx, pose.size(),  elec_dens_atomwise ) );
		energies.set_long_range_container ( lr_type, new_dec );
		get_density_map().is_score_precomputed( false );
	} else {
		//TR << "Checking to see if the density map is up to date with normalization and unweighted score" << std::endl;
		//TR << "Comparing pose.annotated_sequence() (" << pose.annotated_sequence() << ") vs.";
		//TR << " cached (" << pose.data().get< PoseSequence >( pose::datacache::CacheableDataType::POSE_SEQUENCE ).pose_sequence() << ")." << std::endl;
		// calebgeniesse: need to figure out a smarter way to decide when to recompute score
		//                for now, just always recompute
		get_density_map().is_score_precomputed( pose.annotated_sequence() == pose.data().get< PoseSequence >( pose::datacache::CacheableDataType::POSE_SEQUENCE ).pose_sequence() );
	}
	//Pre-calculate the normalization factor and the correlation per
	//atom
	get_density_map().compute_normalization( pose );
	get_density_map().precompute_unweighted_score();

	// if in fact the sequence DID change, now we should set the new seuqence
	pose.data().set( pose::datacache::CacheableDataType::POSE_SEQUENCE, utility::pointer::make_shared< PoseSequence >( pose.annotated_sequence() ) );
}

///////////////////////////////////////////////////////////////////////
///
bool ElecDensAtomwiseEnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const {
	return ( pose.residue( res1 ).aa() == core::chemical::aa_vrt || pose.residue( res2 ).aa() == core::chemical::aa_vrt );
}

///Compute the residue energy
void
ElecDensAtomwiseEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if ( rsd1.aa() != core::chemical::aa_vrt ) {
		if ( rsd2.aa() != core::chemical::aa_vrt ) return;
	} else {
		if ( rsd2.aa() == core::chemical::aa_vrt ) return;
	}

	conformation::Residue const &rsd( rsd1.aa() == core::chemical::aa_vrt ? rsd2 : rsd1 );
	emap[elec_dens_atomwise] = get_density_map().residue_score( rsd );
}

void
ElecDensAtomwiseEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	core::Size const rsd_id = id.rsd();
	core::Size const atm_id = id.atomno();

	// derivatives only defined for (non-VRT) heavyatoms
	if ( pose.residue( rsd_id ).aa() == core::chemical::aa_vrt ) return;

	// if (hydrogen) return
	if ( !pose.residue( rsd_id ).atom_type( atm_id ).is_heavyatom() ) return;

	numeric::xyzVector<core::Real> grad = get_density_map().atom_gradient( pose, rsd_id, atm_id );
	Vector atom_xyz = pose.xyz( id );
	Vector f2 = grad;
	Vector f1( atom_xyz.cross( atom_xyz - f2 ) );
	F1 += weights[ elec_dens_atomwise ] * f1;
	F2 += weights[ elec_dens_atomwise ] * f2;
}

core::Size
ElecDensAtomwiseEnergy::version() const {
	return 1; // Initial versioning
}

} // electron_density_atomwise
} // scoring
} // core

