// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/scoring/FaWaterToBilayerEnergy.cc
/// @brief Implicit Lipid Membrane Model Water-to-Bilayer transfer energy (one-body)
/// @author  Rebecca Alford (ralford3@jhu.edu)

#ifndef INCLUDED_protocols_membrane_scoring_FaWaterToBilayerEnergy_cc
#define INCLUDED_protocols_membrane_scoring_FaWaterToBilayerEnergy_cc

// Unit headers
#include <protocols/membrane/scoring/FaWaterToBilayerEnergy.hh>
#include <protocols/membrane/scoring/FaWaterToBilayerEnergyCreator.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project Headers
#include <protocols/membrane/scoring/MEnvAtomParams.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/Atom.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <utility>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.scoring.FaWaterToBilayerEnergy" );

namespace protocols {
namespace membrane {
namespace scoring {

/// @brief Return a fresh instance of the energy method
core::scoring::methods::EnergyMethodOP
FaWaterToBilayerEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return core::scoring::methods::EnergyMethodOP( new FaWaterToBilayerEnergy() );
}

/// @brief Return Score Types Required for Method
core::scoring::ScoreTypes
FaWaterToBilayerEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::fa_water_to_bilayer );
	return sts;
}

/// @brief Construct Energy Method from Etable
FaWaterToBilayerEnergy::FaWaterToBilayerEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new FaWaterToBilayerEnergyCreator ) ),
	memb_lk_dgrefce_(),
	water_lk_dgrefce_(),
	atypes_list_()
{
	std::string dbfile( "membrane/memb_fa_params_2019.txt" );
	using namespace basic;
	using namespace core;

#ifdef USEMPI
	utility_exit_with_message("fa_water_to_bilayer is not yet threadsafe for MPI mode!");
#endif

#ifdef MULTI_THREADED
	utility_exit_with_message("fa_water_to_bilayer is not yet threadsafe for MPI mode!");
#endif

	utility::io::izstream infile;
	TR << "Reading fa_water_to_bilayer parameters from the database" << std::endl;
	database::open( infile, dbfile );
	if ( !infile.good() ) {
		utility_exit_with_message( "Unable to open database file containing reference energies" );

	}

	std::string line;
	getline( infile, line );
	while ( getline( infile, line ) ) {
		utility::trim(line, "\t\n" );
		if ( line.empty() > 0 ) continue;
		if ( line.find("#",0) == 0 ) continue;

		std::istringstream l( line );
		std::string atype, memb_lk_dgrefce, water_lk_dgrefce;
		l >> atype >> water_lk_dgrefce >> memb_lk_dgrefce;

		if ( l.fail() ) {
			utility_exit_with_message( "bad line: " + line );
		}

		// Add data to maps
		atypes_list_.push_back( atype );
		memb_lk_dgrefce_.push_back( utility::from_string( memb_lk_dgrefce, core::Real(0.0) ) );
		water_lk_dgrefce_.push_back( utility::from_string( water_lk_dgrefce, core::Real(0.0) ) );
	}
}

/// @brief Clone Energy Method
core::scoring::methods::EnergyMethodOP
FaWaterToBilayerEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new FaWaterToBilayerEnergy( *this ) );
}


/// @brief Compute Per-Residue Energies
void
FaWaterToBilayerEnergy::residue_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	core::scoring::EnergyMap & emap
) const {

	using namespace core;
	using namespace core::scoring;

	if ( rsd.name3() == "MEM" || !rsd.is_protein() ) return;
	for ( Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i ) {
		MEnvAtomParamsCOP menv_params = get_menv_params_for_residue( pose, rsd, i );
		emap[ core::scoring::fa_water_to_bilayer ] += eval_fa_wtbe( *menv_params );
	}
}

/// @brief Fianlzie Total Per-Residue Energies
void
FaWaterToBilayerEnergy::finalize_total_energy(
	core::pose::Pose &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {
	emap[ core::scoring::fa_water_to_bilayer ] += 0.0;
}

/// @brief Setup for Computing Derivatives
void
FaWaterToBilayerEnergy::setup_for_derivatives(
	core::pose::Pose &,
	core::scoring::ScoreFunction const & scfxn
) const {
	fa_wtbe_weight_ = scfxn.weights()[ core::scoring::fa_water_to_bilayer ];
}

/// @brief Evaluate Per-Atom Derivatives
void
FaWaterToBilayerEnergy::eval_atom_derivative(
	core::id::AtomID const & atom_id,
	core::pose::Pose const & pose,
	core::kinematics::DomainMap const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const &,
	core::Vector & F1,
	core::Vector & F2
) const {

	using namespace core;
	using namespace core::scoring;

	// Grab residue and atom numbers
	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );

	// Get the actual residue and atom
	conformation::Residue const & rsd1( pose.residue( i ) );
	if ( m > rsd1.nheavyatoms() ) return;
	Vector const heavy_atom_i( rsd1.xyz( m ) );

	// Get the relevant atom type parameters
	MEnvAtomParamsCOP p( get_menv_params_for_residue( pose, rsd1, m ) );

	// Set count pair weight
	Real cp_weight = 1.0;

	// Initialize f1, f2
	Vector f1( 0.0 ), f2( 0.0 );

	Real const deriv = p->dGfreeB() - p->dGfreeW();
	Real dE_dZ_over_r = fa_wtbe_weight_ * deriv * p->hydration_deriv();

	Vector const d_ij = p->memb_coord() - heavy_atom_i;
	Real const d_ij_norm = d_ij.length();
	if ( d_ij_norm == Real(0.0) ) return;

	Real const invd = 1.0 / d_ij_norm;
	f2 = d_ij * invd;
	f1 = p->memb_coord().cross(heavy_atom_i);
	f1 *= invd;

	if ( dE_dZ_over_r != 0.0 ) {
		F1 += dE_dZ_over_r * cp_weight * f1;
		F2 += dE_dZ_over_r * cp_weight * f2;
	}
}

/// @brief Fa_MbenvEnergy is context independent
void
FaWaterToBilayerEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @brief Setup Method for initial scoring
void
FaWaterToBilayerEnergy::setup_for_scoring(
	core::pose::Pose &,
	core::scoring::ScoreFunction const &
) const {}

/// @brief Evaluate Per-Atom Env term
core::Real
FaWaterToBilayerEnergy::eval_fa_wtbe(
	MEnvAtomParams const & p
) const {

	using namespace core;

	Real score( 0.0 );
	score = ( 1 - p.hydration() ) * ( p.dGfreeB() - p.dGfreeW() );
	return score;
}

MEnvAtomParamsCOP
FaWaterToBilayerEnergy::get_menv_params_for_residue(
	core::pose::Pose const & pose,
	core::conformation::Residue const & rsd,
	core::Size atomno
) const {

	using namespace core::pose;
	using namespace core::conformation::membrane;

	// grab implicit lipid information
	ImplicitLipidInfoOP implicit_lipids( pose.conformation().membrane_info()->implicit_lipids() );

	// Get the xyz coordinate
	core::Vector const xyz( rsd.xyz( atomno ) );
	core::Vector const center( pose.conformation().membrane_info()->membrane_center( pose.conformation() ) );
	core::Vector const normal( pose.conformation().membrane_info()->membrane_normal( pose.conformation() ) );

	// Calculate the spatially dependent hydration value for the current conformation
	core::Real zcoord = pose.conformation().membrane_info()->atom_z_position( pose.conformation(), rsd.seqpos(), atomno );
	core::Vector new_xyz( xyz.x(), xyz.y(), zcoord );

	// Calculate projection coordinate
	core::Vector proj_i = center + zcoord * normal;
	core::Vector i_ip = proj_i - xyz;
	core::Vector memb_coord( center - i_ip );

	core::Real hydration = implicit_lipids->f_hydration( new_xyz );
	core::Real hydration_deriv = implicit_lipids->f_hydration_gradient( new_xyz );

	// Get the water-to-bilayer free energy components
	core::Size index( get_atype_index( rsd.atom_type( atomno ).name() ) );
	core::Real dGfreeW = water_lk_dgrefce_[ index ];
	core::Real dGfreeB = memb_lk_dgrefce_[ index ];

	MEnvAtomParamsCOP menv_parameter = MEnvAtomParamsCOP( new MEnvAtomParams( rsd.atom_type(atomno).name(), dGfreeW, dGfreeB, hydration, hydration_deriv, memb_coord ) );
	return menv_parameter;

}

core::Size
FaWaterToBilayerEnergy::get_atype_index( std::string atype_name ) const {
	core::Size index(0);
	for ( core::Size ii = 1; ii <= atypes_list_.size(); ++ii ) {
		if ( atypes_list_[ii] == atype_name ) index = ii;
	}
	debug_assert( index != 0 );
	return index;
}

/// @brief Versioning
core::Size
FaWaterToBilayerEnergy::version() const {
	return 2;
}

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaWaterToBilayerEnergy_cc
