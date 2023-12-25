// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/scoring/ElectricfieldLipidlayer.cc
/// @brief Implicit Lipid Membrane Model electrostatic energy due to the field created by lipid layers(one-body)
/// @author  Rituparna Samanta (rsamant2@jhu.edu)

#ifndef INCLUDED_protocols_membrane_scoring_ElectricfieldLipidlayer_cc
#define INCLUDED_protocols_membrane_scoring_ElectricfieldLipidlayer_cc

// Unit headers
#include <protocols/membrane/scoring/ElectricfieldLipidlayer.hh>
#include <protocols/membrane/scoring/ElectricfieldLipidlayerCreator.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>

// Project Headers
#include <protocols/membrane/scoring/MEnvElectroAtomParams.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/Atom.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// C++ Headers

static basic::Tracer TR( "protocols.membrane.scoring.ElectricfieldLipidlayer" );

namespace protocols {
namespace membrane {
namespace scoring {

/// @brief Return a fresh instance of the energy method
core::scoring::methods::EnergyMethodOP
ElectricfieldLipidlayerCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return core::scoring::methods::EnergyMethodOP( new ElectricfieldLipidlayer );
}

/// @brief Return Score Types Required for Method
core::scoring::ScoreTypes
ElectricfieldLipidlayerCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::f_elec_lipidlayer );
	return sts;
}

/// @brief Construct Energy Method from Etable
ElectricfieldLipidlayer::ElectricfieldLipidlayer() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new ElectricfieldLipidlayerCreator ) )
{}

/// @brief Clone Energy Method
core::scoring::methods::EnergyMethodOP
ElectricfieldLipidlayer::clone() const {
	return core::scoring::methods::EnergyMethodOP( new ElectricfieldLipidlayer );
}


/// @brief Compute Per-Residue Energies
void
ElectricfieldLipidlayer::residue_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	core::scoring::EnergyMap & emap
) const {

	using namespace core;
	using namespace core::scoring;

	if ( rsd.name3() == "MEM" || !rsd.is_protein() ) return;
	for ( core::Size i = 1, i_end = rsd.natoms(); i <= i_end; ++i ) {
		MEnvElectroAtomParamsCOP menv_params = get_menv_params_for_residue( pose, rsd, i );
		emap[ core::scoring::f_elec_lipidlayer ] += eval_felec_lipidlayer( *menv_params );
		// TR << "res: " << rsd.name3() << " atom: " << rsd.atom_type(i).name() << " charge: " << menv_params->charge() << " field: " << menv_params->elec_field() << " score: " << eval_felec_lipidlayer( *menv_params )<<std::endl;

	}


}

/// @brief Setup for Computing Derivatives
void
ElectricfieldLipidlayer::setup_for_derivatives(
	core::pose::Pose &,
	core::scoring::ScoreFunction const & scfxn
) const {
	f_elec_lipidlayer_weight_ = scfxn.weights()[ core::scoring::f_elec_lipidlayer ];
}

/// @brief Evaluate Per-Atom Derivatives
void
ElectricfieldLipidlayer::eval_atom_derivative(
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
	core::Size const i( atom_id.rsd() );
	core::Size const m( atom_id.atomno() );

	// Get the actual residue and atom
	conformation::Residue const & rsd1( pose.residue( i ) );
	if ( m > rsd1.natoms() ) return;
	Vector const atom_i( rsd1.xyz( m ) );

	// Get the relevant atom type parameters
	MEnvElectroAtomParamsCOP p( get_menv_params_for_residue( pose, rsd1, m ) );

	// Set count pair weight
	Real cp_weight = 1.0;

	// Initialize f1, f2
	// Vector f1( 0.0 ), f2( 0.0 );

	F1 += f_elec_lipidlayer_weight_ * p->charge() * p->elec_field_gradient() * cp_weight * p->f1();
	F2 += f_elec_lipidlayer_weight_ * p->charge() * p->elec_field_gradient() * cp_weight * p->f2();
}

/// @brief F_elec_lipidlayer is context independent; doesnt depend on neighbouring entities.
void
ElectricfieldLipidlayer::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @brief Setup Method for initial scoring
void
ElectricfieldLipidlayer::setup_for_scoring(
	core::pose::Pose &,
	core::scoring::ScoreFunction const &
) const {}

/// @brief Evaluate Per-Atom Energy term
core::Real
ElectricfieldLipidlayer::eval_felec_lipidlayer(
	MEnvElectroAtomParams const & p
) const {

	using namespace core;

	Real score( 0.0 );
	score = ( p.charge() ) * ( p.elec_field() ) *(1.00/0.043);
	//1 Kcal/mole = 0.043 eV/molecule
	return score;
}

/// @brief Evaluate Per-Atom Electrostatic term
MEnvElectroAtomParamsCOP
ElectricfieldLipidlayer::get_menv_params_for_residue(
	core::pose::Pose const & pose,
	core::conformation::Residue const & rsd,
	core::Size atomno
) const {

	using namespace core::pose;
	using namespace core::conformation::membrane;

	// // grab implicit lipid information
	ImplicitLipidInfoOP implicit_lipids( pose.conformation().membrane_info()->implicit_lipids() );

	core::conformation::Conformation const & conf( pose.conformation() );
	core::conformation::membrane::MembraneGeometryCOP mp_geometry( conf.membrane_info()->membrane_geometry() );

	if ( !conf.membrane_info()->use_franklin() ) {
		TR.Warning << "Score term optimized with franklin23 transition function" << std::endl;
		TR.Warning << "Recommended to run with -mp:restore_lazaridis_imm_behavior false if using franklin2019.wts or franklin2023.wts" << std::endl;
	}

	core::Real const charge( rsd.atomic_charge( atomno ) );
	core::Real potential_offset( 0.0 ); //adding an offset to have the potential=0 at z->infinity
	potential_offset = implicit_lipids->f_elec_field( 30.0 );
	core::Real elec_field = mp_geometry->fa_elec_lipid( conf, rsd.seqpos(), atomno ) - potential_offset;
	core::Real elec_field_gradient = mp_geometry->fa_elec_lipid_deriv( conf, rsd.seqpos(), atomno );
	core::Vector f1 = mp_geometry->f_transition_f1( conf, rsd.seqpos(), atomno );
	core::Vector f2 = mp_geometry->f_transition_f2( conf, rsd.seqpos(), atomno );


	MEnvElectroAtomParamsCOP melectro_parameter = MEnvElectroAtomParamsCOP( utility::pointer::make_shared<MEnvElectroAtomParams>( rsd.atom_type(atomno).name(), charge, elec_field, elec_field_gradient, f1, f2 ) );
	return melectro_parameter;

}

/// @brief Versioning
core::Size
ElectricfieldLipidlayer::version() const {
	return 1;
}

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaWaterToBilayerEnergy_cc
