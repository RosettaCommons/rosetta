// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/pHEnergy.cc
/// @brief  Energy due to ionization state of the residue at a particular pH (see pHEnergy.hh for a detailed description)
/// @author krishna

// Standard libraries
#include <cmath>

// Unit headers
#include <core/scoring/methods/pHEnergy.hh>
#include <core/scoring/methods/pHEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

#include <basic/options/option.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>


// Option Key includes
#include <basic/options/keys/pH.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

// AMW: oh no this is awful don't let this continue
core::Real pHEnergy::pH_;     //Reference to static variable pH_

/// @details This must return a fresh instance of the pHEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
pHEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new pHEnergy );
}

ScoreTypes
pHEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( e_pH );
	return sts;
}

//ctor
pHEnergy::pHEnergy() :
	parent( methods::EnergyMethodCreatorOP( new pHEnergyCreator ) )
{
	using namespace basic::options;
	pH_ = option[ OptionKeys::pH::value_pH ]();
}


void
pHEnergy::set_pH ( core::Real new_pH_ )
{
	pH_ = new_pH_;
}

/// clone
EnergyMethodOP
pHEnergy::clone() const
{
	return EnergyMethodOP( new pHEnergy );
}


/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
////////////////////////////////////////////////////////////////////////////

// pH score
void
pHEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const {
	using namespace core::chemical;
	core::Real ipKa;
	switch ( rsd.type().aa() ) {
	
	case na_rad:
		ipKa = 6.1;
		break;

	case aa_asp:
		ipKa = 4.0;
		break;

	case aa_glu:
		ipKa = 4.4;
		break;

	case aa_his:
		ipKa = 6.3;
		break;

	case chemical::aa_lys:
		ipKa = 10.4;
		break;

	case chemical::aa_tyr:
		ipKa = 10.0;
		break;

	default :
		return;
	}//end switch

	// ipKa is initialized if control passes to here.
	//protonation probability
	core::Real fprot = 1.0 / ( ( pow( 10.0, ( pH_ - ipKa ) ) ) + 1.0 );
	core::Real score_prot = -0.59 * log( fprot );
	core::Real score_deprot = -0.59 * log( 1.0 - fprot );
	
	core::Real pH_score = 0;
	switch ( rsd.type().aa() ) {
		case na_rad: 
			pH_score = (rsd.type().has_variant_type( PROTONATED_N1_ADENOSINE )) ? score_prot : score_deprot;
			break;
		
		case aa_asp:
		case aa_glu:
		case aa_his:
			pH_score = (rsd.type().has_variant_type( PROTONATED )) ? score_prot : score_deprot;
			break;
			
		case aa_lys:
		case aa_tyr:
			pH_score = (rsd.type().has_variant_type( DEPROTONATED )) ? score_deprot : score_prot;
			break;
			
		default: // should be impossible unless the two switch statements desync.
			break;
	}
	
	emap[ e_pH ] += pH_score; //add pHEnergy to the EmapVector
} // residue_energy


Real
pHEnergy::eval_dof_derivative(
	id::DOF_ID const & ,// dof_id,
	id::TorsionID const & , //tor_id
	pose::Pose const & , // pose
	ScoreFunction const &, //sfxn,
	EnergyMap const & // weights
) const
{
	return 0.0;
}


/// @pHEnergy is context independent
void
pHEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}
core::Size
pHEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}

