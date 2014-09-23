// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/methods/pHEnergy.cc
/// @brief  Energy due to ionization state of the residue at a particular pH (see pHEnergy.hh for a detailed description)
/// @author krishna

// Standard libraries
#include <cmath>

// Unit headers
// AUTO-REMOVED #include <core/scoring/methods/util.hh>
#include <core/scoring/methods/pHEnergy.hh>
#include <core/scoring/methods/pHEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>

// Project headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>


// Option Key includes
#include <basic/options/keys/pH.OptionKeys.gen.hh>

#include <utility/vector1.hh>


// C++

namespace core {
namespace scoring {
namespace methods {

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
) const
{
  core::Real ipKa_, pH_score_, fprot_, score1_, score2_;

  switch ( rsd.type().aa() ) {

    case chemical::aa_asp:
      ipKa_ = 4.0;
      fprot_ = 1.0/((pow (10.0,(pH_ - ipKa_))) + 1.0);                        //protonation probability
      score1_ = -0.59 * log (fprot_);                                        //if residue is protonated
      score2_ = -0.59 * log (1.0 - fprot_);                                  //if residue is not protonated
      pH_score_ = (rsd.type().has_variant_type( chemical::PROTONATED )) ? score1_ : score2_;
      break;

    case chemical::aa_glu:
      ipKa_ = 4.4;
      fprot_ = 1.0/((pow (10.0,(pH_ - ipKa_))) + 1.0);
      score1_ = -0.59 * log (fprot_);
      score2_ = -0.59 * log (1.0 - fprot_);
      pH_score_ = (rsd.type().has_variant_type( chemical::PROTONATED )) ? score1_ : score2_;
      break;

    case chemical::aa_his:
      ipKa_ = 6.3;
      fprot_ = 1.0/((pow (10.0,(pH_ - ipKa_))) + 1.0);
      score1_ = -0.59 * log (fprot_);
      score2_ = -0.59 * log (1.0 - fprot_);
      pH_score_ = (rsd.type().has_variant_type( chemical::PROTONATED )) ? score1_ : score2_;
      break;

    case chemical::aa_lys:
      ipKa_ = 10.4;
      fprot_ = 1.0/((pow (10.0,(pH_ - ipKa_))) + 1.0);
      score1_ = -0.59 * log (fprot_);
      score2_ = -0.59 * log (1.0 - fprot_);
      pH_score_ = (rsd.type().has_variant_type( chemical::DEPROTONATED )) ? score2_ : score1_;
      break;

    case chemical::aa_tyr:
      ipKa_ = 10.0;
      fprot_ = 1.0/((pow (10.0,(pH_ - ipKa_))) + 1.0);
      score1_ = -0.59 * log (fprot_);
      score2_ = -0.59 * log (1.0 - fprot_);
      pH_score_ = (rsd.type().has_variant_type( chemical::DEPROTONATED )) ? score2_ : score1_;
      break;

    default:
      pH_score_ = 0.0;
      break;

  }//end switch

  emap[ e_pH ] += pH_score_; //add pHEnergy to the EmapVector

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

