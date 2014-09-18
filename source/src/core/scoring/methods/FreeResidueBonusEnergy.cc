// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FreeResidueBonusEnergy.cc
/// @brief  Score bonus for residues that form 0 interactions (i.e. bulges) - rough approximation to entropic bonus for flexible residues
/// @author Arvind Kannan


// Unit headers
#include <core/scoring/methods/FreeResidueBonusEnergy.hh>
#include <core/scoring/methods/FreeResidueBonusEnergyCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>

// Utility headers

// C++
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "core.scoring.methods.FreeResidueBonusEnergy" );

/////////////////////////////////////////////////////////////////////////////////////
//
// Created in order to stabilize bulges during SWM calculations
//
/////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the FreeResidueBonusEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
FreeResidueBonusEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new FreeResidueBonusEnergy;
}

ScoreTypes
FreeResidueBonusEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	//sts.push_back( bulge_bonus );
	return sts;
}
	
/// c-tor
FreeResidueBonusEnergy::FreeResidueBonusEnergy() :
	parent( new FreeResidueBonusEnergyCreator )
{
}

/// clone
methods::EnergyMethodOP
FreeResidueBonusEnergy::clone() const
{
	return new FreeResidueBonusEnergy;
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
	
void
FreeResidueBonusEnergy::residue_energy(conformation::Residue const &, pose::Pose const &, EnergyMap &) const
{}


///////////////////////////////////////////////////////////////////////////////
core::Size
FreeResidueBonusEnergy::version() const
{
	return 1;
}


} // methods
} // scoring
} // core
