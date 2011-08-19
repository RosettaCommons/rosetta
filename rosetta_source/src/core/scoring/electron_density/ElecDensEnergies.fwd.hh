// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ElecDensEnergies.fwd.hh
/// @brief  Scoring a structure's fit to electron density
/// @author Frank DiMaio


#ifndef INCLUDED_core_scoring_electron_density_ElecDensEnergies_fwd_hh
#define INCLUDED_core_scoring_electron_density_ElecDensEnergies_fwd_hh
// Utility headers


namespace core {
namespace scoring {
namespace methods {

/// sliding-window fit-to-density
class ElecDensEnergy;

/// ca-only fit-to-density
class ElecDensCenEnergy;

/// all-atom fit-to-density
class ElecDensAllAtomCenEnergy;


}
}
}

#endif
