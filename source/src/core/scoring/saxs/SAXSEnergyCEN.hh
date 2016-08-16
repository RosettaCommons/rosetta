// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/saxs/SAXSEnergy.hh
/// @brief  "Energy" based on a similarity of theoretical SAXS spectrum computed for a pose and the experimental data
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_core_scoring_saxs_SAXSEnergyCEN_hh
#define INCLUDED_core_scoring_saxs_SAXSEnergyCEN_hh

// Package headers
#include <core/scoring/saxs/FormFactor.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorCEN.hh>
#include <core/scoring/saxs/SAXSEnergy.hh>

#include <core/scoring/methods/WholeStructureEnergy.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <string>

#include <core/chemical/ChemicalManager.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace saxs {


class SAXSEnergyCEN : public SAXSEnergy  {
public:

	/// c-tor
	SAXSEnergyCEN() : SAXSEnergy( cen_cfg_file_,
		chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID),saxs_cen_score,methods::EnergyMethodCreatorOP( new SAXSEnergyCreatorCEN )) {}

};


}
}
}

#endif
