// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/energy_methods/CCS_IMMSComplexEnergyCreator.cc
/// @brief Implementation for the class that connects CCS_IMMSComplexEnergy with the ScoringManager
/// @author Akshaya Narayanasamy <akshaya.researcher@gmail.com>

#include <core/energy_methods/CCS_IMMSComplexEnergyCreator.hh>
#include <core/energy_methods/CCS_IMMSEnergy.hh>

#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreType.hh>

#include <utility/vector1.hh>

namespace core {
namespace energy_methods {

/// @details This must return a fresh instance of the CCS_IMMSComplexEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
CCS_IMMSComplexEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const & ) const
{
	return utility::pointer::make_shared< CCS_IMMSComplexEnergy >();
}

core::scoring::ScoreTypes
CCS_IMMSComplexEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back(ccs_imms_complex);
	return sts;
}

} // namespace energy_methods
} // namespace core 