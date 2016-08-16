// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ./protocols/fldsgn/AACompositionEnergyCreator.hh
/// @brief  Declaration for the class that connects AACompositionEnergy with the ScoringManager
/// @author Nobuyasu Koga

#ifndef INCLUDED_protocols_fldsgn_potentials_AACompositionEnergyCreator_HH
#define INCLUDED_protocols_fldsgn_potentials_AACompositionEnergyCreator_HH

#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {

class AACompositionEnergyCreator : public core::scoring::methods::EnergyMethodCreator {
public:

	typedef core::scoring::methods::EnergyMethodOP EnergyMethodOP;
	typedef core::scoring::methods::EnergyMethodOptions EnergyMethodOptions;
	typedef core::scoring::ScoreTypes ScoreTypes;

public:

	/// @brief Instantiate a new AACompositionEnergy
	virtual EnergyMethodOP create_energy_method( EnergyMethodOptions const & ) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual ScoreTypes score_types_for_method() const;

};

} // potentials
} // fldsgn
} // protocols

#endif
