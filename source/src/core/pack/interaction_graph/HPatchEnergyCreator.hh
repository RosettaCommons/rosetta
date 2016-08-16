// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/methods/HPatchEnergyCreator.hh
/// @brief  Declaration for the class that connects HPatchEnergy with the ScoringManager
/// @author Ron Jacak

#ifndef INCLUDED_core_pack_interaction_graph_HPatchEnergyCreator_hh
#define INCLUDED_core_pack_interaction_graph_HPatchEnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace interaction_graph {

class HPatchEnergyCreator : public scoring::methods::EnergyMethodCreator {

public:

	/// @brief Instantiate a new HPatchEnergy
	virtual
	scoring::methods::EnergyMethodOP create_energy_method( scoring::methods::EnergyMethodOptions const & ) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual
	scoring::ScoreTypes score_types_for_method() const;

};

}
}
}

#endif
