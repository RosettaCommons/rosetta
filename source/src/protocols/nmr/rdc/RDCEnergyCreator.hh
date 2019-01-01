// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/rdc/RDCEnergyCreator.hh
/// @brief   Declaration for the class that connects RDCEnergy with the ScoringManager
/// @details last Modified: 08/03/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_rdc_RDCEnergyCreator_hh
#define INCLUDED_protocols_nmr_rdc_RDCEnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nmr {
namespace rdc {

class RDCEnergyCreator : public core::scoring::methods::EnergyMethodCreator
{
public:
	// @brief Instantiate a new RDCEnergy
	virtual core::scoring::methods::EnergyMethodOP
	create_energy_method( core::scoring::methods::EnergyMethodOptions const & ) const;

	// @brief Return the set of score types claimed by the EnergyMethod
	// this EnergyMethodCreator creates in its create_energy_method() function
	virtual core::scoring::ScoreTypes
	score_types_for_method() const;

};

} // namespace rdc
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_rdc_RDCEnergyCreator_hh
