// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/DEEREnergyCreator.hh
/// @brief  Parent class for DEEREnergy structure
/// @details To prevent the energy method from reading from the command line every scoring
///      round and parsing the input file, this method constructs a faux-graph of DEER
///      decay data. The "nodes" are the residues for which at least one set of data
///      is provided, and the "edges" correspond to data sets for two residues.
///      NOTE: There is a graph class in src/utility but it is unnecessarily heavy.
///      Instead this class uses maps
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_energy_methods_DEEREnergyCreator_hh
#define INCLUDED_core_energy_methods_DEEREnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace core {
namespace energy_methods {

class DEEREnergyCreator : public scoring::methods::EnergyMethodCreator
{
public:

	/// @brief  Energy Method creator required by framework (registrator in core/init.cc)
	scoring::methods::EnergyMethodOP
	create_energy_method( scoring::methods::EnergyMethodOptions const & options ) const override;

	/// @brief  Energy Method descriptor to identify it if a weight is provided in the command line
	scoring::ScoreTypes
	score_types_for_method() const override;

};

} // namespace energy_methods
} // namespace core

#endif
