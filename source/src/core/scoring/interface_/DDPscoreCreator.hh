// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/Interface_/DDPscoreCreator.hh
/// @author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)


#ifndef INCLUDED_core_scoring_interface_DDPscoreCreator_hh
#define INCLUDED_core_scoring_interface_DDPscoreCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace interface_ {

class DDPscoreCreator : public methods::EnergyMethodCreator
{
public:
	/// @brief Instantiate a new NVscore
	virtual
	methods::EnergyMethodOP
		create_energy_method(
		methods::EnergyMethodOptions const &
	) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual
	ScoreTypes
	score_types_for_method() const;

};

} // Interface_
} // scoring
} // core

#endif /* INCLUDED_core_scoring_Interface_DDPScoreCreator_HH_ */
