// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/LinearChainbreakEnergyCreator.cc
/// @brief  LinearChainbreakEnergyCreator implementation
/// @author Christopher Miles (cmiles@uw.edu)

#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/types.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/LinearChainbreakEnergy.hh>
#include <core/scoring/methods/LinearChainbreakEnergyCreator.hh>

#include <utility/vector1.hh>


namespace core {
  namespace scoring {
    namespace methods {

      /// @details Creates and initializes a new LinearChainbreakEnergy object with the
      /// specified options
      EnergyMethodOP
      LinearChainbreakEnergyCreator::create_energy_method(const EnergyMethodOptions& opt) const {
				// EnergyMethodOptions is responsible for providing a reasonable default
				// in the event that the caller does not specify a sequence separation
				// constraint. Currently, cst_max_seq_sep() defaults to +inf.
				Size allowable_sequence_separation = opt.cst_max_seq_sep();
				return EnergyMethodOP( new LinearChainbreakEnergy(allowable_sequence_separation) );
      }

      ScoreTypes
      LinearChainbreakEnergyCreator::score_types_for_method() const {
				ScoreTypes types;
				types.push_back(linear_chainbreak);
				types.push_back(overlap_chainbreak);
				return types;
      }

    }  // namespace methods
  }  // namespace scoring
}  // namespace core
