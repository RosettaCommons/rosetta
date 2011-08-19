// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/TriggerUtil.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/TriggerUtil.hh>

// C/C++ headers
#include <iostream>
#include <sstream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>

// Project headers
#include <core/types.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/GenericMonteCarloMover.hh>  // for Trigger def'n

namespace protocols {
namespace nonlocal {

static basic::Tracer TR("protocols.nonlocal.TriggerUtil");

double TriggerUtil::progress(core::Size p, core::Size np) {
  return p / static_cast<double>(np);
}

void TriggerUtil::update_scoring_option(const core::scoring::ScoreType& option,
                                        core::scoring::ScoreFunctionOP scoring,
                                        core::Real setting) {
  scoring->set_weight(option, setting);
}

bool TriggerUtil::ramp_chainbreaks(core::Size stage,
                                   core::Size num_stages,
                                   core::Size cycle,
                                   core::Size num_cycles,
                                   const core::pose::Pose&,
                                   core::scoring::ScoreFunctionOP scoring) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // Modify the score function every N cycles
  if ((cycle % option[OptionKeys::nonlocal::ramp_constraints_cycles]()) != 0)
    return false;

  // Compute the new linear chainbreak weight. The linear_chainbreak term is
  // ramped *up* within a single stage and across multiple stages, as the
  // simulation proceeds.
  double stage_mult = stage / static_cast<double>(num_stages);
  double p = progress(cycle, num_cycles);
  double chainbreak_weight = 1 + 1.5 * p * stage_mult;

  // Retrieve the weight multiplier
  double m = option[jumps::increase_chainbreak]();
  double setting = chainbreak_weight * m;
  update_scoring_option(core::scoring::linear_chainbreak, scoring, setting);
  TR.Debug << "stage: " << stage << " "
           << "cycle: " << cycle << " "
           << "linear_chainbreak weight => " << setting << std::endl;

  // Trigger rescoring
  return true;
}

bool TriggerUtil::ramp_constraints(core::Size stage,
                                   core::Size num_stages,
                                   core::Size cycle,
                                   core::Size num_cycles,
                                   const core::pose::Pose&,
                                   core::scoring::ScoreFunctionOP scoring) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // Modify the score function every N cycles
  if ((cycle % option[OptionKeys::nonlocal::ramp_constraints_cycles]()) != 0)
    return false;

  // Compute the new constraint weight as a function of inter- and intra-
  // stage progress. The atom_pair_constraint term is ramped *up* over the
  // course of a stage and ramped *down* over the course of the simulation.
  double stage_mult = 1 - (stage / (static_cast<double>(num_stages + 1)));
  double p = progress(cycle, num_cycles);
  double constraint_weight = (0.75 + 1.5 * p) * stage_mult;

  // Retrieve the weight multiplier
  double m = option[OptionKeys::constraints::increase_constraints]();
  double setting = constraint_weight * m;
  update_scoring_option(core::scoring::atom_pair_constraint, scoring, setting);
  TR.Debug << "stage: " << stage << " "
           << "cycle: " << cycle << " "
           << "atom_pair_constraint weight => " << setting << std::endl;

  // Trigger rescoring
  return true;
}

bool TriggerUtil::ramp_sequencesep(core::Size stage,
                                   core::Size,
                                   core::Size cycle,
                                   core::Size num_cycles,
                                   const core::pose::Pose& pose,
                                   core::scoring::ScoreFunctionOP scoring) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::scoring::methods::EnergyMethodOptions;

  // Modify the score function every N cycles
  if ((cycle % option[OptionKeys::nonlocal::ramp_constraints_cycles]()) != 0)
    return false;

  // Compute the new maximum sequence separation
  double p = progress(cycle, num_cycles);
  double sequence_sep = 1 + (p * pose.total_residue()) / 20;
  core::Size setting = static_cast<core::Size>(sequence_sep);

  // Update the maximum sequence separation setting, retaining the values
  // of all other energy method options
  EnergyMethodOptions new_options(scoring->energy_method_options());
  new_options.cst_max_seq_sep(setting);

  // Replace the score function
  scoring->set_energy_method_options(new_options);
  TR.Debug << "stage: " << stage << " "
           << "cycle: " << cycle << " "
           << "max_seqsep => " << setting << std::endl;

  // Trigger rescoring
  return true;
}

bool TriggerUtil::write_pose(core::Size,
                             core::Size,
                             core::Size cycle,
                             core::Size,
                             const core::pose::Pose& pose,
                             core::scoring::ScoreFunctionOP) {
  if ((cycle % 100) == 0) {
    std::stringstream ss;
    ss << "nla_cycle_" << cycle << ".pdb";
    core::io::pdb::dump_pdb(pose, ss.str());
  }

  // Don't trigger rescoring
  return false;
}

}  // namespace nonlocal
}  // namespace protocols
