// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/rational_mc/RationalMonteCarlo.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDE_protocols_simple_moves_rational_mc_RationalMonteCarlo_HH
#define INCLUDE_protocols_simple_moves_rational_mc_RationalMonteCarlo_HH

// Unit header
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.fwd.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/function.hpp>
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.hh>

// Package headers
#include <protocols/moves/Mover.hh>

//Auto Headers
#include <protocols/moves/MonteCarlo.fwd.hh>
namespace protocols {
namespace simple_moves {
namespace rational_mc {

/// @brief Trigger API definition
typedef boost::function<void(const core::pose::Pose&)> RationalMonteCarloTrigger;
typedef boost::unordered_map<int, RationalMonteCarloTrigger> Triggers;

/// @class Trial-based Monte Carlo minization primitive. Do not modify this
/// class; almost anything you could possibly want to add is a bad idea.
class RationalMonteCarlo : public protocols::moves::Mover {
  typedef core::pose::Pose Pose;
  typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

 public:
  RationalMonteCarlo(moves::MoverOP mover,
                     ScoreFunctionOP score,
                     unsigned num_trials,
                     double temperature,
                     bool recover_low);

  /// @brief Applies the underlying mover to <pose> the specified number of times
  void apply(Pose& pose);

  /// @brief Returns the number of trials to be attempted in calls to apply()
  unsigned num_trials() const;

  /// @brief Returns true if we should recover the low scoring pose at the end
  /// of apply, false otherwise.
  bool recover_low() const;

  /// @brief Returns this mover's name
  std::string get_name() const;


  // -- Triggers -- //

  /// @brief Registers the specified trigger with this instance. Returns a unique
  /// identifier for referring to this trigger for subsequent operations (e.g. remove).
  int add_trigger(const RationalMonteCarloTrigger& trigger);

  /// @brief Unregisters the specified trigger, warning the user if one is not found.
  void remove_trigger(int trigger_id);


  // -- Analytics -- //

  /// @brief Sets the native pose, enabling a variety of analytics to be collected.
  /// Clears the contents of gdtmss() and rmsds().
  void set_native(const core::pose::Pose& native);

  /// @brief Returns gdtmm as a function of time. Requires a native structure.
  const utility::vector1<double>& gdtmms() const;

  /// @brief Returns rmsd as a function of time. Requires a native structure.
  const utility::vector1<double>& rmsds() const;


 protected:
  /// @brief Executes all triggers attached to this instance. The order of trigger
  /// execution is undefined. Do not assume, depend, or in any way rely upon a
  /// partiular ordering.
  void fire_all_triggers(const Pose& pose);

  /// @brief Built-in trigger for computing analytics given the native structure
  void compute_analytics(const core::pose::Pose& pose);


 private:
  /// @brief Underlying mover
  moves::MoverOP mover_;

  /// @brief Determines whether applications of the base mover should be accepted
  /// by applying the Metropolis criterion to the pose
  moves::MonteCarloOP mc_;

  /// @brief Number of times to execute the underlying mover in calls to apply()
  unsigned num_trials_;

  /// @brief Determines whether the low scoring pose should be recovered at the
  /// conclusion of the apply() method
  bool recover_low_;

  /// @brief Collection of function callbacks
  Triggers triggers_;

  /// @brief Next trigger id to be assigned
  unsigned next_trigger_id_;

  /// @brief Native structure
  Pose native_;

  /// @brief Tracks gdtmm as a function of time if a native was provided
  utility::vector1<double> gdtmms_;

  /// @brief Tracks rmsd as a function of time if a native was provided
  utility::vector1<double> rmsds_;
};

}  // namespace rational_mc
}  // namespace simple_moves
}  // namespace protocols

#endif  // protocols_simple_moves_rational_mc_RationalMonteCarlo_HH_
