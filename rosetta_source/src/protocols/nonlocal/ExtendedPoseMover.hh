// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/ExtendedPoseMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_EXTENDEDPOSEMOVER_HH_
#define PROTOCOLS_NONLOCAL_EXTENDEDPOSEMOVER_HH_

// Unit header
#include <protocols/nonlocal/ExtendedPoseMover.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace nonlocal {

class ExtendedPoseMover : public protocols::moves::Mover {
  typedef std::string string;

 public:
  ExtendedPoseMover(const string& sequence = "",
                    const string& residue_type_set_ = "centroid");

  /// @brief Creates an extended, idealized pose from the sequence and residue
  /// type set specified in the constructor.
  void apply(core::pose::Pose& pose);

  /// @brief Returns the name of this mover
  string get_name() const;

  /// @brief Returns true if this instance is valid (i.e. contains a sequence).
  /// Since RosettaScripts mandates the existence of a no-argument constructor,
  /// we lose the ability to reason about the validity of a particular instance.
  bool valid() const;

  // -- Accessors -- //
  /// @brief Returns the sequence
  const string& sequence() const;

  /// @brief Returns the residue type set
  const string& residue_type_set() const;

  // -- Mutators -- //
  /// @brief Updates the sequence to be used in calls to apply()
  void sequence(const string& sequence);

  /// @brief Updates the residue type set to be used in calls to apply()
  void residue_type_set(const string& residue_type_set);

  // -- RosettaScripts -- //
  /// @brief Creates a replica of this Mover
  protocols::moves::MoverOP clone() const;

  /// @brief Creates a new instance by calling the no-argument constructor
  protocols::moves::MoverOP fresh_instance() const;

  /// @brief Mover-specific parsing required by RosettaScripts
  void parse_my_tag(const utility::tag::TagPtr tag,
                    protocols::moves::DataMap& data,
                    const protocols::filters::Filters_map& filters,
                    const protocols::moves::Movers_map& movers,
                    const core::pose::Pose& pose);

 private:
  string sequence_;
  string residue_type_set_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_EXTENDEDPOSEMOVER_HH_
