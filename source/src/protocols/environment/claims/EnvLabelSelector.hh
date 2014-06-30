// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/environment/claims/EnvLabelSelector.hh
/// @brief  The EnvLabelSelector holds a set of local positions that are converted to seqids at apply.
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_environment_claims_EnvLabelSelector_HH
#define INCLUDED_protocols_environment_claims_EnvLabelSelector_HH

// Unit headers
#include <protocols/environment/claims/EnvLabelSelector.fwd.hh>
#include <core/pack/task/residue_selector/ResidueSelector.hh>

// Package headers
#include <core/environment/LocalPosition.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

namespace protocols {
namespace environment {
namespace claims {

class EnvLabelSelector : public core::pack::task::residue_selector::ResidueSelector {
  typedef core::pack::task::residue_selector::ResidueSubset ResidueSubset;
  typedef core::environment::LocalPositions LocalPositions;

public:
  // derived from base class
  EnvLabelSelector();

  EnvLabelSelector( LocalPositions const& );

  virtual ~EnvLabelSelector();

  virtual void apply( core::pose::Pose const & pose,
                      core::pack::task::residue_selector::ResidueSubset & subset ) const;

  virtual void parse_my_tag(
    utility::tag::TagCOP tag,
    basic::datacache::DataMap & datamap
  );

  void set_local_positions( LocalPositions const& );

  LocalPositions const& local_positions() const{ return positions; }

  virtual
  std::string
  get_name() const;

  static std::string class_name();

private: // data members
  LocalPositions positions;
};

} //namespace claims
} //namespace environment
} //namespace protocols


#endif
