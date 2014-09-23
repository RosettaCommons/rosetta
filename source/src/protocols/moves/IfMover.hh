// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/IfMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_moves_IfMover_hh
#define INCLUDED_protocols_moves_IfMover_hh

// C/C++ headers
#include <string>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

class IfMover : public protocols::moves::Mover {
 public:
  /// @brief No-argument constructor
  IfMover() : protocols::moves::Mover("If") {}

  /// @brief Virtual destructor
  virtual ~IfMover() {};

  protocols::moves::MoverOP clone() const {
    return protocols::moves::MoverOP( new IfMover(*this) );
  }

  protocols::moves::MoverOP fresh_instance() const {
    return protocols::moves::MoverOP( new IfMover() );
  }

  void apply(core::pose::Pose& pose);

  // Required for backwards compatibility.
  // Synonym for `get_additional_output_true_mover().`
  core::pose::PoseOP get_additional_output();

  core::pose::PoseOP get_additional_output_true_mover();
  core::pose::PoseOP get_additional_output_false_mover();

  virtual std::string get_name() const;

  void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

 private:
  protocols::filters::FilterOP filter_;

  /// @brief Invoked when filter evaluates to true
  protocols::moves::MoverOP true_mover_;

  /// @brief Invoked when filter evaluates to false
  protocols::moves::MoverOP false_mover_;
};

} // moves
} // protocols

#endif /*INCLUDED_protocols_ProteinInterfaceDesign_movers_IfMover_HH*/
