// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/TreeBuilderFactory.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/TreeBuilderFactory.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/algorithm/string/case_conv.hpp>

// Utility headers
#include <utility/exit.hh>

// Package headers
#include <protocols/nonlocal/SimpleTreeBuilder.hh>
#include <protocols/nonlocal/StarTreeBuilder.hh>
#include <protocols/nonlocal/TreeBuilder.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

TreeBuilderOP TreeBuilderFactory::get_builder(const std::string& builder_name) {
  // case-insensitive input
  std::string type = builder_name;
  boost::to_lower(type);

  if (type == "star") {
    return TreeBuilderOP( new StarTreeBuilder() );
  } else if (type == "simple") {
    return TreeBuilderOP( new SimpleTreeBuilder() );
  } else {
    utility_exit_with_message("Invalid builder_type: " + type);
  }
  return NULL;
}

}  // namespace nonlocal
}  // namespace protocols
