// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/ExtractSubpose.hh
/// @brief  Extract primary component associated with symdofs and all neighboring components.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_devel_matdes_ExtractSubpose_hh
#define INCLUDED_devel_matdes_ExtractSubpose_hh

// Unit Headers
#include <devel/matdes/ExtractSubpose.fwd.hh>
#include <devel/matdes/ExtractSubposeCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/types.hh>

// Utility Headers

// C++ Headers

namespace devel { 
namespace matdes {

class ExtractSubpose : public protocols::moves::Mover {
public:

  // default constructor
	ExtractSubpose();

	ExtractSubpose(const ExtractSubpose& rval);

	virtual std::string get_name() const { return "ExtractSubpose"; }
  virtual void apply( core::pose::Pose & pose );

  virtual protocols::moves::MoverOP clone() const;
  virtual protocols::moves::MoverOP fresh_instance() const;

  void parse_my_tag(
      utility::tag::TagCOP const tag,
      basic::datacache::DataMap &data,
      protocols::filters::Filters_map const &filters,
      protocols::moves::Movers_map const &movers,
      core::pose::Pose const & pose );

private:

  std::string sym_dof_names_;
	std::string prefix_;
	std::string suffix_;
  core::Real contact_dist_;
	bool extras_;
};

} //namespace matdes 
} //namespace devel

#endif 
