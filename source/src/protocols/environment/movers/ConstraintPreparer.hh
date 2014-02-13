// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/ConstraintPreparer.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_ConstraintPreparer_hh
#define INCLUDED_protocols_environment_ConstraintPreparer_hh

// Unit Headers
#include <protocols/environment/movers/ConstraintPreparer.fwd.hh>

// Package headers
#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/movers/StagePreparer.hh>

// Project headers
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/pose/Pose.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class ConstraintPreparer : public StagePreparer {
typedef StagePreparer Parent;

public:
  ConstraintPreparer();

  void prepare( core::pose::Pose& pose, core::Real progress );

  void parse_my_tag( utility::tag::TagCOP const,
                     basic::datacache::DataMap&,
                     protocols::filters::Filters_map const&,
                     protocols::moves::Movers_map const&,
                     core::pose::Pose const& );

  void cst_file( std::string const& );

  std::string const& cst_file() const { return filename_; }

  void combine_ratio( core::Size const& s ) { combine_ratio_ = s; }
  core::Size const& combine_ratio() const { return combine_ratio_; }

  bool skip_redundant() const { return skip_redundant_; }
  void skip_redundant( bool s ) { skip_redundant_ = s; }

  core::Size skip_redundant_width() const { return skip_redundant_width_; }
  void skip_redundant_width( core::Size const& s ) { skip_redundant_width_ = s; }

  core::Real random_drop_rate() const { return rand_drop_rate_; }
  void random_drop_rate( core::Real const& s ) { rand_drop_rate_ = s; }

  utility::vector1< bool > const& combine_exclude_res() const { return combine_exclude_res_; }
  void combine_exclude_file( std::string const& filename );

  std::string get_name() const;

  virtual
  moves::MoverOP fresh_instance() const { return new ConstraintPreparer(); }

  virtual
  moves::MoverOP clone() const { return new ConstraintPreparer( *this ); }

  //prepares use the prepare method instead of apply
  virtual
  void apply( core::pose::Pose& ) { assert(false); };

  virtual
  claims::EnvClaims yield_claims( core::pose::Pose& ) { };


protected:
  void load_constraints( core::pose::Pose const& );

private:
  core::Size combine_ratio_;
  bool skip_redundant_;
  core::Size skip_redundant_width_;
  core::Real rand_drop_rate_;
  utility::vector1< bool > combine_exclude_res_;

  std::string filename_;
  core::scoring::constraints::ConstraintSetOP constraints_;

}; // end ConstraintPreparer base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_ConstraintPreparer_hh
