// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/PackstatCalculator.hh
/// @brief
/// @author Florian Richter


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_PackstatCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_PackstatCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>

#include <basic/options/option.hh>

#include <utility/vector1.hh>

#include <set>


// option key includes

#include <basic/options/keys/packstat.OptionKeys.gen.hh>


namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

class PackstatCalculator : public core::pose::metrics::StructureDependentCalculator {

public:

  PackstatCalculator(
    core::Size oversample = basic::options::option[basic::options::OptionKeys::packstat::oversample],
    bool remove_nonprotein_res = false
  );


  PackstatCalculator(
    std::set< core::Size > const & special_region,
    core::Size oversample = basic::options::option[basic::options::OptionKeys::packstat::oversample],
    bool remove_nonprotein_res = false
  );


  core::pose::metrics::PoseMetricCalculatorOP clone() const {
    return core::pose::metrics::PoseMetricCalculatorOP( new PackstatCalculator( special_region_, oversample_) ); };

protected:

  virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
  virtual std::string print( std::string const & key ) const;
  virtual void recompute( core::pose::Pose const & this_pose );


private:

  core::Real total_packstat_;
  core::Real special_region_packstat_;
  utility::vector1< core::Real > residue_packstat_;

  core::Size oversample_;
  bool remove_nonprotein_res_;

  std::set< core::Size > special_region_;

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
