// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/ChargeCalculator.hh
/// @brief
/// @author Florian Richter, floric@u.washington.edu, nov 2010


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_ChargeCalculator_HH
#define INCLUDED_protocols_toolbox_pose_metric_calculators_ChargeCalculator_HH

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>

#include <set>

#include <utility/vector1.hh>


namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

class ChargeCalculator : public core::pose::metrics::StructureDependentCalculator {

public:

  ChargeCalculator();


  ChargeCalculator(
    std::set< core::Size > const & special_region
  );

  ~ChargeCalculator();


  core::pose::metrics::PoseMetricCalculatorOP clone() const {
    return core::pose::metrics::PoseMetricCalculatorOP( new ChargeCalculator( special_region_) ); };

protected:

  virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
  virtual std::string print( std::string const & key ) const;
  virtual void recompute( core::pose::Pose const & this_pose );


private:

  core::Real total_charge_;
  core::Size total_pos_charges_;
  core::Size total_neg_charges_;

  core::Real SR_total_charge_;
  core::Size SR_total_pos_charges_;
  core::Size SR_total_neg_charges_;

  std::set< core::Size > special_region_;

};


} // namespace PoseMetricCalculators
} // namespace toolbox
} // namespace protocols

#endif
