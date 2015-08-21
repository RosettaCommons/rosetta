// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_ClashCountCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_ClashCountCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

class ClashCountCalculator : public core::pose::metrics::StructureDependentCalculator {

public:

	ClashCountCalculator( core::Real clash_threshold );

	core::pose::metrics::PoseMetricCalculatorOP clone() const { return core::pose::metrics::PoseMetricCalculatorOP( new ClashCountCalculator( clash_threshold_ ) ); };

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:
	core::Size total_clashes_; //count all clashes above threshold
	core::Real clash_threshold_;
	core::Size bb_clashes_; //count only clashes that do no involve CB, CEN or other sidechain atoms

	// const-ref to scoring database
	core::Real const vdw_scale_factor_;
};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
