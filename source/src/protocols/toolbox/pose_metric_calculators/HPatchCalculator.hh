// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/PoseMetricCalculators/HPatchCalculator.hh
/// @brief A Pose metric which will keep track of the SASA-based hpatch score of a Pose object
/// @author Ron Jacak

#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_HPatchCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_HPatchCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID.fwd.hh>

#include <map>

#include <utility/vector1.hh>

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

class HPatchCalculator : public core::pose::metrics::StructureDependentCalculator {

public:
	HPatchCalculator( bool remove_nonprotein_res = false );
	~HPatchCalculator();

public:
	core::pose::metrics::PoseMetricCalculatorOP clone() const { return core::pose::metrics::PoseMetricCalculatorOP( new HPatchCalculator( remove_nonprotein_res_ ) ); };

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase* valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:
	core::Real total_hpatch_score_;
	bool remove_nonprotein_res_;
	std::map< Size, std::pair< core::Real, core::Real > > patch_scores_;
	std::map< Size, utility::vector1< core::id::AtomID > > atoms_in_patches_;
};

} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
