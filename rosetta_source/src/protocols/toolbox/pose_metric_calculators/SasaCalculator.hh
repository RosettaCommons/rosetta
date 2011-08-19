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
/// @author John Karanicolas


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_SasaCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_SasaCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>

#include <basic/options/option.hh>

#include <utility/vector1.hh>


// option key includes

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>


namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

class SasaCalculator : public core::pose::metrics::StructureDependentCalculator {

public:

	SasaCalculator( core::Real probe_r = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius] ) : probe_radius_(probe_r) {}

	core::pose::metrics::PoseMetricCalculatorOP clone() const { return new SasaCalculator(); };

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	core::Real total_sasa_;
	core::id::AtomID_Map< core::Real > atom_sasa_;
	utility::vector1< core::Real > residue_sasa_;
	core::Real probe_radius_;

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
