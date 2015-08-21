// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/BuriedUnsatisfiedPolarsCalculator.hh
/// @brief
/// @author Florian Richter


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_BuriedUnsatisfiedPolarsCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_BuriedUnsatisfiedPolarsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>

#include <basic/options/option.hh>

#include <utility/vector1.hh>

#include <set>


// option key includes

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>


namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

class BuriedUnsatisfiedPolarsCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:

	BuriedUnsatisfiedPolarsCalculator(
		std::string sasa_calc,
		std::string hbond_calc,
		core::Real burial_cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::atomic_burial_cutoff]
	);


	BuriedUnsatisfiedPolarsCalculator(
		std::string sasa_calc,
		std::string hbond_calc,
		std::set< core::Size > const & special_region,
		core::Real burial_cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::atomic_burial_cutoff]
	);


	core::pose::metrics::PoseMetricCalculatorOP clone() const {
		return core::pose::metrics::PoseMetricCalculatorOP( new BuriedUnsatisfiedPolarsCalculator( name_of_sasa_calc_, name_of_hbond_calc_, burial_sasa_cutoff_) ); };

	std::string const & name_of_hbond_calc() const { return name_of_hbond_calc_; }
	std::string const & name_of_sasa_calc() const { return name_of_sasa_calc_; }

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );


private:

	void assert_calculators();

	static
	core::Size satisfaction_cutoff( std::string atom_type );


	core::Size all_bur_unsat_polars_;
	core::Size special_region_bur_unsat_polars_;
	core::id::AtomID_Map< bool > atom_bur_unsat_;
	utility::vector1< core::Size > residue_bur_unsat_polars_;

	//holds the sasa and atom hbonds calculators necessary for this calculator
	std::string name_of_hbond_calc_, name_of_sasa_calc_;
	core::Real burial_sasa_cutoff_;

	std::set< core::Size > special_region_;

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
