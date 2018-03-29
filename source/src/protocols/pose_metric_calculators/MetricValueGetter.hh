// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Colin Smith


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_MetricValueGetter_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_MetricValueGetter_hh
#include <protocols/pose_metric_calculators/MetricValueGetter.fwd.hh>

#include <core/pose/Pose.fwd.hh>


#include <utility/vector1.hh>
#include <ostream>
#include <basic/MetricValue.fwd.hh>


namespace protocols {
namespace pose_metric_calculators {

class MetricValueGetter {

public:

	MetricValueGetter();

	MetricValueGetter(
		std::string const & calculator,
		std::string const & key,
		basic::MetricValueBaseCOP metric_value_template
	);

	MetricValueGetter(
		MetricValueGetter const & getter
	);

	~MetricValueGetter();

	MetricValueGetter &
	operator = (
		MetricValueGetter const & getter
	);

	std::string const &
	calculator() const;

	void
	calculator(
		std::string const & calculatr
	);

	std::string const &
	key() const;

	void
	key(
		std::string const & key
	);

	basic::MetricValueBaseCOP
	metric_value_template() const;

	void
	metric_value_template(
		basic::MetricValueBaseCOP metric_value_template
	);

	basic::MetricValueBaseOP
	get(
		core::pose::Pose const & pose
	) const;

private:

	std::string calculator_;
	std::string key_;
	basic::MetricValueBaseCOP metric_value_template_;
};


} // namespace pose_metric_calculators
} // namespace protocols

#endif
