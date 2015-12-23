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
/// @author Colin Smith
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>

#include <core/pose/Pose.hh>

#include <utility/vector1.hh>
#include <basic/MetricValue.hh>


namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

MetricValueGetter::MetricValueGetter() :
	calculator_(""),
	key_(""),
	metric_value_template_(/* NULL */)
{}

MetricValueGetter::MetricValueGetter(
	std::string const & calculator,
	std::string const & key,
	basic::MetricValueBaseCOP metric_value_template
) :
	calculator_(calculator),
	key_(key),
	metric_value_template_(metric_value_template)
{}

MetricValueGetter::MetricValueGetter(
	MetricValueGetter const & getter
)
{
	*this = getter;
}

MetricValueGetter::~MetricValueGetter()
{}

MetricValueGetter &
MetricValueGetter::operator = (
	MetricValueGetter const & getter
)
{
	calculator_ = getter.calculator();
	key_ = getter.key();
	metric_value_template_ = getter.metric_value_template();

	return *this;
}

std::string const &
MetricValueGetter::calculator() const
{
	return calculator_;
}

void
MetricValueGetter::calculator(
	std::string const & calculator
)
{
	calculator_ = calculator;
}

std::string const &
MetricValueGetter::key() const
{
	return key_;
}

void
MetricValueGetter::key(
	std::string const & key
)
{
	key_ = key;
}

basic::MetricValueBaseCOP
MetricValueGetter::metric_value_template() const
{
	return metric_value_template_;
}

void
MetricValueGetter::metric_value_template(
	basic::MetricValueBaseCOP metric_value_template
)
{
	metric_value_template_ = metric_value_template;
}

basic::MetricValueBaseOP
MetricValueGetter::get(
	core::pose::Pose const & pose
) const
{
	runtime_assert(metric_value_template_ != 0);
	basic::MetricValueBaseOP new_metric_value(metric_value_template_->clone());
	pose.metric(calculator_, key_, *new_metric_value);
	return new_metric_value;
}

} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols


