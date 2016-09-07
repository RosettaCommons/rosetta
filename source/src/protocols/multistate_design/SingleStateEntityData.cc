// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SingleStateEntityData.cc
/// @brief
/// @author Colin A. Smith

#include <protocols/multistate_design/SingleStateEntityData.hh>

#include <basic/MetricValueIO.hh>


namespace protocols {
namespace multistate_design {

void
SingleStateEntityData::write_checkpoint(
	std::ostream & os
) const
{
	os << "fitness " << fitness_;
	os << " metrics " << metric_value_map_.size();
	for ( auto iter = metric_value_map_.begin(); iter != metric_value_map_.end(); ++iter ) {
		os << " " << iter->first << " ";
		runtime_assert(basic::write_metric_value(os, *(iter->second)));
	}
}

bool
SingleStateEntityData::read_checkpoint(
	std::istream & is
)
{
	std::string word;
	core::Size num_metrics;

	if ( !(is >> word) ) return false;
	if ( word != "fitness" ) return false;
	if ( !(is >> fitness_) ) return false;

	if ( !(is >> word) ) return false;
	if ( word != "metrics" ) return false;
	if ( !(is >> num_metrics) ) return false;

	for ( core::Size i = 1; i <= num_metrics; ++i ) {
		if ( !(is >> word) ) return false;
		basic::MetricValueBaseOP metric_value(basic::read_metric_value(is));
		if ( !metric_value ) return false;
		metric_value_map_[word] = metric_value;
	}

	return true;
}

} // namespace multistate_design
} // namespace protocols
