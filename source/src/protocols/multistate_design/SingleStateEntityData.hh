// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SingleStateEntityData.hh
/// @brief
/// @author Colin A. Smith

#ifndef INCLUDED_protocols_multistate_design_SingleStateEntityData_hh
#define INCLUDED_protocols_multistate_design_SingleStateEntityData_hh

#include <protocols/multistate_design/SingleStateEntityData.fwd.hh>

#include <core/types.hh>
#include <basic/MetricValue.hh>

#include <string>
#include <map>

namespace protocols {
namespace multistate_design {

class SingleStateEntityData {

public:

	typedef std::map<std::string, basic::MetricValueBaseOP > MetricValueMap;

	SingleStateEntityData() : fitness_(0) {}
	virtual ~SingleStateEntityData() {}

	core::Real fitness() const { return fitness_; }
	void fitness(core::Real fitness) { fitness_ = fitness; }

	basic::MetricValueBaseCOP
	metric_value(
		std::string const & name
	) const
	{
		MetricValueMap::const_iterator iter(metric_value_map_.find(name));
		if ( iter != metric_value_map_.end() ) return NULL;
		return iter->second;
	}

	void
	metric_value(
		std::string const & name,
		basic::MetricValueBaseOP metric_value
	)
	{
		metric_value_map_[name] = metric_value;
	}

	MetricValueMap const &
	metric_value_map() const
	{
		return metric_value_map_;
	}

	virtual void write_checkpoint( std::ostream & os ) const;
	virtual bool read_checkpoint( std::istream & is );

private:

	core::Real fitness_;
	MetricValueMap metric_value_map_;
};

} // namespace multistate_design
} // namespace protocols

#endif
