// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/NetChargeFilter.hh
/// @brief definition of filter class NetChargeFilter.
/// @author Dave La (davela@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_NetChargeFilter_hh
#define INCLUDED_protocols_simple_filters_NetChargeFilter_hh

#include <protocols/simple_filters/NetChargeFilter.fwd.hh>


// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/exit.hh>
#include <core/pack/task/TaskFactory.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

class NetChargeFilter : public filters::Filter
{
public:
	NetChargeFilter();
	//NetChargeFilter( core::Size const distance, core::Size const jump_num );
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new NetChargeFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new NetChargeFilter() );
	}

	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );

	~NetChargeFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Size chain_;
	signed int net_charge_max_;
	signed int net_charge_min_;
	core::pack::task::TaskFactoryOP task_factory_; // dflt NULL ; if set, all designable residues are counted for this filter
};

}
}
#endif
