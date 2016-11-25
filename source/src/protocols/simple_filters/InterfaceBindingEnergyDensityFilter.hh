// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/InterfaceBindingEnergyDensityFilter.hh
/// @brief  Filter class for looking at the dGbind/dSASA ratio
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_InterfaceBindingEnergyDensityFilter_hh
#define INCLUDED_protocols_simple_filters_InterfaceBindingEnergyDensityFilter_hh

// Unit headers
#include <protocols/simple_filters/InterfaceBindingEnergyDensityFilter.fwd.hh>

// Package headers
#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.fwd.hh>
#include <protocols/simple_filters/DdgFilter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace simple_filters {

class InterfaceBindingEnergyDensityFilter : public filters::Filter
{
public:
	InterfaceBindingEnergyDensityFilter();
	InterfaceBindingEnergyDensityFilter(
		InterfaceSasaFilterOP sasa_filter,
		DdgFilterOP ddG_filter,
		core::Real threshold
	);

	~InterfaceBindingEnergyDensityFilter() override;

	void set_interface_sasa_filter( InterfaceSasaFilterOP sasa_filter );
	void set_ddG_filter( DdgFilterOP ddG_filter );
	void set_upper_threshold( core::Real threshold );

	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const override;
	filters::FilterOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	InterfaceSasaFilterOP sasa_filter_;
	DdgFilterOP ddG_filter_;

	core::Real upper_threshold_;
};

}
}

#endif
