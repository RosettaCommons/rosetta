// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/constraint_filters/ShowConstraintsFilter.hh
/// @brief iterate over all constraints on the pose and invoke their show methods
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_constraint_filters_ShowConstraintsFilter_hh
#define INCLUDED_protocols_constraint_filters_ShowConstraintsFilter_hh

#include <protocols/constraint_filters/ShowConstraintsFilter.fwd.hh>


// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace constraint_filters {

class ShowConstraintsFilter : public filters::Filter
{
public:
	ShowConstraintsFilter() : filters::Filter( "ShowConstraints" ) {}

	//ShowConstraintsFilter( core::Real const threshold );

	bool apply( core::pose::Pose const & pose ) const override;

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;

	core::Real report_sm( core::pose::Pose const & pose ) const override;

	core::Real compute( core::pose::Pose const & pose ) const;

	filters::FilterOP clone() const override {
		return utility::pointer::make_shared< ShowConstraintsFilter >( *this );
	}

	filters::FilterOP fresh_instance() const override {
		return utility::pointer::make_shared< ShowConstraintsFilter >();
	}

	~ShowConstraintsFilter() override;

	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
};

}
}
#endif
