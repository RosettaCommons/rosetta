// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/constraint_filters/ShowConstraintsFilter.cc
/// @brief iterate over all constraints on the pose and invoke their show methods
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
// Project Headers

#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/constraint_filters/ShowConstraintsFilter.hh>
#include <protocols/constraint_filters/ShowConstraintsFilterCreator.hh>
#include <string>
#include <utility/tag/Tag.fwd.hh>

#include <fstream>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/scoring/constraints/Constraint.hh> // AUTO IWYU For Constraint

namespace protocols {
namespace constraint_filters {

static basic::Tracer TR( "protocols.constraint_filters.ShowConstraintsFilter" );

protocols::filters::FilterOP
ShowConstraintsFilterCreator::create_filter() const { return utility::pointer::make_shared< ShowConstraintsFilter >(); }

std::string
ShowConstraintsFilterCreator::keyname() const { return "ShowConstraints"; }

ShowConstraintsFilter::~ShowConstraintsFilter(){}

void
ShowConstraintsFilter::parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap & )
{
}

bool
ShowConstraintsFilter::apply( core::pose::Pose const & pose ) const {
	core::Real result = compute( pose );
	return result;
}

void
ShowConstraintsFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real result = compute( pose );
	out << "the results is " << result << std::endl;
}

core::Real
ShowConstraintsFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real result = compute( pose );
	return( result );
}

core::Real
ShowConstraintsFilter::compute( core::pose::Pose const & pose ) const {
	core::scoring::constraints::ConstraintCOPs csts = pose.constraint_set()->get_all_constraints();
	TR << "Showing all constraints" << std::endl;
	for ( auto cst : csts ) {
		TR << "==================================" << std::endl;
		cst->show( TR );
		TR << "==================================" << std::endl;
	}
	return 1;
}

std::string ShowConstraintsFilter::name() const {
	return class_name();
}

void ShowConstraintsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "leave 0.5 for boolean behaviour" , "0.5" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "iterate over all constraints in the pose, and invoke their show method", attlist );
}

void ShowConstraintsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ShowConstraintsFilter::provide_xml_schema( xsd );
}

std::string ShowConstraintsFilter::class_name() {
	return "ShowConstraints";
}

}
}
