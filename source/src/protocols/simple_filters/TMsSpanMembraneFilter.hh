// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/TMsSpanMembraneFilter.hh
/// @brief definition of filter class TMsSpanMembraneFilter.
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_simple_filters_TMsSpanMembraneFilter_hh
#define INCLUDED_protocols_simple_filters_TMsSpanMembraneFilter_hh

#include <protocols/simple_filters/TMsSpanMembraneFilter.fwd.hh>


// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <string>

namespace protocols {
namespace simple_filters {

class TMsSpanMembraneFilter : public filters::Filter
{
public:
	TMsSpanMembraneFilter() : filters::Filter( "TMsSpanMembrane" ) {}
	//TMsSpanMembraneFilter( core::Real const threshold );
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const override {
		return utility::pointer::make_shared< TMsSpanMembraneFilter >( *this );
	}
	filters::FilterOP fresh_instance() const override {
		return utility::pointer::make_shared< TMsSpanMembraneFilter >();
	}

	~TMsSpanMembraneFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Real threshold_ = 0.5;
	core::Real min_distance_ = 20.0;
	std::string output_ = "";
	core::Size flank_ = 1;
	core::Size required_dist_ = 10;
};

}
}
#endif
