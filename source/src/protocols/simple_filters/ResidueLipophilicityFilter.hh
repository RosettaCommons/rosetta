// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueLipophilicityFilter.hh
/// @brief definition of filter class ResidueLipophilicityFilter.
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_simple_filters_ResidueLipophilicityFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueLipophilicityFilter_hh

#include <protocols/simple_filters/ResidueLipophilicityFilter.fwd.hh>


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

#include <string>
#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

class ResidueLipophilicityFilter : public filters::Filter
{
public:
	ResidueLipophilicityFilter() : filters::Filter( "ResidueLipophilicity" ) {}
	//ResidueLipophilicityFilter( core::Real const threshold );
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const  override {
		return filters::FilterOP( new ResidueLipophilicityFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override {
		return filters::FilterOP( new ResidueLipophilicityFilter() );
	}

	virtual ~ResidueLipophilicityFilter();
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
	core::Real threshold_;
	std::string output_;
	bool print_splines_;

};

}
}
#endif
