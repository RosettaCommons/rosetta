// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/calc_taskop_filters/OperatorFilter.hh

#ifndef INCLUDED_protocols_simple_filters_OperatorFilter_hh
#define INCLUDED_protocols_simple_filters_OperatorFilter_hh

//unit headers
#include <protocols/calc_taskop_filters/OperatorFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
namespace protocols {
namespace calc_taskop_filters {

enum Operation { SUM, PRODUCT, NORMALIZED_SUM, MAX, MIN, SUBTRACT, ABS, BOOLEAN_OR/*x+y-xy*/,XOR };
/// @brief simply take a list of filters and combine them using the operation above
class Operator : public filters::Filter
{
public:
	Operator();
	~Operator() override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new Operator( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new Operator() );
	}

	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & ) override;
	core::Real compute( core::pose::Pose const & pose ) const;
	utility::vector1< protocols::filters::FilterOP > filters() const;
	void add_filter( protocols::filters::FilterOP f );
	void reset_baseline( core::pose::Pose const & pose, bool const attempt_read_from_checkpoint/* see Sigmoid for details*/ ); /// goes over Sigmoid filters and resets them. Note that this is nonconst, and cannot be called from apply
	core::Real threshold() const{ return threshold_; }
	void threshold( core::Real const t ){ threshold_ = t; }
	Operation operation() const{ return operation_; }
	void operation( Operation const o ){ operation_ = o; }
	void negate( bool const b ){ negate_ = b; }
	bool negate() const{ return negate_; }
	utility::vector1< std::string > relative_pose_names() { return relative_pose_names_; }
	void relative_pose_names( utility::vector1< std::string > const & s ) { relative_pose_names_ = s; }
	bool multi_relative() const { return multi_relative_; }
	void multi_relative( bool const m ){ multi_relative_ = m; }
	void modify_relative_filters_pdb_names();
	bool logarithm() const{ return logarithm_; } //getter
	void logarithm( bool const b ){ logarithm_ = b; } //setter

	void report_subvalues( bool const report ){ report_subvalues_ = report; }
	bool report_subvalues() const { return report_subvalues_; }

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static
	void attributes( utility::tag::AttributeList & attlist );

private:
	utility::vector1< protocols::filters::FilterOP > filters_;
	Operation operation_; // dflt PRODUCT
	core::Real threshold_; // dflt 0
	bool negate_; // dflt false; in optimization, useful to get values between -1 - 0 rather than 0-1
	utility::vector1< std::string > relative_pose_names_; // dflt ""; see below
	bool multi_relative_; //dflt false; if true, searches all of the filters for RelativePoseFilters, replicates them to as many different file names as are listed in relative_pose_names_. Useful in case there are many different states that are all taken into consideration using the same operator
	bool logarithm_; //dflt false; if true, computes the logarithm of the operator's value (10^-9 = -9)
	bool report_subvalues_; //dflt false; report each of the subfilters
};
}
}

#endif
