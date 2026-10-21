// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/RDKitMetricFilter.hh
/// @brief definition of filter class RDKitMetricFilter.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_RDKitMetricFilter_hh
#define INCLUDED_protocols_drug_design_RDKitMetricFilter_hh

//unit headers
#include <protocols/drug_design/RDKitMetricFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <string>

namespace protocols {
namespace drug_design {

class RDKitMetricFilter : public filters::Filter
{
public:
	RDKitMetricFilter():
		Filter( class_name() ),
		lower_threshold_( -9999 ),
		upper_threshold_(  9999 )
	{}

	RDKitMetricFilter( std::string const & residue, std::string const & metric_name, core::Real lower, core::Real upper ):
		Filter( class_name() ),
		residue_( residue ),
		lower_threshold_( lower ),
		upper_threshold_( upper )
	{
		metric( metric_name );
	}

	void residue( std::string const & residue ) { residue_ = residue; }
	std::string const & residue() const { return residue_; }

	void lower( core::Real lower ) { lower_threshold_ = lower; }
	core::Real lower() const { return lower_threshold_; }

	void upper( core::Real upper ) { upper_threshold_ = upper; }
	core::Real upper() const { return upper_threshold_; }

	void metric( std::string const & setting );
	std::string const & metric() const { return metric_; }

	bool apply( core::pose::Pose const & pose ) const override;

	filters::FilterOP clone() const override {
		return filters::FilterOP( new RDKitMetricFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override {
		return filters::FilterOP( new RDKitMetricFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const &pose ) const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string residue_;
	std::string metric_;
	core::Real lower_threshold_;
	core::Real upper_threshold_;
};

}
}

#endif
