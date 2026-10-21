// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/HeteroHeteroBondFilter.hh
/// @brief definition of filter class HeteroHeteroBondFilter.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_HeteroHeteroBondFilter_hh
#define INCLUDED_protocols_drug_design_HeteroHeteroBondFilter_hh

//unit headers
#include <protocols/drug_design/HeteroHeteroBondFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <string>

namespace protocols {
namespace drug_design {

class HeteroHeteroBondFilter : public filters::Filter
{
public:
	HeteroHeteroBondFilter():
		Filter( class_name() ),
		threshold_( 9999 ),
		exclude_ates_( true )
	{}

	HeteroHeteroBondFilter( std::string const & residue, core::Size threshold, bool exclude_ates = true ):
		Filter( class_name() ),
		residue_( residue ),
		threshold_( threshold ),
		exclude_ates_( exclude_ates )
	{}

	void residue( std::string const & residue ) { residue_ = residue; }
	std::string const & residue() const { return residue_; }

	void threshold( core::Size setting ) { threshold_ = setting; }
	core::Size threshold() const { return threshold_; }

	void exclude_ates( bool setting ) { exclude_ates_ = setting; }
	bool exclude_ates() const { return exclude_ates_; }

	bool apply( core::pose::Pose const & pose ) const override;

	filters::FilterOP clone() const override {
		return filters::FilterOP( new HeteroHeteroBondFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override {
		return filters::FilterOP( new HeteroHeteroBondFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Size compute( core::pose::Pose const &pose ) const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string residue_;
	core::Size threshold_;
	bool exclude_ates_;
};

}
}

#endif
