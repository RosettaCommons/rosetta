// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/LigandLocationFilter.hh
/// @brief A filter to detect whether a ligand is inside a pocket
/// @author Yidan Tang (yidan.tang@vanderbilt.edu)

#ifndef INCLUDED_protocols_drug_design_LigandLocationFilter_hh
#define INCLUDED_protocols_drug_design_LigandLocationFilter_hh

//unit headers
#include <protocols/drug_design/LigandLocationFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <numeric/xyzVector.hh>

#include <string>

namespace protocols {
namespace drug_design {

class LigandLocationFilter : public filters::Filter
{
public:
	LigandLocationFilter():
		Filter( class_name() ),
		chain_(),
		radius_(-1),
		center_(0.0)
	{}

	virtual ~LigandLocationFilter() {}

	bool apply( core::pose::Pose const & pose ) const;

	filters::FilterOP clone() const {
		return filters::FilterOP( new LigandLocationFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new LigandLocationFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute_distance( core::pose::Pose const &pose ) const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & );

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string chain_;
	core::Real radius_;
	core::Vector center_;
};

}
}

#endif
