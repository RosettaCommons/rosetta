// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/AtomExistsFilter.hh
/// @brief definition of filter class AtomExistsFilter.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_AtomExistsFilter_hh
#define INCLUDED_protocols_drug_design_AtomExistsFilter_hh

//unit headers
#include <protocols/drug_design/AtomExistsFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <string>

namespace protocols {
namespace drug_design {

class AtomExistsFilter : public filters::Filter
{
public:
	AtomExistsFilter():
		Filter( class_name() )
	{}

	AtomExistsFilter( std::string const & residue, std::string atom_name ):
		Filter( class_name() ),
		residue_( residue ),
		atom_name_( atom_name )
	{}

	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new AtomExistsFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override {
		return filters::FilterOP( new AtomExistsFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const &pose ) const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string residue_;
	std::string atom_name_;
};

}
}

#endif
