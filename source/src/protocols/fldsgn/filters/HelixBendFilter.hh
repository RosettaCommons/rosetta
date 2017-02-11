// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/filters/HelixBendFilter.hh
/// @brief Filter used in 'Principles for designing proteins with cavities formed by curved b-sheets' to control helix geometry.
/// @author Benjamin Basanta (basantab@uw.edu)

#ifndef INCLUDED_protocols_fldsgn_filters_HelixBendFilter_hh
#define INCLUDED_protocols_fldsgn_filters_HelixBendFilter_hh

// Unit headers
#include <protocols/fldsgn/filters/HelixBendFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace fldsgn {
namespace filters {

///@brief Filter used in 'Principles for designing proteins with cavities formed by curved b-sheets' to control helix geometry.
class HelixBendFilter : public protocols::filters::Filter {

public:
	HelixBendFilter();

	// destructor (important for properly forward-declaring smart-pointer members)
	~HelixBendFilter() override;

	// @brief set secondary structure
	void
	secstruct( std::string const & ss );

	/// @brief minimum angle for filtering
	void
	threshold( core::Real const r );

	/// @brief Strand id number
	void
	helix_id( core::Size const r );

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

public:
	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

private:
	mutable std::string secstruct_ ;
	core::Real threshold_ ;
	core::Size helix_id_ ;
	mutable bool filter_status_;
};

} //protocols
} //fldsgn
} //filters

#endif //INCLUDED_protocols_fldsgn_filters_HelixBendFilter_hh
