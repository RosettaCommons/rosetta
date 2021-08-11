// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/filters/HasDisulfideFilter.hh
/// @brief filters based on disulfide presence
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_filters_HasDisulfideFilter_hh
#define INCLUDED_protocols_pose_sewing_filters_HasDisulfideFilter_hh

// Unit headers
#include <protocols/pose_sewing/filters/HasDisulfideFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace pose_sewing {
namespace filters {

///@brief filters based on disulfide presence
class HasDisulfideFilter : public protocols::filters::Filter {

public:
	HasDisulfideFilter();

	// destructor (important for properly forward-declaring smart-pointer members)
	~HasDisulfideFilter() override;

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	void
	set_first_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	void
	set_second_selector( core::select::residue_selector::ResidueSelectorCOP selector );

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
		basic::datacache::DataMap & data
	) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

private:

	core::select::residue_selector::ResidueSelectorCOP first_selector_;
	core::select::residue_selector::ResidueSelectorCOP second_selector_;

};

} //filters
} //pose_sewing
} //protocols

#endif //INCLUDED_protocols_pose_sewing_filters_HasDisulfideFilter_hh
