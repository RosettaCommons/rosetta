// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/LongestContinuousApolarSegmentFilter.hh
/// @brief This filter computes the longest continuous stretch of polar residues within a pose or selection.
/// @author Yang Hsia (yhsia@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_LongestContinuousApolarSegmentFilter_hh
#define INCLUDED_protocols_simple_filters_LongestContinuousApolarSegmentFilter_hh

// Unit headers
#include <protocols/simple_filters/LongestContinuousApolarSegmentFilter.fwd.hh>
#include <protocols/simple_filters/LongestContinuousPolarSegmentFilter.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace simple_filters {

///@brief This filter computes the longest continuous stretch of polar residues within a pose or selection.
class LongestContinuousApolarSegmentFilter : public LongestContinuousPolarSegmentFilter {

public:
	/// @brief Constructor.
	LongestContinuousApolarSegmentFilter();

	/// @brief destructor.
	~LongestContinuousApolarSegmentFilter();

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

	/// @brief Given a residue type, determine whether it's one of the types that this filter
	/// should count.
	/// @details Based on whether the residue type has the APOLAR property.
	bool is_counted( core::chemical::ResidueType const &restype ) const override;

	/// @brief returns type of counted residues (apolar)
	std::string counted_residue_description() const override;
};

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_LongestContinuousApolarSegmentFilter_hh
