// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacySegmentLengthRequirement.hh
///
/// @brief A requirement on the total number of residues for a given SewSegment
/// @author Tim Jacobs, Frank Teets

#ifndef INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacySegmentLengthRequirement_hh
#define INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacySegmentLengthRequirement_hh

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacySegmentLengthRequirement.fwd.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyIntraSegmentRequirement.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Model.hh>

//Utility headers
#include <core/types.hh>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

class LegacySegmentLengthRequirement : public LegacyIntraSegmentRequirement {

public:

	///@brief default constructor
	LegacySegmentLengthRequirement();

	LegacySegmentLengthRequirement(
		core::Size min_length,
		core::Size max_length
	);

	///@brief Do we have the right number of residues in the SewSegment
	bool
	satisfies(
		SewSegment segment
	) const;

	///@brief Can we, through more edge additions, get the right number of residues?
	bool
	violates(
		SewSegment segment
	) const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
	);

	virtual
	void
	show(
		std::ostream & out
	) const;

	static std::string
	class_name();

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );


private:

	core::Size min_length_;
	core::Size max_length_;

};

} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace

#endif
