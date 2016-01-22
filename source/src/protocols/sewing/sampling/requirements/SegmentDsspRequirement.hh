// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SegmentDsspRequirement.hh
///
/// @brief A requirement on the seconday structure of a given SewSegment
/// @author Tim Jacobs, Frank Teets

#ifndef INCLUDED_protocols_sewing_sampling_requirements_SegmentDsspRequirement_hh
#define INCLUDED_protocols_sewing_sampling_requirements_SegmentDsspRequirement_hh

//Unit headers
#include <protocols/sewing/sampling/requirements/SegmentDsspRequirement.fwd.hh>
#include <protocols/sewing/sampling/requirements/IntraSegmentRequirement.hh>

//Package headers
#include <protocols/sewing/conformation/Model.hh>

//Utility headers
#include <set>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

class SegmentDsspRequirement : public IntraSegmentRequirement {

public:

	/// @brief default constructor
	SegmentDsspRequirement();

	/// @brief explicit constructor
	SegmentDsspRequirement(
		std::set<std::string> valid_dssp_codes
	);

	/// @brief Does the segment have a valid DSSP code?
	virtual
	bool
	satisfies(
		SewSegment segment
	) const;

	/// @brief Does the segment have an invalid DSSP code?
	virtual
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


private:

	std::set<std::string> valid_dssp_codes_;

};

} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace

#endif
