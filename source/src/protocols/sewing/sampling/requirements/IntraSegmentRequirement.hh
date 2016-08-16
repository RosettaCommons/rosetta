// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IntraSegmentRequirement.hh
///
/// @brief An interface for all SEWING requirements concerning a single SewSegment
/// @author Frank Teets

#ifndef INCLUDED_protocols_sewing_sampling_requirements_IntraSegmentRequirement_hh
#define INCLUDED_protocols_sewing_sampling_requirements_IntraSegmentRequirement_hh

//Unit headers
#include <protocols/sewing/sampling/requirements/IntraSegmentRequirement.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.fwd.hh>
#include <protocols/sewing/conformation/Model.hh>

//Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

class IntraSegmentRequirement : public utility::pointer::ReferenceCount {

public:

	///@brief does the segment satisfy this requirement?
	virtual
	bool
	satisfies(
		SewSegment segment
	) const = 0;

	///@brief does the segment violate this requirement?
	virtual
	bool
	violates(
		SewSegment segment
	) const = 0;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
	) = 0;

	virtual
	void
	show(
		std::ostream & out
	) const = 0;

private:

	//core::Size segment_index_;

};

} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace

#endif
