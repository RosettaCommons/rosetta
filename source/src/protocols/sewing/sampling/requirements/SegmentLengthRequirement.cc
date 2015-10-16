// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SegmentLengthRequirement.cc
///
/// @brief
/// @author Tim Jacobs, Frank Teets

//Unit headers
#include <protocols/sewing/sampling/requirements/SegmentLengthRequirement.hh>
#include <protocols/sewing/sampling/requirements/SegmentLengthRequirementCreator.hh>

//Utility headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

static basic::Tracer TR("protocols.sewing.sampling.requirements.SegmentLengthRequirement");

//////Creator methods/////
IntraSegmentRequirementOP
SegmentLengthRequirementCreator::create_requirement() const{
	return IntraSegmentRequirementOP( new SegmentLengthRequirement() );
}

std::string
SegmentLengthRequirementCreator::type_name() const {
	return "SegmentLengthRequirement";
}
//////End Creator methods/////

SegmentLengthRequirement::SegmentLengthRequirement():
	min_length_(1),
	max_length_(1)
{}

SegmentLengthRequirement::SegmentLengthRequirement(
	core::Size min_length,
	core::Size max_length
):
	min_length_(min_length),
	max_length_(max_length)
{}

//@brief segments do not change; if the segment does not violate the requirement it must satisfy it
bool
SegmentLengthRequirement::satisfies(
	SewSegment segment
) const {
	return !violates(segment);
}

//@brief check if the segment falls outside the range parameters
bool
SegmentLengthRequirement::violates(
	SewSegment segment
) const {
	core::Size length = segment.residues_.size();
	if ( length < min_length_ ) {
		if ( TR.Debug.visible() ) { TR.Debug << length << " Violation! Too Short!" << std::endl; }
		return true;
	} else if ( length > max_length_ ) {
		if ( TR.Debug.visible() ) { TR.Debug << length << " Violation! Too Long!" << std::endl; }
		return true;
	}
	return false;
}

void
SegmentLengthRequirement::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
) {
	if ( !tag->hasOption("min_length") && !tag->hasOption("max_length") ) {
		utility_exit_with_message("You must provide either max_length, min_length, or both to the SegmentLengthRequirement tag!");
	} else if ( tag->hasOption("min_length") && !tag->hasOption("max_length") ) {
		min_length_ = tag->getOption<core::Size>("min_length");
		max_length_ = -1;/*max size*/
	} else if ( !tag->hasOption("min_length") && tag->hasOption("max_length") ) {
		max_length_ = tag->getOption<core::Size>("max_length");
		min_length_ = 0;
	} else {
		min_length_ = tag->getOption<core::Size>("min_length");
		max_length_ = tag->getOption<core::Size>("max_length");
	}
}

void
SegmentLengthRequirement::show(
	std::ostream & out
) const {
	out << "/////// SegmentLengthRequirement - Segment must contain between "
		<< min_length_ << " and " << max_length_ << " residues" << std::endl;
}

} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace
