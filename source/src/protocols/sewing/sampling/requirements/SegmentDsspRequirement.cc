// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SegmentDsspRequirement.cc
///
/// @brief
/// @author Tim Jacobs, Frank Teets

//Unit headers
#include <protocols/sewing/sampling/requirements/SegmentDsspRequirement.hh>
#include <protocols/sewing/sampling/requirements/SegmentDsspRequirementCreator.hh>

//Utility headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <string>
#include <basic/Tracer.hh>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

static basic::Tracer TR("protocols.sewing.sampling.requirements.SegmentDsspRequirement");

//////Creator methods/////
IntraSegmentRequirementOP
SegmentDsspRequirementCreator::create_requirement() const{
	return IntraSegmentRequirementOP( new SegmentDsspRequirement() );
}

std::string
SegmentDsspRequirementCreator::type_name() const {
	return "SegmentDsspRequirement";
}
//////End Creator methods/////

SegmentDsspRequirement::SegmentDsspRequirement(){}

SegmentDsspRequirement::SegmentDsspRequirement(
	std::set<std::string> valid_dssp_codes
):
	valid_dssp_codes_(valid_dssp_codes)
{}

bool
SegmentDsspRequirement::satisfies(
	SewSegment segment
) const {
	return !violates(segment);
}

bool
SegmentDsspRequirement::violates(
	SewSegment segment
) const {
	if ( valid_dssp_codes_.find(std::string(1, segment.dssp_)) == valid_dssp_codes_.end() ) {
		if ( TR.Debug.visible() ) { TR.Debug << "Violation! " << segment.dssp_<< " is invalid!" << std::endl; }
		return true;
	}
	return false;
}


void
SegmentDsspRequirement::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
) {
	if ( !tag->hasOption("dssp") ) {
		utility_exit_with_message("SegmentDsspRequiremnt must be given a whitespace delimited 'dssp' option");
	}
	valid_dssp_codes_ = utility::split_to_set(tag->getOption<std::string>("dssp"));
}

void
SegmentDsspRequirement::show(
	std::ostream & out
) const {
	out << "/////// SegmentDsspRequirement - Segment must have on of the following DSSP codes: ";
	std::set<std::string>::const_iterator it = valid_dssp_codes_.begin();
	std::set<std::string>::const_iterator it_end = valid_dssp_codes_.end();
	for ( ; it != it_end; ++it ) {
		out << *it << " ";
	}
	out << std::endl;
}

} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace
