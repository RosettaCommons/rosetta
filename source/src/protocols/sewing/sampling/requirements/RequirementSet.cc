// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RequirementSet.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/sampling/requirements/RequirementSet.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.hh>

//Utility headers
#include <basic/Tracer.hh>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

static basic::Tracer TR("protocols.sewing.sampling.requirements.RequirementSet");

RequirementSet::RequirementSet():
	min_segments_(0),
	max_segments_(0)
{}

core::Size
RequirementSet::min_segments() const {
	return min_segments_;
}

void
RequirementSet::min_segments(
	core::Size min_segments
){
	min_segments_ = min_segments;
}

core::Size
RequirementSet::max_segments() const {
	return max_segments_;
}

void
RequirementSet::max_segments(
	core::Size max_segments
){
	max_segments_ = max_segments;
}

void
RequirementSet::add_requirement(
	GlobalRequirementOP requirement
){
	global_requirements_.push_back(requirement);
}

void
RequirementSet::add_requirement(
	core::Size index,
	IntraSegmentRequirementOP requirement
){
	intra_segment_requirements_[index].push_back(requirement);
}

bool
RequirementSet::can_be_added_to(
	AssemblyCOP assembly
) const {
	if ( assembly->segments().size() < max_segments_ ) {
		return true;
	}
	return false;
}

core::Size
RequirementSet::get_max_segments()
const {
	return max_segments_;
}


bool
RequirementSet::satisfies(
	AssemblyCOP assembly
) const {

	if ( assembly->segments().size() < min_segments_ || assembly->segments().size() > max_segments_ ) {
		return false;
	}

	for ( core::Size i=1; i<=global_requirements_.size(); ++i ) {
		if ( !global_requirements_[i]->satisfies(assembly) ) {
			return false;
		}
	}

	utility::vector1<SewSegment> segments = assembly->segments();
	IntraSegmentRequirementsMap::const_iterator intra_it = intra_segment_requirements_.begin();
	IntraSegmentRequirementsMap::const_iterator intra_it_end = intra_segment_requirements_.end();

	for ( ; intra_it != intra_it_end; ++intra_it ) {
		if ( intra_it->first > segments.size() ) {
			return false;
		}
		SewSegment cur_segment = segments[intra_it->first];
		utility::vector1<IntraSegmentRequirementOP> requirements = intra_it->second;
		for ( core::Size i=1; i<=requirements.size(); ++i ) {
			if ( !requirements[i]->satisfies(cur_segment) ) {
				return false;
			}
		}
	}

	return true;
}


///@details This function should return true if any of the requirements contained
///in this set are guaranteed to fail in the case of future edge additions. This is
///implemented by considering a 'hypothetical' assembly in which segments which do not
///yet exist auto-pass all requirements.
bool
RequirementSet::violates(
	AssemblyCOP assembly
) const {

	for ( core::Size i=1; i<=global_requirements_.size(); ++i ) {
		if ( global_requirements_[i]->violates(assembly) ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "Failing global requirements!" << std::endl;
			}
			return true;
		}
	}

	utility::vector1<SewSegment> hypothetical_assembly;
	utility::vector1<SewSegment> segments = assembly->segments();

	core::Size n_hypothetical_segments;
	if ( max_segments_ >= assembly->segments().size() ) {
		n_hypothetical_segments = max_segments_ - assembly->segments().size();
	} else {
		n_hypothetical_segments = 0;
	}

	//A hypothetical segment (model_id: 0, segment_id: 0) that will automatically pass
	//all requirements
	SewSegment hypothetical_segment;

	//Create the hypthetical assembly, which is the true assembly flanked on both the N and C terminus
	//by the number of hypothetical segments that would give an assembly that is the maximum size
	for ( core::Size i=1; i<=n_hypothetical_segments; ++i ) {
		hypothetical_assembly.push_back(hypothetical_segment);
	}
	hypothetical_assembly.insert(hypothetical_assembly.end(), segments.begin(), segments.end());
	for ( core::Size i=1; i<=n_hypothetical_segments; ++i ) {
		hypothetical_assembly.push_back(hypothetical_segment);
	}


	//Go through all windows of lowest_max_segments in our hypothetical assembly. If any
	//of these windows fail to introduce a violation, then return false. If all windows
	//violate then return true.
	bool violates = false;
	for ( core::Size offset = 0; offset <= n_hypothetical_segments; ++offset ) {

		//The current window is not yet violating
		violates = false;
		for ( core::Size i=1; i <= max_segments_; ++i ) {

			core::Size seg_index = i + offset;
			if ( hypothetical_assembly[seg_index] == hypothetical_segment ) {
				continue;
			}

			//Get the requirements for this segment of the hypothetical assembly
			IntraSegmentRequirementsMap::const_iterator intra_reqs_it = intra_segment_requirements_.find(i);

			//If there are no requirements for this index, then it passes
			if ( intra_reqs_it == intra_segment_requirements_.end() ) {
				continue;
			}

			//If this is a real segment that has requirements, check them
			utility::vector1<IntraSegmentRequirementOP> segment_reqs = intra_reqs_it->second;
			for ( core::Size req_i=1; req_i<=segment_reqs.size(); ++req_i ) {
				if ( segment_reqs[req_i]->violates(hypothetical_assembly[seg_index]) ) {
					if ( TR.Debug.visible() ) {
						TR.Debug << "Failing Intra-segment requirement " << req_i<< " on segment " << seg_index << "!" << std::endl;
					}
					violates = true;
					break;
				}
			}

			//If any segment in the window violates, the entire window
			//violates and we can exit early
			if ( violates ) {
				break;
			}
		}

		//We've found a window that has to potential to be
		//valid in the future, and are therefore not violating the requirements
		if ( !violates ) {
			return false;
		}
	}

	//If we were unable to find a window that doesn't violate, then we
	//have an assembly that violates our set
	if ( TR.Debug.visible() ) {
		TR.Debug << "Failing Intra-segment requirements!" << std::endl;
	}
	return true;
}


void
RequirementSet::show(
	std::ostream & out
) const {
	out << "/////////// SEWING RequirementSet /////////////" << std::endl;

	out << "/////// Global Assembly Requirements" << std::endl;
	utility::vector1<GlobalRequirementOP>::const_iterator assem_it = global_requirements_.begin();
	utility::vector1<GlobalRequirementOP>::const_iterator assem_it_end = global_requirements_.end();
	for ( ; assem_it != assem_it_end; ++assem_it ) {
		(*assem_it)->show(out);
	}

	out << "/////// Intra Segment Requirements" << std::endl;
	IntraSegmentRequirementsMap::const_iterator intra_it =
		intra_segment_requirements_.begin();
	IntraSegmentRequirementsMap::const_iterator intra_it_end =
		intra_segment_requirements_.end();
	for ( ; intra_it != intra_it_end; ++intra_it ) {
		out << "/////// Segment " << intra_it->first << ":" << std::endl;
		for ( core::Size i=1; i<=intra_it->second.size(); ++i ) {
			intra_it->second[i]->show(out);
		}
	}

	out << "///////// End SEWING RequirementSet ///////////" << std::endl;
}

} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace
