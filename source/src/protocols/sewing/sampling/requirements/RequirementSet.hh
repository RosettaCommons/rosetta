// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RequirementSet.hh
///
/// @brief A container for all SEWING requirements
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_sampling_requirements_RequirementSet_hh
#define INCLUDED_protocols_sewing_sampling_requirements_RequirementSet_hh

//Unit headers
#include <protocols/sewing/sampling/requirements/RequirementSet.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.fwd.hh>
#include <protocols/sewing/conformation/Model.hh>

#include <protocols/sewing/sampling/requirements/GlobalRequirement.hh>
#include <protocols/sewing/sampling/requirements/IntraSegmentRequirement.hh>

//Utility headers
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

class RequirementSet : public utility::pointer::ReferenceCount {

public:

	typedef std::map< core::Size, utility::vector1<IntraSegmentRequirementOP> > IntraSegmentRequirementsMap;
	//typedef std::map< std::set<core::Size>, utility::vector1<InterSegmentRequirementOP> > InterSegmentRequirementsMap;

	RequirementSet();

	core::Size
	min_segments() const;

	void
	min_segments(core::Size min_segments);

	core::Size
	max_segments() const;

	void
	max_segments(core::Size max_segments);

	void
	add_requirement(
		GlobalRequirementOP requirement
	);

	void
	add_requirement(
		core::Size index,
		IntraSegmentRequirementOP requirement
	);

	// void
	// add_inter_segment_requirement(
	//  std::set<core::Size> segment
	//  InterSegmentRequirementOP inter_segment_requirement
	// );

	///@brief Evaluated if this Assembly satisfies all
	///contained AssemblyRequirements
	bool satisfies(
		AssemblyCOP assembly
	) const;

	///@brief Evaluate if this Assembly violates any contained
	///AssemblyRequirements
	bool violates(
		AssemblyCOP assembly
	) const;

	///@brief Check all Global requirements to see if we
	///can add more edges to this Assembly
	bool can_be_added_to(
		AssemblyCOP assembly
	) const;

	core::Size get_max_segments()
	const;

	void
	show(
		std::ostream & out
	) const;

private:

	core::Size min_segments_;
	core::Size max_segments_;

	//Requirements for the entire Assembly
	utility::vector1<GlobalRequirementOP> global_requirements_;

	//Requirements for individual segments
	IntraSegmentRequirementsMap intra_segment_requirements_;

	//Requirements for interactions between segments
	// utility::vector1< utility::vector <InterSegmentRequirementOP> > inter_segment_requirements_;

};


} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace

#endif
