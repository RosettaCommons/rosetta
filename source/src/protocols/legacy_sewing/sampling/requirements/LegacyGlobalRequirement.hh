// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyGlobalRequirement.hh
///
/// @brief An interface for all LEGACY_SEWING requirements concerning an entire assembly
/// @author Frank Teets

#ifndef INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyGlobalRequirement_hh
#define INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyGlobalRequirement_hh

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalRequirement.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.fwd.hh>

//Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

class LegacyGlobalRequirement : public utility::pointer::ReferenceCount {

public:

	///@brief does the Assembly, as a whole, currently satisfy this requirement?
	//This should return TRUE iff the Assembly as it is now, with no consideration
	//of future Assemblies, satisfies the requirement.
	virtual
	bool
	satisfies(
		AssemblyCOP assembly
	) const = 0;

	///@brief does the Assembly violate this segment? Unlike satisfies, violates
	//implies an irreparable violation, I.E. it should return TRUE iff the Assembly
	//and all possible Assemblies built thereon will violate the Requirement.
	virtual
	bool
	violates(
		AssemblyCOP assembly
	) const = 0;


	///@brief Can we add more edges to the Assembly? Base class implementation
	///returns true, so only implement for requirements that require Assemblies
	///of a specific size
	bool can_be_added_to(
		AssemblyCOP /*assembly*/
	) const {
		return true;
	}

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

};

} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace

#endif
