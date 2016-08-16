// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResidueRetentionRequirement.hh
///
/// @brief A container for all SEWING requirements
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_sampling_requirements_ResidueRetentionRequirement_hh
#define INCLUDED_protocols_sewing_sampling_requirements_ResidueRetentionRequirement_hh

//Unit headers
#include <protocols/sewing/sampling/requirements/ResidueRetentionRequirement.fwd.hh>
#include <protocols/sewing/sampling/requirements/GlobalRequirement.hh>
#include <core/types.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.fwd.hh>

//Utility headers
#include <set>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

class ResidueRetentionRequirement : public GlobalRequirement {

public:

	///@brief default constructor
	ResidueRetentionRequirement();

	ResidueRetentionRequirement(
		int model_id
	);

	void
	model_id(
		int model_id
	);

	void
	required_resnums(
		std::set<core::Size> required_resnums
	);

	void
	add_resnum(
		core::Size resnum
	);

	///@brief Inverse of violated for this requirement
	virtual
	bool
	satisfies(
		AssemblyCOP assembly
	) const;

	///@brief Has the Assembly removed required residues
	///for the specified model
	virtual
	bool
	violates(
		AssemblyCOP assembly
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

	int model_id_;
	std::set<core::Size> required_resnums_;

};


} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace

#endif
