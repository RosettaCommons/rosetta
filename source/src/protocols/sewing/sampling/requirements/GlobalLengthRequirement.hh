// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file GlobalLengthRequirement.hh
///
/// @brief Require that all segments of the specified secondary structure(s) are between the given min and max size
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_sampling_requirements_GlobalLengthRequirement_hh
#define INCLUDED_protocols_sewing_sampling_requirements_GlobalLengthRequirement_hh

//Unit headers
#include <protocols/sewing/sampling/requirements/GlobalLengthRequirement.fwd.hh>
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

class GlobalLengthRequirement : public GlobalRequirement {

public:

	/// @brief default constructor
	GlobalLengthRequirement();

	GlobalLengthRequirement(
		std::set<std::string> valid_dssp_codes,
		core::Size min_length,
		core::Size max_length
	);

	/// @brief Inverse of violated for this requirement
	virtual
	bool
	satisfies(
		AssemblyCOP assembly
	) const;

	/// @brief Has the Assembly removed required residues
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

	std::set<std::string> dssp_codes_;
	core::Size min_length_;
	core::Size max_length_;

};


} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace

#endif
