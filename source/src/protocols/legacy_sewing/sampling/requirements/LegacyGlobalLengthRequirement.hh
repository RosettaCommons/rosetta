// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyGlobalLengthRequirement.hh
///
/// @brief Require that all segments of the specified secondary structure(s) are between the given min and max size
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyGlobalLengthRequirement_hh
#define INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyGlobalLengthRequirement_hh

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalLengthRequirement.fwd.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalRequirement.hh>
#include <core/types.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.fwd.hh>

//Utility headers
#include <set>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

class LegacyGlobalLengthRequirement : public LegacyGlobalRequirement {

public:

	///@brief default constructor
	LegacyGlobalLengthRequirement();

	LegacyGlobalLengthRequirement(
		std::set<std::string> const & valid_dssp_codes,
		core::Size min_length,
		core::Size max_length
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

	static std::string
	class_name();

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );
private:

	std::set<std::string> dssp_codes_;
	core::Size min_length_;
	core::Size max_length_;

};


} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace

#endif
