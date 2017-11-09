// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyAddStartnodeFragments.hh
///
/// @brief A Mover that uses loophash to find fragments that can
/// bridge a gap with minimal modifications to the original pose.
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_legacy_sewing_LegacyAddStartnodeFragments_HH
#define INCLUDED_protocols_legacy_sewing_LegacyAddStartnodeFragments_HH

// Unit Headers
#include <protocols/legacy_sewing/sampling/LegacyAddStartnodeFragments.fwd.hh>
#include <protocols/moves/Mover.hh>

//Protocol headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace legacy_sewing  {

class LegacyAddStartnodeFragments : public protocols::moves::Mover {

public:

	LegacyAddStartnodeFragments();

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;


	void
	apply(
		core::pose::Pose & pose
	) override;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	core::Size start_res_;
	core::Size end_res_;

};

} //legacy_sewing
} //protocols

#endif
