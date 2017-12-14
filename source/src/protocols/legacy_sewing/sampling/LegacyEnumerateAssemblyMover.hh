// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyEnumerateAssemblyMover.hh
///
/// @brief
/// @author Doonam Kim, Tim Jacobs

#ifndef INCLUDED_devel_legacy_sewing_sampling_LegacyEnumerateAssemblyMover_HH
#define INCLUDED_devel_legacy_sewing_sampling_LegacyEnumerateAssemblyMover_HH

// Unit Headers
#include <protocols/legacy_sewing/sampling/LegacyAssemblyMover.hh>
#include <protocols/legacy_sewing/sampling/LegacyEnumerateAssemblyMover.fwd.hh>
#include <protocols/legacy_sewing/conformation/Assembly.hh> // for accessing segments_ and all_segments_
#include <protocols/legacy_sewing/conformation/Model.hh>

//Protocol headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace legacy_sewing  {

class LegacyEnumerateAssemblyMover : public LegacyAssemblyMover {

private:

	typedef LegacyAssemblyMover parent;

public:

	LegacyEnumerateAssemblyMover();

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;


	AssemblyOP
	generate_assembly() override;

	void
	parse_my_tag(
		TagCOP const tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
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
	core::Real min_assembly_score_;
	std::map< int, Model > models;
};

} // legacy_sewing
} // protocols

#endif
