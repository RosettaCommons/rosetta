// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyAppendAssemblyMover.hh
///
/// @brief A Mover that uses loophash to find fragments that can
/// bridge a gap with minimal modifications to the original pose.
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_legacy_sewing_LegacyAppendAssemblyMover_HH
#define INCLUDED_protocols_legacy_sewing_LegacyAppendAssemblyMover_HH

// Unit Headers
#include <protocols/legacy_sewing/sampling/LegacyAppendAssemblyMover.fwd.hh>
#include <protocols/legacy_sewing/sampling/LegacyMonteCarloAssemblyMover.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>

//Protocol headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace legacy_sewing  {

class LegacyAppendAssemblyMover : public LegacyMonteCarloAssemblyMover {

private:

	typedef LegacyMonteCarloAssemblyMover parent;

public:

	LegacyAppendAssemblyMover();

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;


	///@brief The starting node for the Assembly from the append mover
	///will always be the input PDB
	void
	add_starting_model(
		AssemblyOP assembly
	) const override;

	core::pose::Pose
	refine_assembly(
		AssemblyOP & assembly
	) override;

	void
	hash_pdb_model(
		Model const & pdb_model
	);

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

	bool initialized_;

	core::pose::PoseOP partner_pose_;

	std::map<int, Model> models_;
	utility::vector1<core::Size> pose_segment_starts_;
	utility::vector1<core::Size> pose_segment_ends_;
	utility::vector1<core::Size> match_segments_;
};

} //legacy_sewing
} //protocols

#endif
