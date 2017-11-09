// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyAssemblyMover.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_legacy_sewing_sampling_LegacyAssemblyMover_HH
#define INCLUDED_protocols_legacy_sewing_sampling_LegacyAssemblyMover_HH

// Unit Headers
#include <protocols/legacy_sewing/sampling/LegacyAssemblyMover.fwd.hh>

// Package Headers
#include <protocols/legacy_sewing/conformation/Assembly.fwd.hh>
#include <protocols/legacy_sewing/sampling/SewGraph.fwd.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementSet.fwd.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementFactory.fwd.hh>
#include <protocols/legacy_sewing/scoring/LegacyAssemblyScorer.fwd.hh>

// Core Headers
#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>

#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace legacy_sewing  {

class LegacyAssemblyMover : public protocols::moves::Mover {

public:

	LegacyAssemblyMover();

	std::string
	get_name() const;

	virtual
	void
	apply(
		core::pose::Pose & pose
	);

	virtual
	void
	add_starting_model(
		AssemblyOP assembly
	) const;

	virtual
	core::pose::Pose
	get_fullatom_pose(
		AssemblyOP assembly
	) const;

	virtual
	core::pose::Pose
	get_centroid_pose(
		AssemblyOP assembly
	) const;

	bool
	follow_random_edge_from_node(
		AssemblyOP assembly,
		ModelNode const * reference_node
	) const;

	virtual
	AssemblyOP
	generate_assembly() = 0;

	virtual
	core::pose::Pose
	refine_assembly(
		AssemblyOP & assembly
	);

	virtual
	bool
	reinitialize_for_new_input() const;

	void
	append_movie_frame(
		AssemblyOP assembly,
		core::Size cycle
	) const;

	///@brief add several outputs to the score file
	///Pose can't be const due to scoring.
	virtual
	void
	output_stats(
		AssemblyOP const & assembly,
		core::pose::Pose & pose,
		std::string indicator_whether_from_LegacyMonteCarloAssemblyMover_or_LegacyExhaustiveAssemblyMover
	);

	virtual
	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		core::pose::Pose const & pose
	);

	virtual
	void
	parse_requirements(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	virtual
	void
	parse_global_requirements(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	virtual
	void
	parse_intra_segment_requirements(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);
	static std::string
	mover_name();

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static utility::tag::XMLSchemaComplexTypeGeneratorOP
	define_assembly_mover_ct_gen( utility::tag::XMLSchemaDefinition & xsd );

	static std::string
	legacy_assembly_mover_ct_namer( std::string const & );

protected:

	scoring::LegacyAssemblyScoreFunctionOP assembly_scorefxn_;
	SewGraphOP graph_;
	std::string edge_file_;

	//A container for all requirements enforced on the Assembly
	//and the components that are added to the Assembly.
	sampling::requirements::LegacyRequirementSetOP requirement_set_;

	sampling::requirements::LegacyRequirementFactory * requirement_factory_;

	//TODO make this accessible with RosettaScripts
	core::scoring::ScoreFunctionOP cen_scorefxn_;
	core::scoring::ScoreFunctionOP fa_scorefxn_;

	core::Real base_native_bonus_;
	core::Size neighbor_cutoff_;
};

} //legacy_sewing
} //protocols

#endif
