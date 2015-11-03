// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AssemblyMover.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_sewing_sampling_AssemblyMover_HH
#define INCLUDED_protocols_sewing_sampling_AssemblyMover_HH

// Unit Headers
#include <protocols/sewing/sampling/AssemblyMover.fwd.hh>

// Package Headers
#include <protocols/sewing/conformation/Assembly.fwd.hh>
#include <protocols/sewing/sampling/SewGraph.fwd.hh>
#include <protocols/sewing/sampling/requirements/RequirementSet.fwd.hh>
#include <protocols/sewing/sampling/requirements/RequirementFactory.fwd.hh>
#include <protocols/sewing/scoring/AssemblyScorer.fwd.hh>

// Core Headers
#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>

#include <core/pose/Pose.hh>

// Utility Headers

namespace protocols {
namespace sewing  {

class AssemblyMover : public protocols::moves::Mover {

public:

	AssemblyMover();

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
		std::string indicator_whether_from_MonteCarloAssemblyMover_or_ExhaustiveAssemblyMover
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

protected:

	scoring::AssemblyScoreFunctionOP assembly_scorefxn_;
	SewGraphOP graph_;
	std::string edge_file_;

	//A container for all requirements enforced on the Assembly
	//and the components that are added to the Assembly.
	sampling::requirements::RequirementSetOP requirement_set_;

	sampling::requirements::RequirementFactory * requirement_factory_;

	//TODO make this accessible with RosettaScripts
	core::scoring::ScoreFunctionOP cen_scorefxn_;
	core::scoring::ScoreFunctionOP fa_scorefxn_;

	core::Real base_native_bonus_;
	core::Size neighbor_cutoff_;
};

} //sewing
} //protocols

#endif
