// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RepeatAssemblyMover.hh
///
/// @brief An AssemblyMover specific for repeating backbones. This Mover is currently pretty hacky and makes
/// assumes that every node is exactly 3 segments long, and that the assembly is a continuous assembly
/// @author Tim Jacobs


#ifndef INCLUDED_protocol_sewing_sampling_RepeatAssemblyMover_HH
#define INCLUDED_protocol_sewing_sampling_RepeatAssemblyMover_HH

// Unit Headers
#include <protocols/sewing/sampling/RepeatAssemblyMover.fwd.hh>
#include <protocols/sewing/sampling/GivenPathAssemblyMover.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.hh>
#include <protocols/sewing/sampling/SewGraph.hh>

//Protocol headers
#include <core/pose/Pose.hh>

#include <protocols/sewing/scoring/AssemblyScorer.hh>
#include <protocols/sewing/sampling/requirements/RequirementSet.hh>

//Utility headers
#include <numeric/xyzTransform.fwd.hh>

//Unit headers
#include <protocols/sewing/sampling/RepeatAssemblyMover.fwd.hh>


namespace protocols {
namespace sewing  {

class RepeatAssemblyMover : public AssemblyMover {

public:

	typedef AssemblyMover parent;

	RepeatAssemblyMover();

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;


	std::pair< bool, AssemblyOP >
	dfs_cycle_finder(
		ModelNode const * reference_node,
		utility::vector1<ModelNode const *> visited_nodes,
		AssemblyOP assembly
	) const;

	AssemblyOP
	generate_assembly() override;

	///adds single repeat unit data to score file
	void
	output_base_stats(
		AssemblyOP const & assembly
	);

	///adds remaining data to score file
	//void
	//output_stats(
	//        AssemblyOP const & assembly,
	//        core::pose::Pose & pose
	//);

	void
	add_repeats(
		AssemblyOP assembly
	) const;

	utility::vector1< numeric::xyzVector<core::Real> >
	get_segment_coords(
		SewSegment const & segment
	) const;

	void
	transform_segment_coords(
		SewSegment & segment,
		numeric::xyzTransform<core::Real> transformer,
		numeric::xyzVector<float> com,
		core::Size n_transformations
	) const;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		core::pose::Pose const & pose
	) override;

	core::pose::Pose
	refine_assembly(
		AssemblyOP & assembly
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

	core::Size num_repeating_segments_;
	core::Real threshold_score_of_complete_cycle_;
	core::Size num_repeats_;
	core::Size num_trials_;
	scoring::AssemblyScoreFunctionOP clash_scorefxn_;

};

} //sewing
} //protocols

#endif
