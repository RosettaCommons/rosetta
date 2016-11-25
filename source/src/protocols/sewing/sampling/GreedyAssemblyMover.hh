// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file GreedyAssemblyMover.hh
///
/// @brief This Mover generates assemblies by taking the best edge possible for each addition. This process is repeated
/// cycles times and the best overall Assembly that was found is returned
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_sewing_sampling_GreedyAssemblyMover_HH
#define INCLUDED_protocols_sewing_sampling_GreedyAssemblyMover_HH

// Unit Headers
#include <protocols/sewing/sampling/GreedyAssemblyMover.fwd.hh>
#include <protocols/sewing/sampling/AssemblyMover.hh>

//Protocol headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace sewing  {

class GreedyAssemblyMover : public AssemblyMover {

private:

	typedef AssemblyMover parent;

public:

	GreedyAssemblyMover();

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;

	// XRW TEMP  std::string
	// XRW TEMP  get_name() const;

	AssemblyOP
	generate_assembly() override;

	//void
	//follow_random_edge_from_node(
	// AssemblyOP assembly,
	// ModelNode const * node
	//);

	void
	parse_my_tag(
		TagCOP tag,
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


	// virtual
	// void
	// output_stats(
	//  AssemblyOP const & assembly,
	//  core::pose::Pose & pose
	// );

private:

	//Best complete assembly seen throughout the simulation
	AssemblyOP best_complete_assembly_;
	core::Real best_score_;

	core::Size cycles_;
	core::Size max_edges_per_node_;

};

} //sewing
} //protocols

#endif
