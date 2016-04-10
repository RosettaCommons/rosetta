//G vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;

	std::string
	get_name() const;

	virtual
	AssemblyOP
	generate_assembly();

	//void
	//follow_random_edge_from_node(
	//	AssemblyOP assembly,
	//	ModelNode const * node
	//);

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		core::pose::Pose const & pose
	);

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
