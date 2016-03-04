// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MonteCarloAssemblyMover.hh
///
/// @brief Assembly substructures by MonteCarlo way
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_sewing_sampling_MonteCarloAssemblyMover_HH
#define INCLUDED_protocols_sewing_sampling_MonteCarloAssemblyMover_HH

// Unit Headers
#include <protocols/sewing/sampling/MonteCarloAssemblyMover.fwd.hh>
#include <protocols/sewing/sampling/AssemblyMover.hh>

//Protocol headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace sewing  {

enum move_type {
	ADD_EDGE = 1,
	DELETE_EDGE,
	SWITCH_EDGE
};

class MonteCarloAssemblyMover : public AssemblyMover {

private:

	typedef AssemblyMover parent;

public:

	MonteCarloAssemblyMover();

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;

	std::string
	get_name() const;

	virtual
	AssemblyOP
	generate_assembly();

	///@brief Add a new edge to the Assembly. If this fails for any reason then
	///revert to the pre-operation state and return false.
	bool
	add_edge(
		AssemblyOP const & pre_op_assembly,
		utility::vector1<AssemblyOP> const & pre_op_assembly_list,
		AssemblyOP & assembly,
		utility::vector1<AssemblyOP> & assembly_list
	) const;

	///@brief Remove the most recently added edge from the Assembly
	///return void since this operation can never fail
	void
	delete_edge(
		AssemblyOP & assembly,
		utility::vector1<AssemblyOP> & assembly_list
	) const;

	///@brief Replace the most recently added node (try a different edge,
	///or start with a new node). Implemented as a delete followed by an
	///add and returns the result of the add operation. A return of false
	///reverts to the pre-operation state.
	bool
	switch_edge(
		AssemblyOP const & pre_op_assembly,
		utility::vector1<AssemblyOP> const & pre_op_assembly_list,
		AssemblyOP & assembly,
		utility::vector1<AssemblyOP> & assembly_list
	) const;

	void
	boltzman(
		utility::vector1<AssemblyOP> const & pre_op_assembly_list,
		utility::vector1<AssemblyOP> & assembly_list,
		AssemblyOP const & pre_op_assembly,
		AssemblyOP & working_assembly,
		core::Size cur_cycle,
		AssemblyOP & best_complete_assembly,
		core::Real & best_score
	) const;

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


	core::Size cycles_;
	core::Real max_temperature_;
	core::Real min_temperature_;
	core::Real min_assembly_score_;

	bool use_best_assembly_score_;

	core::Real add_probability_;
	core::Real delete_probability_;
	core::Real switch_probability_;
	bool remove_non_contiguous_assembly_;
};

} //sewing
} //protocols

#endif
