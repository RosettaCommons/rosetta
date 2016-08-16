// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AppendAssemblyMover.hh
///
/// @brief A Mover that uses loophash to find fragments that can
/// bridge a gap with minimal modifications to the original pose.
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_sewing_AppendAssemblyMover_HH
#define INCLUDED_protocols_sewing_AppendAssemblyMover_HH

// Unit Headers
#include <protocols/sewing/sampling/AppendAssemblyMover.fwd.hh>
#include <protocols/sewing/sampling/MonteCarloAssemblyMover.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.hh>

//Protocol headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace sewing  {

class AppendAssemblyMover : public MonteCarloAssemblyMover {

private:

	typedef MonteCarloAssemblyMover parent;

public:

	AppendAssemblyMover();

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;

	std::string
	get_name() const;

	///@brief The starting node for the Assembly from the append mover
	///will always be the input PDB
	virtual
	void
	add_starting_model(
		AssemblyOP assembly
	) const;

	core::pose::Pose
	refine_assembly(
		AssemblyOP & assembly
	);

	void
	hash_pdb_model(
		Model const & pdb_model
	);

	virtual
	void
	apply(
		core::pose::Pose & pose
	);

	virtual
	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

private:

	bool initialized_;

	core::pose::PoseOP partner_pose_;

	std::map<int, Model> models_;
	utility::vector1<core::Size> pose_segment_starts_;
	utility::vector1<core::Size> pose_segment_ends_;
	utility::vector1<core::Size> match_segments_;
};

} //sewing
} //protocols

#endif
