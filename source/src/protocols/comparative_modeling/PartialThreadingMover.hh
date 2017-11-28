// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/PartialThreadingMover.hh
/// @brief
/// @author James Thompson
/// @author fpd some updates

// libRosetta headers

#ifndef INCLUDED_protocols_comparative_modeling_PartialThreadingMover_HH
#define INCLUDED_protocols_comparative_modeling_PartialThreadingMover_HH

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

class PartialThreadingMover : public protocols::moves::Mover {

public:
	/// @brief align describes the association between the query and template
	/// sequences, template_pose is the conformation from which to build a
	/// threading model.
	PartialThreadingMover(
		core::sequence::SequenceAlignment const & align,
		core::pose::Pose const & template_pose,
		bool const skip_repack = false
	);

	~PartialThreadingMover() override = default;

	/// @brief Threads the given Pose onto the template_pose with the
	/// SequenceAlignment provided.
	void apply( core::pose::Pose & query_pose ) override;

	std::string get_name() const override;

private:
	core::pose::Pose template_pose_;
	core::sequence::SequenceAlignment align_;
	bool skip_repack_;

}; // PartialThreadingMover

} // comparative_modeling
} // protocols

#endif
