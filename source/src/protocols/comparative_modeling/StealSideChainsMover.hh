// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief simple mover for stealing side chains from one pose and sticking them
/// on another pose.
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_StealSideChainsMover_hh
#define INCLUDED_protocols_comparative_modeling_StealSideChainsMover_hh

#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.hh>

#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

class StealSideChainsMover : public protocols::moves::Mover {

public:

	StealSideChainsMover(
		core::pose::Pose const & source,
		core::id::SequenceMapping map
	);

	/// maps from pose to source_
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	core::pose::Pose const & source_;
	core::id::SequenceMapping map_;
};

} // comparative_modeling
} // protocols

#endif
