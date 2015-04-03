// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AlignChainMover.hh
/// @brief switch the chain order

#ifndef INCLUDED_protocols_simple_moves_AlignChainMover_hh
#define INCLUDED_protocols_simple_moves_AlignChainMover_hh

#include <protocols/simple_moves/AlignChainMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// C++ Headers
namespace protocols {
namespace simple_moves {

/// @brief Align a chain in the working pose to a chain. CA superposition
/// @details  
///  source_chain: the chain number in the working pose. 0 means the entire pose.
///  target_chain: the chain number in the target pose. 0 means the entire pose.
///  **target and source chains must have the same length of residues.
class AlignChainMover : public moves::Mover {
public:
	AlignChainMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	core::pose::PoseOP pose() const;
	void pose( core::pose::PoseOP pose );

	core::Size source_chain() const { return source_chain_; }
	void source_chain( core::Size const s ){ source_chain_ = s; }

	core::Size target_chain() const { return target_chain_; }
	void target_chain( core::Size const s ){ target_chain_ = s; }
private:
	core::pose::PoseOP pose_; //dflt NULL;
	core::Size source_chain_, target_chain_ ; //dflt 0; which chains to align from source (current pose) on target (pose from disk). 0 means all chains
};


} // simple_moves
} // protocols

#endif
