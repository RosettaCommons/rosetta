// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DeleteChainMover.hh
/// @brief Delete the chain specified by chain number

#ifndef INCLUDED_protocols_simple_moves_DeleteChainMover_hh
#define INCLUDED_protocols_simple_moves_DeleteChainMover_hh

#include <protocols/simple_moves/DeleteChainMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

// C++ Headers
namespace protocols {
namespace simple_moves {

class DeleteChainMover : public moves::Mover {
public:

	DeleteChainMover();

	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	void chain_num(
		core::Size chain_num
	);

	core::Size chain_num();

	virtual void apply( core::pose::Pose & pose );

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

private:
	core::Size chain_num_;
};


} // simple_moves
} // protocols

#endif
