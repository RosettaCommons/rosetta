// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ForceDisulfidesMover.hh
/// @brief
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_simple_moves_CutChainMover_hh
#define INCLUDED_protocols_simple_moves_CutChainMover_hh

// Unit headers
#include <protocols/simple_moves/CutChainMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {


/// @brief simple mover that sperates Fv from Fl into two seperate chains
class CutChainMover : public moves::Mover {

public:
	CutChainMover();
	~CutChainMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;


	core::Size
	chain_cut( core::pose::Pose & pose);

	core::Size
	chain_cut( core::pose::Pose & pose, core::Size start_res,core::Size end_res);

	void
	create_subpose(core::pose::Pose & pose );

	void
	foldTree (core::pose::Pose & pose);

	core::Real bond_length() const;
	core::Size chain_id() const;

	void bond_length(core::Real const);
	void chain_id(core::Size const);

private:
	core::Real bond_length_;
	core::Size chain_id_;

};


} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_CutChainMover_HH
