// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Carl Walkey ( cwalkey at uw dot edu )

#ifndef INCLUDED_protocols_simple_moves_AddResidueLabelMover_hh
#define INCLUDED_protocols_simple_moves_AddResidueLabelMover_hh

#include <protocols/simple_moves/AddResidueLabelMover.fwd.hh>

#include <protocols/moves/Mover.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>


namespace protocols {
namespace simple_moves {

class AddResidueLabelMover : public protocols::moves::Mover {
public:
	AddResidueLabelMover();
	AddResidueLabelMover(core::select::residue_selector::ResidueSelectorCOP, std::string);
	~AddResidueLabelMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

private:
	core::select::residue_selector::ResidueSelectorCOP selector_;
	std::string label_;

};

} // pilot
} // apps

#endif
