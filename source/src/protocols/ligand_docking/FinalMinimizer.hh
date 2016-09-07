// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_FinalMinimizer_hh
#define INCLUDED_protocols_ligand_docking_FinalMinimizer_hh

// Unit Headers
#include <protocols/ligand_docking/FinalMinimizer.fwd.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>

// Package Headers
//// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers

#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

/// @brief
class FinalMinimizer : public protocols::moves::Mover
{
public:
	FinalMinimizer();
	FinalMinimizer(core::scoring::ScoreFunctionOP score_fxn, MoveMapBuilderOP movemap_builder);
	~FinalMinimizer() override;
	FinalMinimizer(FinalMinimizer const & that);

	void apply( core::pose::Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	MoveMapBuilderOP movemap_builder_;

	protocols::simple_moves::MinMoverOP const
	get_final_min_mover(core::pose::Pose const & pose) const;
};

} //namespace ligand_docking
} //namespace protocols

#endif
