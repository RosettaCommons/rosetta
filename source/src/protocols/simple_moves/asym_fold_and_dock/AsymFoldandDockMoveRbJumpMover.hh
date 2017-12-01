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
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_simple_moves_asym_fold_and_dock_AsymFoldandDockMoveRbJumpMover_hh
#define INCLUDED_protocols_simple_moves_asym_fold_and_dock_AsymFoldandDockMoveRbJumpMover_hh

// Unit headers
#include <protocols/simple_moves/asym_fold_and_dock/AsymFoldandDockMoveRbJumpMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/datacache/DataMap.fwd.hh>


// Utility Headers

namespace protocols {
namespace simple_moves {
namespace asym_fold_and_dock {
///////////////////////////////////////////////////////////////////////////////
class AsymFoldandDockMoveRbJumpMover : public moves::Mover
{
public:

	// default constructor
	AsymFoldandDockMoveRbJumpMover();
	AsymFoldandDockMoveRbJumpMover( core::Size chain_start );

	~AsymFoldandDockMoveRbJumpMover(){}

	void apply( core::pose::Pose & pose ) override;
	void find_new_jump_residue( core::pose::Pose & pose );
	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;
	//virtual std::string get_name() const override;
	std::string get_name() const override;

	static
	std::string mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Size chain_start_;

};

}
} // asym_fold_and_dock
} // rosetta
#endif
