// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#include <protocols/simple_moves/ClearConstraintsMover.hh>
#include <protocols/simple_moves/ClearConstraintsMoverCreator.hh>

#include <core/pose/Pose.hh>

namespace protocols {
namespace simple_moves {

ClearConstraintsMover::ClearConstraintsMover(){}
ClearConstraintsMover::~ClearConstraintsMover(){}

void ClearConstraintsMover::apply( core::pose::Pose & pose )
{
	pose.remove_constraints();
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
ClearConstraintsMover::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
}
	
moves::MoverOP ClearConstraintsMover::clone() const { return moves::MoverOP( new ClearConstraintsMover( *this ) ); }
moves::MoverOP ClearConstraintsMover::fresh_instance() const { return moves::MoverOP( new ClearConstraintsMover ); }

protocols::moves::MoverOP
ClearConstraintsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ClearConstraintsMover );
}

std::string
ClearConstraintsMoverCreator::keyname() const
{
	return ClearConstraintsMoverCreator::mover_name();
}

std::string
ClearConstraintsMoverCreator::mover_name()
{
	return "ClearConstraintsMover";
}

std::string
ClearConstraintsMover::get_name() const {
	return "ClearConstraintsMover";
}
	
} // moves
} // protocols
