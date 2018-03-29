// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_constraint_movers_ClearConstraintsMover_hh
#define INCLUDED_protocols_constraint_movers_ClearConstraintsMover_hh

#include <protocols/moves/Mover.hh>

namespace protocols {
namespace constraint_movers {

class ClearConstraintsMover : public moves::Mover {
public:
	ClearConstraintsMover();
	~ClearConstraintsMover() override;

	void apply( Pose & ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

} // moves
} // protocols

#endif
