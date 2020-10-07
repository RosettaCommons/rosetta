// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddSapMathConstraintMover.hh
/// @brief Mover that adds the SapMathConstraint to the pose
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_AddSapMathConstraintMover_hh
#define INCLUDED_protocols_simple_moves_AddSapMathConstraintMover_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/AddSapMathConstraintMover.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapMathConstraint.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>


// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {


/// @brief Adds variant types to selected residues
class AddSapMathConstraintMover : public protocols::moves::Mover
{
public:
	AddSapMathConstraintMover();

	AddSapMathConstraintMover( AddSapMathConstraintMover const & other );

	AddSapMathConstraintMover & operator=( AddSapMathConstraintMover const & ot );

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &) override;


	void add_constraint( core::Real weight, std::string const & name );

	void set_upper_bound( core::Real upper );

	void set_lower_bound( core::Real lower );

	void set_penalty_per_unit( core::Real penalty );


	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	core::pack::guidance_scoreterms::sap::SapMathConstraint cst_;

};

} // moves
} // protocols


#endif
