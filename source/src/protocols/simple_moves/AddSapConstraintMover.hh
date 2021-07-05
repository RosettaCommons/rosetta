// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddSapConstraintMover.hh
/// @brief Mover that adds the SapConstraint to the pose
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_AddSapConstraintMover_hh
#define INCLUDED_protocols_simple_moves_AddSapConstraintMover_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/AddSapConstraintMover.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.fwd.hh>


// Project headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>


// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers


namespace protocols {
namespace simple_moves {



/// @brief Adds variant types to selected residues
class AddSapConstraintMover : public protocols::moves::Mover
{
public:
	AddSapConstraintMover();

	AddSapConstraintMover( AddSapConstraintMover const & other );

	AddSapConstraintMover & operator=( AddSapConstraintMover const & ot );

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &) override;

	/// @brief Set the ResidueSelector used by this mover.
	void set_score_selector( core::select::residue_selector::ResidueSelectorCOP const & selector_in );

	void set_sap_calculate_selector( core::select::residue_selector::ResidueSelectorCOP const & selector_in );

	void set_sasa_selector( core::select::residue_selector::ResidueSelectorCOP const & selector_in );

	void set_sap_goal( core::Real goal );

	void set_sap_lb_goal( core::Real goal );

	void set_penalty_per_sap( core::Real penalty );

	void set_packing_correction( core::Real correction );

	void set_speed( std::string const & speed );

	void set_full_accuracy_when_scoring( bool full_accuracy );

	void set_name( std::string const & name );


	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	core::pack::guidance_scoreterms::sap::SapConstraintOptionsOP options_;

};

} // moves
} // protocols


#endif
