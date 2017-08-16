// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/rigid/RollMover.hh
/// @brief
/// @author
#ifndef INCLUDED_protocols_rigid_RollMover_hh
#define INCLUDED_protocols_rigid_RollMover_hh
// Unit Headers
#include <protocols/rigid/RollMover.fwd.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

namespace protocols {
namespace rigid {

/// @details
class RollMover : public protocols::moves::Mover {

public:

	/// @brief
	RollMover();

	RollMover(
		core::Size start_res,
		core::Size stop_res,
		core::Real min_angle,
		core::Real max_angle,
		numeric::xyzVector< core::Real > axis,
		numeric::xyzVector< core::Real > translate
	);

	~RollMover() override;

	void apply( core::pose::Pose & pose ) override;
	void set_min_max_angles( core::Real min_angle, core::Real max_angle );

	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP clone() const override;


	void
	parse_my_tag(
		TagCOP /*tag*/,
		basic::datacache::DataMap & /*data*/,
		Filters_map const & /*filters*/,
		moves::Movers_map const & /*movers*/,
		Pose const & /*pose*/) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:

	core::Size start_res_;
	core::Size stop_res_;
	core::Real angle_;
	core::Real min_angle_;
	core::Real max_angle_;
	numeric::xyzVector< core::Real > axis_;
	numeric::xyzVector< core::Real > translate_;
	bool random_roll_;
	core::Real random_roll_angle_;
	core::Real random_roll_trans_;
};//end RollMover


}//namespace rigid
}//namespace protocols

#endif
