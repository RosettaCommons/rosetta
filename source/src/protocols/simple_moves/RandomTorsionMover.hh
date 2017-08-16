// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/moves/RandomTorsionMover.hh
/// @brief RandomTorsionMover class declaration file
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_simple_moves_RandomTorsionMover_HH
#define INCLUDED_protocols_simple_moves_RandomTorsionMover_HH

// Unit Headers
#include <protocols/simple_moves/RandomTorsionMover.fwd.hh>

// protocols headers
#include <protocols/moves/Mover.hh>

// core headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

namespace protocols {
namespace simple_moves {

/// @details Simple class that randomly purturbs a random torsion selected from a movemap
class RandomTorsionMover : public protocols::moves::Mover {
public:

	// default ctor
	RandomTorsionMover();

	// ctor
	RandomTorsionMover( core::kinematics::MoveMapOP move_map, core::Real max_angle, core::Size num_moves );

	// cctor
	RandomTorsionMover( RandomTorsionMover const & other );

	// dtor
	~RandomTorsionMover() override;

	// mover interface
	void apply( core::pose::Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	// mover specific
	void setup_torsion_list( core::pose::Pose & pose );

private:
	// the movemap
	core::kinematics::MoveMapOP move_map_;

	// list of torsion ids allowed to be moved
	utility::vector1< core::id::TorsionID > torsion_id_list_;

	// max angle change
	core::Real max_angle_;

	// numer of moves to perform
	core::Size num_moves_;
};

} //namespace simple_moves
} //namespace protocols

#endif // INCLUDED_protocols_moves_RandomTorsionMover_HH
