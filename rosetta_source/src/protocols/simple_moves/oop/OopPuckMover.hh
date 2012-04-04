// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/simple_moves/oop/OopPuckMover.hh
/// @brief
/// @author
#ifndef INCLUDED_protocols_simple_moves_oop_OopPuckMover_hh
#define INCLUDED_protocols_simple_moves_oop_OopPuckMover_hh
// Unit Headers
#include <protocols/simple_moves/oop/OopPuckMover.fwd.hh>
#include <protocols/simple_moves/oop/OopMover.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
//#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace oop {

class OopPuckPlusMover : public OopMover {

public:

	OopPuckPlusMover( core::Size oop_seq_position );
	virtual ~OopPuckPlusMover();
	virtual std::string get_name() const;

};//end OopPuckPlusMover

class OopPuckMinusMover : public OopMover {

public:

	OopPuckMinusMover( core::Size oop_seq_position );
	virtual ~OopPuckMinusMover();
	virtual std::string get_name() const;

};//end OopPuckMinusMover


/*
///@details
class OopPuckMover : public protocols::moves::Mover {

public:

	///@brief
	OopPuckMover(
	);

	OopPuckMover( utility::vector1< core::Size > oop_seq_positions, bool oop_puck_plus, bool oop_puck_minus, bool random, bool patch );

	OopPuckMover( utility::vector1< core::Size > oop_seq_positions, bool oop_puck_plus, bool oop_puck_minus, bool random, bool patch, bool oop_puck_small, core::Real max_small_angle );

	OopPuckMover( utility::vector1< core::Size > oop_seq_positions );

	virtual ~OopPuckMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:

	bool oop_puck_plus_;
	bool oop_puck_minus_;
	bool random_;
	utility::vector1< core::Size > const oop_seq_positions_;
	utility::vector1< std::string > available_moves_;
	bool patch_;
	bool oop_puck_small_;
	core::Real max_small_angle_;

};//end OopPuckMover
*/


}//namespace oop
}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_oop_OopPuckMover_hh
