// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/ncbb/oop/OopPuckMover.hh
/// @brief
/// @author
#ifndef INCLUDED_protocols_simple_moves_oop_OopPuckMover_hh
#define INCLUDED_protocols_simple_moves_oop_OopPuckMover_hh
// Unit Headers
#include <protocols/ncbb/oop/OopPuckMover.fwd.hh>
#include <protocols/ncbb/oop/OopMover.hh>
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

//kdrew: defines puck mover class for D chiral oop residues in half-chair conformation
class OopDPuckPlusMover : public OopMover {

public:

	OopDPuckPlusMover( core::Size oop_seq_position );
	virtual ~OopDPuckPlusMover();
	virtual std::string get_name() const;

};//end OopDPuckPlusMover

//kdrew: defines puck mover class for D chiral oop residues in boat conformation
class OopDPuckMinusMover : public OopMover {

public:

	OopDPuckMinusMover( core::Size oop_seq_position );
	virtual ~OopDPuckMinusMover();
	virtual std::string get_name() const;

};//end OopDPuckMinusMover


}//namespace oop
}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_oop_OopPuckMover_hh
