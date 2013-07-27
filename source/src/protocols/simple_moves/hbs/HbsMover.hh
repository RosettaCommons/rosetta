// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/simple_moves/hbs/HbsMover.hh
/// @brief
/// @author Kevin Drew, kdrew@nyu.edu
#ifndef INCLUDED_protocols_simple_moves_hbs_HbsMover_hh
#define INCLUDED_protocols_simple_moves_hbs_HbsMover_hh
// Unit Headers
#include <protocols/simple_moves/hbs/HbsMover.fwd.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
//#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace hbs {

///@details
class HbsMover : public protocols::moves::Mover {

public:

	///@brief
	HbsMover( core::Size seq_position );
	HbsMover( core::Size seq_position, core::Real phi_angle, core::Real psi_angle );

	virtual ~HbsMover();
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	
	virtual void set_phi( core::Real angle ) { phi_angle_ = angle; }
	virtual void set_psi( core::Real angle ) { psi_angle_ = angle; }

  core::Size seq_pos_;
private:

	core::Real phi_angle_;
	core::Real psi_angle_;

};//end HbsMover


}//namespace hbs
}//namespace simple_moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_hbs_HbsMover_hh
