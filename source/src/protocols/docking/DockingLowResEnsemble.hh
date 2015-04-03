// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/docking/DockinLowResEnsemble
/// @brief low resolution mode for ensemble docking
/// @details
/// @author Daisuke Kuroda


#ifndef INCLUDED_protocols_docking_DockingLowRes_Ensemble_hh
#define INCLUDED_protocols_docking_DockingLowRes_Ensemble_hh

#include <protocols/docking/types.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockingLowResEnsemble.fwd.hh>
#include <protocols/docking/DockingEnsemble.fwd.hh>

// Package headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/docking/ConformerSwitchMover.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>

#include <string>

// option key includes


#include <utility/vector1.hh>
#include <iostream>

namespace protocols {
namespace docking {

class DockingLowResEnsemble : public DockingLowRes
{
	//typedef core::Real Real;
public:

	/// @brief Default constructor
	DockingLowResEnsemble();

	// destructor
	~DockingLowResEnsemble();

	/// @brief Constructor with two arguments.  The first is scorefunction to be used for docking, the second is the
	///		DockJumps.
	DockingLowResEnsemble(
        core::scoring::ScoreFunctionCOP scorefxn,
        DockJumps const movable_jumps
    );
	
    virtual protocols::moves::MoverOP clone() const;
    
    virtual void show( std::ostream & out=std::cout ) const;
    
	// option setters
	void set_ensemble1( DockingEnsembleOP ensemble1 );
	void set_ensemble2( DockingEnsembleOP ensemble2 );
    
protected:
	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	virtual void finalize_setup( core::pose::Pose & pose); // Comment out by DK

private:
    protocols::docking::ConformerSwitchMoverOP ensemble1_mover_;
	protocols::docking::ConformerSwitchMoverOP ensemble2_mover_;
};

} // docking
} // protocols

#endif
