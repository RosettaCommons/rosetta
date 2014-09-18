// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockingLowResEnsemble
/// @brief protocols that are specific to low resolution ensemble docking
/// @detailed
/// @author Daisuke Kuroda

#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockingLowResEnsemble.hh>
#include <protocols/docking/DockingEnsemble.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>

// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/docking/ConformerSwitchMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>

// Utility Headers
#include <utility/tools/make_vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <string>

//Utility Headers
// AUTO-REMOVED #include <numeric/conversions.hh>

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.docking.DockingLowResEnsemble" );

using namespace core;

namespace protocols {
namespace docking {

// default constructor
DockingLowResEnsemble::DockingLowResEnsemble() : DockingLowRes()
{
	type( "DockingLowResEnsemble" );
}
    
// constructor with arguments    
DockingLowResEnsemble::DockingLowResEnsemble(
    core::scoring::ScoreFunctionCOP scorefxn,
    DockJumps const movable_jumps
) : DockingLowRes( scorefxn, movable_jumps )
{
    type( "DockingLowResEnsemble" );
}

//destructor
DockingLowResEnsemble::~DockingLowResEnsemble() {}

    
protocols::moves::MoverOP DockingLowResEnsemble::clone() const {
	return new DockingLowResEnsemble(*this);
}

 
void DockingLowResEnsemble::set_ensemble1( DockingEnsembleOP ensemble1 )
{
    ensemble1_mover_ = new protocols::docking::ConformerSwitchMover( ensemble1 );
}

void DockingLowResEnsemble::set_ensemble2( DockingEnsembleOP ensemble2 )
{
    ensemble2_mover_ = new protocols::docking::ConformerSwitchMover( ensemble2 );
}

void DockingLowResEnsemble::finalize_setup( core::pose::Pose & pose){
	using namespace moves;

    DockingLowRes::finalize_setup( pose );
    
    // add if ens_mover1 ...
    // inner_move_sequence
	docking_lowres_protocol_->add_mover( ensemble1_mover_ );
	docking_lowres_protocol_->add_mover( ensemble2_mover_ );        
}

void DockingLowResEnsemble::show( std::ostream & out ) const
{
    DockingLowRes::show( out );
    
    using namespace ObjexxFCL::format;
    std::string line_marker = "///";
    out << line_marker << " Ensemble 1: " << ( ( ensemble1_mover_ ) ? ( "on" ) : ( "off " ) );
    out << space( 59 ) << line_marker << std::endl;
    out << line_marker << " Ensemble 2: " << ( ( ensemble2_mover_ ) ? ( "on" ) : ( "off " ) );
    out << space( 59 ) << line_marker << std::endl;
}
    
} // namespace docking
} // namespace protocols
