// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  Contains currently: Relax Baseclass, ClassicRelax Stage 1,2,3, ClassicRelax
///
///
/// @author Mike Tyka


#ifndef INCLUDED_protocols_relax_util_hh
#define INCLUDED_protocols_relax_util_hh
// AUTO-REMOVED #include <core/scoring/MembraneTopology.fwd.hh> //pba

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/relax/RelaxProtocolBase.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers


namespace protocols {
namespace relax {

//pba  read in membrane topology
void setup_membrane_topology( core::pose::Pose & pose, std::string spanfile );

void relax_pose( core::pose::Pose& pose, core::scoring::ScoreFunctionOP scorefxn, std::string const& tag );

RelaxProtocolBaseOP generate_relax_from_cmd( bool NULL_if_no_cmd = false );

}
} // protocols

#endif
