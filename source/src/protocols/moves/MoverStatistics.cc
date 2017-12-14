// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Keeps track of the acceptance_rate of a Mover
/// @author Monica Berrondo August 16 2007


// Package headers
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverStatistics.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <string>
#include <sstream>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


static basic::Tracer TR( "protocols.moves.TrialMover" );

//static basic::Tracer trDebug("protocols.moves.TrialMover", basic::t_debug );
//MY_TRACERS("protocols.moves.TrialMover")


using namespace ObjexxFCL::format;

namespace protocols {
namespace moves {

/// @details Auto-generated virtual destructor
MoverStatistics::~MoverStatistics() = default;

using Real = core::Real;
using basic::Error;
using basic::Warning;

void MoverStatistics::print( MonteCarloOP mc, std::string const & type )
{
	//clear_score();
	//return;

	//T("protocols.moves.TrialMover.energies") << "trialE ";
	if ( TR.Trace.visible() ) { //change from Debug --> Trace since it produces output every step!
		std::ostringstream outstring;

		for ( double & ii : score_ ) {
			outstring << F( 9, 3, ii ) << " ";
		}
		outstring << F( 9, 3, mc->last_accepted_score() ) << "  " << F( 9, 3, mc->lowest_score() ) << "  " << type << "   ";

		switch( mc->mc_accepted() ) {
		case 0 : outstring << "reject             "; break;
		case 1 : outstring << "thermal accept     "; break;
		case 2 : outstring << "downhill accept    "; break;
		case 3 : outstring << "lowest accept      "; break;
		}
		//TR.Debug << outstring.str() << "\n";
		TR.Trace << outstring.str() << std::endl;
	}
	clear_score();
	return;
}


}

}


