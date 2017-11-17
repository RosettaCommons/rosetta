// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       apps/pilot/membrane_scorefxn.cc
///
/// @brief      Report Membrane Per-Residue Scores
/// @details    Report scores for membrane framework supported centroid energy terms.
///             Tested terms include MPEnv, MPPair, MPCBeta, MPTermini and MPTMProj.
///             MPNonHelix is currently not included in the integration test.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (5/2/14)

// App headers
#include <devel/init.hh>

// Options System
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/jd2/JobDistributor.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
#include <cstdlib>
#include <algorithm>

static basic::Tracer TR( "membrane_scorefxn" );

using namespace protocols;

/// @brief Membrane Relax Mover: Create a Membrane Pose and Apply Relax Protocol
class MembraneSfxnMover : public moves::Mover {

public:

	MembraneSfxnMover() :
		moves::Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MembraneSfxnMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace protocols::membrane;
		using namespace core::scoring;

		// Create a membrane pose
		AddMembraneMoverOP mp( new AddMembraneMover() );
		mp->apply(pose);

		// Create new scoring function from membrane score weights
		ScoreFunctionOP sfxn( new ScoreFunction() );

		sfxn->set_weight( MPEnv, 1.0);
		sfxn->set_weight( MPPair, 1.0);
		sfxn->set_weight( MPCbeta, 1.0);
		sfxn->set_weight( MPTermini, 1.0);
		sfxn->set_weight( MPTMProj, 1.0);
		sfxn->set_weight( Mpair, 1.0 );
		sfxn->set_weight( Menv, 1.0 );
		sfxn->set_weight( Mcbeta, 1.0 );
		// sfxn->set_weight( Menv_termini, 1.0 );
		// sfxn->set_weight( Menv_tm_proj, 1.0 );

		sfxn->score(pose);

		// Show Pose energies
		pose.energies().show(std::cout);
	}
};

typedef utility::pointer::shared_ptr< MembraneSfxnMover > MembraneSfxnMoverOP;

/// @brief Main Function
int main( int argc, char* argv[] )
{
	try {

		using namespace core::scoring;
		using namespace protocols::membrane;
		using namespace protocols::jd2;

		using namespace basic::options;

		// Initialize Options System, RG, and All Factory_Registrators
		devel::init(argc, argv);

		TR << "Membrane Centroid Scoring Function ScoreFunction Fingerprint Test" << std::endl;
		TR << "@ralford - updated 5/2/14" << std::endl;
		TR << "=======================================================" << std::endl;

		// Initialize Membrane protein From Membrane Mover

		MembraneSfxnMoverOP mp( new MembraneSfxnMover() );
		JobDistributor::get_instance()->go(mp);

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
