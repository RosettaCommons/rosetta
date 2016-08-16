// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector1.hh>
#include <numeric/random/random.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/sasa.hh>

#include <utility/vector1.hh>

#include <string>
#include <iostream>

int
main( int argc, char* argv [] ) {

	try {

	// options, random initialization
	devel::init( argc, argv );

	core::scoring::disulfides::FullatomDisulfidePotentialOP potent( new core::scoring::disulfides::FullatomDisulfidePotential );

	potent->print_score_functions();

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
