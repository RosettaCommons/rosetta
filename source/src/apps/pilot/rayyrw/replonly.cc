// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
//	@author Will Sheffler
//	@author Ray Wang


// libRosetta headers
//#include <core/options/option.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <utility/excn/Exceptions.hh>
#include <iostream>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
//score structures
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
//
#include <protocols/simple_moves/RepulsiveOnlyMover.hh>
#include <protocols/moves/Mover.hh>

int
main (int argc, char *argv[])
{
	try{
		using namespace core;

		devel::init(argc, argv);
		protocols::moves::MoverOP mover = new protocols::simple_moves::RepulsiveOnlyMover();
		scoring::ScoreFunctionOP score = scoring::getScoreFunction();

		pose::Pose pose;
	//changing fa_standard to centroid
		chemical::ResidueTypeSetCAP centroid_residue_set = chemical::ChemicalManager::get_instance()->residue_type_set(chemical::CENTROID );

	//Will changed the sequence of the argument here by checking out src/core/io/pdb/pose_io.hh
		core::import_pose::pose_from_pdb( pose, *centroid_residue_set, basic::options::option[basic::options::OptionKeys::in::file::s][1]) ;


		score->show(pose);
		pose.dump_scored_pdb("before.pdb",*score);
		mover->apply(pose);
		score->show(pose);
		pose.dump_scored_pdb("after.pdb",*score);


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
		return 0;

}
