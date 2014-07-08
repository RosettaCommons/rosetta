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
//

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <iostream>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/excn/Exceptions.hh>
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
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/CompositionMover.hh>

//
#include <protocols/relax/FastRelax.hh>
//

// arguments for program:
// fast_relax_replonly.linuxgccrelease -in:file:s start.pdb -out:file:silent relax_replonly.out -out:file:silent_struct_type binary -nstruct 20
int
main ( int argc, char *argv[] ) {
	try{
	using namespace core;
	using namespace protocols;
	using namespace protocols::jd2;
	using namespace protocols::moves;

	devel::init(argc, argv);

	scoring::ScoreFunctionOP score = scoring::get_score_function();

	CompositionMoverOP container( new CompositionMover );
	container->add_mover( new protocols::simple_moves::RepulsiveOnlyMover() );
	container->add_mover( new protocols::relax::FastRelax(score) );

	JobDistributor::get_instance()->go(container);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

