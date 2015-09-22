// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <utility/io/mpistream.hh>
#include <core/kinematics/MoveMap.hh>

//Protocol Headers
#include <protocols/pockets/PocketGrid.hh>
//#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace core::pose::datacache;
using namespace core::id;
using namespace protocols::simple_moves;
using namespace protocols::rigid;


OPT_KEY( String, eggshell1_fname )
OPT_KEY( String, eggshell2_fname )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.karen_pocket_compare.main" );

/// General testing code
int main( int argc, char * argv [] ) {

	try{

		NEW_OPT( eggshell1_fname, "eggshell1_fname", "fname" );
		NEW_OPT( eggshell2_fname, "eggshell2_fname", "fname" );

		//initializes Rosetta functions
		devel::init(argc, argv);

		std::string const fname1 ( option[ eggshell1_fname ] );
		std::string const fname2 ( option[ eggshell2_fname ] );

		TR << "Starting eggshell compare" << std::endl;

		protocols::pockets::EggshellGrid eggshell1( fname1 );
		// eggshell1.load_eggshell( fname1 ) ;
		protocols::pockets::EggshellGrid eggshell2( fname2 );
		// eggshell2.load_eggshell( fname2 ) ;

		// call function to compare template and comparison EggshellGrids, report score
		core::Real d1 = eggshell1.get_eggshell_distance( eggshell2 );
		TR << "Distance is: " << d1 << std::endl;

		//utility::io::ozstream fout;
		//fout.open("distance.txt", std::ios::out);
		//fout << d1 << std::endl;
		//fout.close();
		//  fout.clear();

		TR << "Done!" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}

