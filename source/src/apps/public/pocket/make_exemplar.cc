// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This File Is Made Available Under The Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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
#include <core/pose/util.hh>
#include <core/pose/PDB_Info.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
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
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

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
using namespace protocols::rigid;
using namespace protocols::simple_moves;


OPT_KEY( String, central_relax_pdb_num )

static basic::Tracer TR( "apps.public.make_exemplar.main", basic::t_debug );

/// General testing code
int main( int argc, char * argv [] ) {

	try{

	NEW_OPT( central_relax_pdb_num, "target residue(s)", "-1");

	TR << "Calling init" << std::endl;
	//initializes Rosetta functions
	devel::init(argc, argv);

	TR << "done" << std::endl;
	//allows output when running the program
	//sets input residue numbers to resid and resid_c
	std::string const resid_c = option[central_relax_pdb_num];
	TR << "Starting pocket compare" << std::endl;

	// create pose for comparison pose from pdb
	std::string const comparison_pdb_name ( basic::options::start_file() );
	pose::Pose comparison_pose;
	core::import_pose::pose_from_pdb( comparison_pose, comparison_pdb_name );
	TR << "set pdb"<< "    Number of residues: " << comparison_pose.total_residue() << std::endl;

	std::string tag = "";
        if (!option[ OptionKeys::out::output_tag ]().empty()){
          tag = "." + option[ OptionKeys::out::output_tag ]();
        }

	std::vector< conformation::ResidueOP > residues = protocols::pockets::PocketGrid::getRelaxResidues(comparison_pose, resid_c);
	// call function to make a grid around a target residue (seqpos)
	protocols::pockets::PocketGrid comparison_pg( residues );
	//call function to define the pocket
	comparison_pg.zeroAngle();
 	comparison_pg.autoexpanding_pocket_eval( residues, comparison_pose ) ;
	//output the pocket in a pdb
	//comparison_pg.dumpGridToFile();
	TR << "Pocket defined" << std::endl;
std::cout << "Pocket score (unweighted) is: " << comparison_pg.netTargetPocketVolume() << std::endl;
	// dump to file, with name based on input pdb
	//std::stringstream out_fname;
	//out_fname << comparison_pdb_name << tag << ".pocket.pdb";
	//std::cout<<out_fname.str() <<std::endl;
	//comparison_pg.dumpTargetPocketsToPDB( out_fname.str() );
  std::stringstream out_exfname;
  out_exfname << comparison_pdb_name << tag << ".exemplar.pdb";
  //std::cout<<out_exfname.str() <<std::endl;
	comparison_pg.dumpExemplarToFile( out_exfname.str() );	

	//std::stringstream out_pfname;
	//out_pfname << comparison_pdb_name << tag << ".pocket";
  //std::cout<<out_pfname.str() <<std::endl;
	//comparison_pg.dumpTargetPocketsToFile( out_pfname.str() );

	//utility::io::ozstream fout;
	//fout.open("distance.txt", std::ios::out);
	//fout << d1 << std::endl;
	//fout.close();
	//		fout.clear();

	TR << "Done!" << std::endl;
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
	return 0;

}

