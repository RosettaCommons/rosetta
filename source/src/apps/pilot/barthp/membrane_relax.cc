// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief


// libRosetta headers

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <protocols/loops/ccd_closure.hh>
#include <protocols/relax_protocols.hh>

#include <core/types.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/MembraneTopology.hh>

#include <basic/prof.hh> // profiling

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>

//#include <core/chemical/residue_io.hh>

#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
//#include <core/id/AtomID_Map.Pose.hh>

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>

#include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef


//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;



using namespace core;
using namespace protocols;

using utility::vector1;

using io::pdb::dump_pdb;





///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	devel::init(argc, argv);

	// is this also done inside devel::init?
	//numeric::random::RandomGenerator::initializeRandomGenerators(
	//	 111111, numeric::random::_RND_TestRun_, "ran3");

	using namespace protocols::moves;
	using namespace scoring;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
  using namespace core::chemical;
  using namespace core::io::silent;

  // setup residue types
  ResidueTypeSetCAP rsd_set;
  if ( option[ in::file::fullatom ]() )
     rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  // configure silent-file data object
  ScoreFunctionOP scorefxn = get_score_function();

  // configure silent-file data object
  core::io::silent::SilentFileData sfd;
  pose::Pose pose;
  std::string infile  = *(option[ in::file::silent ]().begin());
  std::string const spanfile = option[ in::file::spanfile ]();
  std::string outfile = option[ out::file::silent ]();
  utility::io::ozstream output;
  output.open( outfile.c_str() );
  std::cout << "spanfile: " << spanfile << "\n";
  //  core::scoring::MembraneTopologyOP topology=new core:scoring:MembraneTopology;

  core::scoring::MembraneTopologyOP topologyOP = new core::scoring::MembraneTopology;
  pose.data().set( MEMBRANE_TOPOLOGY, topologyOP );
  //  core::scoring::MembraneTopology & topology=*( static_cast< core::scoring::MembraneTopology * >( pose.data().get_ptr( basic::MEMBRANE_TOPOLOGY )() ));
  core::scoring::MembraneTopology & topology=*( static_cast< core::scoring::MembraneTopology * >( pose.data().get_ptr( MEMBRANE_TOPOLOGY )() ));
  topology.initialize(spanfile);

  SilentFileData silent_file_data;
  std::string silent_file( "fa_memb.out" );

  std::cout << "READING start.pdb" << std::endl;
	core::import_pose::pose_from_pdb( pose, *rsd_set, "start.pdb" ); // default is standard fullatom residue_set

	kinematics::MoveMap mm;
	// setup moving dofs
	//for ( int i=30; i<= 35; ++i ) {
		mm.set_bb ( true );
		mm.set_chi( true );
	//}

	using namespace optimization;
	/*using pose::Pose;
	using id::AtomID;
	using id::DOF_ID;
	using id::PHI;
	using id::THETA;
	using id::D;*/

  pose.data().set( MEMBRANE_TOPOLOGY, topologyOP );

	//scoring::ScoreFunction scorefxn;
	//scorefxn.set_weight( scoring::envsmooth, 1.0 );
	//scorefxn.set_weight( scoring::fa_elec , 1.0 );
	//scorefxn.set_weight( scoring::fa_rep , 1.0 );
  /*scorefxn->set_weight( fa_atr, 0.80 );
  scorefxn->set_weight( fa_rep, 0.44 );
  scorefxn->set_weight( fa_sol, 0.00 );
  scorefxn->set_weight( fa_pair, 0.49 );
  scorefxn->set_weight( fa_dun, 0.00 );
  scorefxn->set_weight( rama, 0.2 );
  scorefxn->set_weight( omega, 0.5 );
  scorefxn->set_weight( hbond_lr_bb, 0.00 );
  scorefxn->set_weight( hbond_sr_bb, 0.00 );
  scorefxn->set_weight( hbond_bb_sc, 0.00 );
  scorefxn->set_weight( hbond_sc   , 1.10 );
  scorefxn->set_weight( fa_mbenv   , 2.00 );
  */

	(*scorefxn)(pose);

  output << "SCORE: TOTAL\t" ;
  scorefxn->show_line_headers(output);
  output << "SCORE: ";
  scorefxn->show_line(output,pose);

	scorefxn->show(std::cout, pose);

	// setup the options

  MinimizerOptions options( "dfpmin_armijo_nonmonotone", 0.001, true , /*was 0.000001*/
                          false , false );

	AtomTreeMinimizer minimizer;
	std::cout << "MINTEST: " << "\n";
	std::cout << "start score: " << (*scorefxn)( pose ) << "\n";
  int starttime = time(NULL);
  int startclock = clock();
	minimizer.run( pose, mm, *scorefxn, options );
  int endtime = time(NULL);
  int endclock = clock();
  std::cout << "Time: " << (endtime - starttime) << "   " << (endclock - startclock) << std::endl;
	pose.dump_pdb( "toto.pdb" );
  scorefxn->show(std::cout, pose);
	std::cout << "end score: " << (*scorefxn)( pose ) << "\n";


//	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );
//
////	// aiming for standard packer weights
//	scorefxn->set_weight( envsmooth, 0.80 );
//
//	(*scorefxn)(*pose);
//	scorefxn->show(std::cout, *pose);
//
////	using namespace optimization;
////	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );
////	scorefxn->set_weight( fa_atr, 0.80 );
////	scorefxn->set_weight( fa_rep, 0.44 );
////	scorefxn->set_weight( fa_sol, 0.65 );
////	scorefxn->set_weight( fa_pair, 0.49 );
////	scorefxn->set_weight( fa_dun, 0.56 );
////	scorefxn->set_weight( rama, 0.2 );
////	scorefxn->set_weight( omega, 0.5 );
////	scorefxn->set_weight( hbond_lr_bb, 1.17 );
////	scorefxn->set_weight( hbond_sr_bb, 1.17 );
////	scorefxn->set_weight( hbond_bb_sc, 1.17 );
////	scorefxn->set_weight( hbond_sc   , 1.10 );
//
//
//
//	kinematics::MoveMap mm;
//	mm.set_bb (  true );
//	mm.set_chi(  true );
//	AtomTreeMinimizer minimizer;
//	std::cout << "MINTEST: p_aa_pp" << "\n";
//	std::cout << "start score: " << (*scorefxn)( *pose ) << "\n";
//	int starttime = time(NULL);
//	int startclock = clock();
//	minimizer.run( *pose, mm, *scorefxn, options );
//	int endtime = time(NULL);
//	int endclock = clock();
//	std::cout << "Time: " << (endtime - starttime) << "   " << (endclock - startclock) << std::endl;
//
//	(*scorefxn)(*pose);
//	scorefxn->show(std::cout, *pose);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
