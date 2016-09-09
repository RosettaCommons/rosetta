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

// libRosetta headers
//#include <basic/options/option.hh>

#include <core/types.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/id/AtomID.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/VariantType.hh>

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/TenANeighborGraph.hh>


#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>


//#include <core/pack/rotamer_set/OptEData.hh>

#include <utility/graph/Graph.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>

#include <basic/basic.hh>

#include <basic/database/open.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>

#include <protocols/dna/util.hh>
//#include <protocols/dna/classes.hh>

#include <protocols/loops/Loops.hh>

#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/IRCollection.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

//Auto Headers
#include <core/import_pose/import_pose.hh>


//silly using/typedef

using namespace core;
using namespace pose;
using namespace chemical;
using namespace scoring;
using namespace optimization;

using utility::vector1;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void
place_waters_and_minimize( Pose & pose )
{

  ScoreFunction scorefxn;
  scorefxn.set_weight( fa_atr, 0.80 );
  scorefxn.set_weight( fa_rep, 0.44 );
  scorefxn.set_weight( fa_sol, 0.65 );
//  scorefxn.set_weight( gb_elec, 0.20 );

  scorefxn.set_weight( h2o_intra, 0.01 );
  scorefxn.set_weight( h2o_hbond, 1.0 );

  scorefxn.set_weight( hbond_sr_bb, 1.0 );
  scorefxn.set_weight( hbond_lr_bb, 1.0 );
  scorefxn.set_weight( hbond_bb_sc, 1.0 );
  scorefxn.set_weight( hbond_sc, 1.0 );

//  scorefxn.set_weight( pro_close, 1.0 );

	Energy score_orig = scorefxn( pose );

	std::cout << "Score before pack " << score_orig << std::endl;

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	task->set_bump_check( true );
//	pack::pack_rotamers( pose, scorefxn, task);

	Energy end_score = scorefxn( pose );
	std::cout << "Score after pack " << end_score << std::endl;

    kinematics::MoveMap mm;
    mm.set_bb( false );
    mm.set_chi( true );

  	{
    // setup the options
    MinimizerOptions options( "lbfgs_armijo_nonmonotone", 1.0e-3, true /*use_nblist*/, false /*deriv_check*/ );

		AtomTreeMinimizer minimizer;
//    minimizer.run( pose, mm, scorefxn, options );
//    dump_pdb( pose, "post_minimization.pdb" );
    }

	Energy min_score = scorefxn( pose );
	std::cout << "Score after minimization " << min_score << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void
build_asn_motifs()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace io::pdb;


	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace optimization;

	Pose pose;
	core::import_pose::pose_from_file( pose, "2GTX_asn_ala3.pdb" , core::import_pose::PDB_file);

  std::string weights( "soft_rep_design" );
  ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( weights ) );

	protocols::motifs::MotifLibrary motif_lib;

	motif_lib.add_from_file( "asn_combined_library" );

	// We want to build lots of motifs around ASN

	utility::vector1< core::Size > target_positions;
	target_positions.push_back( pose.pdb_info()->pdb2pose( 'I', 1 ) );

	protocols::motifs::IRCollection inv_rotamers( pose, motif_lib, target_positions );

	// Use the Loop facility to define flexible segments
	protocols::loops::Loops flexible_regions;
	flexible_regions.add_loop( pose.pdb_info()->pdb2pose( 'A', 57 ), pose.pdb_info()->pdb2pose( 'A', 71 ) );

	// Mutate the entire loop to alanines

	inv_rotamers.incorporate_motifs( pose, flexible_regions );

	return;

}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		//using namespace core;
		devel::init( argc, argv );

		build_asn_motifs();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
