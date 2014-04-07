// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief
/// @author jk+dj

#include <iostream>
#include <iomanip>

// Protocol Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>

// Core Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <core/id/AtomID_Map.hh>
#include <protocols/pockets/PocketConstraint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/Energies.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

#include <basic/options/option_macros.hh>

// Numeric Headers
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( Integer, num_angles )

/// General testing code
int
main( int argc, char * argv [] )
{
	try {
  NEW_OPT ( num_angles, "Number of different pose angles to measure score at", 1);

	devel::init(argc, argv);
  int angles = option[ num_angles ];
  if (angles <1){
    fprintf (stderr, "Error: invalid number of angles.  Must be greather than 0\n");
    return -1;
  }

	std::string const output_tag = option[ OptionKeys::out::output_tag ]();
	pose::Pose input_pose;

	//read in pdb file from command line
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( input_pose, input_pdb_name );


	scoring::ScoreFunctionOP scorefxn( getScoreFunction() );

	(*scorefxn)(input_pose);
  core::Real const starting_total_score = input_pose.energies().total_energies()[ total_score ];
	std::cout << "Total score at start without constraint is: " << starting_total_score << std::endl;


  std::string resid(option[ OptionKeys::pocket_grid::central_relax_pdb_num ]);


  std::vector< conformation::ResidueOP > residues = protocols::pockets::PocketGrid::getRelaxResidues(input_pose, resid);

/*
  int  central_relax_pdb_number;
  char chain = ' ';
  std::size_t fpos( resid.find(':') );
  if ( fpos != std::string::npos ) {
    central_relax_pdb_number = ObjexxFCL::int_of( resid.substr(0,fpos) );
    if (fpos != resid.size()-1 ) {
      chain = resid[ fpos+1 ];
    }
  } else {
    central_relax_pdb_number = ObjexxFCL::int_of( resid );
  }

  int seqpos = 0;
  for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j ) {
    if ( input_pose.pdb_info()->number(j) == central_relax_pdb_number ) {
      //seqpos_ = j;
      if (chain != ' '){
        if ( input_pose.pdb_info()->chain(j) == chain ) {
          seqpos = j;
        }
      }else{
        seqpos = j;
      }
    }
  }

  if ( seqpos == 0 ) {
*/
  if ( residues.size() == 0 ) {
    std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
    exit(1);
  }

  core::Real constraint_pocket_score = 0;
  core::Real largest_pocket_score = 0;

  protocols::pockets::PocketGrid pg( residues );
  for (int i=0; i<angles; ++i){
    if (i>0){
/*      core::Real x,y,z;
      x = (int) (numeric::random::uniform() *89 +1);
      y = (int) (numeric::random::uniform() *89 +1);
      z = (int) (numeric::random::uniform() *89 +1);
      numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_degrees(x) );
      numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_degrees(y) );
      numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_degrees(z) );
      core::Vector v(0,0,0);
      input_pose.apply_transform_Rx_plus_v(x_rot_mat, v);
      input_pose.apply_transform_Rx_plus_v(y_rot_mat, v);
      input_pose.apply_transform_Rx_plus_v(z_rot_mat, v);
*/
      pg.randomAngle();
    }else{
      pg.zeroAngle();
    }
    //protocols::pockets::PocketGrid pg( input_pose.conformation().residue(seqpos) );
    //if ( pg.autoexpanding_pocket_eval( input_pose.conformation().residue(seqpos), input_pose ) ){
    if ( pg.autoexpanding_pocket_eval( residues, input_pose ) ){
      constraint_pocket_score += pg.netTargetPocketVolume();
      largest_pocket_score += pg.netTargetPocketVolume();
      //std::cout<<pg.netTargetPocketVolume()<<std::endl;
    }
	//else std::cout<<"0\n";
    if (i==0){
      if (option[ OptionKeys::pocket_grid::pocket_dump_pdbs ]()){
        pg.dumpGridToFile();
      }
      if (option[ OptionKeys::pocket_grid::pocket_dump_exemplars ]()){
        pg.dumpExemplarToFile();
      }
    }
  }
  constraint_pocket_score /= angles;
  largest_pocket_score /= angles;
	std::cout << "Pocket score (unweighted) is: " << constraint_pocket_score << std::endl;
	std::cout << "Largest pocket score (unweighted) is: " << largest_pocket_score << std::endl;
  } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
	return 0;
}



