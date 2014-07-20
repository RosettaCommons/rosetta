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
#include <core/pose/PDB_Info.hh>
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

// RMS alignmnet headers
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <core/pose/util.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( Integer, num_angles )
OPT_KEY( String, template_pdb_name )
OPT_KEY( String, contact_list )

//set to store pdb info keys
std::set <std::string> interface;

bool
is_interface_heavyatom(
        core::pose::Pose const & pose,
         core::pose::Pose const & ,//pose2,
        core::Size resno,
        core::Size atomno
)
{
        // ws get residue "key" for set
        std::ostringstream residuestream;
        residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
        std::string res_id = residuestream.str();

        core::conformation::Residue const & rsd = pose.residue(resno);
        if ( interface.count( res_id ) > 0 ) return rsd.is_protein() && !rsd.atom_is_hydrogen(atomno);

        return false;
}

Real
iface_pdb_superimpose_pose(
        pose::Pose & mod_pose,
        pose::Pose const & ref_pose
)
{
  id::AtomID_Map< id::AtomID > atom_map;
  core::pose::initialize_atomid_map( atom_map, mod_pose, id::BOGUS_ATOM_ID );
  for ( Size ii = 1; ii <= mod_pose.total_residue(); ++ii ) {
    if ( ! mod_pose.residue(ii).has("CA") ) continue;
    if ( ! mod_pose.residue(ii).is_protein() ) continue;
    for ( Size jj = 1; jj <= ref_pose.total_residue(); ++jj ) {
      if ( ! ref_pose.residue(jj).has("CA") ) continue;
      if ( ! ref_pose.residue(jj).is_protein() ) continue;
      if ( mod_pose.pdb_info()->chain(ii) != ref_pose.pdb_info()->chain(jj)) continue;
      if ( mod_pose.pdb_info()->number(ii) != ref_pose.pdb_info()->number(jj)) continue;
      if ( is_interface_heavyatom ( ref_pose, mod_pose, jj, 2) ){
        id::AtomID const id1( mod_pose.residue(ii).atom_index("CA"), ii );
        id::AtomID const id2( ref_pose.residue(jj).atom_index("CA"), jj );
        atom_map.set( id1, id2 );
      }

      break;
    }

  }
  return superimpose_pose( mod_pose, ref_pose, atom_map );
}


/// General testing code
int
main( int argc, char * argv [] )
{
	try {
  	NEW_OPT ( num_angles, "Number of different pose angles to measure score at", 1);
        NEW_OPT( template_pdb_name, "template pdb", "template.pdb" );
        NEW_OPT ( contact_list, "File name for optional list of contact residues to check","");

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

	//sets template input pdb name
        std::string const template_fname ( option[ template_pdb_name ] );

        //sets pdb as a Rosetta pose
        pose::Pose template_pose;
        core::import_pose::pose_from_pdb( template_pose, template_fname );

  	std::string const cfilename = option[ contact_list ];
  	if ( cfilename != "" ){
    		std::ifstream ifs(cfilename.c_str(), std::ifstream::in);
    		if (!ifs.is_open()){
      			std::cout<< "Error opening contact list file "<<cfilename<<std::endl;
      			return -100;
    		}
    		//ifb.open (cfilename,std::ios::in);
    		//std::ostream ios(&ifb);
    		std::string intres;
    		while (ifs.good()){
      			ifs >> intres;
      			interface.insert(intres);
    		}

    		iface_pdb_superimpose_pose( input_pose, template_pose);
  	}else{


    		// align comparison pose to template pose
        	//      protocols::simple_moves::SuperimposeMoverOP sp_mover = new protocols::simple_moves::SuperimposeMover( template_pose );
        	protocols::simple_moves::SuperimposeMoverOP sp_mover = new protocols::simple_moves::SuperimposeMover();
        	sp_mover->set_reference_pose( template_pose, 1, template_pose.total_residue() );
        	sp_mover->set_target_range( 1, template_pose.total_residue() );
        	sp_mover->apply( input_pose );
  	}


	scoring::ScoreFunctionOP scorefxn = get_score_function();

	(*scorefxn)(input_pose);
  core::Real const starting_total_score = input_pose.energies().total_energies()[ total_score ];
	std::cout << "Total score at start without constraint is: " << starting_total_score << std::endl;


  std::string resid(option[ OptionKeys::pocket_grid::central_relax_pdb_num ]);



  std::vector< conformation::ResidueOP > tmpresidues = protocols::pockets::PocketGrid::getRelaxResidues(input_pose, resid);
  if ( tmpresidues.size() == 0 ) {
    std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
    exit(1);
  }

  for (int i=0; i<angles; ++i){
    core::Real x,y,z;
    if (i>0){
      x = (int) (numeric::random::uniform() *89 +1);
      y = (int) (numeric::random::uniform() *89 +1);
      z = (int) (numeric::random::uniform() *89 +1);
    }else{
      x=0;
      y=0;
      z=0;
    }

    std::stringstream rotpdbtag;
    std::stringstream tag;
    std::stringstream tag2;
    rotpdbtag << input_pdb_name << ".alignedto." << template_fname << ".pdb";
    tag << input_pdb_name << ".rotation_experiment." << x <<"-"<<y<<"-"<<z<< ".rotated";
    tag2 << input_pdb_name << ".rotation_experiment." << x <<"-"<<y<<"-"<<z<< ".reset";

    input_pose.dump_pdb(rotpdbtag.str());

    pose::Pose rotated_pose(input_pose);
    numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_degrees(x) );
    numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_degrees(y) );
    numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_degrees(z) );
    numeric::xyzMatrix<core::Real> x_inv_rot_mat( numeric::x_rotation_matrix_degrees(-1. * x) );
    numeric::xyzMatrix<core::Real> y_inv_rot_mat( numeric::y_rotation_matrix_degrees(-1. * y) );
    numeric::xyzMatrix<core::Real> z_inv_rot_mat( numeric::z_rotation_matrix_degrees(-1. * z) );
    core::Vector v(0,0,0);
    rotated_pose.apply_transform_Rx_plus_v(x_rot_mat, v);
    rotated_pose.apply_transform_Rx_plus_v(y_rot_mat, v);
    rotated_pose.apply_transform_Rx_plus_v(z_rot_mat, v);

    std::stringstream outpdb;
    outpdb << tag.str() << ".pdb";
    rotated_pose.dump_pdb(outpdb.str());

    std::vector< conformation::ResidueOP > residues = protocols::pockets::PocketGrid::getRelaxResidues(rotated_pose, resid);
    protocols::pockets::PocketGrid pg( residues );


    pg.autoexpanding_pocket_eval( residues, rotated_pose );


    //pg.dumpGridToFile();

    protocols::pockets::EggshellGrid comparison_eggshell_grid( pg );

    // dump to files
    std::stringstream outpdb_pfname;
    outpdb_pfname << tag.str()<< ".pocket.pdb";
    pg.dumpTargetPocketsToPDB( outpdb_pfname.str() );
    std::stringstream out_pfname;
    out_pfname << tag.str()<< ".pocket";
    pg.dumpTargetPocketsToFile( out_pfname.str() );
    std::stringstream out_fname;
    out_fname << tag.str()<< ".eggshell";
    comparison_eggshell_grid.dump_eggshell( out_fname.str() );
    std::stringstream outpdb_fname;
    outpdb_fname << tag.str()<< ".eggshell.pdb";
    comparison_eggshell_grid.write_eggshell_to_pdb( outpdb_fname.str() );

    std::stringstream outpdb_pfname2;
    outpdb_pfname2 << tag2.str()<< ".pocket.pdb";
    pg.dumpTargetPocketsToPDB( outpdb_pfname2.str(), z_inv_rot_mat, y_inv_rot_mat, x_inv_rot_mat );
    std::stringstream out_pfname2;
    out_pfname2 << tag2.str()<< ".pocket";
    pg.dumpTargetPocketsToFile( out_pfname2.str(), z_inv_rot_mat, y_inv_rot_mat, x_inv_rot_mat );
    std::stringstream out_fname2;
    out_fname2 << tag2.str()<< ".eggshell";
    comparison_eggshell_grid.dump_eggshell( out_fname2.str(), z_inv_rot_mat, y_inv_rot_mat, x_inv_rot_mat );
    std::stringstream outpdb_fname2;
    outpdb_fname2 << tag2.str()<< ".eggshell.pdb";
    comparison_eggshell_grid.write_eggshell_to_pdb( outpdb_fname2.str(), z_inv_rot_mat, y_inv_rot_mat, x_inv_rot_mat );


  }

  } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
	return 0;
}



