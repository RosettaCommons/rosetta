// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk + dj

#include <iostream>
#include <iomanip>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/scoring/rms_util.hh>

// Numeric Headers
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( String, contact_list )
OPT_KEY( Integer, num_angles )

static thread_local basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );
bool is_interface_residue (core::pose::Pose const & pose, core::Size resno);

//set to store pdb info keys
std::set <std::string> interface;
//stores resid of the ligand residue

void define_interface( core::pose::Pose & input_pose ) {
  core::Size lig_res_num = 0;
	for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j ) {
		if (!input_pose.residue(j).is_protein()){
			lig_res_num = j;
			break;
		}
		//TR<<input_pose.residue(j).type().name();
    //TR<<" ";
	}
  //TR << "\n";
	if (lig_res_num == 0){
		TR << "No ligand given in reference PDB structure.  Cannot identify interface."<<std::endl;
		exit (1);
	}

	//TR <<"sele ";

        scoring::ScoreFunctionOP scorefxn( get_score_function() );
        (*scorefxn)(input_pose);

	EnergyGraph & energy_graph(input_pose.energies().energy_graph());
	for ( graph::Graph::EdgeListIter
                                iru  = energy_graph.get_node( lig_res_num )->lower_edge_list_begin(),
                                irue = energy_graph.get_node( lig_res_num )->lower_edge_list_end();
                                iru != irue; ++iru ) {
		EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
		Size const j( edge->get_first_node_ind() );

		// the pair energies cached in the link
		EnergyMap const & emap( edge->fill_energy_map());
		Real const attr( emap[ fa_atr ] );
		//TR<<"\n"<<j<<": "<<attr<<"\n";
		if (attr < -.2){
			// create string id to store in set
			std::ostringstream residuestream;

			TR << "resi "<< input_pose.pdb_info()->number(j)<<" or ";

	        	residuestream << input_pose.pdb_info()->chain(j) << input_pose.pdb_info()->number(j);
			std::string res_id = residuestream.str();
			interface.insert(res_id);


		}
	}


	TR << std::endl;
	TR << lig_res_num<< std::endl;


	input_pose.delete_polymer_residue(lig_res_num);



}

/// General testing code
int
main( int argc, char * argv [] )
{
	try {

  NEW_OPT ( contact_list, "File name for optional list of contact residues to check","");
   NEW_OPT ( num_angles, "Number of different pose angles to measure score at", 1);
  devel::init(argc, argv);
  pose::Pose input_pose;
  int angles = option[ num_angles ];
  if (angles <1){
    fprintf (stderr, "Error: invalid number of angles.  Must be greather than 0\n");
    return -1;
  }


  //read in pdb file from command line
  std::string const input_pdb_name ( basic::options::start_file() );
  core::import_pose::pose_from_pdb( input_pose, input_pdb_name );


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


  }else{
    define_interface(input_pose);
  }

  std::filebuf fb,fb2;
  std::stringstream filename, filename2;
  if (!option[ OptionKeys::out::output_tag ]().empty()){
    filename<<option[ OptionKeys::out::output_tag ]()<<".pscore";
    filename2<<option[ OptionKeys::out::output_tag ]()<<".lpscore";
  }else{
    filename<<basic::options::start_file()<<".pscore";
    filename2<<basic::options::start_file()<<".lpscore";
  }
  fb.open (filename.str().c_str(),std::ios::out);
  fb2.open (filename2.str().c_str(),std::ios::out);
  std::ostream os(&fb);
  std::ostream os2(&fb2);

  //Loop over all residues and get score
  for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j )
  {
    if (is_interface_residue(input_pose, j)){
      core::Real constraint_pocket_score = 0;
      core::Real largest_pocket_score = 0;

      for (int i=0; i<angles; ++i){
        if (i>0){
          core::Real x,y,z;
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
        }

        protocols::pockets::PocketGrid pg( input_pose.conformation().residue(j) );

        if (pg.autoexpanding_pocket_eval( input_pose.conformation().residue(j), input_pose ) ){
          constraint_pocket_score += pg.netTargetPocketVolume();
          largest_pocket_score += pg.largestTargetPocketVolume();
        }
      }
      constraint_pocket_score /= angles;
      largest_pocket_score /= angles;
      os << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<<"\t" << constraint_pocket_score << std::endl;
      os2 << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<<"\t" << largest_pocket_score << std::endl;
      std::cout << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<<"\t" << constraint_pocket_score << std::endl;
      std::cout << input_pose.pdb_info()->chain(j)<<input_pose.pdb_info()->number(j)<<"\tLargest: " << largest_pocket_score << std::endl;

    }

  }

  fb.close();
  fb2.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;

}

bool is_interface_residue (core::pose::Pose const & pose, core::Size resno){
  std::ostringstream residuestream;
  residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
  std::string res_id = residuestream.str();
  //mjo commenting out 'rsd' because it is unused and causes a warning
	//core::conformation::Residue const & rsd = pose.residue(resno);
  if ( interface.count( res_id ) > 0 ) return true;
  return false;
}

