// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

// Protocol Headers
#include <devel/init.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <basic/options/option_macros.hh>
#include <protocols/pockets/FingerprintMultifunc.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>

#include <basic/options/option_macros.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;


OPT_KEY( String, central_relax_pdb_num )
OPT_KEY( String, input_protein_file )
OPT_KEY( String, input_ligand_file )
OPT_KEY( String, known_ligand_file )
OPT_KEY( Integer, num_poses )
OPT_KEY( Integer, num_angles )
OPT_KEY( Integer, num_runs )
OPT_KEY( Integer, num_particles )
OPT_KEY( Boolean, cheat )
OPT_KEY( Boolean, trim_pocket )
OPT_KEY( Real, steric_weight )
OPT_KEY( Real, missing_point_weight )
OPT_KEY( Real, extra_point_weight )
OPT_KEY( Real, origin_cutoff )

int main( int argc, char * argv [] ) {

  NEW_OPT( central_relax_pdb_num, "target residue", "-1" );
  NEW_OPT( input_protein_file, "protein file name", "protein.pdb" );
  NEW_OPT( input_ligand_file, "ligand file name", "ligand.pdb" );
  NEW_OPT( known_ligand_file, "known ligand file name", "known_ligand.pdb" );
	NEW_OPT( num_poses, "No. of poses to search for initail stochastic search", 1000 );
	NEW_OPT( num_angles, "Number of different pose angles to measure score at", 1);
	NEW_OPT( num_runs, "no. of runs for PSO", 100 );
  NEW_OPT( num_particles, "no. of particles for PSO", 100 );
	NEW_OPT( cheat, "move pocket CoM over Ligand MCoM", false );
	NEW_OPT( trim_pocket, "trim the non-plaid pocket using a known ligand", false );
  NEW_OPT( steric_weight, "steric weight for PSO", 5.0 );
  NEW_OPT( missing_point_weight, "missing point weight", 20.0 );
  NEW_OPT( extra_point_weight, "extra point weight", 20.0 );
  NEW_OPT( origin_cutoff, "value for setting minimum and maximum origin cut off", 7.0 );
	//NEW_OPT( angle_increment, "angle increment", 20 );

	devel::init(argc, argv);

	std::string const input_protein = option[ input_protein_file ];
	std::string const input_ligand = option[ input_ligand_file ];
	std::string const known_ligand = option[ known_ligand_file ];
	//	int num_pose_search  = option[ num_poses ];
	std::string const resid = option[ central_relax_pdb_num ];
  int angles = option[ num_angles ];
  //int particle_size = option[ num_particles ];
  //int run_size = option[ num_runs ];
	core::Real const steric_wt = option[ steric_weight ];
	core::Real const missing_pt_wt = option[ missing_point_weight ];
	core::Real const extra_pt_wt = option[ extra_point_weight ];
	//core::Real const origin_space = option[ origin_cutoff ];
	//int const ang_inc  = option[ angle_increment ];

	protocols::pockets::NonPlaidFingerprint npf;

	pose::Pose protein_pose;
	core::import_pose::pose_from_pdb( protein_pose, input_protein );

	int  central_relax_residue_number;
  char chain = ' ';
	std::size_t fpos( resid.find(':') );
  if ( fpos != std::string::npos ) {
    central_relax_residue_number = ObjexxFCL::int_of( resid.substr(0,fpos) );
    if (fpos != resid.size()-1 ) {
      chain = resid[ fpos+1 ];
    }
  } else {
    central_relax_residue_number = ObjexxFCL::int_of( resid );
  }
  int seqpos = 0;
  for ( int j = 1, resnum = protein_pose.total_residue(); j <= resnum; ++j ) {
    if ( protein_pose.pdb_info()->number(j) == central_relax_residue_number ) {
      //seqpos_ = j;
      if (chain != ' '){
        if ( protein_pose.pdb_info()->chain(j) == chain ) {
          seqpos = j;
        }
      }else{
        seqpos = j;
      }
    }
  }
  if ( seqpos == 0 ) {
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
    exit(1);
  }

  utility::vector1<core::Real> original_pocket_angle_transform(3, 0.);
	if(angles <1){
		fprintf (stderr, "Error: invalid number of angles.  Must be greather than 0\n");
    return -1;
  }else if (angles > 1)
		{
			core::Real best_vol(0), curr_vol(1);
			for (int i=0; i<angles; ++i){
				core::pose::Pose temp_pose;
				temp_pose = protein_pose;
				core::Real x = ( numeric::random::uniform() * numeric::constants::r::pi_2 ) + 0.0001;
				core::Real y = ( numeric::random::uniform() * numeric::constants::r::pi_2 ) + 0.0001;
				core::Real z = ( numeric::random::uniform() * numeric::constants::r::pi_2 ) + 0.0001;
				numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_radians(x) );
				numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_radians(y) );
				numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_radians(z) );
				numeric::xyzMatrix<core::Real> tot_rot_mat = z_rot_mat * y_rot_mat * x_rot_mat;
				core::Vector v(0,0,0);
				temp_pose.apply_transform_Rx_plus_v(tot_rot_mat, v);
				protocols::pockets::PocketGrid	pg( temp_pose.conformation().residue(seqpos) );
				pg.autoexpanding_pocket_eval( temp_pose.conformation().residue(seqpos), temp_pose ) ;
				curr_vol = pg.netTargetPocketVolume();
				std::cout<<"curr_volume "<<curr_vol<<std::endl;
				if(curr_vol > best_vol){
					best_vol = curr_vol;
					original_pocket_angle_transform[1] = x;
					original_pocket_angle_transform[2] = y;
					original_pocket_angle_transform[3] = z;
				}
			}

			numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_radians( original_pocket_angle_transform[1] ) );
			numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_radians( original_pocket_angle_transform[2] ) );
			numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_radians( original_pocket_angle_transform[3] ) );
			numeric::xyzMatrix<core::Real> bestxyz_rot_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
			core::Vector v(0,0,0);
			protein_pose.apply_transform_Rx_plus_v(bestxyz_rot_mat, v);

			core::pose::Pose best_pose;
			best_pose = protein_pose;

			protocols::pockets::PocketGrid	pg( best_pose.conformation().residue(seqpos) );

			pg.autoexpanding_pocket_eval( best_pose.conformation().residue(seqpos), best_pose ) ;
			std::cout<<"best_volume: "<<pg.netTargetPocketVolume()<<std::endl;
			//pg.dumpGridToFile();
			npf.setup_from_PocketGrid( best_pose, pg );

		}else if (angles == 1){
		protocols::pockets::PocketGrid	pg( protein_pose.conformation().residue(seqpos) );
		pg.autoexpanding_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose ) ;
		//pg.dumpGridToFile();
		npf.setup_from_PocketGrid( protein_pose, pg );
	}

	if (option[ trim_pocket ]()){
		//calc lig_COM and move pock_COM to lig_com of known ligand
		pose::Pose known_ligand_pose;
		core::import_pose::pose_from_pdb( known_ligand_pose, known_ligand );
		core::Size lig_res_num = 0;
		for ( int j = 1, resnum = known_ligand_pose.total_residue(); j <= resnum; ++j ) {
			if (!known_ligand_pose.residue(j).is_protein()){
				lig_res_num = j;
				break;
			}
		}
		if (lig_res_num == 0){
			std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
			exit(1);
		}
		numeric::xyzVector<core::Real> known_ligand_CoM(0.);
		conformation::Residue const & curr_rsd = known_ligand_pose.conformation().residue(lig_res_num);
		for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
			known_ligand_CoM.x() += curr_rsd.atom(i).xyz()(1);
			known_ligand_CoM.y() += curr_rsd.atom(i).xyz()(2);
			known_ligand_CoM.z() += curr_rsd.atom(i).xyz()(3);
		}
		known_ligand_CoM /= curr_rsd.nheavyatoms();
		npf.CHEAT_CoM( known_ligand_CoM );
		npf.trim_based_on_known_ligand(known_ligand_pose);
	}

  int dot_index1 = input_protein.rfind(".", input_protein.size());
  assert(dot_index1 != -1 && "No dot found in filename");
	std::string protein_name = input_protein.substr(0,dot_index1);
  int dot_index2 = input_ligand.rfind(".", input_ligand.size());
  assert(dot_index2 != -1 && "No dot found in filename");
	std::string ligand_name = input_ligand.substr(0,dot_index2);
	std::string tag = ligand_name + "_" + protein_name + "_" + resid;

	pose::Pose small_mol_pose;
	core::import_pose::pose_from_pdb( small_mol_pose, input_ligand );
	core::pose::Pose original_pose = small_mol_pose;
	numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_radians(original_pocket_angle_transform[1] ) );
	numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_radians(original_pocket_angle_transform[2] ) );
	numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_radians(original_pocket_angle_transform[3] ) );
	numeric::xyzMatrix<core::Real> bestxyz_rot_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
	core::Vector v(0,0,0);
	small_mol_pose.apply_transform_Rx_plus_v(bestxyz_rot_mat, v);

	if (option[ cheat ]()){
		//Calculate CoM of Ligand, move Pocket COM to COM of input_ligand
		core::Size lig_res_num = 0;
  for ( int j = 1, resnum = small_mol_pose.total_residue(); j <= resnum; ++j ) {
    if (!small_mol_pose.residue(j).is_protein()){
      lig_res_num = j;
      break;
    }
  }
  if (lig_res_num == 0){
		std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
    exit(1);
  }
	numeric::xyzVector<core::Real> input_ligand_CoM(0.);
	conformation::Residue const & curr_rsd = small_mol_pose.conformation().residue(lig_res_num);
  for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
    input_ligand_CoM.x() += curr_rsd.atom(i).xyz()(1);
    input_ligand_CoM.y() += curr_rsd.atom(i).xyz()(2);
    input_ligand_CoM.z() += curr_rsd.atom(i).xyz()(3);
  }
  input_ligand_CoM /= curr_rsd.nheavyatoms();
	npf.CHEAT_CoM( input_ligand_CoM );
	}

	protocols::pockets::PlaidFingerprint pf( small_mol_pose, npf );

	numeric::xyzVector<core::Real> pocket_CoM = npf.CoM();
	numeric::xyzVector<core::Real> smallmol_CoM = pf.CoM();

	std::cout << "DARC_SCORE : " <<  pf.fp_compare( npf, missing_pt_wt, steric_wt, extra_pt_wt ) <<std::endl;

	return 0;

}
