// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

//GPU enabling is not default
//To test how many threads are fastest for your computer,
//use -gpu:threads 1024 (or other number) on the command line

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
#include <basic/options/option_macros.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>


using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( String, central_relax_pdb_num )
OPT_KEY( String, input_protein_file )
OPT_KEY( String, template_protein )
OPT_KEY( String, input_ligand_file )
OPT_KEY( String, known_ligand_file )
OPT_KEY( Integer, num_angles )
OPT_KEY( Boolean, cheat )
OPT_KEY( Boolean, trim_pocket )
OPT_KEY( Boolean, resize_adt_grid )
OPT_KEY( Boolean, print_pocket )
OPT_KEY( Boolean, print_fingerprints )
OPT_KEY( Boolean, print_output_complex )
OPT_KEY( Real, origin_cutoff )
OPT_KEY( Real, gc_x )
OPT_KEY( Real, gc_y )
OPT_KEY( Real, gc_z )
OPT_KEY( Real, gd_x )
OPT_KEY( Real, gd_y )
OPT_KEY( Real, gd_z )
OPT_KEY( Real, gs )

int main( int argc, char * argv [] ) {
	try{

  NEW_OPT( central_relax_pdb_num, "target residue", "-1" );
  NEW_OPT( input_protein_file, "protein file name", "protein.pdb" );
	NEW_OPT( template_protein, "template protein file name", "template.pdb" );
  NEW_OPT( input_ligand_file, "ligand file name", "ligand.pdb" );
  NEW_OPT( known_ligand_file, "known ligand file name", "known_ligand.pdb" );
	NEW_OPT( num_angles, "Number of different pose angles to measure score at", 1);
	NEW_OPT( cheat, "move pocket CoM over Ligand MCoM", false );
	NEW_OPT( trim_pocket, "trim the non-plaid pocket using a known ligand", false );
	NEW_OPT( resize_adt_grid, "resize grid based on user entered AUTODOCK grid values", false );
	NEW_OPT( print_pocket, "call 'dumpgridtofile()' to print pocket", false );
	NEW_OPT( print_fingerprints, "print fingerprint points into pdb files", false );
	//NEW_OPT( angle_increment, "angle increment", 20 );
  NEW_OPT( gc_x, "gid center : X ", 1.0 );
  NEW_OPT( gc_y, "gid center : Y ", 1.0 );
  NEW_OPT( gc_z, "gid center : Z ", 1.0 );
  NEW_OPT( gd_x, "gid dimension : X ", 20.0 );
  NEW_OPT( gd_y, "gid dimension : Y ", 20.0 );
  NEW_OPT( gd_z, "gid dimension : Z ", 20.0 );

	devel::init(argc, argv);
	std::string const input_protein = option[ input_protein_file ];
	std::string const template_pdb = option[ template_protein ];
	std::string const input_ligand = option[ input_ligand_file ];
	std::string const known_ligand = option[ known_ligand_file ];
	std::string const resid = option[ central_relax_pdb_num ];
  int angles = option[ num_angles ];
	//int const ang_inc  = option[ angle_increment ];
	core::Real const grid_cen_x = option[ gc_x ];
	core::Real const grid_cen_y = option[ gc_y ];
	core::Real const grid_cen_z = option[ gc_z ];
	core::Real grid_dim_x = option[ gd_x ];
	core::Real grid_dim_y = option[ gd_y ];
	core::Real grid_dim_z = option[ gd_z ];

  //create 'tag' for eggshell output filename
  int dot_index1 = input_protein.rfind(".", input_protein.size());
  assert(dot_index1 != -1 && "No dot found in filename");
	std::string protein_name = input_protein.substr(0,dot_index1);
	std::string eggshell_pdb_tag = "eggshell_" + protein_name + "_" + resid + ".pdb";
	std::string eggshell_triplet_tag = "eggshell_" + protein_name + "_" + resid + ".txt";

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

	if (option[ resize_adt_grid ]()){
		using namespace basic::options;
		core::Real const spacing = option[ OptionKeys::pocket_grid::pocket_grid_spacing ]();
		//Modify adt grid point values to use in DAVID's PocketGrid code
		grid_dim_x = (grid_dim_x/2) * spacing;
		grid_dim_y = (grid_dim_y/2) * spacing;
		grid_dim_z = (grid_dim_z/2) * spacing;
		protocols::pockets::PocketGrid	pg( grid_cen_x, grid_cen_y, grid_cen_z, grid_dim_x, grid_dim_y, grid_dim_z );
		numeric::xyzVector<core::Real> grid_center (0.);
		grid_center.x() = grid_cen_x;
		grid_center.y() = grid_cen_y;
		grid_center.z() = grid_cen_z;
		pg.DARC_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose, grid_center ) ;
		npf.setup_from_PocketGrid( protein_pose, pg );

	}	else {
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
				npf.setup_from_PocketGrid( best_pose, pg );
			}

		else if (angles == 1){
			protocols::pockets::PocketGrid	pg( protein_pose.conformation().residue(seqpos) );
			pg.autoexpanding_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose ) ;
			npf.setup_from_PocketGrid( protein_pose, pg );
		}
	}

	//print the eggshell pdb file
	npf.write_eggshell_to_pdb_file(eggshell_pdb_tag);
	npf.print_to_file(eggshell_triplet_tag);
	std::cout<< "Written eggshell to pdb file : "<< eggshell_pdb_tag << std::endl;
	std::cout<< "Written eggshell to triplet file: "<< eggshell_triplet_tag << std::endl;
	std::cout<< "DONE!"<< std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
	return 0;
}
