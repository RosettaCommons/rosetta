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
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <protocols/pockets/DarcParticleSwarmMinimizer.hh>
#include <core/optimization/Minimizer.hh>
#include <basic/options/option_macros.hh>
#include <protocols/pockets/FingerprintMultifunc.hh>

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
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <utility/options/StringOption.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( String, optimization_method )
OPT_KEY( String, central_relax_pdb_num )
OPT_KEY( String, input_protein_file )
OPT_KEY( String, template_protein )
OPT_KEY( String, input_ligand_file )
OPT_KEY( String, known_ligand_file )
OPT_KEY( Integer, num_poses )
OPT_KEY( Integer, num_angles )
OPT_KEY( Integer, num_runs )
OPT_KEY( Integer, num_particles )
OPT_KEY( Boolean, cheat )
OPT_KEY( Boolean, trim_pocket )
OPT_KEY( Boolean, resize_adt_grid )
OPT_KEY( Boolean, score_only )
OPT_KEY( Boolean, print_pocket )
OPT_KEY( Boolean, print_fingerprints )
OPT_KEY( Boolean, print_output_complex )
OPT_KEY( Real, steric_weight )
OPT_KEY( Real, missing_point_weight )
OPT_KEY( Real, extra_point_weight )
OPT_KEY( Real, origin_cutoff )
OPT_KEY( Real, gc_x )
OPT_KEY( Real, gc_y )
OPT_KEY( Real, gc_z )
OPT_KEY( Real, gd_x )
OPT_KEY( Real, gd_y )
OPT_KEY( Real, gd_z )
OPT_KEY( Real, gs )
OPT_KEY( String, eggshell_triplet )
OPT_KEY( String, eggshell_pdb )
OPT_KEY( String, bound_ligand )
OPT_KEY( Integer, add_grid_size )
OPT_KEY( Integer, add_ext_grid_size )

int main( int argc, char * argv [] ) {
	try{

	NEW_OPT( optimization_method, "optimization_method", "PSO" );
  NEW_OPT( central_relax_pdb_num, "target residue", "-1" );
  NEW_OPT( input_protein_file, "protein file name", "protein.pdb" );
	NEW_OPT( template_protein, "template protein file name", "template.pdb" );
  NEW_OPT( input_ligand_file, "ligand file name", "ligand.pdb" );
  NEW_OPT( known_ligand_file, "known ligand file name", "known_ligand.pdb" );
	NEW_OPT( num_poses, "No. of poses to search for initail stochastic search", 100 );
	NEW_OPT( num_angles, "Number of different pose angles to measure score at", 1);
	NEW_OPT( num_runs, "no. of runs for PSO", 200 );
  NEW_OPT( num_particles, "no. of particles for PSO", 200 );
	NEW_OPT( cheat, "move pocket CoM over Ligand MCoM", false );
	NEW_OPT( trim_pocket, "trim the non-plaid pocket using a known ligand", false );
	NEW_OPT( resize_adt_grid, "resize grid based on user entered AUTODOCK grid values", false );
	NEW_OPT( score_only, "calculate DARC score without docking", false );
	NEW_OPT( print_pocket, "call 'dumpgridtofile()' to print pocket", false );
  NEW_OPT( print_output_complex, "print DARC output ligand model with protein as PDB file", true );
	NEW_OPT( print_fingerprints, "print fingerprint points into pdb files", false );
  NEW_OPT( missing_point_weight, "missing point weight", 21.616 );
  NEW_OPT( extra_point_weight, "extra point weight", 9.58689 );
  NEW_OPT( steric_weight, "steric weight for PSO", 1.43091 );
  NEW_OPT( origin_cutoff, "value for setting minimum and maximum origin cut off", 5.0 );
	//NEW_OPT( angle_increment, "angle increment", 20 );
  NEW_OPT( gc_x, "gid center : X ", 1.0 );
  NEW_OPT( gc_y, "gid center : Y ", 1.0 );
  NEW_OPT( gc_z, "gid center : Z ", 1.0 );
  NEW_OPT( gd_x, "gid dimension : X ", 20.0 );
  NEW_OPT( gd_y, "gid dimension : Y ", 20.0 );
  NEW_OPT( gd_z, "gid dimension : Z ", 20.0 );
	NEW_OPT( eggshell_pdb, "use input eggshell pdb file insted of generating new eggshell", "" );
	NEW_OPT( eggshell_triplet, "use input eggshell triplet file insted of generating new eggshell", "" );
  NEW_OPT( bound_ligand, "use bound ligand to set the grid for generating eggshell", "" );
	NEW_OPT( add_grid_size, "add grid dimension along x,y,z axis", 2 );
	NEW_OPT( add_ext_grid_size, "add extra grid dimension along x,y,z axis", 2 );

	devel::init(argc, argv);
	std::string const optimization_type = option[ optimization_method ];
	std::string const input_protein = option[ input_protein_file ];
	std::string const template_pdb = option[ template_protein ];
	std::string const input_ligand = option[ input_ligand_file ];
	std::string const known_ligand = option[ known_ligand_file ];
	std::string const bound_ligand_file = option[ bound_ligand ];
	int num_pose_search  = option[ num_poses ];
	std::string const resid = option[ central_relax_pdb_num ];
  int angles = option[ num_angles ];
  int particle_size = option[ num_particles ];
  int run_size = option[ num_runs ];
	core::Real const steric_wt = option[ steric_weight ];
	core::Real const missing_pt_wt = option[ missing_point_weight ];
	core::Real const extra_pt_wt = option[ extra_point_weight ];
	core::Real const origin_space = option[ origin_cutoff ];
	//int const ang_inc  = option[ angle_increment ];
	core::Real const grid_cen_x = option[ gc_x ];
	core::Real const grid_cen_y = option[ gc_y ];
	core::Real const grid_cen_z = option[ gc_z ];
	core::Real grid_dim_x = option[ gd_x ];
	core::Real grid_dim_y = option[ gd_y ];
	core::Real grid_dim_z = option[ gd_z ];
	std::string const input_eggshell_pdb = option[ eggshell_pdb ];
	std::string const input_eggshell_triplet = option[ eggshell_triplet ];
	int add_grid_dim  = option[ add_grid_size ];
	int add_ext_grid_dim  = option[ add_ext_grid_size ];

	using namespace basic::options;
	core::Real const spacing = option[ OptionKeys::pocket_grid::pocket_grid_spacing ]();

	protocols::pockets::NonPlaidFingerprint npf;
	pose::Pose protein_pose;
	core::import_pose::pose_from_pdb( protein_pose, input_protein );
	utility::vector1<core::Real> original_pocket_angle_transform(3, 0.);

	//find sequence position for the target residue
	//(not needed if using eggshell pdb or triplet file to setup pocket)
	int seqpos = 0;
	if ( (input_eggshell_triplet.empty()) && (input_eggshell_pdb.empty()) ) {
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
		for ( int j = 1, resnum = protein_pose.total_residue(); j <= resnum; ++j ) {
    if ( protein_pose.pdb_info()->number(j) == central_relax_residue_number ) {
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
	}

	//use input eggshell triplet text file to setup pocket
	if (!input_eggshell_triplet.empty()){
		npf.setup_from_eggshell_triplet_file( input_eggshell_triplet );
	}

	//use input eggshell PDB file to setup pocket
	else if (!input_eggshell_pdb.empty()){
		//NOTE:: currently not working for trimmed pocket
		npf.setup_from_eggshell_pdb_file( input_eggshell_pdb );
	}

	//use autodock grid values for grid size and setup PocketGrid
	else if (option[ resize_adt_grid ]()){
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
	}

	//use bound ligand to get grid size and setup PocketGrid
	else if (!bound_ligand_file.empty()){
		pose::Pose bound_ligand_pose;
		core::import_pose::pose_from_pdb( bound_ligand_pose, bound_ligand_file );
		core::Size lig_res_num = 0;
		for ( int j = 1, resnum = bound_ligand_pose.total_residue(); j <= resnum; ++j ) {
			if (!bound_ligand_pose.residue(j).is_protein()){
				lig_res_num = j;
				break;
			}
		}
		if (lig_res_num == 0){
			std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
			exit(1);
		}
		numeric::xyzVector<core::Real> input_ligand_CoM(0.);
		conformation::Residue const & curr_rsd = bound_ligand_pose.conformation().residue(lig_res_num);
		core::Real minx(999.), miny(999.), minz(999.), maxx(-999.), maxy(-999.), maxz(-999.);
		for(Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i) {
			if (curr_rsd.atom(i).xyz()(1) > maxx){maxx = curr_rsd.atom(i).xyz()(1);}
			if (curr_rsd.atom(i).xyz()(1) < minx){minx = curr_rsd.atom(i).xyz()(1);}
			if (curr_rsd.atom(i).xyz()(2) > maxy){maxy = curr_rsd.atom(i).xyz()(2);}
			if (curr_rsd.atom(i).xyz()(2) < miny){miny = curr_rsd.atom(i).xyz()(2);}
			if (curr_rsd.atom(i).xyz()(3) > maxz){maxz = curr_rsd.atom(i).xyz()(3);}
			if (curr_rsd.atom(i).xyz()(3) < minz){minz = curr_rsd.atom(i).xyz()(3);}
		}
		core::Real x_from_grd_cen, y_from_grd_cen, z_from_grd_cen, x_halfwidth, y_halfwidth, z_halfwidth;
		core::Real const cen_x = (maxx + minx)/2;
		core::Real const cen_y = (maxy + miny)/2;
		core::Real const cen_z = (maxz + minz)/2;
		x_halfwidth = std::abs(maxx - minx)/2;
		y_halfwidth = std::abs(maxy - miny)/2;
		z_halfwidth = std::abs(maxz - minz)/2;
		numeric::xyzVector<core::Real> const grid_center(cen_x,cen_y,cen_z);
		x_from_grd_cen = x_halfwidth + add_grid_dim;
		y_from_grd_cen = y_halfwidth + add_grid_dim;
		z_from_grd_cen = z_halfwidth + add_grid_dim;

		protocols::pockets::PocketGrid	pg( cen_x, cen_y, cen_z, x_from_grd_cen, y_from_grd_cen, z_from_grd_cen );
		pg.DARC_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose, grid_center ) ;
		x_from_grd_cen = x_halfwidth + add_grid_dim + add_ext_grid_dim;
		y_from_grd_cen = y_halfwidth + add_grid_dim + add_ext_grid_dim;
		z_from_grd_cen = z_halfwidth + add_grid_dim + add_ext_grid_dim;
		protocols::pockets::PocketGrid	ext_grd( cen_x, cen_y, cen_z, x_from_grd_cen, y_from_grd_cen, z_from_grd_cen );
		ext_grd.DARC_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose, grid_center ) ;
		npf.setup_from_PocketGrid( protein_pose, pg, ext_grd );
	}

	//use default grid size centered around target residue to setup pocketgrid
	else if (angles <1){
		//no. grid rotation cant be < 1
		fprintf (stderr, "Error: invalid number of angles.  Must be greather than 0\n");
		return -1;
	}
	//rotate grid to choose best pocket volume to choose largest pocket
	else if (angles > 1){
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
	//no grid rotation
	else if (angles == 1){
		protocols::pockets::PocketGrid	pg( protein_pose.conformation().residue(seqpos) );
		pg.autoexpanding_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose ) ;
		npf.setup_from_PocketGrid( protein_pose, pg );
	}

	//trim pocket based on known ligand
	if (option[ trim_pocket ]()){
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
		//calc lig_COM and move pock_COM to lig_com of known ligand
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

	pose::Pose small_mol_pose;
	core::import_pose::pose_from_pdb( small_mol_pose, input_ligand );
	core::pose::Pose original_pose = small_mol_pose;

	numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_radians(original_pocket_angle_transform[1] ) );
	numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_radians(original_pocket_angle_transform[2] ) );
	numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_radians(original_pocket_angle_transform[3] ) );
	numeric::xyzMatrix<core::Real> bestxyz_rot_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
	core::Vector v(0,0,0);
	small_mol_pose.apply_transform_Rx_plus_v(bestxyz_rot_mat, v);

	//CHEAT! Calculate CoM of Ligand, move Pocket COM to COM of input_ligand
	if (option[ cheat ]()){
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
	}//END CHEAT!

	//setup ligand fingerprint
	protocols::pockets::PlaidFingerprint pf( small_mol_pose, npf );

	//create 'tag' for output filenames
	int dot_index1 = input_protein.rfind(".", input_protein.size());
  assert(dot_index1 != -1 && "No dot found in filename");
	std::string protein_name = input_protein.substr(0,dot_index1);
  int dot_index2 = input_ligand.rfind(".", input_ligand.size());
  assert(dot_index2 != -1 && "No dot found in filename");
	std::string ligand_name = input_ligand.substr(0,dot_index2);
	std::string tag = ligand_name + "_" + protein_name + "_" + resid;
	std::string eggshell_pdb_tag = "eggshell_" + protein_name + "_" + resid + ".pdb";
	std::string eggshell_triplet_tag = "eggshell_" + protein_name + "_" + resid + ".txt";

	//print the eggshell pdb file
	using namespace basic::options;
	if (option[ OptionKeys::fingerprint::print_eggshell ]()){
		npf.write_eggshell_to_pdb_file(eggshell_pdb_tag);
		npf.print_to_file(eggshell_triplet_tag);
	}

	//print unaligned fingerprint files
	if (option[ print_fingerprints ]()){
		std::string np_output_filename = "npf_" + tag + ".txt";
		std::string np_output_pdbname = "npf_" + tag + ".pdb";
		npf.print_to_file(np_output_filename);
		npf.print_to_pdb(np_output_pdbname);
		std::string p_output_filename = "pf_" + tag + ".txt";
		std::string p_output_pdbname = "pf_" + tag + ".pdb";
		pf.print_to_file(p_output_filename);
		pf.print_to_pdb(p_output_pdbname);
		std::cout << "SCORE : unaligned  : " <<  pf.fp_compare( npf, missing_pt_wt, steric_wt, extra_pt_wt ) <<std::endl;
		std::string pose_name = "unal_pose_" + tag + ".pdb";
		std::string fp_name = "unal_fp_" + tag + ".pdb";
		std::cout<< "JK ERROR - this code is not yet conformer-enabled, fix it in the app by removing the 0 in dump_oriented_pose_and_fp_to_pdb below..." << std::endl;
		exit(1);
		pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, 0., 0., 0., original_pocket_angle_transform, 0 );
	}

	if (option[ score_only ]()){
		std::cout << "SCORE : unaligned  : " <<  pf.fp_compare( npf, missing_pt_wt, steric_wt, extra_pt_wt ) <<std::endl;
	}
	else if (optimization_type == "DFP"){
		core::Real best_score = std::numeric_limits<core::Real>::max();
		//core::Real dfp_score;
		utility::vector1<core::Real> out_vars(6);
		std::cout<< "JK ERROR - this code is not yet conformer-enabled, fix it in the app by removing the 1 in FingerprintMultifunc constructor below..." << std::endl;
		exit(1);
		protocols::pockets::FingerprintMultifunc fpm(npf, pf, missing_pt_wt, steric_wt, extra_pt_wt, 1);
		//core::optimization::MinimizerOptions options("dfpmin", 0.0000001, false, false, false);
		//core::optimization::MinimizerOptions options("dfpmin_armijo", 0.000001, false, false, false);
		core::optimization::MinimizerOptions options("dfpmin_armijo_nonmonotone", 0.000001, false, false, false);
		core::optimization::Minimizer dfp(fpm, options);
		//generate multiple stating vars
		for (int j = 0; j < num_pose_search; ++j ){
			utility::vector1<core::Real> p_opt(6, 0.);
			//			core::Real curr_CoM_offset_x = (int) (numeric::random::uniform() * origin_space);
			//			core::Real curr_CoM_offset_y = (int) (numeric::random::uniform() * origin_space);
			//			core::Real curr_CoM_offset_z = (int) (numeric::random::uniform() * origin_space);
			core::Real curr_angle1 = ( (int) (numeric::random::uniform() * 359.999) ) * numeric::constants::r::pi_over_180;
			core::Real curr_angle2 = ( (int) (numeric::random::uniform() * 359.999) ) * numeric::constants::r::pi_over_180;
			core::Real curr_angle3 = ( (int) (numeric::random::uniform() * 359.999) ) * numeric::constants::r::pi_over_180;
			//p_opt[1] = curr_CoM_offset_x;
			//p_opt[2] = curr_CoM_offset_y;
			//p_opt[3] = curr_CoM_offset_z;
			p_opt[4] = curr_angle1;
			p_opt[5] = curr_angle2;
			p_opt[6] = curr_angle3;
			core::Real curr_score = dfp.run(p_opt);
      std::cout<<"curr_score "<< curr_score <<std::endl;
			if ( curr_score < best_score ) {
				best_score = curr_score;
				std::cout<<"best_score "<< curr_score <<std::endl;
			}
		}
		std::cout<<"DARC_DFP_score : " << best_score <<std::endl;
		//dfp_score = dfp.run(p_opt);
		//fpm.dump(out_vars);
	}

	else if (optimization_type == "PSO"){
		utility::vector1<core::Real> p_min(6);
		p_min[1] = origin_space * -1;
		p_min[2] = origin_space * -1;
		p_min[3] = origin_space * -1;
		p_min[4] = 0.;
		p_min[5] = 0.;
		p_min[6] = 0.;
		utility::vector1<core::Real> p_max(6);
		p_max[1] = origin_space;
		p_max[2] = origin_space;
		p_max[3] = origin_space;
		p_max[4] = numeric::constants::r::pi_2;
		p_max[5] = numeric::constants::r::pi_2;
		p_max[6] = numeric::constants::r::pi_2;

#ifdef USEOPENCL
    npf.gpu_setup( missing_pt_wt, steric_wt, extra_pt_wt, particle_size, pf );
#endif
		ParticleOPs particles;
		std::cout<< "JK ERROR - this code is not yet conformer-enabled, fix it in the app by removing the 1 in FingerprintMultifunc constructor below..." << std::endl;
		exit(1);
		protocols::pockets::FingerprintMultifunc fpm(npf, pf, missing_pt_wt, steric_wt, extra_pt_wt, 1);
		//protocols::pockets::DarcParticleSwarmMinimizer pso( npf, pf, missing_pt_wt, steric_wt, extra_pt_wt, p_min, p_max);
		core::optimization::ParticleSwarmMinimizer pso(p_min, p_max);
		particles = pso.run(run_size, fpm, particle_size);
		ParticleOP p = particles[1];
		core::optimization::Particle parti(*p);
		core::Real fit_best = -(parti.fitness_pbest());
		utility::vector1<core::Real> best_vars(6);
		best_vars = parti.pbest();
		std::string complex_filename = "DARC_" + tag + ".pdb";
		std::cout<<"BEST FITNESS:"<<"	"<<complex_filename<<"	"<<fit_best<<std::endl;
		//fpm.dump(vars);
		//std::string header = "pso-optimization";
		//pso.print_particles(particles, header);
		//std::cout<<"BEST dofs ["<<tag<<"]:  [  "<<best_vars[1]<<",   "<<best_vars[2]<<",   "<<best_vars[3]<<",   "<<best_vars[4]<<",   "<<best_vars[5]<<",   "<<best_vars[6]<<"   ]"<<std::endl;
		numeric::xyzVector<core::Real> optimized_origin(0.);
		optimized_origin.x() = best_vars[1];
		optimized_origin.y() = best_vars[2];
		optimized_origin.z() = best_vars[3];
		//print PSO optimized fingerprints
		if (option[ print_fingerprints ]()){
			std::string pose_name = "pso_pose_" + tag + ".pdb";
			std::string fp_name = "pso_fp_" + tag + ".pdb";
			std::cout<< "JK ERROR - this code is not yet conformer-enabled, fix it in the app by removing the 0 in dump_oriented_pose_and_fp_to_pdb below..." << std::endl;
			exit(1);
			pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, best_vars[4], best_vars[5], best_vars[6], original_pocket_angle_transform, optimized_origin, 0 );
		}//END printing PSO fingerprints

		//Calculate RMSD
		std::cout<< "JK ERROR - this code is not yet conformer-enabled, fix it in the app by removing the zero in the call to get_oriented_pose below..." << std::endl;
		exit(1);
		core::pose::Pose oriented_pose = pf.get_oriented_pose(npf, best_vars[4], best_vars[5], best_vars[6], original_pocket_angle_transform, optimized_origin, 0 );
		std::string pso_pose_name = "LIGAND_" + tag + ".pdb";
		core::Real rmsd_value = pf.rmsd(original_pose, oriented_pose);
		std::cout<<"RMSD ["<<tag<<"]:"<<rmsd_value<<std::endl;
		std::string rmsd_file = "rmsd.out";
    ofstream RMSDfile(rmsd_file.c_str(), ofstream::app);
    RMSDfile <<rmsd_value<<"\n";
    RMSDfile.close();
		//print protein-ligand complex into a pdb file
		if (option[ print_output_complex ]()){
			protein_pose.append_residue_by_jump(oriented_pose.residue( 1 ), protein_pose.total_residue(),"", "",  true);
			protein_pose.dump_pdb(complex_filename);
		}
	}else {
		std::cout<<"ERROR! : Wrong optimization_method "<<std::endl;
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
	return 0;
}
