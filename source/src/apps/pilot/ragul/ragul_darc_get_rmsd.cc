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
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/after_opts.hh>

#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/random/random.hh>

#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>

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
OPT_KEY( Integer, num_poses )
OPT_KEY( Integer, num_angles )
OPT_KEY( Integer, num_runs )
OPT_KEY( Integer, num_particles )
OPT_KEY( Boolean, cheat )
OPT_KEY( Boolean, trim_pocket )
OPT_KEY( Boolean, print_pocket )
OPT_KEY( Boolean, print_fingerprints )
OPT_KEY( Boolean, print_output_complex )
OPT_KEY( Real, steric_weight )
OPT_KEY( Real, missing_point_weight )
OPT_KEY( Real, extra_point_weight )
OPT_KEY( Real, origin_cutoff )

int main( int argc, char * argv [] ) {
	try{

  NEW_OPT( central_relax_pdb_num, "target residue", "-1" );
  NEW_OPT( input_protein_file, "protein file name", "protein.pdb" );
	NEW_OPT( template_protein, "template protein file name", "template.pdb" );
  NEW_OPT( input_ligand_file, "ligand file name", "ligand.pdb" );
  NEW_OPT( known_ligand_file, "known ligand file name", "known_ligand.pdb" );
	NEW_OPT( num_poses, "No. of poses to search for initail stochastic search", 1000 );
	NEW_OPT( num_angles, "Number of different pose angles to measure score at", 1);
	NEW_OPT( num_runs, "no. of runs for PSO", 200 );
  NEW_OPT( num_particles, "no. of particles for PSO", 200 );
	NEW_OPT( cheat, "move pocket CoM over Ligand MCoM", false );
	NEW_OPT( trim_pocket, "trim the non-plaid pocket using a known ligand", false );
	NEW_OPT( print_pocket, "call 'dumpgridtofile()' to print pocket", false );
	NEW_OPT( print_fingerprints, "print fingerprint points into pdb files", false );
	NEW_OPT( print_output_complex, "print DARC output model protein-ligand complex as PDB file", false );
  NEW_OPT( steric_weight, "steric weight for PSO", 5.0 );
  NEW_OPT( missing_point_weight, "missing point weight", 20.0 );
  NEW_OPT( extra_point_weight, "extra point weight", 20.0 );
  NEW_OPT( origin_cutoff, "value for setting minimum and maximum origin cut off", 9.0 );
	//NEW_OPT( angle_increment, "angle increment", 20 );

	devel::init(argc, argv);

	std::string const input_protein = option[ input_protein_file ];
	std::string const template_pdb = option[ template_protein ];
	std::string const input_ligand = option[ input_ligand_file ];
	std::string const known_ligand = option[ known_ligand_file ];
	//	int num_pose_search  = option[ num_poses ];
	std::string const resid = option[ central_relax_pdb_num ];
  int angles = option[ num_angles ];
  int particle_size = option[ num_particles ];
  int run_size = option[ num_runs ];
	core::Real const steric_wt = option[ steric_weight ];
	core::Real const missing_pt_wt = option[ missing_point_weight ];
	core::Real const extra_pt_wt = option[ extra_point_weight ];
	core::Real const origin_space = option[ origin_cutoff ];
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
			npf.setup_from_PocketGrid( best_pose, pg );
			//print the pocket pdb file
			if (option[ print_pocket ]()){
				pg.dumpGridToFile();
			}//END printing pocket
		}
	else if (angles == 1){
		protocols::pockets::PocketGrid	pg( protein_pose.conformation().residue(seqpos) );
		pg.autoexpanding_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose ) ;
		npf.setup_from_PocketGrid( protein_pose, pg );
		//print the pocket pdb file
		if (option[ print_pocket ]()){
			pg.dumpGridToFile();
		}//END printing pocket
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


	protocols::pockets::PlaidFingerprint pf( small_mol_pose, npf );

	numeric::xyzVector<core::Real> pocket_CoM = npf.CoM();
	numeric::xyzVector<core::Real> smallmol_CoM = pf.CoM();

	//create 'tag' for output filenames
	int dot_index1 = input_protein.rfind(".", input_protein.size());
  assert(dot_index1 != -1 && "No dot found in filename");
	std::string protein_name = input_protein.substr(0,dot_index1);
  int dot_index2 = input_ligand.rfind(".", input_ligand.size());
  assert(dot_index2 != -1 && "No dot found in filename");
	std::string ligand_name = input_ligand.substr(0,dot_index2);
	std::string tag = ligand_name + "_" + protein_name + "_" + resid;

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
	std::cout<< "JK WARNING this code is not yet conformer-enabled, it will arbitrarily use the first conformer. Fix it in the app by removing the 0 in dump_oriented_pose_and_fp_to_pdb below..." << std::endl;
	pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, 0., 0., 0., original_pocket_angle_transform, 0 );
	}//END printing unaligned fingerprints

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


	ParticleOPs particles;
	std::cout<< "JK WARNING this code is not yet conformer-enabled, it will arbitrarily use the first conformer. Fix it in the app by removing the 1 in FingerprintMultifunc constructor below..." << std::endl;
	protocols::pockets::FingerprintMultifunc fpm(npf, pf, missing_pt_wt, steric_wt, extra_pt_wt, 1);
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
		std::cout<< "JK WARNING this code is not yet conformer-enabled, fix it in the app by removing the zero in the call to dump_oriented_pose_and_fp_to_pdb below..." << std::endl;
		pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, best_vars[4], best_vars[5], best_vars[6], original_pocket_angle_transform, optimized_origin, 0 );
	}//END printing PSO fingerprints

	//calculate RMSD and print into a file
	std::cout<< "JK WARNING this code is not yet conformer-enabled, fix it in the app by removing the zero in the call to get_oriented_pose below..." << std::endl;
	core::pose::Pose oriented_pose = pf.get_oriented_pose(npf, best_vars[4], best_vars[5], best_vars[6], original_pocket_angle_transform, optimized_origin, 0 );
	std::string pso_pose_name = "LIGAND_" + tag + ".pdb";
	oriented_pose.dump_pdb(pso_pose_name);
	core::Real rmsd_value = pf.rmsd(original_pose, oriented_pose);
	std::string rmsd_file = "rmsd.out";
	ofstream RMSDfile(rmsd_file.c_str(), ofstream::app);
	RMSDfile <<rmsd_value<<"\n";
	RMSDfile.close();


	//print protein-ligand complex into a pdb file
	if (option[ print_output_complex ]()){
		std::string Plineread;
		std::string Llineread;
		ofstream PLfile(complex_filename.c_str());
		//ifstream Pfile(complex_filename.c_str());
		ifstream Pfile(input_protein.c_str());
		if (!Pfile) {
			std::cout<< "Can't open Protein-pose file " << input_protein << std::endl;
			exit(1);
		}
		while (std::getline(Pfile, Plineread)) {
			if (Plineread[0] == 'E' && Plineread[1] == 'N' && Plineread[2] == 'D') continue;
			if (Plineread[0] == 'T' && Plineread[1] == 'E' && Plineread[2] == 'R') continue;
			PLfile << Plineread<<"\n";
		}
		PLfile <<"TER\n";
		ifstream Lfile(pso_pose_name.c_str());
		if (!Lfile) {
			std::cout<< "Can't open Ligand-pose file " << pso_pose_name << std::endl;
			exit(1);
		}
		while (std::getline(Lfile, Llineread)) {
			if (Llineread[0] == 'E' && Llineread[1] == 'N' && Llineread[2] == 'D') continue;
			if (Llineread[0] == 'T' && Llineread[1] == 'E' && Llineread[2] == 'R') continue;
			PLfile << Llineread<<"\n";
		}
		PLfile <<"END\n";
		Pfile.close();
		Lfile.close();
		PLfile.close();
	}	//END printing DARC_complex

	//delete the optimized ligand pose pdb file
	if(	remove(pso_pose_name.c_str()) != 0) perror( "Error deleting Ligand_PSO_pose_pdb file" );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
	return 0;

}
