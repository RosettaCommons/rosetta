// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>

#include <protocols/moves/ScoreMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>


#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
// Utility Headers
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
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

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/moves/SuperimposeMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

using namespace core;
using namespace std;
using namespace core::id;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::optimization;
using namespace core::pose::datacache;
using namespace core::pose::metrics;
using namespace core::scoring::constraints;
using namespace conformation;
using namespace protocols::moves;

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
OPT_KEY( Boolean, pre_align )
OPT_KEY( Boolean, trim_pocket )
OPT_KEY( Boolean, print_pocket )
OPT_KEY( Boolean, print_fingerprints )
OPT_KEY( Boolean, print_output_complex )
OPT_KEY( Real, steric_weight )
OPT_KEY( Real, missing_point_weight )
OPT_KEY( Real, extra_point_weight )
OPT_KEY( Real, origin_cutoff )
OPT_KEY( Real, cst_force_constant )

int main( int argc, char * argv [] ) {
	try{

  NEW_OPT( central_relax_pdb_num, "target residue", "-1" );
  NEW_OPT( input_protein_file, "protein file name", "protein.pdb" );
	NEW_OPT( template_protein, "template protein file name", "template.pdb" );
  NEW_OPT( input_ligand_file, "ligand file name", "ligand.pdb" );
  NEW_OPT( known_ligand_file, "known ligand file name", "known_ligand.pdb" );
	NEW_OPT( num_poses, "No. of poses to search for initail stochastic search", 1000 );
	NEW_OPT( num_angles, "Number of different pose angles to measure score at", 1);
	NEW_OPT( num_runs, "no. of runs for PSO", 100 );
  NEW_OPT( num_particles, "no. of particles for PSO", 100 );
	NEW_OPT( cheat, "move pocket CoM over Ligand MCoM", false );
	NEW_OPT( pre_align, "align protein to template before DARC", false );
	NEW_OPT( trim_pocket, "trim the non-plaid pocket using a known ligand", false );
	NEW_OPT( print_pocket, "call 'dumpgridtofile()' to print pocket", false );
	NEW_OPT( print_fingerprints, "print fingerprint points into pdb files", false );
	NEW_OPT( print_output_complex, "print DARC output model protein-ligand complex as PDB file", true );
  NEW_OPT( steric_weight, "steric weight for PSO", 5.0 );
  NEW_OPT( missing_point_weight, "missing point weight", 20.0 );
  NEW_OPT( extra_point_weight, "extra point weight", 20.0 );
  NEW_OPT( origin_cutoff, "value for setting minimum and maximum origin cut off", 7.0 );
	NEW_OPT( cst_force_constant, "coordinate constraint force constant", 0.5 );

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
	core::import_pose::pose_from_file( protein_pose, input_protein , core::import_pose::PDB_file);

	if (option[ pre_align ]()){
		pose::Pose template_pose;
		core::import_pose::pose_from_file( template_pose, template_pdb , core::import_pose::PDB_file);
		protocols::moves::SuperimposeMoverOP sp_mover = new protocols::moves::SuperimposeMover( template_pose );
		sp_mover->apply( protein_pose );
		std::string aligned_pose_name = "aligned_pose.pdb";
		protein_pose.dump_pdb(aligned_pose_name);
	}

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
  for ( int j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
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
		core::import_pose::pose_from_file( known_ligand_pose, known_ligand , core::import_pose::PDB_file);
		core::Size lig_res_num = 0;
		for ( int j = 1, resnum = known_ligand_pose.size(); j <= resnum; ++j ) {
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
	core::import_pose::pose_from_file( small_mol_pose, input_ligand , core::import_pose::PDB_file);
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
		for ( int j = 1, resnum = small_mol_pose.size(); j <= resnum; ++j ) {
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
	pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, 0., 0., 0., original_pocket_angle_transform );
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
		pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, best_vars[4], best_vars[5], best_vars[6], original_pocket_angle_transform, optimized_origin );
	}//END printing PSO fingerprints

	if (option[ print_output_complex ]()){
		//print protein-ligand complex into a pdb file
		core::pose::Pose oriented_pose = pf.get_oriented_pose(npf, best_vars[4], best_vars[5], best_vars[6], original_pocket_angle_transform, optimized_origin );
		std::string pso_pose_name = "LIGAND_" + tag + ".pdb";
		oriented_pose.dump_pdb(pso_pose_name);
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
		//delete the optimized ligand pose pdb file
	if(	remove(pso_pose_name.c_str()) != 0) perror( "Error deleting Ligand_PSO_pose_pdb file" );
	}	//END printing DARC_complex


	//STARTING MINIMIZATION

	std::string mini_complex_pdb = "mini_" + tag + ".pdb";

	//setup scorefxn
	scoring::ScoreFunctionOP scorefxn(get_score_function());
	scoring::ScoreFunctionOP repack_scorefxn(get_score_function());

	//setup the bound pose
	pose::Pose bound_pose;
	core::import_pose::pose_from_file( bound_pose, complex_filename , core::import_pose::PDB_file);

	//Apply constraint
  if ( bound_pose.residue( bound_pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		bound_pose.append_residue_by_jump
      ( *ResidueFactory::create_residue( bound_pose.residue(1).residue_type_set().name_map( "VRT" ) ),
        bound_pose.size()/2 );
  }

  Size nres = bound_pose.size();
  Real const coord_sdev( option[ OptionKeys::cst_force_constant ] );
  // default is 0.5 (from idealize) -- maybe too small
  for ( Size i = 1; i<= nres; ++i ) {
    if ( (Size)i==(Size)bound_pose.fold_tree().root() ) continue;
    Residue const & nat_i_rsd( bound_pose.residue(i) );
    for ( Size ii = 1; ii<= nat_i_rsd.nheavyatoms(); ++ii ) {
      bound_pose.add_constraint( new CoordinateConstraint(
																													AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
																													new HarmonicFunc( 0.0, coord_sdev ) ) );
    }
  }

  scorefxn->set_weight( coordinate_constraint, 0.5 );

  // setting degrees of freedom which can move during minimization - everything
	kinematics::MoveMap mm_all;
  mm_all.set_chi( true );
  mm_all.set_bb( true );
  mm_all.set_jump( true );

	// minimize protein
  AtomTreeMinimizer minimizer;
  MinimizerOptions min_options( "lbfgs_armijo_nonmonotone", 0.00001, true, false );

  minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
  minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
  (*scorefxn)(bound_pose);
	bound_pose.remove_constraints();
  (*scorefxn)(bound_pose);

  // Setup packer task for repacking
	pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( bound_pose ));
  base_packer_task->set_bump_check( false );
  base_packer_task->initialize_from_command_line();
  base_packer_task->or_include_current( true );

  for ( Size ii = 1; ii <= bound_pose.size(); ++ii ) {
    base_packer_task->nonconst_residue_task(ii).restrict_to_repacking();
  }

	// First repack
	pack::pack_rotamers( bound_pose, *repack_scorefxn, base_packer_task );


  // iterate over minimizing and repacking
  for ( Size iter = 1; iter <= 5; ++iter ) {
    minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
    minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
    minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
		pack::pack_rotamers( bound_pose, *repack_scorefxn, base_packer_task );
	}
	// final minimization
	minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
	minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
	minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
	minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
  (*scorefxn)(bound_pose);
  bound_pose.dump_scored_pdb( mini_complex_pdb, *scorefxn );

	//create ligand file from the minimized protein-ligand complex

	std::string MPLlineread;
	std::string Mini_Lig_File = "mini_" + ligand_name + ".pdb";
	ofstream MLfile(Mini_Lig_File.c_str());
	ifstream MPLfile(mini_complex_pdb.c_str());
	while (std::getline(MPLfile, MPLlineread)) {
		if (MPLlineread[0] == 'H' && MPLlineread[1] == 'E' && MPLlineread[2] == 'T' && MPLlineread[3] == 'A' && MPLlineread[4] == 'T' && MPLlineread[5] == 'M'){
			MLfile << MPLlineread<<"\n";
		}
	}
		MLfile <<"END\n";
	// close files
		MLfile.close();
		MPLfile.close();

		//calculate RMSD and print into a file
		pose::Pose minimized_ligand_pose;
		core::import_pose::pose_from_file( minimized_ligand_pose, Mini_Lig_File , core::import_pose::PDB_file);
		core::Real rmsd_value = pf.rmsd(original_pose, minimized_ligand_pose);
		std::string rmsd_file = "rmsd.out";
		ofstream RMSDfile(rmsd_file.c_str(), ofstream::app);
		RMSDfile <<pf.rmsd(original_pose, minimized_ligand_pose)<<"\n";
		RMSDfile.close();

    } catch (utility::excn::Exception const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
	return 0;

}
