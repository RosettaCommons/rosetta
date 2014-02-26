// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

// Protocol Headers
#include <devel/init.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/FingerprintMultifunc.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <core/import_pose/import_pose.hh>

//reqd minimization headers
#include <protocols/simple_moves/ScoreMover.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <utility/file/file_sys_util.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/kinematics/FoldTree.hh>

using namespace std;
using namespace core;
using namespace basic::options;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

//reqd minimization namespace
using namespace core::pose::datacache;
using namespace core::optimization;
using namespace core::pose::metrics;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::id;
using namespace conformation;
using namespace protocols::simple_moves;
using namespace protocols::rigid;

OPT_KEY( String, protein )
OPT_KEY( StringVector, ligand )
OPT_KEY( String, ray_file )
OPT_KEY( Integer, num_runs )
OPT_KEY( Integer, num_particles )
OPT_KEY( Real, steric_weight )
OPT_KEY( Real, missing_point_weight )
OPT_KEY( Real, extra_point_weight )
OPT_KEY( Real, origin_cutoff )
OPT_KEY( Boolean, print_output_complex )
OPT_KEY( Boolean, search_conformers )
OPT_KEY( Boolean, darc_score_only )
OPT_KEY( Boolean, use_ligand_filename )
OPT_KEY( String, ligand_pdb_list)
//reqd minimization options
OPT_KEY( Boolean, minimize_output_complex )
OPT_KEY( Real, cst_force_constant )

static basic::Tracer TR( "apps.pilot.ragul_run_darc_with_input_eggshell.main" );

int main( int argc, char * argv [] ) {

  try {
  NEW_OPT( protein, "protein file name", "protein.pdb" );
  NEW_OPT( ligand, "input ligand(s)", "" );
  NEW_OPT( ray_file, "input eggshell(ray) triplet file name", "" );
  NEW_OPT( num_runs, "no. of runs for PSO", 100 );
  NEW_OPT( num_particles, "no. of particles for PSO", 100 );
  NEW_OPT( origin_cutoff, "value for setting minimum and maximum origin cut off", 7.0 );
  NEW_OPT( missing_point_weight, "missing point weight", 21.6 );
  NEW_OPT( extra_point_weight, "extra point weight", 9.5 );
  NEW_OPT( steric_weight, "steric weight for PSO", 1.4 );
  NEW_OPT( print_output_complex, "print DARC output ligand model with protein as PDB file", true );
  NEW_OPT( minimize_output_complex, "minimize the best scoring DARC output model", false );
  NEW_OPT( cst_force_constant, "coordinate constraint force constant", 0.5 );
  NEW_OPT( search_conformers, "optimize conformer during docking", true );
  NEW_OPT( darc_score_only, "DARC score only (no pso docking)", false );
  NEW_OPT( use_ligand_filename, "append ligand file name to output files, instead of ligand code", false );
	NEW_OPT( ligand_pdb_list, "List of pdb files of the ligand structures", "" );

	using utility::file::FileName;
	using utility::file::file_exists;

  devel::init(argc, argv);

  std::string const input_protein = option[ protein ];
  std::string const input_ligand = "temp";//change this
  std::string const input_eggshell_triplet = option[ ray_file ];
  int particle_size = option[ num_particles ];
  int run_size = option[ num_runs ];
  core::Real const steric_wt = option[ steric_weight ];
  core::Real const missing_pt_wt = option[ missing_point_weight ];
  core::Real const extra_pt_wt = option[ extra_point_weight ];
  core::Real const origin_space = option[ origin_cutoff ];

  protocols::pockets::NonPlaidFingerprint npf;
  pose::Pose protein_pose;
  core::import_pose::pose_from_pdb( protein_pose, input_protein );
  pose::Pose bound_pose = protein_pose;

  //use input eggshell triplet file to setup nonplaid fingerprint (pocket)
  if (!input_eggshell_triplet.empty()){
		std::cout<<"Reading eggshell"<<std::endl;
    npf.setup_from_eggshell_triplet_file( input_eggshell_triplet );
  }
  else if(input_eggshell_triplet.empty()){
    std::cout<<"ERROR! : no input eggshell file to setup pocket"<<std::endl;
    exit(1);
  }

	// Input ligands can be read in by several ways
	// For more than one ligand, a list of ligand pdb filenames can be stored in a text file and given as input using the flag -ligand_pdb_list
	// use the flag -extra_res_batch_path for input params file
	// For single ligand, just use the flag -ligand Ex: -ligand a.pdb
	// The -ligand flag can also be used for more than one ligand; Ex: -ligand a.pdb b.pdb (or) -ligand *.pdb
	utility::vector1< pose::Pose > ligand_poses;
	// read ligand(s) specified by the flag 'ligand'
	if (option[ ligand ].user()){
		utility::vector1<string> input_ligand_list = option[ ligand ]();
		for (core::Size f=1; f <= input_ligand_list.size(); f++) {
			std::string const input_ligand_name = input_ligand_list[f];
			std::cout<<"Reading ligand " << input_ligand_name <<std::endl;
			core::pose::Pose ligand_pose;
			core::import_pose::pose_from_pdb( ligand_pose, input_ligand_name );
			ligand_poses.push_back( ligand_pose );
		}
	}
	// read list file. open the file specified by the flag 'ligand_pdb_list' and read in all the lines in it
	else if (option[ ligand_pdb_list ].user()){
		std::vector< FileName > input_ligand_pdb_file_names;
		std::string ligand_pdb_list_file_name( option[ ligand_pdb_list ].value() );
		std::ifstream ligand_data( ligand_pdb_list_file_name.c_str() );
		std::string ligand_line;
		if ( !ligand_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + ligand_pdb_list_file_name + '\n' );
		}
		while( getline( ligand_data, ligand_line ) ) {
			input_ligand_pdb_file_names.push_back( FileName( ligand_line ) );
		}
		ligand_data.close();
		// iterate over FileName vector and read in the input LIGAND PDB files
		std::vector< FileName >::iterator ligand_pdb( input_ligand_pdb_file_names.begin() ), ligand_last_pdb(input_ligand_pdb_file_names.end());
		while ( ligand_pdb != ligand_last_pdb ) {
			// check to make sure the file exists
			if ( !file_exists( *ligand_pdb ) ) {
				utility_exit_with_message( "Ligand pdb " + std::string(*ligand_pdb) + " not found! skipping" );
			}
			std::cout<< "Reading in pose " << *ligand_pdb << std::endl;
			core::pose::Pose ligand_pose;
			core::import_pose::pose_from_pdb( ligand_pose, *ligand_pdb );
			ligand_poses.push_back( ligand_pose );
			ligand_pdb++;
		}
	}
  else if( !(option[ ligand ].user()) && !(option[ ligand_pdb_list ].user()) ){
		std::cout<<"ERROR! : no input ligand for DARC docking"<<std::endl;
		exit(1);
	}

	//open file 'darc_score.sc' to print all darc scores
	std::string outfname;
	if (!option[ OptionKeys::out::output_tag ]().empty()){
		outfname = "darc_score." + option[ OptionKeys::out::output_tag ]() + ".sc";
  }else{
    outfname = "darc_score.sc";
  }

	std::ofstream darc_score_file;
	darc_score_file.open(outfname.c_str(), std::ofstream::app);

	//loop through all ligand and get darc score
	for (core::Size lig_num = 1; lig_num <= ligand_poses.size(); ++lig_num) {
		core::pose::Pose small_mol_pose = ligand_poses[lig_num];
		Size const nconformers = small_mol_pose.total_residue();
		//create 'tag' for output filenames
		int pfounddir = input_protein.find_last_of("/\\");
		int pfounddot = input_protein.find_last_of(".");
		std::string protein_name = input_protein.substr((pfounddir+1),(pfounddot-(pfounddir+1)));//get the basename of the protein pdb file
		std::string ligand_code = small_mol_pose.residue( 1 ).name3();//get the 3 letter code of the ligand
		if (option[ use_ligand_filename ]) {
		  std::string ligand_file_name = small_mol_pose.pdb_info()->name();
		  int lfounddir = ligand_file_name.find_last_of("/\\");
		  int lfounddot = ligand_file_name.find_last_of(".");
		  std::string ligand_name = ligand_file_name.substr((lfounddir+1),(lfounddot-(lfounddir+1)));//get the basename of the ligand pdb file
		  ligand_code = ligand_name;
		}
		std::string tag = protein_name + "_" + ligand_code;
		std::string pso_pose_name = "LIGAND_" + tag + ".pdb";
		std::string mini_pso_pose_name = "mini_LIGAND_" + tag + ".pdb";
		std::string darc_complex_filename = "DARC_" + tag + ".pdb";
		std::string mini_complex_filename = "mini_" + tag + ".pdb";

		//set max and min ligand translation & rotation space
		utility::vector1<core::Real> p_min(7);
		p_min[1] = origin_space * -1;
		p_min[2] = origin_space * -1;
		p_min[3] = origin_space * -1;
		p_min[4] = 0.;
		p_min[5] = 0.;
		p_min[6] = 0.;
		p_min[7] = 0.;
		utility::vector1<core::Real> p_max(7);
		p_max[1] = origin_space;
		p_max[2] = origin_space;
		p_max[3] = origin_space;
		p_max[4] = numeric::constants::r::pi_2;
		p_max[5] = numeric::constants::r::pi_2;
		p_max[6] = numeric::constants::r::pi_2;
		p_max[7] = nconformers - 0.00001;		// note: conformer is indexed starting at 0

		if(option[ darc_score_only ]()){
			utility::vector1< core::pose::PoseOP > rot_poses(small_mol_pose.total_residue());
			for (Size ii = 1; ii <= small_mol_pose.total_residue(); ++ii){
				rot_poses[ii] = new core::pose::Pose(small_mol_pose, ii, ii);
				protocols::pockets::PlaidFingerprint conf_pf( *rot_poses[ii], npf );
				std::cout << "DARC SCORE (unaligned) :" <<  conf_pf.fp_compare( npf, missing_pt_wt, steric_wt, extra_pt_wt ) <<std::endl;
			}
			return 0;
		}

		//initialize particles, best DARC score, best optimized values
		std::cout<<"Starting PSO"<<std::endl;
		ParticleOPs particles;
		core::Real best_DARC_score = 9999.;
		utility::vector1<core::Real> best_vars(7);

		protocols::pockets::PlaidFingerprint pf(small_mol_pose, npf);
		core::pose::Pose oriented_pose;
		core::Size best_conformer; // note: indexed from zero

		//setup GPU
#ifdef USEOPENCL
		npf.gpu_setup( missing_pt_wt, steric_wt, extra_pt_wt, particle_size, pf );
#endif

		if (option[ search_conformers ]()){
			// search conformers during the PSO
			protocols::pockets::FingerprintMultifunc fpm(npf, pf, missing_pt_wt, steric_wt, extra_pt_wt, nconformers);
			protocols::pockets::DarcParticleSwarmMinimizer pso( npf, pf, missing_pt_wt, steric_wt, extra_pt_wt, p_min, p_max);
			particles = pso.run(run_size, fpm, particle_size);
			ParticleOP p = particles[1];
			core::optimization::Particle parti(*p);
			best_DARC_score = -(parti.fitness_pbest());
			best_vars = parti.pbest();
			// note: conformer is indexed starting at 0
			best_conformer=((core::Size)(floor(best_vars[7])) % nconformers);
			numeric::xyzVector<core::Real> optimized_origin(best_vars[1], best_vars[2], best_vars[3]);
			oriented_pose = pf.get_oriented_pose(npf, best_vars[4], best_vars[5], best_vars[6], optimized_origin, best_conformer );
		} else {
			// run DARC for each conformer separately, then take the best score at the end
			p_max[7] = 0.00001;
			utility::vector1< core::pose::PoseOP > rot_poses(small_mol_pose.total_residue());
			for (Size ii = 1; ii <= small_mol_pose.total_residue(); ++ii){
				rot_poses[ii] = new core::pose::Pose(small_mol_pose, ii, ii);
				protocols::pockets::PlaidFingerprint conf_pf( *rot_poses[ii], npf );
				//std::cout << "DARC SCORE (unaligned) :" <<  conf_pf.fp_compare( npf, missing_pt_wt, steric_wt, extra_pt_wt ) <<std::endl;
				protocols::pockets::FingerprintMultifunc fpm(npf, conf_pf, missing_pt_wt, steric_wt, extra_pt_wt, 1);
				protocols::pockets::DarcParticleSwarmMinimizer pso( npf, conf_pf, missing_pt_wt, steric_wt, extra_pt_wt, p_min, p_max);
				//core::optimization::ParticleSwarmMinimizer pso(p_min, p_max);
				particles = pso.run(run_size, fpm, particle_size);
				ParticleOP p = particles[1];
				core::optimization::Particle parti(*p);
				core::Real fit_best = -(parti.fitness_pbest());
				std::cout<< "SCORES : " << tag <<"	"<< "Conformer " << ii << " has DARC Score : " << fit_best <<std::endl;
				if(fit_best < best_DARC_score){
					best_DARC_score = fit_best;
					best_vars = parti.pbest();
					pf = conf_pf;
					numeric::xyzVector<core::Real> optimized_origin(best_vars[1], best_vars[2], best_vars[3]);
					oriented_pose = pf.get_oriented_pose(npf, best_vars[4], best_vars[5], best_vars[6], optimized_origin, 0 );
					best_conformer=ii-1; // note: index starts at 0
				}
			}
		}
		//Append the best ligand conformer pose to protein pose and print(optional) the complex PDB file
		bound_pose.append_residue_by_jump(oriented_pose.residue( 1 ), protein_pose.total_residue(),"", "",  true);
		pose::Pose pre_min_darc_pose = bound_pose;
		if (option[ print_output_complex ]()){
		  bound_pose.dump_pdb(darc_complex_filename);
		  oriented_pose.dump_pdb(pso_pose_name);
		}
		//print the best DARC score (without minimization)
		if (!option[ minimize_output_complex ]()){
		  // note: conformer is indexed starting at 0
		  darc_score_file << tag << "	" << best_DARC_score << "\n";
		  std::cout<< "SCORES : " << tag <<"	"<< "Conformer " << best_conformer << " has DARC Score : " << best_DARC_score <<std::endl;
		}

		//minimize the protein-ligand complex
		else if (option[ minimize_output_complex ]()){
		  //setup scorefxn
		  scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(core::scoring::PRE_TALARIS_2013_STANDARD_WTS, core::scoring::SCORE12_PATCH) );
		  scoring::ScoreFunctionOP repack_scorefxn( ScoreFunctionFactory::create_score_function(core::scoring::PRE_TALARIS_2013_STANDARD_WTS, core::scoring::SCORE12_PATCH) );

		  //Register calculators
			std::string sasa_calc_name = "sasa";
			std::string hbond_calc_name = "hbond";
			std::string packstat_calc_name = "packstat";
			std::string burunsat_calc_name = "burunsat";
			core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculatorLegacy;
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( sasa_calc_name, sasa_calculator );
			core::pose::metrics::PoseMetricCalculatorOP hb_calc = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( hbond_calc_name, hb_calc );
			core::pose::metrics::PoseMetricCalculatorOP packstat_calc =	new protocols::toolbox::pose_metric_calculators::PackstatCalculator();
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( packstat_calc_name, packstat_calc );
			core::pose::metrics::PoseMetricCalculatorOP burunsat_calc = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, hbond_calc_name);
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( burunsat_calc_name, burunsat_calc );

			// add constraint
			int nres = bound_pose.total_residue();
			Real coord_sdev( option[ OptionKeys::cst_force_constant ] );
			//take reciprocal and sqrt to pass as force constant
			coord_sdev = sqrt(1/coord_sdev);
			//std::cout<<" coord sdev "<< coord_sdev <<std::endl;
			Real cst_weight = 1;

			ConstraintSetOP cst_set( new ConstraintSet() );
			core::scoring::func::HarmonicFuncOP spring = new core::scoring::func::HarmonicFunc( 0 /*mean*/, coord_sdev /*std-dev*/);
			conformation::Conformation const & conformation( bound_pose.conformation() );

			// jk we need an anchor in order to use CoordinateConstraint !!!
			Size const my_anchor = 1;
			core::kinematics::FoldTree fold_tree=bound_pose.fold_tree();
			core::kinematics::FoldTree rerooted_fold_tree = fold_tree;
			rerooted_fold_tree.reorder( my_anchor );
			bound_pose.fold_tree( rerooted_fold_tree);

			for (int i=1; i <= nres; i++){
				//Residue const  & reside = pose.residue( i );
				Residue const & nat_i_rsd( bound_pose.residue(i) );
				for ( Size ii = 1; ii<= nat_i_rsd.nheavyatoms(); ++ii ) {
					AtomID CAi ( ii, i );
					cst_set->add_constraint
						(  new CoordinateConstraint
							 ( CAi, AtomID(1,my_anchor), conformation.xyz( CAi ), spring )
							 );
				}
			}
			bound_pose.constraint_set( cst_set );
			scorefxn->set_weight( coordinate_constraint, cst_weight );

			TR << "Starting minimization...." << std::endl;
			//AtomTreeMinimizer minimizer;
			AtomTreeMinimizer minimizer;
			MinimizerOptions min_options( "dfpmin", 0.00001, true, false );
			kinematics::MoveMap mm_all;
			mm_all.set_chi( true );
			mm_all.set_bb( true );
			mm_all.set_jump( true );
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			(*scorefxn)(bound_pose);
			TR << "Post minimization 1 constrained score: " << bound_pose.energies().total_energies()[ total_score ] << std::endl;
			//bound_pose.dump_scored_pdb( "1.pdb", *scorefxn );

			// Setup packer task for repacking
			pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( bound_pose ));
			base_packer_task->set_bump_check( false );
			base_packer_task->initialize_from_command_line();
			base_packer_task->or_include_current( true );
			for ( Size ii = 1; ii <= bound_pose.total_residue(); ++ii ) {
				base_packer_task->nonconst_residue_task(ii).restrict_to_repacking();
			}
			// First repack
			pack::pack_rotamers( bound_pose, *repack_scorefxn, base_packer_task );
			// Report Scores
			(*scorefxn)(bound_pose);
			//bound_pose.dump_scored_pdb( "repacked_once.pdb", *scorefxn );
			TR << "Score after repacking once: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

			// iterate over minimizing and repacking
			for ( Size iter = 1; iter <= 5; ++iter ) {
				cst_weight = cst_weight/2;
				scorefxn->set_weight( coordinate_constraint, cst_weight );
				minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
				minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
				minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
				(*scorefxn)(bound_pose);
				TR << "Current score after minimizing: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
				pack::pack_rotamers( bound_pose, *repack_scorefxn, base_packer_task );
				(*scorefxn)(bound_pose);
				TR << "Current score after repacking: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
			}

			// final minimization
			bound_pose.remove_constraints();
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			(*scorefxn)(bound_pose);

			if (option[ print_output_complex ]){
				//align minimized pose to the original docked pose and dump pdb complex and ligand
				protocols::simple_moves::SuperimposeMoverOP sp_mover = new protocols::simple_moves::SuperimposeMover();
				sp_mover->set_reference_pose( pre_min_darc_pose, 1, (pre_min_darc_pose.total_residue()-1) );
				sp_mover->set_target_range( 1, (bound_pose.total_residue()-1) );
				sp_mover->apply( bound_pose );
				//spilt into chains and print minimized ligand pose
				utility::vector1< core::pose::PoseOP > chains = bound_pose.split_by_chain();
				Size chain_num = chains.size();
				core::pose::Pose minimized_lig_pose = *chains[chain_num];
				minimized_lig_pose.dump_pdb(mini_pso_pose_name);
				bound_pose.dump_scored_pdb( mini_complex_filename, *scorefxn );
			}

			//  TR << "Final score: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
			//  TR << "Successfully finished minimizing input." << std::endl;

		//setup the unbound pose
		core::pose::Pose unbound_pose = bound_pose;
		core::Real const unbound_dist = 80.;
		Size const rb_jump = 1; // use the first jump as the one between partners
		protocols::rigid::RigidBodyTransMover trans_mover( unbound_pose, rb_jump );
		trans_mover.trans_axis( trans_mover.trans_axis() );
		trans_mover.step_size(unbound_dist);
		trans_mover.apply( unbound_pose );
		(*scorefxn)(unbound_pose);

		// define containers for metrics for total complex
		basic::MetricValue<Real> tot_sasa_mval;
		basic::MetricValue<Size> tot_hb_mval;
		basic::MetricValue<Real> tot_packstat_mval;
		basic::MetricValue<Size> tot_unsat_mval;

		// calculate and store total metrics for bound and unbound poses
		core::Real bound_energy = 0.0, unbound_energy = 0.0, Interface_Energy = 0.0;
		core::Size  bound_hb = 0,   unbound_hb = 0, Interface_HB = 0;
		core::Real bound_packstat = 0.0, unbound_packstat = 0.0, Total_packstats = 0.0;
		core::Size  bound_unsat = 0, unbound_unsat = 0, Interface_unsat = 0;

		//calculate interface Energy
		Interface_Energy = bound_energy - unbound_energy;

		//delta sasa calculation
		bound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
		unbound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);

		//interface hb calculation
		bound_pose.metric(hbond_calc_name,"all_Hbonds", tot_hb_mval);
		bound_hb = tot_hb_mval.value();
		unbound_pose.metric(hbond_calc_name,"all_Hbonds", tot_hb_mval);
		unbound_hb = tot_hb_mval.value();
		Interface_HB = bound_hb - unbound_hb;

		//packstat calculation
		bound_pose.metric(packstat_calc_name,"total_packstat", tot_packstat_mval);
		bound_packstat = tot_packstat_mval.value();
		unbound_pose.metric(packstat_calc_name,"total_packstat", tot_packstat_mval);
		unbound_packstat = tot_packstat_mval.value();
		Total_packstats = bound_packstat - unbound_packstat;

		//unsat polar calculation
		bound_pose.metric(burunsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
		bound_unsat = tot_unsat_mval.value();
		unbound_pose.metric(burunsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
		unbound_unsat = tot_unsat_mval.value();
		Interface_unsat = bound_unsat - unbound_unsat;

		darc_score_file<<"DARC_scores:"<<"	"<<tag<<"    "<<best_DARC_score<<"    "<<bound_energy<<"    "<<Interface_Energy<<"    "<<Interface_HB<<"    "<<Total_packstats<<"    "<<Interface_unsat<<"\n";
		std::cout<<"SCORES_DESC:"<<"	"<<"TAG"<<"    "<<"DARC_score"<<"    "<<"Total_Energy"<<"    "<<"Interface_Energy"<<"    "<<"Interface_HB"<<"    "<<"Total_packstats"<<"    "<<"Interface_unsat"<<std::endl;
		std::cout<<"DARC_scores:"<<"	"<<tag<<"    "<<best_DARC_score<<"    "<<bound_energy<<"    "<<Interface_Energy<<"    "<<Interface_HB<<"    "<<Total_packstats<<"    "<<Interface_unsat<<std::endl;
		}
	}
	darc_score_file.close();
        } catch ( utility::excn::EXCN_Base const & e ) {
                std::cout << "caught exception " << e.msg() << std::endl;
        }
	return 0;

}
