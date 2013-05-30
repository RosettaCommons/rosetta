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
#include <protocols/pockets/PocketGrid.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/FingerprintMultifunc.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <core/import_pose/import_pose.hh>

//reqd minimization headers
#include <protocols/simple_moves/ScoreMover.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>

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
OPT_KEY( String, ligand )
OPT_KEY( String, eggshell_triplet )
OPT_KEY( Integer, num_runs )
OPT_KEY( Integer, num_particles )
OPT_KEY( Real, steric_weight )
OPT_KEY( Real, missing_point_weight )
OPT_KEY( Real, extra_point_weight )
OPT_KEY( Real, origin_cutoff )
OPT_KEY( Boolean, print_output_complex )

//reqd minimization options
OPT_KEY( Boolean, minimize_output_complex )
OPT_KEY( Real, cst_force_constant )

static basic::Tracer TR( "apps.pilot.ragul_run_darc_with_input_eggshell.main" );

int main( int argc, char * argv [] ) {

  NEW_OPT( protein, "protein file name", "protein.pdb" );
  NEW_OPT( ligand, "ligand file name", "ligand.pdb" );
  NEW_OPT( eggshell_triplet, "input eggshell triplet file name", "" );
  NEW_OPT( num_runs, "no. of runs for PSO", 100 );
  NEW_OPT( num_particles, "no. of particles for PSO", 100 );
  NEW_OPT( origin_cutoff, "value for setting minimum and maximum origin cut off", 7.0 );
  NEW_OPT( missing_point_weight, "missing point weight", 21.6 );
  NEW_OPT( extra_point_weight, "extra point weight", 9.5 );
  NEW_OPT( steric_weight, "steric weight for PSO", 1.4 );
	NEW_OPT( print_output_complex, "print DARC output ligand model with protein as PDB file", true );
	NEW_OPT( minimize_output_complex, "minimize the best scoring DARC output model", false );
	NEW_OPT( cst_force_constant, "coordinate constraint force constant", 0.5 );

  devel::init(argc, argv);

  std::string const input_protein = option[ protein ];
  std::string const input_ligand = option[ ligand ];
  std::string const input_eggshell_triplet = option[ eggshell_triplet ];
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
		npf.setup_from_eggshell_triplet_file( input_eggshell_triplet );
	}
	else if(input_eggshell_triplet.empty()){
		std::cout<<"ERROR! : no input eggshell file to setup pocket"<<std::endl;
			exit(1);
	}

	//ligand pose can be a single pdb file or multi-conformer pdb file
  pose::Pose small_mol_pose;
  core::import_pose::pose_from_pdb( small_mol_pose, input_ligand );

	//create 'tag' for output filenames
	int dot_index1 = input_protein.rfind(".", input_protein.size());
	assert(dot_index1 != -1 && "No dot found in filename");
	std::string protein_name = input_protein.substr(0,dot_index1);
	int dot_index2 = input_ligand.rfind(".", input_ligand.size());
	assert(dot_index2 != -1 && "No dot found in filename");
	std::string ligand_name = input_ligand.substr(0,dot_index2);
	std::string tag = ligand_name + "_" + protein_name;
	std::string pso_pose_name = "LIGAND_" + tag + ".pdb";
	std::string darc_complex_filename = "DARC_" + tag + ".pdb";
	std::string mini_complex_filename = "mini_" + tag + ".pdb";

	//set max and min ligand translation & rotation space
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

	//initialize particles, best DARC score, best optimized values
	ParticleOPs particles;
	core::Real best_DARC_score = 9999.;
	utility::vector1<core::Real> best_vars(6);

	//Find the residue number for the ligand
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

	//if ligand is a multi-conformer file, run DARC for each conformer and find the best score
	protocols::pockets::PlaidFingerprint pf(small_mol_pose, npf);

  //setup GPU
#ifdef USEOPENCL
	npf.gpu_setup( missing_pt_wt, steric_wt, extra_pt_wt, particle_size, pf );
#endif

	utility::vector1< core::pose::PoseOP > rot_poses(small_mol_pose.total_residue());
	for (Size ii = 1; ii <= small_mol_pose.total_residue(); ++ii){
		rot_poses[ii] = new core::pose::Pose(small_mol_pose, ii, ii);
		protocols::pockets::PlaidFingerprint conf_pf( *rot_poses[ii], npf );
		std::cout<< "JK ERROR - this code is not yet conformer-enabled, fix it in the app by removing the 1 in FingerprintMultifunc constructor below..." << std::endl;
		exit(1);
		protocols::pockets::FingerprintMultifunc fpm(npf, conf_pf, missing_pt_wt, steric_wt, extra_pt_wt, 1);
		protocols::pockets::DarcParticleSwarmMinimizer pso( npf, pf, missing_pt_wt, steric_wt, extra_pt_wt, p_min, p_max);
		//core::optimization::ParticleSwarmMinimizer pso(p_min, p_max);
		particles = pso.run(run_size, fpm, particle_size);
		ParticleOP p = particles[1];
		core::optimization::Particle parti(*p);
		core::Real fit_best = -(parti.fitness_pbest());
		if(fit_best < best_DARC_score){
			best_DARC_score = fit_best;
			best_vars = parti.pbest();
			pf = conf_pf;
		}
	}

	//Append the best ligand conformer pose to protein pose and print(optional) the complex PDB file
	numeric::xyzVector<core::Real> optimized_origin(best_vars[1], best_vars[2], best_vars[3]);
	std::cout<< "JK ERROR - this code is not yet conformer-enabled, fix it in the app by removing the zero in the call to get_oriented_pose below..." << std::endl;
	exit(1);
	core::pose::Pose oriented_pose = pf.get_oriented_pose(npf, best_vars[4], best_vars[5], best_vars[6], optimized_origin, 0 );
	bound_pose.append_residue_by_jump(oriented_pose.residue( 1 ), protein_pose.total_residue(),"", "",  true);
	if (option[ print_output_complex ]()){
		bound_pose.dump_pdb(darc_complex_filename);
	}

	//print the best DARC score (without minimization)
	if (!option[ minimize_output_complex ]()){
		std::cout<< "SCORES : " << tag <<"\n"<< "DARC Score : " << best_DARC_score <<std::endl;
	}

	//minimize the protein-ligand complex
	else if (option[ minimize_output_complex ]()){

		//setup scorefxn
		scoring::ScoreFunctionOP scorefxn(getScoreFunction());
		scoring::ScoreFunctionOP repack_scorefxn(getScoreFunction());

		//Register calculators
		std::string sasa_calc_name = "sasa";
		std::string hbond_calc_name = "hbond";
		std::string packstat_calc_name = "packstat";
		std::string burunsat_calc_name = "burunsat";
		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculator;
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( sasa_calc_name, sasa_calculator );
		core::pose::metrics::PoseMetricCalculatorOP hb_calc = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( hbond_calc_name, hb_calc );
		core::pose::metrics::PoseMetricCalculatorOP packstat_calc =	new protocols::toolbox::pose_metric_calculators::PackstatCalculator();
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( packstat_calc_name, packstat_calc );
		core::pose::metrics::PoseMetricCalculatorOP burunsat_calc = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, hbond_calc_name);
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( burunsat_calc_name, burunsat_calc );

		// add constraint
		Size nres = bound_pose.total_residue();
		Real const coord_sdev( option[ OptionKeys::cst_force_constant ] );
		// default is 0.5 (from idealize) -- maybe too small
		for ( Size i = 1; i< nres; ++i ) {
			if ( (Size)i==(Size)bound_pose.fold_tree().root() ) continue;
			Residue const & nat_i_rsd( bound_pose.residue(i) );
			for ( Size ii = 1; ii<= nat_i_rsd.nheavyatoms(); ++ii ) {
				bound_pose.add_constraint( new CoordinateConstraint(
																														AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
																														new HarmonicFunc( 0.0, coord_sdev ) ) );
			}
		}

		// score Initial bound pose
		scorefxn->set_weight( coordinate_constraint, 0.5 );
		(*scorefxn)(bound_pose);
		TR << "Initial score: " << bound_pose.energies().total_energies()[ total_score ] << std::endl;

		// setting degrees of freedom which can move during minimization - everything
		kinematics::MoveMap mm_all;
		mm_all.set_chi( true );
		mm_all.set_bb( true );
		mm_all.set_jump( true );

		// start minimizing protein
		TR << "Starting minimization...." << std::endl;
		AtomTreeMinimizer minimizer;
		MinimizerOptions min_options( "dfpmin", 0.00001, true, false );
		minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
		minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
		(*scorefxn)(bound_pose);
		TR << "Post minimization 1 constrained score: " << bound_pose.energies().total_energies()[ total_score ] << std::endl;
		bound_pose.remove_constraints();
		(*scorefxn)(bound_pose);
		TR << "Post minimization 1 UNconstrained score: " << bound_pose.energies().total_energies()[ total_score ] << std::endl;

		// Setup packer task for repacking
		pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( bound_pose ));
		base_packer_task->set_bump_check( false );
		base_packer_task->initialize_from_command_line();
		base_packer_task->or_include_current( true );
		for ( Size ii = 1; ii <= bound_pose.total_residue(); ++ii ) {
			base_packer_task->nonconst_residue_task(ii).restrict_to_repacking();
		}

		// First repack and report score
		pack::pack_rotamers( bound_pose, *repack_scorefxn, base_packer_task );
		(*scorefxn)(bound_pose);
		TR << "Score after repacking once: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;

		// iterate over minimizing and repacking
		for ( Size iter = 1; iter <= 5; ++iter ) {
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
			(*scorefxn)(bound_pose);
			//    TR << "Current score after minimizing: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
			pack::pack_rotamers( bound_pose, *repack_scorefxn, base_packer_task );
			(*scorefxn)(bound_pose);
			//    TR << "Current score after repacking: " << bound_pose.energies().total_energies()[ total_score ] << std::endl << std::endl;
		}

		// final minimization
		minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
		minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
		minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
		minimizer.run( bound_pose, mm_all, *scorefxn, min_options );
		(*scorefxn)(bound_pose);
		if (option[ print_output_complex ]){
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
		core::Real bound_sasa = 0.0, unbound_sasa = 0.0, Total_BSA = 0.0;
		core::Size  bound_hb = 0,   unbound_hb = 0, Interface_HB = 0;
		core::Real bound_packstat = 0.0, unbound_packstat = 0.0, Total_packstats = 0.0;
		core::Size  bound_unsat = 0, unbound_unsat = 0, Interface_unsat = 0;

		//calculate interface Energy
		bound_energy = bound_pose.energies().total_energy();
		unbound_energy = unbound_pose.energies().total_energy();
		Interface_Energy = bound_energy - unbound_energy;

		//delta sasa calculation
		bound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
		bound_sasa = tot_sasa_mval.value();
		unbound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
		unbound_sasa = tot_sasa_mval.value();
		Total_BSA = unbound_sasa - bound_sasa;

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

		std::cout<< "SCORES : " << tag <<"\n"<< "DARC Score : " << best_DARC_score <<"\n"<< "Total Energy : " << bound_energy <<"\n"<< " Interface Energy : "<< Interface_Energy << std::endl;

	}

	return 0;

}
