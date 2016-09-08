// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <core/pack/task/ResfileReader.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <ios>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <basic/Tracer.hh>


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.hh>

#include <devel/init.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

//protocols
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
//#include <protocols/looprelax/looprelax_main.hh>
#include <protocols/comparative_modeling/ConstraintRemodelMover.hh>

#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
//#include "james_util.hh" //for calculation of burial
#include <basic/Tracer.hh>
#include <time.h>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//C++ filechek
#include <sys/stat.h>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace std;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace ddg;
using namespace id;
using namespace protocols;
using namespace moves;


ScoreFunction&
reduce_fa_rep(float fraction_fa_rep, ScoreFunction & s){
	s.set_weight( core::scoring::score_type_from_name("fa_rep"),
								s.get_weight(core::scoring::score_type_from_name("fa_rep"))*fraction_fa_rep);
	return s;
}

void
minimize_with_constraints(pose::Pose & p, ScoreFunction & s,std::string output_tag){
	core::optimization::AtomTreeMinimizer min_struc;
	float const minimizer_tol = 0.0000001;
	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", minimizer_tol, true /*use_nb_list*/,
																								false /*deriv_check_in*/, false /*deriv_check_verbose_in*/);
	options.nblist_auto_update( true );
	options.max_iter(5000); //otherwise, they don't seem to converge
	core::kinematics::MoveMap mm;
	if(!basic::options::option[sc_min_only]()){
		mm.set_bb(true);
	}else{
		mm.set_bb(false);
	}
	mm.set_chi(true);


//	s.show(std::cout, p);
	std::string out_pdb_prefix = basic::options::option[OptionKeys::ddg::out_pdb_prefix ]();

	if(basic::options::option[OptionKeys::ddg::ramp_repulsive]()){
		//set scorefxn fa_rep to 1/10 of original weight and then minimize
		ScoreFunction one_tenth_orig(s);
		reduce_fa_rep(0.1,one_tenth_orig);
		//min_struc.run(p,mm,s,options);
		min_struc.run(p,mm,one_tenth_orig,options);
//		std::cout << "one tenth repulsive fa_rep score-function" << std::endl;
//		one_tenth_orig.show(std::cout, p);

		//then set scorefxn fa_rep to 1/3 of original weight and then minimize
		ScoreFunction one_third_orig(s);
		reduce_fa_rep(0.33,one_third_orig);
		min_struc.run(p,mm,one_third_orig,options);
//		std::cout << "one third repulsive fa_rep score-function" << std::endl;
//		one_third_orig.show(std::cout, p);
		//then set scorefxn fa_rep to original weight and then minimize
	}
	min_struc.run(p,mm,s,options);
//	s.show(std::cout, p);

	p.dump_pdb(output_tag+"_0001.pdb");
}


void
setup_ca_constraints(pose::Pose & pose, ScoreFunction & s, float const CA_cutoff, float const cst_tol){
	//create constraints for all residues
	//type: HARMONIC
	//static float const CA_cutoff(9.0);
	int nres = pose.size();
	for(int i = 1; i <= nres; i++){
		Vector const CA_i( pose.residue(i).xyz(" CA "));
		for(int j = 1; j <=nres; j++){
			Vector const CA_j(pose.residue(j).xyz(" CA "));
			Real const CA_dist = (CA_i - CA_j).length();
			if(CA_dist < CA_cutoff){
				//verbose_output << "c-alpha constraints added to residues " << i << " and " << j << std::endl;
				ConstraintCOP cst(new AtomPairConstraint( AtomID(pose.residue(i).atom_index(" CA "),i),AtomID(pose.residue(j).atom_index(" CA "),j),new HarmonicFunc(CA_dist, cst_tol)));
				pose.add_constraint(cst);
			}
		}
	}

	s.set_weight(atom_pair_constraint, basic::options::option[OptionKeys::ddg::constraint_weight]());

}


void
setup_movers(simple_moves::SmallMoverOP small, simple_moves::ShearMoverOP shear,
						 core::Real small_H_angle_max, core::Real small_E_angle_max, core::Real small_L_angle_max,
						 core::Real shear_H_angle_max, core::Real shear_E_angle_max, core::Real shear_L_angle_max){
	//smaller moves to boost acceptance rate
	small->angle_max( 'H', small_H_angle_max ); //def 0.5
	small->angle_max( 'E', small_E_angle_max ); //def 0.5
 	small->angle_max( 'L', small_L_angle_max ); //def 0.75

	shear->angle_max( 'H', shear_H_angle_max); //def 2.0
	shear->angle_max( 'E', shear_E_angle_max); //def 2.0
	shear->angle_max( 'L', shear_L_angle_max); //def 3.0
}


double
run_mc(pose::Pose & p, ScoreFunction & s,
			 core::Real temperature, int numstruct,
			 std::string output_tag, bool output_pdbs){

	using namespace protocols;
	using namespace moves;
	using namespace core::scoring::constraints;

	core::Size nmoves = (core::Size)p.size()/4; //number of moves for each move type
	//	std::cout << "nmoves " << nmoves << " will be performed for each move type at each monte carlo step" << std::endl;
	core::Size nrounds = 1000;

	core::kinematics::MoveMapOP mm= new core::kinematics::MoveMap();
	mm->set_bb(true);

	//vector to hold CA-rmsd from start to assess whether the level of variability is correct
	utility::vector1<core::Real> CA_rmsd_from_start;

	//small mover
	simple_moves::SmallMoverOP small_mover(new simple_moves::SmallMover( mm, temperature, nmoves)) ;
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover(mm, temperature, nmoves));

	setup_movers(small_mover,shear_mover,
							 0.2,0.2,0.4,
							 1.6,1.6,2.0);

	moves::RandomMoverOP apply_random_move( new moves::RandomMover());
	apply_random_move->add_mover( small_mover, .5);
	apply_random_move->add_mover( shear_mover, .5);

	//separate set of movers for low temperature (reduce angle-max)
	simple_moves::SmallMoverOP small_mover_low(new simple_moves::SmallMover( mm, (temperature*0.25), nmoves)) ;
	simple_moves::ShearMoverOP shear_mover_low( new simple_moves::ShearMover(mm, (temperature*0.25), nmoves));

	setup_movers(small_mover_low,shear_mover_low,
							 0.1,0.1,0.2,
							 1.0,1.0,1.5);

	moves::RandomMoverOP apply_random_move_low( new moves::RandomMover());
	apply_random_move_low->add_mover( small_mover_low, .5);
	apply_random_move_low->add_mover( shear_mover_low, .5);

	pose::Pose init_pose(p);
	s(init_pose);
	//s.show(std::cout,init_pose);


	MonteCarloOP mc(new moves::MonteCarlo(init_pose,s,temperature));
	for(int ns =1;ns<=numstruct;ns++){

		std::ostringstream curr;
		curr << ns;

		struct stat stFileInfo;
		int file_stat = stat((basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+"lowest."+curr.str()+".pdb").c_str(),
												 &stFileInfo);
		if(file_stat != 0 || !output_pdbs){ //file doesn't exist or file can exist and output_pdbs set to false

			time_t time_per_decoy = time(NULL);
			mc->reset_counters();
			mc->reset(init_pose);
			mc->set_temperature(temperature);
			p = init_pose;

			moves::TrialMoverOP tm( new moves::TrialMover(apply_random_move,mc));
			RepeatMoverOP full_cycle(new moves::RepeatMover( tm, nrounds ));
			full_cycle->apply( p );
			//		std::cout << "ended trial at high temp" << std::endl;
			//		mc->show_scores();
			//		mc->show_counters();
			//std::cout << "monte carlo scores at high temp" << std::endl;
			pose::Pose last_accepted = mc->last_accepted_pose();


			if(output_pdbs){
				core::io::pdb::dump_pdb(last_accepted, basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+"last_accepted_high_temp."+curr.str()+".pdb");
			}

			mc->reset_counters();

			mc->set_lowest_score_pose(last_accepted); //high energy mark for annealing back to native
			core::Real low_temp = temperature * 0.25;
			mc->set_temperature(low_temp);
			//std::cout << "double-checking" << std::endl;
			//mc->show_scores();
			//mc->show_counters();
			//std::cout << "end double-checking" << std::endl;

			//		std::cout << "original temperature was " << temperature << " and new low temperature is " << low_temp << std::endl;
			moves::TrialMoverOP ltm( new moves::TrialMover(apply_random_move_low,mc));
			RepeatMoverOP full_cycle_2(new moves::RepeatMover( ltm, nrounds ));
			full_cycle_2->apply( p );
			//		std::cout << "end low temp trial, showing monte carlo stats" << std::endl;
			//		mc->show_scores();
			//		mc->show_counters();
			//		std::cout << "end show monte-carlo stats" << std::endl;

			if(output_pdbs){
				core::io::pdb::dump_pdb(mc->last_accepted_pose(), basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+"last_accepted."+curr.str()+".pdb");
				core::io::pdb::dump_pdb(mc->lowest_score_pose(), basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+"lowest."+curr.str()+".pdb");
			}
			pose::Pose lowest_score_pose = mc->lowest_score_pose();
			//		std::cout << "CA rmsd from start being stored: " << core::scoring::CA_rmsd(lowest_score_pose,init_pose) << std::endl;
			if(basic::options::option[OptionKeys::ddg::min_with_cst]()){
				//std::cout << "min score" << std::endl;
				ScoreFunctionOP s = ScoreFunctionFactory::create_score_function( option[OptionKeys::ddg::min_cst_weights]);

				//new constraints based on latest structure
				setup_ca_constraints(lowest_score_pose, (*s), 9.0, 0.5);
				ConstraintSetCOP cs = lowest_score_pose.constraint_set();
				cs->show_definition(std::cout,lowest_score_pose);

				minimize_with_constraints(lowest_score_pose, (*s),basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+(output_tag + "." + curr.str()));
				//			std::cout << "CA rmsd from start being stored: " << core::scoring::CA_rmsd(lowest_score_pose,init_pose) << std::endl;
				CA_rmsd_from_start.push_back(core::scoring::CA_rmsd(lowest_score_pose,init_pose));
				//std::cout << "min score" << std::endl;
				time_t time_per_decoy_finish = time(NULL);
				//			std::cout << "time to finish decoy " << (time_per_decoy_finish-time_per_decoy) << std::endl;
			}else{
				//			std::cout << "CA rmsd from start being stored: " << core::scoring::CA_rmsd(lowest_score_pose,init_pose) << std::endl;
				CA_rmsd_from_start.push_back(core::scoring::CA_rmsd(lowest_score_pose,init_pose));
				time_t time_per_decoy_finish = time(NULL);
				//			std::cout << "time to finish decoy " << (time_per_decoy_finish-time_per_decoy) << std::endl;
			}
		} //if(file_stat != 0){
		else{
			std::cout << "file:  " << (basic::options::option[OptionKeys::ddg::last_accepted_pose_dir]()+"lowest."+curr.str()+".pdb") << " already exists, skipping to next iteration" << std::endl;
		}
	}
	double ca_sum=0;
	for(int i=1; i<=CA_rmsd_from_start.size();i++){
		ca_sum += CA_rmsd_from_start[i];
	}
	return ca_sum/CA_rmsd_from_start.size();
}


core::Real
set_temp_based_on_ens_diversity(pose::Pose & p,ScoreFunction & s, core::Real avg_ca_rmsd_target){
	//acceptance rate should be around 40-50% and temperature should yield ensemble diversity of avg_ca_rmsd_target
	//first alter acceptance rate to be around 40-50% with the scorefunction provided, and then alter the temperature to match the ca-rmsd-target.
	int numstruct = 5;
	core::Real avg_rms = -1;
	core::Real temperature = 10;

	while((avg_rms) > (avg_ca_rmsd_target+0.05) || (avg_rms) < (avg_ca_rmsd_target-0.05)){
		avg_rms = run_mc(p, s,temperature,numstruct,"NO_OUTPUT",false);
		if(avg_rms > (avg_ca_rmsd_target+0.05)){
			temperature -= 0.75;
	//		std::cout << "[TEST]: avg rmsd is " << avg_rms << " which is greater than target: " << avg_ca_rmsd_target << " so temperature is decreased to: " << temperature << std::endl;
		}
		else if(avg_rms < (avg_ca_rmsd_target-0.05)){
			temperature += 0.75;
//			std::cout << "[TEST]: avg rmsd is " << avg_rms << " which is less than target: " << avg_ca_rmsd_target << " so temperature is increased to: " << temperature << std::endl;
		}
	}
//	std::cout << "[TEST]: average ensemble variation: " << avg_rms << " target is: " << avg_ca_rmsd_target << " and temperature is: " << temperature << std::endl;
	return temperature;
}


void
create_ensemble(pose::Pose & p, ScoreFunction & s, std::string output_tag){

	using namespace protocols;
	using namespace moves;

	core::Real temperature=0;

	if(basic::options::option[OptionKeys::ddg::ens_variation].user()){
		//set temperature based on what kind of diversity you want  to see in your ensemble.
		core::Real avg_ca_rmsd = basic::options::option[OptionKeys::ddg::ens_variation]();
//		std::cout << "setting temperature based on ensemble variation" << std::endl;
		temperature = set_temp_based_on_ens_diversity(p,s,avg_ca_rmsd);
	}else if(basic::options::option[OptionKeys::ddg::temperature].user()){
		temperature = basic::options::option[OptionKeys::ddg::temperature]();
	}
	double avg_ca_rmsd = run_mc(p,s,temperature,basic::options::option[out::nstruct](),output_tag,true);
//	std::cout << "average production run CA-rmsd to start: " << avg_ca_rmsd << std::endl;
}


int
main( int argc, char* argv [] )
{
    try {
	using namespace core;
	using namespace core::pose;
	using namespace utility;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	// options, random initialization. MAKE SURE THIS COMES BEFORE OPTIONS
	devel::init( argc, argv );
	//tolerance for constraints; defaults to 2.0
	Real cst_tol = basic::options::option[ OptionKeys::ddg::harmonic_ca_tether ]();

	//	basic::options::option[run::dont_randomize_missing_coords](true);
	//gone?
	pose::Pose pose;

	ScoreFunctionOP s = get_score_function();
	ScoreFunctionOP scorefxn( new ScoreFunction());

	scorefxn->set_weight(core::scoring::score_type_from_name("rama"), 4.0);
	scorefxn->set_weight(core::scoring::score_type_from_name("omega"), 1.0);
	scorefxn->set_weight(core::scoring::score_type_from_name("fa_dun"), 1.0);
	scorefxn->set_weight(core::scoring::score_type_from_name("p_aa_pp"), 2.0);


	vector1<file::FileName> files;
	if(basic::options::option[in::file::s].user()){
		std::cout << "using -s option" << std::endl;
		files=option[in::file::s]();
	}else if( option[in::file::l].user()){
		std::cout << "using -l option " << std::endl;

		utility::vector1<file::FileName> list = basic::options::option[ in::file::l ]();
		for(unsigned int h=1;h<=list.size();h++){
			utility::io::izstream pdbs(list[h]);
			std::string fname;
			while(pdbs >> fname){
				files.push_back(fname);
			}
		}
	}
	for(unsigned int f=1; f<=files.size();f++){
//		std::cout << "examining file: " << files[f] << std::endl;
		core::import_pose::pose_from_file(pose, files[f], core::import_pose::PDB_file);
		setup_ca_constraints(pose,(*scorefxn),9.0,cst_tol);
		ConstraintSetCOP cs = pose.constraint_set();
		cs->show_definition(std::cout,pose);
		//	scorefxn->show(std::cout,pose);
		//create constraints for all residues
		//type: HARMONIC

		//then minimize
		std::string output = pose.pdb_info()->name();
		std::string no_pdb = output.erase(output.find(".pdb",0));
//		std::cout << "output pdb prefix is: " << (no_pdb.erase(0,(no_pdb.find_last_of("/")+1))) << std::endl;
		create_ensemble(pose, *scorefxn , (no_pdb.erase(0,(no_pdb.find_last_of("/")+1))));

	}
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
