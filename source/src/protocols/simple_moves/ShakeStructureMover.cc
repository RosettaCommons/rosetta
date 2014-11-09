// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Liz Kellogg

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueMatcher.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>

// AUTO-REMOVED #include <core/init/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <core/pack/task/ResfileReader.hh>

// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <sstream>
// AUTO-REMOVED #include <ios>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <ObjexxFCL/format.hh>

// C++ headers
#include <cstdlib>
#include <string>

#include <basic/Tracer.hh>


#include <core/chemical/ResidueTypeSet.fwd.hh>

// AUTO-REMOVED #include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/Constraints.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>


// AUTO-REMOVED #include <basic/options/option.hh>

// AUTO-REMOVED #include <basic/basic.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/ProteinSilentStruct.hh>

// AUTO-REMOVED #include <core/pack/rotamer_trials.hh>
// AUTO-REMOVED #include <core/io/silent/silent.fwd.hh>

//protocols
#include <protocols/simple_moves/ShakeStructureMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
// AUTO-REMOVED #include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/Mover.hh>

// AUTO-REMOVED #include <utility/file/FileName.hh>
#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <time.h>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers

#include <utility/vector0.hh>



namespace protocols {
namespace simple_moves {


ShakeStructureMover::ShakeStructureMover() :
	protocols::moves::Mover("ShakeStructureMover"),
	mc_temp(0),
	ramp_fa_rep(false),
	min_cst(false),
	scorefxn(/* 0 */),
	ensemble_ca_rmsd(-1),
	ensemble_ca_rmsd_tolerance(0.75f),
	is_properly_initialized(false),
	harmonic_ca_cst_std_dev(2.0f),
	scorefunction_initialized(false),
	sc_min_only(false),
	nrounds(1000),
	cst_weight(1.0f),
	skip_low_temp_phase(true),
	min_scorefxn(/* 0 */),
	testing_phase(false)
{}

ShakeStructureMover::ShakeStructureMover(core::scoring::ScoreFunctionOP s) :
	protocols::moves::Mover("ShakeStructureMover"),
	mc_temp(0),
	ramp_fa_rep(false),
	min_cst(false),
	scorefxn(s),
	ensemble_ca_rmsd(-1),
	ensemble_ca_rmsd_tolerance(0.75),
	is_properly_initialized(false),
	harmonic_ca_cst_std_dev(2.0),
	scorefunction_initialized(true),
	sc_min_only(false),
	nrounds(1000),
	cst_weight(1.0),
	skip_low_temp_phase(true),
	min_scorefxn(/* 0 */),
	testing_phase(false)
{}

ShakeStructureMover::ShakeStructureMover(
	core::scoring::ScoreFunctionOP s,
	core::Real temperature
) :
	protocols::moves::Mover("ShakeStructureMover"),
	mc_temp(temperature),
	ramp_fa_rep(false),
	min_cst(false),
	scorefxn(s),
	ensemble_ca_rmsd(-1),
	ensemble_ca_rmsd_tolerance(0.75),
	is_properly_initialized(false),
	harmonic_ca_cst_std_dev(2.0),
	scorefunction_initialized(true),
	sc_min_only(false),
	nrounds(1000),
	cst_weight(1.0),
	skip_low_temp_phase(true),
	min_scorefxn(/* 0 */),
	testing_phase(false)
{}

ShakeStructureMover::ShakeStructureMover(
	core::scoring::ScoreFunctionOP s,
	core::Real ens_diversity, core::Real ens_div_tolerance
):
	protocols::moves::Mover("ShakeStructureMover"),
	mc_temp(0),
	ramp_fa_rep(false),
	min_cst(false),
	scorefxn(s),
	ensemble_ca_rmsd(ens_diversity),
	ensemble_ca_rmsd_tolerance(ens_div_tolerance),
	is_properly_initialized(false),
	harmonic_ca_cst_std_dev(2.0),
	scorefunction_initialized(true),
	sc_min_only(false),
	nrounds(1000),
	cst_weight(1.0),
	skip_low_temp_phase(true),
	min_scorefxn(/* 0 */),
	testing_phase(false)
{}

ShakeStructureMover::~ShakeStructureMover(){}

void
ShakeStructureMover::set_skip_low_temp_phase(bool truefalse){
	skip_low_temp_phase = truefalse;
}

void
ShakeStructureMover::set_mc_temperature(core::Real temp){
	mc_temp=temp;
}

void
ShakeStructureMover::set_nrounds( int new_nrounds ){
	nrounds=new_nrounds;
}


void
ShakeStructureMover::set_ramp_fa_rep(bool truefalse){
	ramp_fa_rep=truefalse;
}

void
ShakeStructureMover::set_minimize_with_cst(bool truefalse){
	min_cst=truefalse;
}

void
ShakeStructureMover::set_scorefunction(core::scoring::ScoreFunctionOP s){
	scorefxn = s;
	scorefunction_initialized=true;
}

void
ShakeStructureMover::set_ensemble_diversity(core::Real ca_rmsd){
	ensemble_ca_rmsd=ca_rmsd;
	ensemble_ca_rmsd_tolerance=((ca_rmsd*0.25)+0.1);
	mc_temp=-1;
}


void
ShakeStructureMover::set_rmsd_target_tolerance(core::Real tol){
	ensemble_ca_rmsd_tolerance=tol;
}

void
ShakeStructureMover::set_sc_min(bool truefalse){
	sc_min_only=truefalse;
	min_cst=true;
}

void
ShakeStructureMover::set_testing_phase( bool truefalse ){
	testing_phase = truefalse;
}

void
ShakeStructureMover::set_mc_temp( core::Real temperature ){
	mc_temp = temperature;
}

void
ShakeStructureMover::set_is_properly_initialized( bool truefalse ){
	is_properly_initialized = truefalse;
}

void
ShakeStructureMover::set_min_scorefunction(core::scoring::ScoreFunctionOP scfxn ){
	min_scorefxn = scfxn;
}

//accessors
core::Real
ShakeStructureMover::get_mc_temperature(){
	return mc_temp;
}

bool
ShakeStructureMover::get_ramp_fa_rep(){
	return ramp_fa_rep;
}

bool
ShakeStructureMover::get_minimize_with_cst(){
	return min_cst;

}

core::scoring::ScoreFunctionOP
ShakeStructureMover::get_scorefunction(){
	return scorefxn;
}

core::Real
ShakeStructureMover::get_ensemble_diversity(){
	return ensemble_ca_rmsd;
}

core::Real
ShakeStructureMover::get_rmsd_target_tolerance(){
	return ensemble_ca_rmsd_tolerance;
}

bool
ShakeStructureMover::get_sc_min(){
	return sc_min_only;
}

core::Real
ShakeStructureMover::get_harmonic_ca_cst_std_dev(){
	return harmonic_ca_cst_std_dev;
}

core::Real
ShakeStructureMover::get_ensemble_ca_rmsd(){
	return ensemble_ca_rmsd;
}

bool
ShakeStructureMover::get_skip_low_temp_phase(){
	return skip_low_temp_phase;
}

bool
ShakeStructureMover::get_min_cst(){
	return min_cst;
}

bool
ShakeStructureMover::get_testing_phase(){
	return testing_phase;
}

bool
ShakeStructureMover::get_scorefunction_initialized(){
	return scorefunction_initialized;
}

core::scoring::ScoreFunctionOP
ShakeStructureMover::get_min_scorefunction(){
	return min_scorefxn;
}

core::Size
ShakeStructureMover::get_nrounds(){
	return nrounds;
}

//run-time


void
ShakeStructureMover::apply(core::pose::Pose & pose){

	using namespace core;
	using namespace core::pose;
	using namespace utility;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	if(!is_properly_initialized){
		setup_for_run(pose);
	}

	create_ensemble(pose, *scorefxn);
}


std::string
ShakeStructureMover::get_name() const {
	return "ShakeStructureMover";
}

void
ShakeStructureMover::reduce_fa_rep(float fraction_fa_rep, core::scoring::ScoreFunction & s){
	s.set_weight(
		core::scoring::score_type_from_name("fa_rep"),
		s.get_weight(core::scoring::score_type_from_name("fa_rep"))*fraction_fa_rep);
}

void
ShakeStructureMover::setup_for_run(core::pose::Pose & p){
	if(!scorefunction_initialized){
		core::scoring::ScoreFunctionOP s( new core::scoring::ScoreFunction() );
		//std::cout << "[DEBUG]: scorefunction being initialized"<< std::endl;
		s->set_weight(core::scoring::score_type_from_name("rama"), 4.0);
		s->set_weight(core::scoring::score_type_from_name("omega"), 1.0);
		s->set_weight(core::scoring::score_type_from_name("fa_dun"), 1.0);
		s->set_weight(core::scoring::score_type_from_name("p_aa_pp"), 2.0);
		scorefxn = s;
	}
	setup_ca_constraints(p,(*scorefxn),9.0,harmonic_ca_cst_std_dev);

	if(min_cst || sc_min_only){
		min_scorefxn = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	}

	if(mc_temp <= 0 && ensemble_ca_rmsd > 0){
		//set mc_temp based on ensemble_ca_rmsd
		testing_phase=true;
		mc_temp = set_temp_based_on_ens_diversity(p,(*scorefxn));
		is_properly_initialized = true;
	}else if(mc_temp > 0 && ensemble_ca_rmsd < 0){
		//ready to go!
		is_properly_initialized = true;
	}else{
		//what to do, what to do?
		is_properly_initialized = false;
	}


}

void
ShakeStructureMover::minimize_with_constraints(core::pose::Pose & p,
															  core::scoring::ScoreFunction & s){

	core::optimization::AtomTreeMinimizer min_struc;
	float const minimizer_tol = 0.0000001;
	core::optimization::MinimizerOptions options(
		"dfpmin_armijo_nonmonotone",
		minimizer_tol,
		true /*use_nb_list*/,
		false /*deriv_check_in*/,
		false /*deriv_check_verbose_in*/);
	options.nblist_auto_update( true );
	//      options.max_iter(5000);
	core::kinematics::MoveMap mm;
	if(!sc_min_only){
		mm.set_bb(true);
	}else{
		mm.set_bb(false);
	}
	mm.set_chi(true);

	if(ramp_fa_rep){
		core::scoring::ScoreFunctionOP one_tenth_orig(s.clone());
		reduce_fa_rep(0.1,*one_tenth_orig);
		min_struc.run(p,mm,*one_tenth_orig,options);
		core::scoring::ScoreFunctionOP one_third_orig(s.clone());
		reduce_fa_rep(0.33,*one_third_orig);
		min_struc.run(p,mm,*one_third_orig,options);
	}
	min_struc.run(p,mm,s,options);

}

void
ShakeStructureMover::setup_ca_constraints(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction & s,
	float const CA_cutoff,
	float const cst_tol
){

	using namespace core::scoring::constraints;
	using namespace core::id;
	using namespace core;
	//create constraints for all residues
	//type: HARMONIC
	//static float const CA_cutoff(9.0);
	int nres = pose.total_residue();
	for(int i = 1; i <= nres; i++){
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;
		Vector const CA_i( pose.residue(i).xyz(" CA "));
		for(int j = 1; j <=nres; j++){
			if ( pose.residue(j).aa() == core::chemical::aa_vrt ) continue;
			Vector const CA_j(pose.residue(j).xyz(" CA "));
			Real const CA_dist = (CA_i - CA_j).length();
			if(CA_dist < CA_cutoff){
				ConstraintCOP cst( ConstraintOP( new AtomPairConstraint( AtomID(pose.residue(i).atom_index(" CA "),i),AtomID(pose.residue(j).atom_index(" CA "),j), core::scoring::func::FuncOP(new core::scoring::func::HarmonicFunc(CA_dist, cst_tol)) ) ) );
				pose.add_constraint(cst);
			}
		}
	}

	s.set_weight(core::scoring::atom_pair_constraint, cst_weight);

}

void
ShakeStructureMover::setup_movers(
	simple_moves::SmallMoverOP gsmall, simple_moves::ShearMoverOP gshear,
	core::Real small_H_angle_max, core::Real small_E_angle_max, core::Real small_L_angle_max,
	core::Real shear_H_angle_max, core::Real shear_E_angle_max, core::Real shear_L_angle_max
){
	//smaller moves to boost acceptance rate
	gsmall->angle_max( 'H', small_H_angle_max ); //def 0.5
	gsmall->angle_max( 'E', small_E_angle_max ); //def 0.5
	gsmall->angle_max( 'L', small_L_angle_max ); //def 0.75

	gshear->angle_max( 'H', shear_H_angle_max); //def 2.0
	gshear->angle_max( 'E', shear_E_angle_max); //def 2.0
	gshear->angle_max( 'L', shear_L_angle_max); //def 3.0
}



void
ShakeStructureMover::run_mc(
	core::pose::Pose & p, core::scoring::ScoreFunction & s,
	core::Real temperature
){

	using namespace protocols;
	using namespace moves;

	core::pose::Pose init(p); //for comparison at the end

	core::Size nmoves = (core::Size)p.total_residue()/4; //number of moves for each move type

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
	mm->set_bb(true);

	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( mm, temperature, nmoves) ) ;
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover(mm, temperature, nmoves) );

	setup_movers(
		small_mover,shear_mover,
		0.2,0.2,0.4,
		1.6,1.6,2.0);

	moves::RandomMoverOP apply_random_move( new moves::RandomMover() );
	apply_random_move->add_mover( small_mover, .5);
	apply_random_move->add_mover( shear_mover, .5);

	simple_moves::SmallMoverOP small_mover_low( new simple_moves::SmallMover( mm, (temperature*0.25), nmoves) ) ;
	simple_moves::ShearMoverOP shear_mover_low( new simple_moves::ShearMover(mm, (temperature*0.25), nmoves) );

	setup_movers(
		small_mover_low,shear_mover_low,
		0.1,0.1,0.2,
		1.0,1.0,1.5);

	moves::RandomMoverOP apply_random_move_low( new moves::RandomMover() );
	apply_random_move_low->add_mover( small_mover_low, .5);
	apply_random_move_low->add_mover( shear_mover_low, .5);

	protocols::moves::MonteCarloOP mc( new moves::MonteCarlo(p,s,temperature) );
	//time_t time_per_decoy = time(NULL);

	//			mc->reset_counters();
	//			mc->reset(p);
	mc->set_temperature(temperature);

	moves::TrialMoverOP tm( new moves::TrialMover(apply_random_move,mc) );
	RepeatMoverOP full_cycle( new moves::RepeatMover( tm, nrounds ) );
	full_cycle->apply( p );
	//			mc->show_counters();

	if(!skip_low_temp_phase){
		mc->reset_counters();

		mc->set_lowest_score_pose(mc->last_accepted_pose()); //high energy mark for annealing back to native
		core::Real low_temp = temperature * 0.25;
		mc->set_temperature(low_temp);

		moves::TrialMoverOP ltm( new moves::TrialMover(apply_random_move_low,mc) );
		RepeatMoverOP full_cycle_2( new moves::RepeatMover( ltm, nrounds ) );
		full_cycle_2->apply( p );
		//mc->show_counters();
		p = mc->lowest_score_pose();
	} else {
		p = mc->last_accepted_pose();
	}

	if(min_cst && !testing_phase){
		//remove constraints
		p.remove_constraints((p.constraint_set())->get_all_constraints());

		//std::cout << "minimizing with constraints" << std::endl;
		//temporarily hardcode score12 into minimization.
		setup_ca_constraints(p, (*min_scorefxn), 9.0, 0.5);
		minimize_with_constraints(p, (*min_scorefxn));

		//std::cout << " CA rmsd of current: " << core::scoring::CA_rmsd(p,init) << std::endl;
		//time_t time_per_decoy_finish = time(NULL);
		//std::cout << "time to finish decoy " << (time_per_decoy_finish-time_per_decoy) << std::endl;
	}else{
		//time_t time_per_decoy_finish = time(NULL);
		//std::cout << " CA rmsd of current: " << core::scoring::CA_rmsd(p,init) << std::endl;
		//std::cout << "time to finish decoy " << (time_per_decoy_finish-time_per_decoy) << std::endl;
	}
}

core::Real
ShakeStructureMover::set_temp_based_on_ens_diversity(
	core::pose::Pose & p,
	core::scoring::ScoreFunction & s
){

	int numstruct = 10;
	core::Real avg_ca_rmsd = -1;
	core::Real temperature = 10;
	core::Real increment = 0.75;

	core::pose::Pose init(p);
	core::pose::Pose test(p);
	while((avg_ca_rmsd) > (ensemble_ca_rmsd+ensemble_ca_rmsd_tolerance) || (avg_ca_rmsd) < (ensemble_ca_rmsd-ensemble_ca_rmsd_tolerance)){
		double avg_ca_rmsd=0;
		for(int i =0;i<numstruct;i++){
			test = init;
			run_mc(test, s,temperature);
			//std::cout << "[DEBUG]: rmsd in set temp " << core::scoring::CA_rmsd(test,init) << std::endl;
			avg_ca_rmsd+=core::scoring::CA_rmsd(test,init);
		}
		avg_ca_rmsd = avg_ca_rmsd /numstruct;
		if(avg_ca_rmsd > (ensemble_ca_rmsd+ensemble_ca_rmsd_tolerance)){
			temperature -= increment;
			//std::cout <<  "avg rmsd is " << avg_ca_rmsd << " which is greater than upper bound: " << ensemble_ca_rmsd+ensemble_ca_rmsd_tolerance << " so temp is decreased to " << temperature << std::endl;
		}
		else if(avg_ca_rmsd < (ensemble_ca_rmsd-ensemble_ca_rmsd_tolerance)){
			temperature += increment;
			//std::cout <<  "avg rmsd is " << avg_ca_rmsd << " which is less than lower bound: " << ensemble_ca_rmsd-ensemble_ca_rmsd_tolerance << " so temp is increased to " << temperature << std::endl;
		}
	}
	testing_phase=false;
	return temperature;
}

void
ShakeStructureMover::create_ensemble(core::pose::Pose & p, core::scoring::ScoreFunction & s){
	using namespace protocols;
	using namespace moves;

	run_mc(p,s,mc_temp);
}

} //simple_moves
} //protocols
