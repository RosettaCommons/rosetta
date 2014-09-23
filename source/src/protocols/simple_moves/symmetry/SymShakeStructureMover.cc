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
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueMatcher.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>

// AUTO-REMOVED #include <core/init/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>
// AUTO-REMOVED #include <core/pack/task/ResfileReader.hh>

// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <sstream>
// AUTO-REMOVED #include <ios>
// AUTO-REMOVED #include <utility/io/izstream.hh>

// C++ headers
#include <cstdlib>
#include <string>

#include <basic/Tracer.hh>


// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.fwd.hh>

// AUTO-REMOVED #include <core/scoring/sasa.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/Constraints.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
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
#include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/io/silent/silent.fwd.hh>

//protocols
#include <protocols/simple_moves/symmetry/SymShakeStructureMover.hh>
#include <protocols/simple_moves/ShakeStructureMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
// AUTO-REMOVED #include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/Mover.hh>

// AUTO-REMOVED #include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <time.h>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>



namespace protocols {
namespace simple_moves {
namespace symmetry {

    SymShakeStructureMover::SymShakeStructureMover():
      ShakeStructureMover()
    {}

    SymShakeStructureMover::SymShakeStructureMover(core::scoring::ScoreFunctionOP s):
      ShakeStructureMover(s)
    {}

    SymShakeStructureMover::SymShakeStructureMover(core::scoring::ScoreFunctionOP s,
                                             core::Real temperature):
      ShakeStructureMover(s, temperature)
    {}

    SymShakeStructureMover::SymShakeStructureMover(core::scoring::ScoreFunctionOP s,
                                             core::Real ens_diversity, core::Real ens_div_tolerance):
			ShakeStructureMover(s, ens_diversity,ens_div_tolerance )
    {}

    SymShakeStructureMover::~SymShakeStructureMover(){}

//    void
//    SymShakeStructureMover::set_scorefunction(core::scoring::ScoreFunction & s){
//      scorefxn = s;
//      scorefunction_initialized=true;
//    }

	std::string
	SymShakeStructureMover::get_name() const {
		return "SymShakeStructureMover";
	}

		void
		SymShakeStructureMover::reduce_fa_rep(float fraction_fa_rep, core::scoring::ScoreFunction & s){
      s.set_weight( core::scoring::score_type_from_name("fa_rep"),
                    s.get_weight(core::scoring::score_type_from_name("fa_rep"))*fraction_fa_rep);
    }

    void
    SymShakeStructureMover::setup_for_run(core::pose::Pose & p){
      if(!get_scorefunction_initialized() ){
        core::scoring::ScoreFunctionOP s( new core::scoring::symmetry::SymmetricScoreFunction() );
        //std::cout << "[DEBUG]: scorefunction being initialized"<< std::endl;
        s->set_weight(core::scoring::score_type_from_name("rama"), 4.0);
        s->set_weight(core::scoring::score_type_from_name("omega"), 1.0);
        s->set_weight(core::scoring::score_type_from_name("fa_dun"), 1.0);
        s->set_weight(core::scoring::score_type_from_name("p_aa_pp"), 2.0);
        set_scorefunction(s);
      }
      setup_ca_constraints(p,(*get_scorefunction() ),9.0,get_harmonic_ca_cst_std_dev() );

			if(get_min_cst() || get_sc_min() ){
				set_min_scorefunction( core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );
			}

      if(get_mc_temperature() <= 0 && get_ensemble_ca_rmsd() > 0){
        //set mc_temp based on ensemble_ca_rmsd
				set_testing_phase(true);
        set_mc_temp( set_temp_based_on_ens_diversity(p,(*get_scorefunction() )) );
        set_is_properly_initialized(true);
      }else if(get_mc_temperature() > 0 && get_ensemble_ca_rmsd() < 0){
        //ready to go!
        set_is_properly_initialized(true);
      }else{
        //what to do, what to do?
        set_is_properly_initialized(false);
      }


    }

    void
    SymShakeStructureMover::minimize_with_constraints(core::pose::Pose & p,
                                                   core::scoring::ScoreFunction & s){

      core::optimization::symmetry::SymAtomTreeMinimizer min_struc;
      float const minimizer_tol = 0.0000001;
      core::optimization::MinimizerOptions options( "dfpmin_armijo_nonmonotone",
                                                    minimizer_tol,
                                                    true /*use_nb_list*/,
                                                    false /*deriv_check_in*/,
                                                    false /*deriv_check_verbose_in*/);
      options.nblist_auto_update( true );
			//      options.max_iter(5000);
      core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
      if(!get_sc_min() ){
        mm->set_bb(true);
      }else{
        mm->set_bb(false);
      }
      mm->set_chi(true);

			core::pose::symmetry::make_symmetric_movemap( p, *mm );

      if( get_ramp_fa_rep() ){
        core::scoring::symmetry::SymmetricScoreFunctionOP one_tenth_orig( core::scoring::symmetry::symmetrize_scorefunction(s) );
        reduce_fa_rep(0.1,*one_tenth_orig);
        min_struc.run(p,*mm,*one_tenth_orig,options);
        core::scoring::symmetry::SymmetricScoreFunctionOP one_third_orig( core::scoring::symmetry::symmetrize_scorefunction(s) );
        reduce_fa_rep(0.33,*one_third_orig);
        min_struc.run(p,*mm,*one_third_orig,options);
      }
			min_struc.run(p,*mm,s,options);

    }

		void
    SymShakeStructureMover::run_mc(core::pose::Pose & p, core::scoring::ScoreFunction & s,
                                core::Real temperature){

      using namespace protocols;
      using namespace simple_moves;

			core::pose::Pose init(p); //for comparison at the end

      core::Size nmoves = (core::Size)p.total_residue()/4; //number of moves for each move type

      core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
      mm->set_bb(true);

			// make symmetric movemap
			core::pose::symmetry::make_symmetric_movemap( p, *mm );

      simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( mm, temperature, nmoves) ) ;
      simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover(mm, temperature, nmoves) );

      setup_movers(small_mover,shear_mover,
                   0.2,0.2,0.4,
                   1.6,1.6,2.0);

      moves::RandomMoverOP apply_random_move( new moves::RandomMover() );
      apply_random_move->add_mover( small_mover, .5);
      apply_random_move->add_mover( shear_mover, .5);

      simple_moves::SmallMoverOP small_mover_low( new simple_moves::SmallMover( mm, (temperature*0.25), nmoves) ) ;
      simple_moves::ShearMoverOP shear_mover_low( new simple_moves::ShearMover(mm, (temperature*0.25), nmoves) );

      setup_movers(small_mover_low,shear_mover_low,
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
			protocols::moves::RepeatMoverOP full_cycle( new moves::RepeatMover( tm, get_nrounds() ) );
			full_cycle->apply( p );
//			mc->show_counters();

			if(!get_skip_low_temp_phase() ){
				mc->reset_counters();

				mc->set_lowest_score_pose(mc->last_accepted_pose()); //high energy mark for annealing back to native
				core::Real low_temp = temperature * 0.25;
				mc->set_temperature(low_temp);

				moves::TrialMoverOP ltm( new moves::TrialMover(apply_random_move_low,mc) );
				protocols::moves::RepeatMoverOP full_cycle_2( new moves::RepeatMover( ltm, get_nrounds() ) );
				full_cycle_2->apply( p );
				//mc->show_counters();
				p = mc->lowest_score_pose();
			}
			else{
				p = mc->last_accepted_pose();
			}

			if(get_min_cst() && !get_testing_phase() ){
				//remove constraints
				p.remove_constraints((p.constraint_set())->get_all_constraints());

				//std::cout << "minimizing with constraints" << std::endl;
				//temporarily hardcode score12 into minimization.
				setup_ca_constraints(p, (*get_min_scorefunction() ), 9.0, 0.5);
				minimize_with_constraints(p, (*get_min_scorefunction() ));

				//std::cout << " CA rmsd of current: " << core::scoring::CA_rmsd(p,init) << std::endl;
				//time_t time_per_decoy_finish = time(NULL);
				//std::cout << "time to finish decoy " << (time_per_decoy_finish-time_per_decoy) << std::endl;
			}else{
				//time_t time_per_decoy_finish = time(NULL);
				//std::cout << " CA rmsd of current: " << core::scoring::CA_rmsd(p,init) << std::endl;
				//std::cout << "time to finish decoy " << (time_per_decoy_finish-time_per_decoy) << std::endl;
			}
		}

} // namespace symmetry
} // namespace simple_moves
} // namespace protocols

