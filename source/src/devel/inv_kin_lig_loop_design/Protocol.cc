// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/InvKinLigLoopDesign/Protocol.cc
///
/// @brief
/// @author

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>


#include <devel/inv_kin_lig_loop_design/Protocol.hh>
#include <devel/inv_kin_lig_loop_design/Mover.hh>
#include <devel/inv_kin_lig_loop_design/JumpManager.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/RamaCheck.hh>
#include <protocols/loops/loops_main.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/pack/pack_rotamers.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#define REPORTCYCLES1(c0,m0)                    " Cycles:\t" << c0 << "/" << m0
#define REPORTCYCLES2(c0,m0,c1,m1)              " Cycles:\t" << c0 << "/" << m0 << "\t" << c1 << "/" << m1
#define REPORTCYCLES3(c0,m0,c1,m1,c2,m2)        " Cycles:\t" << c0 << "/" << m0 << "\t" << c1 << "/" << m1 << "\t" << c2 << "/" << m2
#define REPORTCYCLES4(c0,m0,c1,m1,c2,m2,c3,m3)  " Cycles:\t" << c0 << "/" << m0 << "\t" << c1 << "/" << m1 << "\t" << c2 << "/" << m2 << "\t" << c3 << "/" << m3


namespace devel {

namespace inv_kin_lig_loop_design {


namespace {

bool contains(vector<Loop> const& loops, int seqpos ) {
	 for ( auto const & loop : loops ) {
		if ( loop.lo <= seqpos && seqpos <= loop.hi ) {
			return true;
		}
	}
	return false;
}

bool is_anchor(vector<Loop> const& loops, int seqpos ) {
	 for ( auto const & loop : loops ) {
		if ( loop.to == seqpos ) {
			return true;
		}
	}
	return false;
}

core::kinematics::MoveMap get_move_map(int num_res, vector<Loop> const& loops) {
	core::kinematics::MoveMap rval;
	for ( int i = 1; i <= num_res; ++i ) {
		rval.set_bb( i, contains(loops,i) );
		rval.set_chi( i, is_anchor(loops,i) );

		//cout << "get_move_map: " << i << " " << contains(loops,i) << " " << is_anchor(loops,i) << endl;

	}
	rval.set_jump( false );
	return rval;
}

}

//     Protocol::Protocol() {
//     }

// const double Protocol::DEFAULT_MIN_ACCEPTANCE_RATE = 5e-2;
//     const double Protocol::DEFAULT_MIN_WEIGHT_REP      = 1e-2;

void Protocol::setParams(TagCOP tag0, vector<Loop> const& loops_) {

	loops = loops_;

	max_cycles               = tag0->getTag("param")->getOption<int>("max_cycles",10);
	max_cycles_start         = tag0->getTag("param")->getOption<int>("max_cycles_start",10);
	max_cycles_start_3mer    = tag0->getTag("param")->getOption<int>("max_cycles_start_3mer",1000);
	max_cycles_start_3mer_anchor = tag0->getTag("param")->getOption<int>("max_cycles_start_3mer_anchor",100);
	max_cycles_start_3mer_1mer   = tag0->getTag("param")->getOption<int>("max_cycles_start_3mer_1mer",100);
	max_cycles_start_1mer    = tag0->getTag("param")->getOption<int>("max_cycles_start_1mer",100);
	max_cycles_ramp          = tag0->getTag("param")->getOption<int>("max_cycles_ramp",10);
	max_cycles_ramp_sm       = tag0->getTag("param")->getOption<int>("max_cycles_ramp_sm",1);
	max_cycles_ramp_sm_min   = tag0->getTag("param")->getOption<int>("max_cycles_ramp_sm_min",10);
	max_cycles_design        = tag0->getTag("param")->getOption<int>("max_cycles_design",10);
	max_cycles_design_sm     = tag0->getTag("param")->getOption<int>("max_cycles_design_sm",10);
	max_cycles_design_sm_min = tag0->getTag("param")->getOption<int>("max_cycles_design_sm_min",10);

} // Protocol::setParams

void Protocol::phase_lores() {

	//setLoRes();

	core::optimization::AtomTreeMinimizer atm;

	double temp = 1;

	core::pose::PoseOP pose_init( new core::pose::Pose );
	(*pose_init) = (*pose);

	protocols::moves::MonteCarlo min_start(*pose, *score_fxn_lores, temp);

	// StageI / start - find an acceptable starting configuration for subsequent ljrep optimization

	JumpManager jm;
	Mover mover(pose.get());

	for ( int cycles_start = 0; cycles_start < max_cycles_start; ++cycles_start ) {
		cout << ">>> StageI / start - " << REPORTCYCLES1(cycles_start,max_cycles_start) << endl;

		(*pose) = (*pose_init);

		protocols::moves::MonteCarlo min_start_3mer(*pose,*score_fxn_lores,temp);

		cout << ">>> StageI / start / 3mer - inserting 3mers " << endl;

		for ( int cycles_start_3mer = 0; cycles_start_3mer < max_cycles_start_3mer /* && get_acceptance_rate(min_start_3mer,min_acceptance_rate_start_3mer)*/ ; ++cycles_start_3mer ) {

			//     if( getVerbose() ) {
			cout << "start_3mer: " << REPORTCYCLES2(cycles_start,max_cycles_start,cycles_start_3mer,max_cycles_start_3mer) << endl;
			//     }

			(*pose) = (*pose_init);

			vector<Loop> loops_shuffled = loops; // need a list of all loops not just anchored loops
			numeric::random::random_permutation(loops_shuffled.begin(),loops_shuffled.end(),numeric::random::rg());

			//g.setUseNb(false);

			// start with a completely random configuration
			for ( auto & loop : loops_shuffled ) { // !!! first generate a completely random
				mover.randomFragments( loop.lo, loop.hi, 1 );
				mover.randomFragments( loop.lo, loop.hi, 3 );

				if ( loop.to != 0 ) {

					mover.randomRotamer( loop.to );

					if ( loop.tag->hasTag("hbond") ) {
						jm.set_random_hbond_jump( *pose, loop );
					}

				} // cycles_start_3mer_anchor

			}

			protocols::moves::MonteCarlo min_start_3mer_inner(*pose,*score_fxn_lores,temp);

			min_start_3mer_inner.boltzmann( *pose );
			min_start_3mer_inner.recover_low( *pose );

			for ( auto & loop : loops_shuffled ) {

				int const a = loop.lo;
				int const b = loop.hi;

				if ( loop.to != 0 ) {

					for ( int cycles_start_3mer_anchor = 0; cycles_start_3mer_anchor < max_cycles_start_3mer_anchor; ++cycles_start_3mer_anchor ) {

						mover.randomRotamer( loop.to );
						min_start_3mer_inner.boltzmann( *pose );

					}

				}

				min_start_3mer_inner.recover_low( *pose );

				for ( int cycles_start_3mer_1mer = 0; cycles_start_3mer_1mer < max_cycles_start_3mer_1mer; ++cycles_start_3mer_1mer ) {

					cout << "start_3mer_1mer: " << REPORTCYCLES3(cycles_start,max_cycles_start,cycles_start_3mer,max_cycles_start_3mer,cycles_start_3mer_1mer,max_cycles_start_3mer_1mer) << endl;

					mover.insertFragment(a,b,1);
					atm.run(*pose,*move_map,*score_fxn_lores,core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",1,true,false,false) );

					min_start_3mer_inner.boltzmann(*pose);

					pose->energies().show_total_headers( cout ); cout << endl;
					pose->energies().show_totals( cout ); cout << endl;

				} // cycles_start_3mer_anchor

				min_start_3mer_inner.recover_low(*pose); // ???

			} // k

			min_start_3mer.boltzmann(*pose);

		} // cycles_start_3mer

		min_start_3mer.recover_low(*pose);

		// Stage I / start / 1mer - optimize the given conformation by 1mer & ccd moves
		cout << ">>> StageI / start / 1mer - inserting 1mers" << endl;

		protocols::moves::MonteCarlo min_start_1mer(*pose,*score_fxn_lores,temp);

		// Store the initial FoldTree so it can be set back to the pose after loop closure
		core::kinematics::FoldTree orig_ft = pose->fold_tree();

		for ( int cycles_start_1mer = 0; cycles_start_1mer < max_cycles_start_1mer /*&& get_acceptance_rate(min_start_1mer,DEFAULT_MIN_ACCEPTANCE_RATE)*/ ; ++cycles_start_1mer ) {

			cout << "Stage I / start / 1mer " << REPORTCYCLES2(cycles_start,max_cycles_start,cycles_start_1mer,max_cycles_start_1mer) << endl;

			pose->energies().show_total_headers( cout ); cout << endl;
			pose->energies().show_totals( cout ); cout << endl;

			const int k_rand_loop = numeric::random::rg().random_range(0,INT_MAX) % loops.size(); // this is biased towards shorter loops

			const Loop& loop = loops[k_rand_loop];

			mover.smallMoves(loop.lo,loop.hi);
			//mover.insertFragment(a,b,1);

			atm.run(*pose,*move_map,*score_fxn_lores,core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.01,true,false,false) );

			if ( loop.to ) {
				protocols::loops::loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover;
				ccd_loop_closure_mover.movemap( move_map );
				ccd_loop_closure_mover.tolerance( 0.1 );
				ccd_loop_closure_mover.rama()->max_rama_score_increase( 5.0 );
				ccd_loop_closure_mover.max_total_torsion_delta_per_residue( 5.0, 5.0, 5.0 );

				// Set the loop to the CCD mover and then generate an appropriate FoldTree for loop closure
				ccd_loop_closure_mover.loop ( protocols::loops::Loop( loop.lo-1, loop.to-1, loop.lo ) );
				protocols::loops::set_single_loop_fold_tree( *pose, ccd_loop_closure_mover.loop() );
				protocols::loops::add_cutpoint_variants( *pose );
				ccd_loop_closure_mover.apply( *pose );
				protocols::loops::remove_cutpoint_variants( *pose );

				// Set the loop to the CCD mover and then generate an appropriate FoldTree for loop closure
				// The current FoldTree will be thrown asunder
				ccd_loop_closure_mover.loop ( protocols::loops::Loop( loop.to+1, loop.hi+1, loop.hi ) );
				protocols::loops::set_single_loop_fold_tree( *pose, ccd_loop_closure_mover.loop() );
				protocols::loops::add_cutpoint_variants( *pose );
				ccd_loop_closure_mover.apply( *pose );
				protocols::loops::remove_cutpoint_variants( *pose );

				// Set the pose's FoldTree back to the incoming FoldTree
				pose->fold_tree( orig_ft );

				(*score_fxn_lores)( *pose );
				pose->energies().show_total_headers( cout ); cout << endl;
				pose->energies().show_totals( cout ); cout << endl;

			} else {
				assert( false );
			}

			atm.run(*pose,*move_map,*score_fxn_lores,core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.01,true,false,false) );

			pose->energies().show_total_headers( cout ); cout << endl;
			pose->energies().show_totals( cout ); cout << endl;

			min_start_1mer.boltzmann(*pose);
			//min_start_1mer.recover_low(*pose);

		} // cycles_start_1mer

		min_start_1mer.recover_low(*pose);
		min_start.boltzmann(*pose );

	} // cycles_start

	min_start.recover_low(*pose);

}

void Protocol::phase_hires() {
	// let's do this in a random order
	vector<Loop> loops_shuffled = loops;
	numeric::random::random_permutation(loops_shuffled.begin(),loops_shuffled.end(),numeric::random::rg());

	Mover mover(pose.get());
	core::optimization::AtomTreeMinimizer atm;

	double temp = 1;

	const double MIN_WEIGHT_REP = 0.001;
	const double MAX_WEIGHT_REP = score_fxn_hires->get_weight( core::scoring::fa_rep );

	// Store the initial FoldTree so it can be set back to the pose after loop closure
	core::kinematics::FoldTree orig_ft = pose->fold_tree();

	// LJRep Ramp up
	for ( int cycles_ramp = 0; cycles_ramp <= max_cycles_ramp; ++cycles_ramp ) {

		cout << ">>> Stage II / start / ramp - " << REPORTCYCLES1(cycles_ramp,max_cycles_ramp) << endl;

		double f;
		if ( cycles_ramp == 0 ) {
			f = MIN_WEIGHT_REP;
		} else {
			f = MAX_WEIGHT_REP*(static_cast<double>(cycles_ramp + 1) / max_cycles_ramp);
		}

		score_fxn_hires->set_weight( core::scoring::fa_rep, f );

		protocols::moves::MonteCarlo min_ramp(*pose, *score_fxn_hires, temp);

		min_ramp.boltzmann( *pose );

		for ( Size k = 0; k < loops.size(); ++k ) {
			const Loop& loop = loops_shuffled[k];

			int const a = loop.lo;
			int const b = loop.hi;

			cout << "hi: " << a << " " << b << endl;

			for ( int cycles_ramp_sm = 0; cycles_ramp_sm < max_cycles_ramp_sm; ++cycles_ramp_sm ) {

				mover.smallMoves(a,b);

				for ( int cycles_ramp_sm_min = 0; cycles_ramp_sm_min < max_cycles_ramp_sm_min ; ++ cycles_ramp_sm_min ) {

					if ( loop.to ) {
						protocols::loops::loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover;
						ccd_loop_closure_mover.movemap( move_map );
						ccd_loop_closure_mover.tolerance( 0.1 );
						ccd_loop_closure_mover.rama()->max_rama_score_increase( 5.0 );
						ccd_loop_closure_mover.max_total_torsion_delta_per_residue( 5.0, 5.0, 5.0 );

						// Set the loop to the CCD mover and then generate an appropriate FoldTree for loop closure
						ccd_loop_closure_mover.loop ( protocols::loops::Loop( loop.lo-1, loop.to-1, loop.lo ) );
						protocols::loops::set_single_loop_fold_tree( *pose, ccd_loop_closure_mover.loop() );
						protocols::loops::add_cutpoint_variants( *pose );

						ccd_loop_closure_mover.apply( *pose );
						protocols::loops::remove_cutpoint_variants( *pose );

						cout << "dev=" << ccd_loop_closure_mover.deviation() <<
							" tdelta=" << ccd_loop_closure_mover.torsion_delta() <<
							" rdelta=" << ccd_loop_closure_mover.rama_delta() << endl;

						// Set the loop to the CCD mover and then generate an appropriate FoldTree for loop closure
						// The current FoldTree will be thrown asunder
						ccd_loop_closure_mover.loop ( protocols::loops::Loop( loop.to+1, loop.hi+1, loop.hi ) );
						protocols::loops::set_single_loop_fold_tree( *pose, ccd_loop_closure_mover.loop() );
						protocols::loops::add_cutpoint_variants( *pose );

						ccd_loop_closure_mover.apply( *pose );
						protocols::loops::remove_cutpoint_variants( *pose );

						// Set the pose's FoldTree back to the incoming FoldTree
						pose->fold_tree( orig_ft );

						cout << "dev=" << ccd_loop_closure_mover.deviation() <<
							" tdelta=" << ccd_loop_closure_mover.torsion_delta() <<
							" rdelta=" << ccd_loop_closure_mover.rama_delta() << endl;
					} else {
						assert( false );
					}

					atm.run(*pose,*move_map,*score_fxn_lores,core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",1,true,false,false) );

					min_ramp.boltzmann(*pose); // !!! don't recover here

				} // cycles_ramp_min

				min_ramp.boltzmann(*pose);

			} // cycles_ramp_small_moves

		} // cycles_ramp_loop

		min_ramp.recover_low(*pose);

	} // cycles_rep

	// stage 3 - design

	utility::vector1<bool> residues_allowed_to_be_packed( pose->n_residue(), false);

	core::pack::task::PackerTaskOP design_task( core::pack::task::TaskFactory::create_packer_task( *pose ) );
	//design_task->initialize_from_command_line(); // .read_resfile().or_include_current( true );

	for ( auto & loop : loops ) {
		for ( int i = loop.lo; i <= loop.hi ; ++i ) {
			if ( i != loop.to ) {
				residues_allowed_to_be_packed[i] = true;
			}
		}
	}

	design_task->restrict_to_residues( residues_allowed_to_be_packed );

	for ( int cycles_design = 0; cycles_design < max_cycles_design; ++cycles_design ) {

		cout << "Protocol::execute - DESIGN - starting cycle #" << cycles_design << " / " << max_cycles_design << endl;

		cout << "calling pack_rotamers..." << endl;
		core::pack::pack_rotamers( *pose, *score_fxn_hires, design_task );
		cout << "done calling pack_rotamers" << endl;

		protocols::moves::MonteCarlo min_design_sm(*pose, *score_fxn_hires, temp);

		for ( Size k = 0; k < loops_shuffled.size(); ++k ) {

			cout << "Protocol::execute - DESIGN - small moves + min after design cycle# " << cycles_design << " small moves + insertions into loop#" << k << endl;

			Loop& loop = loops_shuffled[k];

			int const a = loop.lo;
			int const b = loop.hi;

			protocols::moves::MonteCarlo min_design_sm(*pose,*score_fxn_hires,temp);

			for ( int cycles_design_sm = 0; cycles_design_sm < max_cycles_design_sm; ++cycles_design_sm ) {

				mover.smallMoves(a,b);

				for ( int cycles_design_sm_min = 0; cycles_design_sm_min < max_cycles_design_sm_min; ++cycles_design_sm_min ) {

					if ( loop.to ) {
						protocols::loops::loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover;
						ccd_loop_closure_mover.movemap( move_map );
						ccd_loop_closure_mover.tolerance( 0.1 );
						ccd_loop_closure_mover.rama()->max_rama_score_increase( 5.0 );
						ccd_loop_closure_mover.max_total_torsion_delta_per_residue( 5.0, 5.0, 5.0 );

						// Set the loop to the CCD mover and then generate an appropriate FoldTree for loop closure
						ccd_loop_closure_mover.loop ( protocols::loops::Loop( loop.lo-1, loop.to-1, loop.lo ) );
						protocols::loops::set_single_loop_fold_tree( *pose, ccd_loop_closure_mover.loop() );
						protocols::loops::add_cutpoint_variants( *pose );

						ccd_loop_closure_mover.apply( *pose );
						protocols::loops::remove_cutpoint_variants( *pose );


						// Set the loop to the CCD mover and then generate an appropriate FoldTree for loop closure
						// The current FoldTree will be thrown asunder
						ccd_loop_closure_mover.loop ( protocols::loops::Loop( loop.to+1, loop.hi+1, loop.hi ) );
						protocols::loops::set_single_loop_fold_tree( *pose, ccd_loop_closure_mover.loop() );
						protocols::loops::add_cutpoint_variants( *pose );

						ccd_loop_closure_mover.apply( *pose );
						protocols::loops::remove_cutpoint_variants( *pose );

						// Set the pose's FoldTree back to the incoming FoldTree
						pose->fold_tree( orig_ft );

					} else {
						assert( false );
					}

					atm.run(*pose,*move_map,*score_fxn_lores,core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",1,true,false,false) );

					min_design_sm.boltzmann(*pose);

				} // cycles_design_sm_min

			} // cycles_design_sm

		} // k

		min_design_sm.recover_low(*pose);

	} // cycles_design

}


void Protocol::apply(core::pose::PoseOP pose_) {
	pose = pose_;
	move_map = core::kinematics::MoveMapOP( new core::kinematics::MoveMap( get_move_map(pose->n_residue(),loops) ) );

	const double CHAINBREAK_WEIGHT = 10;
	const double RAMA_WEIGHT       = 10;
	const double DUNBRACK_WEIGHT   = 1;
	const double LORES_REP_WEIGHT  = 0.001;

	score_fxn_lores = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
	score_fxn_lores->set_weight( core::scoring::chainbreak, CHAINBREAK_WEIGHT );
	score_fxn_lores->set_weight( core::scoring::rama,       RAMA_WEIGHT       );
	score_fxn_lores->set_weight( core::scoring::fa_dun,     DUNBRACK_WEIGHT   );
	score_fxn_lores->set_weight( core::scoring::fa_rep,     LORES_REP_WEIGHT  );

	//score_fxn_lores = core::scoring::ScoreFunctionFactory::create_score_function( "murphp_lores" );
	//score_fxn_lores = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::SOFT_REP_WTS );

	score_fxn_hires = core::scoring::get_score_function();
	score_fxn_hires->set_weight( core::scoring::chainbreak, CHAINBREAK_WEIGHT );
	score_fxn_hires->set_weight( core::scoring::rama,       RAMA_WEIGHT       );

	//for( int i = 0; i <

	cout << "Protocol::apply: running lo_res phase" << endl;

	phase_lores();

	cout << "Protocol::apply: running hi_res_phase" << endl;

	phase_hires();


	// Score new structure.  Cached energies (including *residue* energies)
	// must be up-to-date in order to get sensible output.  If you remove these
	// lines, you *must* insert equivalent logic at the end of all apply() methods
	// (or at least for all movers that might be passed to this function).
	(*score_fxn_hires)( *pose );
	//score_fxn_hires->accumulate_residue_total_energies( *pose );


	/*

	core::scoring::Energies const& energies = pose->energies();

	cout << "Protocol::apply::getting score_types:" << endl;

	core::scoring::ScoreTypes score_types = score_fxn_hires->get_nonzero_weighted_scoretypes();


	cout << "Protocol::apply::ScoreTypes: " << score_types.size() << ": ";
	for( Size i = 1; i <= score_types.size(); ++i ) {
	cout << score_types[i] << ",";
	}
	cout << endl;


	cout << "Protocol::apply::Energies:" << endl;

	for( Size i = 1; i <= pose->n_residue(); ++i ) {

	//core::scoring::EnergyMap emap = energies();

	core::scoring::EMapVector const& emap_1b = energies.onebody_energies(i);

	cout << "i " << i << " emap_1b:\t";
	for( Size k = 1; k <= score_types.size(); ++k ) {
	if( score_types[k] == core::scoring::chainbreak ) {
	cout << "cb: ";
	cout << emap_1b[ score_types[k] ] << "\t";
	}
	}
	cout << " ";


	core::scoring::EMapVector const& emap_2b = energies.residue_total_energies(i);
	cout << " emap_2b:\t";
	for( Size k = 1; k <= score_types.size(); ++k ) {
	if( score_types[k] == core::scoring::chainbreak ) {
	cout << emap_2b[ score_types[k] ] << "\t";
	}
	}

	cout << endl;

	}


	for( Size k = 1; k <= score_types.size(); ++k ) {
	cout << score_types[k] << " " << energies.total_energies()[ score_types[k] ] << endl;
	}

	//       for( int i = 1; i <= pose->n_residues(); ++i ) {
	// for( int j = 1; j <= pose->n_residues(); ++j ) {
	// }
	//       }

	//cout << "done with pose->energies" << endl;

	*/

}


}

}
