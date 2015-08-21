// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopMover_Perturb_CCD.cc
/// @brief kinematic loop closure main protocols
/// @author Chu Wang
/// @author Mike Tyka

//// Unit Headers
#include <protocols/loops/loops_main.hh>
//#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_CCDCreator.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/conformation/Residue.hh>

// Rosetta Headers
#include <core/chemical/VariantType.hh>

#include <core/conformation/util.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/fragment/FragSet.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <core/pack/task/TaskFactory.fwd.hh>

#include <core/pose/symmetry/util.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

#include <protocols/loops/util.hh>
#include <basic/Tracer.hh> // tracer output

//Utility and numeric Headers
#include <numeric/random/random.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>

// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>

//Auto Headers


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

namespace protocols {
namespace loops {
namespace loop_mover {
namespace perturb {

///////////////////////////////////////////////////////////////////////////////
using namespace core;

static thread_local basic::Tracer TR( "protocols.loops.loop_mover.perturb.LoopMover_Perturb_CCD" );

//constructors
LoopMover_Perturb_CCD::LoopMover_Perturb_CCD() :
	IndependentLoopMover()
{
	set_scorefxn( get_cen_scorefxn() );

	protocols::moves::Mover::type("LoopMover_Perturb_CCD");
	set_default_settings();
}


LoopMover_Perturb_CCD::LoopMover_Perturb_CCD(
	protocols::loops::LoopsOP loops_in
) : IndependentLoopMover( loops_in )
{
	set_scorefxn( get_cen_scorefxn() );
	protocols::moves::Mover::type("LoopMover_Perturb_CCD");
	set_default_settings();
}


LoopMover_Perturb_CCD::LoopMover_Perturb_CCD(
	protocols::loops::LoopsOP loops_in,
	core::scoring::ScoreFunctionOP scorefxn
) : IndependentLoopMover( loops_in )
{
	if ( scorefxn ) {
		set_scorefxn( scorefxn );
	} else {
		set_scorefxn( get_cen_scorefxn() );
	}

	protocols::moves::Mover::type("LoopMover_Perturb_CCD");
	set_default_settings();
}


LoopMover_Perturb_CCD::LoopMover_Perturb_CCD(
	protocols::loops::LoopsOP loops_in,
	core::scoring::ScoreFunctionOP scorefxn,
	core::fragment::FragSetOP fragset
) : IndependentLoopMover( loops_in )
{
	if ( scorefxn ) {
		set_scorefxn( scorefxn );
	} else {
		set_scorefxn( get_cen_scorefxn() );
	}

	protocols::moves::Mover::type("LoopMover_Perturb_CCD");
	set_default_settings();

	add_fragments( fragset );
}

//destructor
LoopMover_Perturb_CCD::~LoopMover_Perturb_CCD() {}

std::string
LoopMover_Perturb_CCD::get_name() const {
	return "LoopMover_Perturb_CCD";
}

void
LoopMover_Perturb_CCD::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Loops:\n" << get_loops();
	output << "Scorefunction: " << get_scorefxn()->get_name() << std::endl;
}

//clone
protocols::moves::MoverOP LoopMover_Perturb_CCD::clone() const {
	return protocols::moves::MoverOP( new LoopMover_Perturb_CCD(*this) );
}

loop_mover::LoopResult LoopMover_Perturb_CCD::model_loop(
	core::pose::Pose & pose,
	protocols::loops::Loop const & loop
) {
	using namespace scoring;
	using namespace optimization;
	using namespace basic::options;
	using namespace core::kinematics;
	using namespace protocols::simple_moves;

	bool const verbose( true );
	bool const local_debug( false );

	// store starting fold tree and cut pose
	kinematics::FoldTree f_orig=pose.fold_tree();
	set_single_loop_fold_tree( pose, loop );

	/// prepare fragment movers
	MoveMapOP movemap( new MoveMap() );
	movemap->set_bb_true_range( loop.start(), loop.stop() );

	Loops one_loop_loops;
	one_loop_loops.add_loop( loop );

	std::vector< FragmentMoverOP > fragmover;
	for ( std::vector< core::fragment::FragSetOP >::const_iterator
			it = frag_libs_.begin(), it_end = frag_libs_.end();
			it != it_end; ++it ) {
		ClassicFragmentMoverOP cfm( new ClassicFragmentMover( *it, movemap ) );
		cfm->set_check_ss( false );
		cfm->enable_end_bias_check( false );
		fragmover.push_back( cfm );
	}


	core::pose::Pose native_pose;
	if ( get_native_pose() ) {
		native_pose = *get_native_pose();
	} else {
		native_pose = pose;
	}

	Size const loop_begin( loop.start() ), loop_end( loop.stop() ), loop_cut( loop.cut() );
	Size const loop_size( loop_end - loop_begin + 1 );
	runtime_assert( loop.is_terminal( pose ) || pose.fold_tree().is_cutpoint( loop_cut ) );

	// see if we should skip large frags. Don't if that would mean no fragment insertions:
	// mt temporary removed this.
	////  bool skip_long_frags( true ); // PBHACK DONT CHECKIN if false
	////  {
	////   bool found_any_frags( false );
	////   frag_libs;
	////   for ( std::vector< core::fragment::ConstantLengthFragSetOP >::const_reverse_iterator
	////       it = frag_libs.rbegin(), it_end = frag_libs.rend(); it != it_end; ++it ) {
	////    Size const frag_size(  );
	////    if ( loop_size >= frag_size * 2 + 1 ) found_any_frags = true;
	////   }
	////   if ( !found_any_frags ) skip_long_frags = false;
	////  }

	tr() << "perturb_one_loop_with_ccd: " << loop_begin << ' ' << loop_size << std::endl;

	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );

	if ( !loop.is_terminal( pose ) ) {
		// set cutpoint variant for chainbreak scoring.
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, loop_cut );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, loop_cut+1 );
	}

	pose::Pose start_pose;
	start_pose = pose;

	if ( local_debug ) {
		std::ofstream out("score.tmp_input_cen");
		out << "scoring before cen_perturb: " << ( *scorefxn() )(pose) << std::endl;
		/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies(pose);
		scorefxn()->show( out );
		out << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
		tr() << "before cen_perturb: "
			<< pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
		out << pose.energies();
	}

	kinematics::MoveMap mm_one_loop;
	utility::vector1<bool> allow_sc_move( pose.total_residue(), false );
	loops_set_move_map( one_loop_loops, allow_sc_move, mm_one_loop);

	tr() << "loop rmsd before initial fragment perturbation:" << loop_rmsd( pose, native_pose, one_loop_loops ) << std::endl;
	if ( loop.is_extended() ) {
		if ( local_debug ) {
			pose.dump_pdb("tmp_cen_preidl.pdb");
			( *scorefxn() )(pose);
			tr() << "preidl: "
				<< pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
		}

		// avoid potential problems at the ends:
		if ( loop_begin > 1 && !pose.fold_tree().is_cutpoint( loop_begin-1 ) &&
				pose.residue( loop_begin ).is_bonded( loop_begin-1 ) ) {
			conformation::insert_ideal_bonds_at_polymer_junction( loop_begin-1, pose.conformation() );
			pose.set_omega( loop_begin - 1, init_omega );
		}
		if ( loop_end < pose.total_residue() && !pose.fold_tree().is_cutpoint( loop_end ) &&
				pose.residue( loop_end ).is_bonded( loop_end+1 ) ) {
			conformation::insert_ideal_bonds_at_polymer_junction( loop_end, pose.conformation() );
		}

		if ( local_debug ) {
			pose.dump_pdb("tmp_cen_endidl.pdb");
			( *scorefxn() )(pose);
			tr() << "endidl: "
				<< pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
		}

		for ( Size i= loop_begin; i<= loop_end; ++i ) {
			conformation::idealize_position( i, pose.conformation() );
		}

		if ( local_debug ) {
			pose.dump_pdb("tmp_cen_postidl.pdb");
			( *scorefxn() )(pose);
			tr() << "postidl: "
				<< pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
		}


		for ( Size i = loop_begin; i <= loop_end; ++ i ) {
			pose.set_phi  ( i, init_phi   );
			pose.set_psi  ( i, init_psi   );
			pose.set_omega( i, init_omega );
		}

		if ( local_debug ) {
			pose.dump_pdb("tmp_cen_init.pdb");
			( *scorefxn() )(pose);
			tr() << "extended: "
				<< pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
		}

		for ( core::Size i = loop.start(); i <= loop.stop(); ++i ) {
			for ( std::vector< FragmentMoverOP >::const_iterator
					it = fragmover.begin(),it_end = fragmover.end(); it != it_end; ++it ) {
				(*it)->apply( pose );
			}
		}


		if ( local_debug ) {
			pose.dump_pdb("tmp_cen_ran.pdb");
			( *scorefxn() )(pose);
			tr() << "random frags: "
				<< pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
		}
	}
	tr() << "loop rmsd after initial fragment perturbation:" << loop_rmsd( pose, native_pose, one_loop_loops ) << std::endl;

	// scheduler
	bool const fast = option[OptionKeys::loops::fast];
	int outer_cycles = option[ OptionKeys::loops::perturb_outer_cycles ]();
	int inner_cycles( fast ? std::min( Size(250), loop_size*5 ) : std::min( Size(1000), loop_size*20 ) );

	if ( option[OptionKeys::loops::debug] ) inner_cycles = 10;

	// Monte Carlo
	float const init_temp( 2.0 ), final_temp( 1.0 );
	float const gamma = std::pow( (final_temp/init_temp), (1.0f/(outer_cycles*inner_cycles)) );
	float temperature = init_temp;
	if ( local_debug ) { // hacking
		( *scorefxn() )(pose);
		tr() << "before mc ctor: "
			<< pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
	}

	protocols::moves::MonteCarlo mc( pose, *scorefxn(), temperature);

	// minimizer
	AtomTreeMinimizerOP minimizer;
	float const dummy_tol( 0.001 ); // linmin sets tol internally
	bool const use_nblist( false ), deriv_check( false ); // true ); // false );
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		// minimizer = dynamic_cast<AtomTreeMinimizer*> (new core::optimization::symmetry::SymAtomTreeMinimizer);
		minimizer = AtomTreeMinimizerOP( new core::optimization::symmetry::SymAtomTreeMinimizer );
	} else {
		minimizer = AtomTreeMinimizerOP( new core::optimization::AtomTreeMinimizer );
	}

	MinimizerOptions options( "linmin", dummy_tol, use_nblist, deriv_check);

	mc.show_scores();

	for ( int i=1; i<=outer_cycles; ++i ) {
		// recover low
		if ( local_debug ) { // debug
			( *scorefxn() )( pose );
			tr() << "befor rLOW: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) <<
				" rmsd: " << F(9,3,loop_rmsd( pose, native_pose, one_loop_loops )) << std::endl;
		}

		mc.recover_low(pose);

		if ( local_debug ) { // debug
			( *scorefxn() )( pose );
			tr() << "after rLOW: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) <<
				" rmsd: " << F(9,3,loop_rmsd( pose, native_pose, one_loop_loops )) << std::endl;
		}

		for ( int j=1; j<=inner_cycles; ++j ) {
			// change temperature
			temperature *= gamma;
			mc.set_temperature( temperature );
			// use rbegin and rend to iterate from larger size to smaller size
			//for ( std::map< Size, protocols::frags::TorsionFragmentLibraryOP >::const_reverse_iterator
			//    it = frag_libs.rbegin(), it_end = frag_libs.rend(); it != it_end; ++it ) {
			for ( std::vector< FragmentMoverOP >::const_iterator
					it = fragmover.begin(),it_end = fragmover.end(); it != it_end; ++it ) {

				// skip if the loop is too short
				//Size const frag_size( it->first );

				//if ( loop_size < frag_size || ( skip_long_frags && loop_size < frag_size * 2 + 1 ) ) continue;
				// monte carlo fragment ccd minimization
				if ( verbose ) {
					tr() << "perturb out/in/frag: "
						<< i << "/" << outer_cycles << " "
						<< j << "/" << inner_cycles << std::endl;
				}

				if ( local_debug ) { // debug
					( *scorefxn() )( pose );
					tr() << "befor frag: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) <<
						" rmsd: " << F(9,3,loop_rmsd( pose, native_pose, one_loop_loops )) << std::endl;
				}

				//insert_fragment( loop_begin, loop_end, pose, it->second );
				(*it)->apply( pose );

				if ( local_debug ) { // debug
					( *scorefxn() )( pose );
					tr() << "after frag: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) <<
						" rmsd: " << F(9,3,loop_rmsd( pose, native_pose, one_loop_loops )) << std::endl;
				}

				if ( !loop.is_terminal( pose ) ) ccd_close_loops( pose, one_loop_loops, mm_one_loop);

				if ( local_debug ) { // debug
					( *scorefxn() )( pose );
					tr() << "after ccd:  " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) <<
						" rmsd: " << F(9,3,loop_rmsd( pose, native_pose, one_loop_loops )) << std::endl;
				}

				minimizer->run( pose, mm_one_loop, *scorefxn(), options );

				if ( local_debug ) { // debug
					( *scorefxn() )( pose );
					tr() << "after min:  " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) <<
						" rmsd: " << F(9,3,loop_rmsd( pose, native_pose, one_loop_loops )) << std::endl;
				}

				std::string move_type = "frag_ccd_min";
				mc.boltzmann( pose, move_type );
				mc.show_scores();
			} // for each fragment size
		} // inner_cycles
	} // outer_cycles
	pose = mc.lowest_score_pose();

	if ( local_debug ) {
		std::ofstream out("score.tmp_perturb_cen");
		out << "scoring after cen_perturb: " << ( *scorefxn() )(pose) << std::endl;
		/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies(pose);
		scorefxn()->show( out );
		out << pose.energies();
	}

	// return to original fold tree
	pose.fold_tree( f_orig );
	return loop_mover::Success;

}

protocols::loops::LoopsCOP LoopMover_Perturb_CCD::get_loops() const {
	return loops();
}

core::scoring::ScoreFunctionOP LoopMover_Perturb_CCD::get_scorefxn() const {
	return scorefxn();
}

basic::Tracer & LoopMover_Perturb_CCD::tr() const
{
	return TR;
}

LoopMover_Perturb_CCDCreator::~LoopMover_Perturb_CCDCreator() {}

moves::MoverOP LoopMover_Perturb_CCDCreator::create_mover() const {
	return moves::MoverOP( new LoopMover_Perturb_CCD() );
}

std::string LoopMover_Perturb_CCDCreator::keyname() const {
	return "LoopMover_Perturb_CCD";
}

std::ostream &operator<< ( std::ostream &os, LoopMover_Perturb_CCD const &mover )
{
	mover.show(os);
	return os;
}

} // namespace perturb
} // namespace loop_mover
} // namespace loops
} // namespace protocols
