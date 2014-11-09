// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_mover/perturb/LoopMover_Perturb_QuickCCD.cc
/// @brief kinematic loop closure main protocols
/// @author Mike Tyka
/// @author James Thompson

//// Unit Headers
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_QuickCCD.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_QuickCCDCreator.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/conformation/Residue.hh>
//// Rosetta Headers
#include <core/chemical/VariantType.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
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

#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <basic/Tracer.hh> // tracer output
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>

//Utility Headers
#include <numeric/random/random.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>
#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <utility/vector1.hh>

//Auto Headers


namespace protocols {
namespace loops {
namespace loop_mover {
namespace perturb {

static thread_local basic::Tracer TR( "protocols.loops.loop_mover.perturb.LoopMover_Perturb_QuickCCD" );

///////////////////////////////////////////////////////////////////////////////
using namespace core;


//constructors

LoopMover_Perturb_QuickCCD::LoopMover_Perturb_QuickCCD() :
	IndependentLoopMover()
{
	set_scorefxn( get_cen_scorefxn() );

	protocols::moves::Mover::type("LoopMover_Perturb_QuickCCD");
	set_default_settings();
}

LoopMover_Perturb_QuickCCD::LoopMover_Perturb_QuickCCD(
	protocols::loops::LoopsOP loops_in
) : IndependentLoopMover( loops_in )
{
	set_scorefxn( get_cen_scorefxn() );

	protocols::moves::Mover::type("LoopMover_Perturb_QuickCCD");
	set_default_settings();
}



LoopMover_Perturb_QuickCCD::LoopMover_Perturb_QuickCCD(
	protocols::loops::LoopsOP loops_in,
	core::scoring::ScoreFunctionOP scorefxn
) : IndependentLoopMover( loops_in )
{
	if ( scorefxn ) {
		set_scorefxn(  scorefxn );
	} else {
		set_scorefxn( get_cen_scorefxn() );
	}

	protocols::moves::Mover::type("LoopMover_Perturb_QuickCCD");
	set_default_settings();
}

LoopMover_Perturb_QuickCCD::LoopMover_Perturb_QuickCCD(
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

	protocols::moves::Mover::type("LoopMover_Perturb_QuickCCD");
	set_default_settings();
	add_fragments( fragset );
}

//destructors
LoopMover_Perturb_QuickCCD::~LoopMover_Perturb_QuickCCD(){}

std::string
LoopMover_Perturb_QuickCCD::get_name() const {
	return "LoopMover_Perturb_QuickCCD";
}

//clone
protocols::moves::MoverOP LoopMover_Perturb_QuickCCD::clone() const {
		return protocols::moves::MoverOP( new LoopMover_Perturb_QuickCCD( *this ) );
}

LoopResult LoopMover_Perturb_QuickCCD::model_loop(
	core::pose::Pose & pose,
	protocols::loops::Loop const & loop
) {
	using namespace kinematics;
	using namespace scoring;
	using namespace optimization;
	using namespace protocols::simple_moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace numeric::random;

	tr().Info << "***** CCD CLOSURE *****" << std::endl;
	bool debug = option[ basic::options::OptionKeys::loops::debug ]();

	core::Size const nres =  pose.total_residue();

	// store starting fold tree and cut pose
 	kinematics::FoldTree f_orig=pose.fold_tree();
	set_single_loop_fold_tree( pose, loop );

	// generate movemap after fold_tree is set
	kinematics::MoveMapOP mm_one_loop( new kinematics::MoveMap() );
	kinematics::MoveMapOP mm_one_loop_symm( new kinematics::MoveMap() );
	set_move_map_for_centroid_loop( loop, *mm_one_loop );
	set_move_map_for_centroid_loop( loop, *mm_one_loop_symm );
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		core::pose::symmetry::make_symmetric_movemap( pose, *mm_one_loop_symm );
	}

	int const loop_size( loop.stop() - loop.start() + 1 );
	int cycles2 = 2;
	// Minimum of 15 cycles
	// Maximum of 5*15 cycles;
	int base_cycles( std::max( 15, static_cast<int>( 5*std::min(loop_size,15) )));
	int cycles3 = base_cycles;
	tr().Info << "Number of cycles: cycles2 and cycles3 "
		<< cycles2 << " " << cycles3
		<< std::endl;

	//bool chainbreak_present =  ( loop.start() != 1 && loop.stop() != nres );
	// special case ... vrt res at last position
	//chainbreak_present &= (loop.stop() != nres-1 || pose.residue( nres ).aa() != core::chemical::aa_vrt );
	bool chainbreak_present = ( loop.start() != 1 && loop.stop() != nres &&
	                            !pose.residue( loop.start() ).is_lower_terminus() &&
	                            !pose.residue( loop.stop() ).is_upper_terminus() );

	// set loop.cut() variant residue for chainbreak score if chanbreak is present
	if ( chainbreak_present ) {
		core::pose::add_variant_type_to_pose_residue(
			pose, chemical::CUTPOINT_LOWER, loop.cut()
		);
		core::pose::add_variant_type_to_pose_residue(
			pose, chemical::CUTPOINT_UPPER, loop.cut()+1
		);
	}


	( *scorefxn() )(pose);
	core::pose::Pose start_pose = pose;

	// either extend or at least idealize the loop (just in case).
	if ( loop.is_extended() ) set_extended_torsions( pose, loop );
	else                      idealize_loop( pose, loop );

	// omega needs to be moveable for FragmentMovers, so use a separate movemap
	kinematics::MoveMapOP frag_mover_movemap( new kinematics::MoveMap() );
	frag_mover_movemap->set_bb_true_range( loop.start(), loop.stop() );
	//fpd need to leave movemap for fragment insertion nonsymmetric
	//    (so that insertions will happen when scoring subunit is not 1st)

	utility::vector1< FragmentMoverOP > fragmover;
	for ( utility::vector1< core::fragment::FragSetOP >::const_iterator
				it = frag_libs().begin(), it_end = frag_libs().end();
				it != it_end; it++ ) {
		ClassicFragmentMoverOP cfm( new ClassicFragmentMover( *it, frag_mover_movemap ) );
		cfm->set_check_ss( false );
		cfm->enable_end_bias_check( false );
		fragmover.push_back( cfm );
	}

	core::Real m_Temperature_ = 2.0;
	moves::MonteCarloOP mc_( new moves::MonteCarlo( pose, *scorefxn(), m_Temperature_ ) );

	if ( randomize_loop_ ) {
		// insert random fragment as many times as the loop is long (not quite
		// the exact same as the old code)
		for ( Size i = loop.start(); i <= loop.stop(); ++i ) {
			for ( utility::vector1< FragmentMoverOP >::const_iterator
					it = fragmover.begin(),it_end = fragmover.end(); it != it_end; it++
			) {
				(*it)->apply( pose );
			}
		}
	}

	// Set up Minimizer object
	AtomTreeMinimizerOP mzr;
	float const dummy_tol( 0.001 ); // linmin sets tol internally
	bool const use_nblist( true ), deriv_check( false ); // true ); // false );
	MinimizerOptions options( "linmin", dummy_tol, use_nblist, deriv_check);
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		mzr = AtomTreeMinimizerOP( new core::optimization::symmetry::SymAtomTreeMinimizer );
	} else {
		mzr = AtomTreeMinimizerOP( new core::optimization::AtomTreeMinimizer );
	}


	// Set up MonteCarlo Object
	core::Real const init_temp = 2.0;
	core::Real temperature = init_temp;
 	mc_->reset( pose );
	mc_->set_temperature( temperature );

  int   starttime    = time(NULL);
  int   frag_count   = 0;
	scorefxn()->show_line_headers( tr().Info );
	tr().Info << std::endl;

	core::Real const final_temp( 1.0 );
	core::Real const gamma = std::pow(
		(final_temp/init_temp), (1.0/(cycles2*cycles3))
	);

	bool 	has_constraints         = pose.constraint_set()->has_constraints();
	float final_constraint_weight = option[ basic::options::OptionKeys::constraints::cst_weight ](); // this is really stupid.

	if ( chainbreak_present ) scorefxn()->set_weight( chainbreak, 1.0 );

	// for symmetric case, upweight the chainbreak by # of monomers
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::conformation::symmetry::SymmetricConformation const & symm_conf (
					dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()) );
		core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		scorefxn()->set_weight( chainbreak, symm_info->subunits() );
	}

	for( int c2 = 1; c2 <= cycles2; ++c2 ) {
		mc_->recover_low( pose );

		// ramp up constraints
		if ( has_constraints ) {
			if( c2 != cycles2 ) {
				scorefxn()->set_weight( core::scoring::atom_pair_constraint, final_constraint_weight*float(c2)/float( cycles2 ) );
			} else {
				scorefxn()->set_weight( core::scoring::atom_pair_constraint, final_constraint_weight * 0.2 );
			}
		}

		( *scorefxn() )(pose);
		if ( tr().visible() ) { scorefxn()->show_line( tr().Info , pose ); }
		tr().Info << std::endl;
		mc_->score_function( *scorefxn() );

		for ( int c3 = 1; c3 <= cycles3; ++c3 ) {
			temperature *= gamma;
			mc_->set_temperature( temperature );
			for ( utility::vector1< FragmentMoverOP >::const_iterator
					it = fragmover.begin(),it_end = fragmover.end();
					it != it_end; it++
			) {
				(*it)->apply( pose );
				if ( chainbreak_present ) {
					fast_ccd_close_loops( pose, loop,  *mm_one_loop );
				}

				mzr->run( pose, *mm_one_loop_symm, *scorefxn(), options );
				mc_->boltzmann( pose, "QuickCCD" );
				frag_count++;
			}

			if ( debug ) {
				( *scorefxn() )(pose);
				scorefxn()->show_line( tr().Info , pose );
				tr().Info << std::endl;
			}
		} // for c3 in cycles3
	} // for c2 in cycles2

	int looptime = time(NULL) - starttime;
	tr() << "FragCount: " << frag_count << std::endl;
	tr() << "Looptime " << looptime << std::endl;

	pose = mc_->lowest_score_pose();
	scorefxn()->show(  tr().Info , pose );
	tr().Info << "-------------------------" << std::endl;
	mc_->show_counters();

	// restore original CB wt if symmetric
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		scorefxn()->set_weight( chainbreak, 1.0 );
	}

	// Check chain break !
	if ( chainbreak_present ) {
		using namespace core::scoring;
		( *scorefxn() )( pose );
		core::Real chain_break_score = std::max(
			static_cast<float> (pose.energies().total_energies()[ chainbreak ]),
			static_cast<float> (pose.energies().total_energies()[ linear_chainbreak ])
		);

		core::Real chain_break_tol
			= option[ basic::options::OptionKeys::loops::chain_break_tol ]();
		tr().Info << "Chainbreak: " << chain_break_score << " Max: "
			<< chain_break_tol << std::endl;

		//remove cutpoints before all return statements because the code is now determining loop failure without cutpoints.
		loops::remove_cutpoint_variants( pose );
		pose.fold_tree( f_orig );

		if ( chain_break_score > (chain_break_tol*10) )	return ExtendFailure;   // if we have a really bad chain break, extend loop definitions
		if ( chain_break_score > chain_break_tol )      return Failure;         // if we only have a slight chainbreak problem, try again
	} // if ( chainbreak_present )

	// return to original fold tree
	loops::remove_cutpoint_variants( pose );
	pose.fold_tree( f_orig );
	return Success;
}

basic::Tracer & LoopMover_Perturb_QuickCCD::tr() const
{
    return TR;
}

LoopMover_Perturb_QuickCCDCreator::~LoopMover_Perturb_QuickCCDCreator() {}

moves::MoverOP LoopMover_Perturb_QuickCCDCreator::create_mover() const {
  return moves::MoverOP( new LoopMover_Perturb_QuickCCD() );
}

std::string LoopMover_Perturb_QuickCCDCreator::keyname() const {
  return "LoopMover_Perturb_QuickCCD";
}


//////////////////////////////////////////////////////////////////////////////////
/// @details  CCD close the loop [loop_begin,loop_end].
/// Wraps protocols::loops::loop_closure::ccd::CCDLoopClosureMover.apply() using most of its default options.
/// rama scores are not checked, however, and the secondary structure is "fixed" afterward.
/// @remark   This is a misnomer; it actually closes a single loop only. ~Labonte
void fast_ccd_close_loops(
	core::pose::Pose & pose,
	Loop const & loop,
	kinematics::MoveMap & mm
)
{
	loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover(
			loop, kinematics::MoveMapCOP( kinematics::MoveMapOP( new kinematics::MoveMap( mm ) ) ) );
	ccd_loop_closure_mover.check_rama_scores( false );
	ccd_loop_closure_mover.apply( pose );

	// fix secondary structure??
	for (core::Size i=loop.start(); i<=loop.stop(); ++i) {
		char ss_i = pose.conformation().secstruct( i );
		if ( ss_i != 'L' && ss_i != 'H' && ss_i != 'E')
			pose.set_secstruct( i , 'L' );
	}
}

} // namespace perturb
} // namespace loop_mover
} // namespace loops
} // namespace protocols
