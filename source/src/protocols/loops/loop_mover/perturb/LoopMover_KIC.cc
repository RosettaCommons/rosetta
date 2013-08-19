// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopMover_Perturb_KIC.cc
/// @brief kinematic loop closure main protocols
/// @author Chu Wang
/// @author Daniel J. Mandell
/// @author Mike Tyka
/// @author James Thompson
/// @author Roland A. Pache
/// @author Amelie Stein, amelie.stein@ucsf.edu, Oct 2012 -- next-generation KIC


#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KICCreator.hh>

//// Unit Headers
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverStatus.hh>
#include <core/conformation/Residue.hh>

//
//// Rosetta Headers
#include <core/chemical/VariantType.hh>

#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/optimization/AtomTreeMinimizer.hh>
//#include <core/optimization/CartesianMinimizer.hh> // AS Feb 6, 2013 -- actually we might not need this and can also remove the other two minimizer headers, if the MinMover can take care of it all?
//#include <core/optimization/MinimizerOptions.hh>
#include <protocols/simple_moves/MinMover.hh> // AS FEb 6, 2013 -- rerouting minimization to include cartmin
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

//#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
// AUTO-REMOVED #include <core/pack/rotamer_trials.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>

#include <basic/Tracer.hh>

//Utility Headers
#include <numeric/random/random.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>

//Auto Headers


//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


namespace protocols {
namespace loops {
namespace loop_mover {
namespace perturb {

///////////////////////////////////////////////////////////////////////////////
using namespace core;

static numeric::random::RandomGenerator RG(42444);
static basic::Tracer TR("protocols.loops.loop_mover.perturb.LoopMover_Perturb_KIC");

LoopMover_Perturb_KIC::LoopMover_Perturb_KIC() :
	IndependentLoopMover()
{
	set_scorefxn( get_cen_scorefxn() );
	loop_mover::loops_set_chainbreak_weight( scorefxn(), 1 );

	protocols::moves::Mover::type("LoopMover_Perturb_KIC");
	set_default_settings();
}


LoopMover_Perturb_KIC::LoopMover_Perturb_KIC(
	protocols::loops::LoopsOP  loops_in
) : IndependentLoopMover( loops_in )
{
	set_scorefxn( get_cen_scorefxn() );
	loop_mover::loops_set_chainbreak_weight( scorefxn(), 1 );

	protocols::moves::Mover::type("LoopMover_Perturb_KIC");
	set_default_settings();

}

LoopMover_Perturb_KIC::LoopMover_Perturb_KIC(
	protocols::loops::LoopsOP  loops_in,
	core::scoring::ScoreFunctionOP  scorefxn
) : IndependentLoopMover( loops_in )
{
	if( scorefxn ){
		set_scorefxn( scorefxn );
	}else{
		set_scorefxn( get_cen_scorefxn() );
 		loop_mover::loops_set_chainbreak_weight( scorefxn(), 1 );
	}
	protocols::moves::Mover::type("LoopMover_Perturb_KIC");
	set_default_settings();
}

//destructor
LoopMover_Perturb_KIC::~LoopMover_Perturb_KIC(){}

//clone
protocols::moves::MoverOP LoopMover_Perturb_KIC::clone() const {
	return new LoopMover_Perturb_KIC(*this);
}

void LoopMover_Perturb_KIC::set_default_settings()
{
	if ( basic::options::option[ basic::options::OptionKeys::loops::strict_loops ].user() ) {
		set_strict_loops( basic::options::option[ basic::options::OptionKeys::loops::strict_loops ]() );
	}
	else set_strict_loops(true); // obey loop definitions in kinematic mode
	max_seglen_ = basic::options::option[basic::options::OptionKeys::loops::kic_max_seglen];
	recover_low_ = ( ! basic::options::option[basic::options::OptionKeys::loops::kic_recover_last] );
	max_kic_build_attempts_ = basic::options::option[basic::options::OptionKeys::loops::max_kic_build_attempts];
	remodel_kic_attempts_ = basic::options::option[basic::options::OptionKeys::loops::remodel_kic_attempts];
}

void  LoopMover_Perturb_KIC::set_extended_torsions(
												   core::pose::Pose & ,
												   Loop const &
												   )
{
	// do nothing for now, overriding LoopMover::set_extended_torsions()
}


/// @detailed
/// Uses kinematic_mover to remodel a protein segment. If the 'extended' flag in the loop
/// definition for the segment is set to '1', will idealize all bond lengths, bond angles, and phi,
/// psi, and omega torsions before modeling. This stage is carried out entirely with a centroid
/// representation. Applies to only one loop, given as an argument.
loop_mover::LoopResult LoopMover_Perturb_KIC::model_loop(
	core::pose::Pose & pose,
  protocols::loops::Loop const & loop
){
	static int cur_struct=0; // for movie output
	// Dont allow loops < 3 residues.
	if( (loop.stop() - loop.start() < 2 )){
		tr().Error << "[WARNING] KinematicMover cannot handle loops smaller than 3 residues. Doing nothing. " << std::endl;
		return loop_mover::CriticalFailure;
	}

	// Objects representing one loop
	Loops one_loop_loops;
	one_loop_loops.add_loop( loop );

	using namespace scoring;
	using namespace optimization;
	using namespace basic::options;

	//bool const verbose( true );
	bool const local_debug( false ); 
	bool const local_movie( false );

	core::pose::Pose native_pose;
	if( get_native_pose() ){
		native_pose = *get_native_pose();
	}else{
		native_pose = pose;
	}

	Size const loop_begin( loop.start() ), loop_end( loop.stop() ), loop_cut( loop.cut() );
	Size const loop_size( loop_end - loop_begin + 1 );
	runtime_assert( loop.is_terminal( pose ) || pose.fold_tree().is_cutpoint( loop_cut ) );
	std::ofstream loop_outfile; // for movie

	tr() << "perturb_one_loop_with_KIC: " << loop_begin << ' ' << loop_size << std::endl;

	// set cutpoint variant for chainbreak scoring.
	core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, loop_cut );
	core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, loop_cut+1 );

	if (local_debug) {
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
	utility::vector1<bool> allow_sc_move_one_loop( pose.total_residue(), false );
	loops_set_move_map( one_loop_loops, allow_sc_move_one_loop, mm_one_loop);
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		core::pose::symmetry::make_symmetric_movemap( pose, mm_one_loop );
	}


	// scheduler
	bool const fast = option[OptionKeys::loops::fast];
	int outer_cycles( 3 );
	if ( option[ OptionKeys::loops::outer_cycles ].user() ) {
		outer_cycles = option[ OptionKeys::loops::outer_cycles ]();
	}
	if ( option[ OptionKeys::run::test_cycles ]() ) {
		outer_cycles = 3;
	}
	int inner_cycles( fast ? std::min( Size(250), loop_size*5 ) : std::min( Size(1000), loop_size*20 ) );
	if ( option[ OptionKeys::loops::max_inner_cycles ].user() ) {
		inner_cycles = option[ OptionKeys::loops::max_inner_cycles ]();
	}
	if ( option[ OptionKeys::run::test_cycles ]() ) {
		inner_cycles = 3;
	}

	// Monte Carlo vars
	float const init_temp( option[ OptionKeys::loops::remodel_init_temp ]() );
	float const	final_temp( option[ OptionKeys::loops::remodel_final_temp ]() );
	float const gamma = std::pow( (final_temp/init_temp), (1.0f/(outer_cycles*inner_cycles)) );
	float temperature = init_temp;
	if ( local_debug ) { // hacking
		( *scorefxn() )(pose);
		tr() << "before mc ctor: "
			<< pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) << std::endl;
	}

	// minimizer
	/*
	AtomTreeMinimizerOP minimizer;
	float const dummy_tol( 0.001 ); // linmin sets tol internally
	bool const use_nblist( false ), deriv_check( false ); // true ); // false );
	MinimizerOptions options( "linmin", dummy_tol, use_nblist, deriv_check);
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		minimizer = new core::optimization::symmetry::SymAtomTreeMinimizer;
	} else {
		minimizer = new core::optimization::AtomTreeMinimizer;
	}
	*/
	// AS Feb 6 2013: rewriting the minimizer section to use the MinMover, which allows seamless integration of cartesian minimization
	protocols::simple_moves::MinMoverOP min_mover;
	MoveMapOP mm_one_loop_OP = new core::kinematics::MoveMap( mm_one_loop ); // it appears we need a move map OP to do this... 

	const std::string min_type = "linmin";
	core::Real dummy_tol( 0.001 ); // linmin sets tol internally -- MinMover doesn't accept const input... 
	bool use_nblist( false ), deriv_check( false ), use_cartmin ( option[ OptionKeys::loops::kic_with_cartmin ]() ); // true ); // false );
	if ( use_cartmin ) runtime_assert( scorefxn()->get_weight( core::scoring::cart_bonded ) > 1e-3 ); // AS -- actually I'm not sure if this makes any sense in centroid... ask Frank?
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		min_mover = new simple_moves::symmetry::SymMinMover( mm_one_loop_OP, scorefxn(), min_type, dummy_tol, use_nblist, deriv_check );
		//min_mover = new simple_moves::symmetry::SymMinMover(); 
	} else {
		min_mover = new protocols::simple_moves::MinMover( mm_one_loop_OP, scorefxn(), min_type, dummy_tol, use_nblist, deriv_check );
		//min_mover = new protocols::simple_moves::MinMover(); // version above doesn't work, but setting all one by one does.. try again with scorefxn() w/o the star though
	}
	/*
	min_mover->movemap( mm_one_loop_OP );
	min_mover->score_function( scorefxn() );
	min_mover->min_type( min_type );
	min_mover->tolerance( dummy_tol );
	min_mover->nb_list( use_nblist );
	min_mover->deriv_check( deriv_check );
	 */
	min_mover->cartesian( use_cartmin );
	

	// show temps
	tr() << "remodel init temp: " << init_temp << std::endl;
	tr() << "remodel final temp: " << final_temp << std::endl;

	// perform the initial perturbation
	// setup the kinematic mover
	
	loop_closure::kinematic_closure::KinematicMover myKinematicMover;
	
	//tr() << "kinematic mover generated, now setting up perturber... " << std::endl;
	
	// AS: select from different perturbers implemented for NGK
	// torsion-restricted sampling
	if ( basic::options::option[ basic::options::OptionKeys::loops::restrict_kic_sampling_to_torsion_string ].user() || basic::options::option[ basic::options::OptionKeys::loops::derive_torsion_string_from_native_pose ]() ) {
		std::string torsion_bins = basic::options::option[ basic::options::OptionKeys::loops::restrict_kic_sampling_to_torsion_string ]();
		
		// derive torsion string from native/input pose, if requested -- warning: this overwrites the externally provided one
		if ( basic::options::option[ basic::options::OptionKeys::loops::derive_torsion_string_from_native_pose ]() )
			torsion_bins = torsion_features_string( native_pose ); // does this work at this point in the code? The loop needs to be set up!
		
		// check if the user-provided string has the correct length
		runtime_assert(torsion_bins.size() == loop_end - loop_begin + 1); // this is where torsion-restricted sampling with multiple loops currently fails -- however, proper access to the torsion bins in the perturber would also be implemented
		
		loop_closure::kinematic_closure::TorsionRestrictedKinematicPerturberOP perturber = new loop_closure::kinematic_closure::TorsionRestrictedKinematicPerturber ( &myKinematicMover, torsion_bins );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() ); // to make the code cleaner this should really be a generic function, even though some classes may not use it
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_rama2b ]() && basic::options::option[ basic::options::OptionKeys::loops::taboo_sampling ]() ) {  // TabooSampling with rama2b (neighbor-dependent phi/psi lookup)
		loop_closure::kinematic_closure::NeighborDependentTabooSamplingKinematicPerturberOP 
		perturber =
        new loop_closure::kinematic_closure::NeighborDependentTabooSamplingKinematicPerturber( &myKinematicMover );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::taboo_sampling ] ) { // TabooSampling (std rama)
		loop_closure::kinematic_closure::TabooSamplingKinematicPerturberOP perturber = new loop_closure::kinematic_closure::TabooSamplingKinematicPerturber ( &myKinematicMover );
		
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() ); 
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_rama2b ]() ) { // rama2b
		loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturberOP 
		perturber =
        new loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturber( &myKinematicMover );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );
	} else { // default behavior [for now] -- std KIC
		loop_closure::kinematic_closure::TorsionSamplingKinematicPerturberOP 
		perturber =
		new loop_closure::kinematic_closure::TorsionSamplingKinematicPerturber( &myKinematicMover );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );
	}
	
	//tr() << "perturber set!" << std::endl;
	
	myKinematicMover.set_vary_bondangles( true );
	//myKinematicMover.set_vary_bondangles( false );  // trying without varying angles
	
	myKinematicMover.set_sample_nonpivot_torsions( option[ OptionKeys::loops::nonpivot_torsion_sampling ]());
	myKinematicMover.set_rama_check( true );
	
	myKinematicMover.set_loop_begin_and_end( loop_begin, loop_end ); // AS -- for restricted torsion bin sampling, the mover needs to know about the start of the defined loop, not just the segment that is sampled in a given move
	
	
	// for Taboo Sampling we need to update the sequence here every time, in case it changes (e.g. when modeling multiple loops)
	// this setup would also allow us to combine taboo sampling with other types of sampling, e.g. to increase diversity after several rounds of standard/random KIC sampling
	utility::vector1< core::chemical::AA > loop_sequence;
	loop_sequence.resize(0); // make sure it is empty
	for (core::Size cur_res = loop_begin; cur_res <= loop_end; cur_res++) { 
		loop_sequence.push_back(pose.aa(cur_res)); 
	}
	myKinematicMover.update_sequence( loop_sequence ); // should only have an effect on the TabooSamplingKinematicPerturber
	
	
	
	Size kic_start, kic_middle, kic_end; // three pivot residues for kinematic loop closure
	kic_start = loop_begin;
	kic_end = loop_end;
	Size middle_offset = (kic_end - kic_start) / 2; // need to ensure this isn't a proline
	kic_middle = kic_start + middle_offset;
	tr() << "kinematic initial perturb with start_res: "  << kic_start << "  middle res: " << kic_middle << "  end_res: "
	<< kic_end << std::endl;
	myKinematicMover.set_pivots(kic_start, kic_middle, kic_end);
	myKinematicMover.set_temperature(temperature);
	
	std::string torsion_features = torsion_features_string( pose );
	tr() << "loop rmsd before initial kinematic perturbation:" << loop_rmsd( pose, native_pose, one_loop_loops ) /* << " -- torsion bins: " << torsion_features */ << std::endl;

	if (loop.is_extended() ) {
		myKinematicMover.set_idealize_loop_first( true ); // start without any native angles or lengths
		core::Size nits=0;

		core::Real previous_bump_overlap_factor=myKinematicMover.get_bump_overlap_factor();//RAP
		if ( !option[ OptionKeys::loops::kic_bump_overlap_factor ].user() ) { // AS: allow external changes to the bump overlap factor
			//RAP: temporarily lower the bump_overlap_factor to 0.4 to speed-up initial loop closure after centroid radii fix
			myKinematicMover.set_bump_overlap_factor(0.4);//RAP
		}
		while (nits < max_kic_build_attempts_) {
			tr() << "Attempting loop building: " << nits << " ... " << std::endl;
			myKinematicMover.apply( pose );
			if (myKinematicMover.last_move_succeeded()) {
				set_last_move_status(protocols::moves::MS_SUCCESS);
				tr() << "initial kinematic perturbation complete" << std::endl;
				myKinematicMover.set_idealize_loop_first( false ); // now the loop is idealized
				
				// AS Jan 4, 2013: adding option to leave centroid mode directly after successful initial closure
				if ( option[ OptionKeys::loops::kic_leave_centroid_after_initial_closure ]() ) {
					return loop_mover::Success;
				}
				
				break;
			}
			nits++;
		}
		if ( !option[ OptionKeys::loops::kic_bump_overlap_factor ].user() ) { 
			//RAP: restore bump_overlap_factor to original value
			myKinematicMover.set_bump_overlap_factor(previous_bump_overlap_factor);//RAP
		}
		if (!myKinematicMover.last_move_succeeded()) {
			tr().Error << "[WARNING] Failed to build loop with kinematic Mover during initial kinematic perturbation after " << nits << " trials: " << loop << std::endl;
			set_last_move_status(protocols::moves::FAIL_RETRY);
			//pose.fold_tree( f_orig ); // DJM: doing above in LoopRelaxMover now
			return loop_mover::CriticalFailure;
		}
		( *scorefxn() )(pose);
		//minimizer->run( pose, mm_one_loop, *scorefxn(), options );
		
		if ( local_debug ) // AS
			tr() << " minimization next " << std::endl;
		
		min_mover->score_function( *scorefxn() ); // at the moment this doesn't change, but if we ever implemented ramp it would
		min_mover->apply( pose ); // the other options were already registered in setup
		
		
		tr() << "loop rmsd after initial kinematic perturbation:" << loop_rmsd( pose, native_pose, one_loop_loops ) << std::endl;

	}
	else {
		tr() << "not performing initial kinematic perturbation" << std::endl;
		if (option[ OptionKeys::loops::vicinity_sampling ]()) {
			// AS Oct 3 2012: replace TorsionRestrictedKinematicPerturber by VicinitySamplingKinematicPerturber
			loop_closure::kinematic_closure::VicinitySamplingKinematicPerturberOP v_perturber =
			new loop_closure::kinematic_closure::VicinitySamplingKinematicPerturber( &myKinematicMover );
			v_perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
			v_perturber->set_degree_vicinity( option[ OptionKeys::loops::vicinity_degree ]() );
			myKinematicMover.set_perturber( v_perturber );
			
		}
	}

	if (local_movie) {
		std::string outname_base = option[ OptionKeys::loops::output_pdb ]().name();
		std::string outname_prefix = option[ OptionKeys::out::prefix ];
		std::string outname = outname_prefix + outname_base + "_centroid_movie_" +
			right_string_of(cur_struct,4,'0') + ".pdb";
		loop_outfile.open(outname.c_str(), std::ios::out | std::ios::binary);
		loop_outfile << "MODEL" << std::endl;
		utility::vector1<Size> indices(loop_end - loop_begin + 3);
		for (Size i=loop_begin-1, j=1; i<=loop_end+1; i++, j++) {
			indices[j]=i;
		}
		//pose.dump_pdb(loop_outfile, indices, "init_perturb");
		loop_outfile << "ENDMDL" << std::endl;
	}

	// Monte Carlo object
	protocols::moves::MonteCarlo mc( pose, *scorefxn(), temperature);
	mc.show_scores();

	for( int i=1; i<=outer_cycles; ++i ) {
		if ( local_debug) { // debug
			( *scorefxn() )( pose );
			tr() << "befor rLOW: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) <<
				" rmsd: " << F(9,3,loop_rmsd( pose, native_pose, one_loop_loops )) << std::endl;
		}

		// recover low
		if ( recover_low_ ) {
			mc.recover_low(pose);
		}

		if ( local_debug) { // debug
			( *scorefxn() )( pose );
			tr() << "after rLOW: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() ) <<
				" rmsd: " << F(9,3,loop_rmsd( pose, native_pose, one_loop_loops )) << std::endl;
		}

		for( int j=1; j<=inner_cycles; ++j ) {
			// change temperature
			temperature *= gamma;
			mc.set_temperature( temperature );
			core::Size nits=0;
			while (nits < remodel_kic_attempts_) {
				nits++;
				if (option[ OptionKeys::loops::always_remodel_full_loop ]()) { // warning: for long loops this could be extremely inefficient
					kic_start = loop_begin;
					kic_end = loop_end;
					kic_middle = RG.random_range(kic_start+2, kic_end-2);
				} else { // randomly selected sub-segments
					// AS Oct 2012: in the previous KIC implementation there was a directional bias here, as the start pivot is always selected first, which means that the C-terminal part of the loop is sampled more than the N-terminal part
					if ( option[ OptionKeys::loops::legacy_kic ]() || j % 2 == 0 ) {
						kic_start = RG.random_range(loop_begin,loop_end-2);
						// choose a random end residue so the length is >= 3, <= min(loop_end, start+maxlen)
						kic_end = RG.random_range(kic_start+2, std::min((kic_start+max_seglen_ - 1), loop_end));
						Size middle_offset = (kic_end - kic_start) / 2;
						kic_middle = kic_start + middle_offset;
					} else {
						//tr() << " -- selection starting from the end of the loop -- " << std::endl;
						kic_end = RG.random_range(loop_begin+2,loop_end);
						kic_start = RG.random_range(std::max((kic_end - std::min(max_seglen_, kic_end) + 1), loop_begin), kic_end-2);
						Size middle_offset = (kic_end - kic_start) / 2;
						kic_middle = kic_start + middle_offset;
					}
				}
				myKinematicMover.set_pivots(kic_start, kic_middle, kic_end);
				myKinematicMover.set_temperature(temperature);
				myKinematicMover.apply( pose );
				if (myKinematicMover.last_move_succeeded()) {
					break;
				}
			}
			if (myKinematicMover.last_move_succeeded()) {
				//fpd symmetrize 'mm_one_loop'
				if ( core::pose::symmetry::is_symmetric( pose ) )  {
					core::pose::symmetry::make_symmetric_movemap( pose, mm_one_loop );
				}
				( *scorefxn() )(pose);
				//minimizer->run( pose, mm_one_loop, *scorefxn(), options );
				
				if ( local_debug ) // AS
					tr() << " minimization next " << std::endl;

				
				min_mover->score_function( *scorefxn() );
				min_mover->apply( pose ); // the other options were already registered in setup

				std::string move_type = "kinematic_perturb";
				bool accepted = mc.boltzmann( pose, move_type );
				if (accepted) {
					tr() << "new centroid perturb rmsd: " << loop_rmsd( pose, native_pose, one_loop_loops ) << std::endl;
					if (local_movie) {
						loop_outfile << "MODEL" << std::endl;
						utility::vector1<Size> indices(loop_end - loop_begin + 3);
						for (Size i=loop_begin-1, j=1; i<=loop_end+1; i++, j++) {
							indices[j]=i;
						}
						pose.dump_pdb(loop_outfile, indices, "init_perturb");
						loop_outfile << "ENDMDL" << std::endl;
					}
					//tr << "chainbreak score: " << pose.energies().total_energies()[ core::scoring::chainbreak ] << std::endl;
				}
				//mc.show_scores();
			}else{
				tr().Error << "[WARNING] Failed to build loop with kinematic Mover after " << nits << " trials: " << loop << std::endl;
				// return to original fold tree
				//pose.fold_tree( f_orig ); // DJM: doing above in LoopRelaxMover now
				return loop_mover::Failure;
			}
		} // inner_cycles
	} // outer_cycles
	
	if ( recover_low_ ) {
		pose = mc.lowest_score_pose();
	}
	else {
		pose = mc.last_accepted_pose();
	}
	if (local_movie) {
		// this assumes there is only one loop.
		Size begin_loop=loop.start();
		Size end_loop=loop.stop();
		loop_outfile << "MODEL" << std::endl;
		utility::vector1<Size> indices(end_loop - begin_loop + 3);
		for (Size i=begin_loop-1, j=1; i<=end_loop+1; i++, j++) {
			indices[j]=i;
		}
		pose.dump_pdb(loop_outfile, indices, "final_perturb");
		loop_outfile << "ENDMDL" << std::endl;
	}
	if (local_debug) {
		std::ofstream out("score.tmp_perturb_cen");
		out << "scoring after cen_perturb: " << ( *scorefxn() )(pose) << std::endl;
		/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies(pose);
		scorefxn()->show( out );
		out << pose.energies();
	}

	// return to original fold tree
	//pose.fold_tree( f_orig ); // DJM: doing above in LoopRelaxMover now

	return loop_mover::Success;
}

std::string
LoopMover_Perturb_KIC::get_name() const {
	return "LoopMover_Perturb_KIC";
}

void
LoopMover_Perturb_KIC::show(std::ostream & output) const
{
	Mover::show(output);
	output <<   "Maximum KIC segment length:             " << max_seglen_ <<
				"\nRecover the lowest energy conformation: " << (recover_low_ ? "True" : "False") <<
				"\nMaximum KIC build attempts:             " << max_kic_build_attempts_ <<
				"\nRemodel KIC attempts:                   " << remodel_kic_attempts_;
}

std::ostream &operator<< ( std::ostream &os, LoopMover_Perturb_KIC const &mover )
{
	mover.show(os);
	return os;
}

basic::Tracer & LoopMover_Perturb_KIC::tr() const
{
    return TR;
}

LoopMover_Perturb_KICCreator::~LoopMover_Perturb_KICCreator() {}

moves::MoverOP LoopMover_Perturb_KICCreator::create_mover() const {
  return new LoopMover_Perturb_KIC();
}

std::string LoopMover_Perturb_KICCreator::keyname() const {
  return "LoopMover_Perturb_KIC";
}

} // namespace perturb
} // namespace loop_mover
} // namespace loops
} // namespace protocols
