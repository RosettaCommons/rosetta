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
/// @author Amelie Stein, amelie.stein@ucsf.edu, Oct 2012 -- next-generation KIC


#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KICCreator.hh>

//// Unit Headers
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverStatus.hh>
#include <core/conformation/Residue.hh>

// added by Rebecca Alford - fixes the build?
#include <core/conformation/ppo_torsion_bin.hh>


//// Rosetta Headers
#include <core/chemical/VariantType.hh>

#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/optimization/AtomTreeMinimizer.hh>
//#include <core/optimization/MinimizerOptions.hh>
#include <protocols/simple_moves/MinMover.hh> // AS FEb 6, 2013 -- rerouting minimization to include cartmin

#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
// #include <core/scoring/Energies.hh> // AS April 25, 2013 -- needed to get chainbreak energy
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/symmetry/util.hh>

//#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>

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
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <fstream>

//Auto Headers

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

///////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace core::conformation;

static thread_local basic::Tracer TR( "protocols.loops.loop_mover.refine.LoopMover_Refine_KIC" );

LoopMover_Refine_KIC::LoopMover_Refine_KIC() : LoopMover()
{
	init( get_fa_scorefxn() );
}
LoopMover_Refine_KIC::LoopMover_Refine_KIC(
	protocols::loops::LoopsOP  loops_in
) : LoopMover( loops_in )
{
	init( get_fa_scorefxn() );
}

LoopMover_Refine_KIC::LoopMover_Refine_KIC(
	protocols::loops::LoopsOP const loops_in,
	core::scoring::ScoreFunctionCOP  scorefxn
) : LoopMover( loops_in )
{
	init( scorefxn );
}

//destructor
LoopMover_Refine_KIC::~LoopMover_Refine_KIC(){}

//clone
protocols::moves::MoverOP LoopMover_Refine_KIC::clone() const {
	return protocols::moves::MoverOP( new LoopMover_Refine_KIC(*this) );
}

void LoopMover_Refine_KIC::init( core::scoring::ScoreFunctionCOP  scorefxn )
{
	set_scorefxn( scorefxn->clone() );
	protocols::moves::Mover::type("LoopMover_Refine_KIC");
	set_default_settings();
}

void
LoopMover_Refine_KIC::set_default_settings(){
	using namespace core;
	using namespace basic::options;

	fix_natsc_ = option[OptionKeys::loops::fix_natsc];
	redesign_loop_ = false;
	flank_residue_min_ = false; // by JQX
	neighbor_dist_ = option[ OptionKeys::loops::neighbor_dist ]();
	max_seglen_ = option[ OptionKeys::loops::kic_max_seglen ];
	recover_low_ = ( ! option[ OptionKeys::loops::kic_recover_last ] );
	min_after_repack_ = option[ OptionKeys::loops::kic_min_after_repack ];
	optimize_only_kic_region_sidechains_after_move_ =
		option[ OptionKeys::loops::optimize_only_kic_region_sidechains_after_move ];

}

void LoopMover_Refine_KIC::set_task_factory( core::pack::task::TaskFactoryOP value ){ task_factory = value; }
bool LoopMover_Refine_KIC::get_task_factory(){ return task_factory != 0; }


/// detailed
/// Refines the polypeptide segments in the LoopMover_Refine_KIC::loops_ class variable from their current conformations.
/// Uses Rosetta's all-atom respresentation and high-resolution scoring function (with an upweighted chain
/// break term). At the beginning of this stage, unless the flag -loops:fix_natsc has been set, all residues within
/// the neighbor distance of a loop (defined by -loops:neighbor_dist) are repacked and then subject to rotamer trials.
/// The backbones of all loop residues, and the side-chains of all loop residues and neighbors are then subject to
/// DFPmin. If -loops:fix_natsc is set, only the loop residues (and not the neighbors) will be subject to repacking,
/// rotamer trials, and minimization. Consequently, if this stage has been preceeded by LoopMover_Perturb_KIC:model_loop,
/// and the -loops:fix_natsc flag is omitted, the side-chains surrounding the loop will be optimized for the
/// perturbed loop conformation, rather than the starting conformation that preceded the call to model_loop.
void LoopMover_Refine_KIC::apply(
	core::pose::Pose & pose
){
	// AMW: cppcheck flags, but I'm not changing this scope because it's static
	static int cur_struct=0; // for movie output

	/// must be called once the Pose has become available.
	resolve_loop_indices( pose );

	using namespace core;
	using namespace optimization;
	using namespace scoring;
	using namespace basic::options;
	using namespace loop_closure::kinematic_closure;

	// Set the native pose if provided
	core::pose::Pose native_pose;
	if ( get_native_pose() ) {
		native_pose = *get_native_pose();
	} else {
		native_pose = pose;
	}

	// 'verbose' should be an option, or better yet Tracer's output filtering capabilities should be used instead
	bool const verbose( true );
	bool const local_debug( option[ OptionKeys::loops::debug ].user() );
	bool const local_movie( false );
	std::ofstream loop_outfile; // for movie

	// prepare for movie output if requested
	if ( local_movie ) {
		std::string outname_base = option[ OptionKeys::loops::output_pdb ]().name();
		std::string outname_prefix = option[ OptionKeys::out::prefix ];
		std::string outname = outname_prefix + outname_base + "_fullatom_movie_" +
			right_string_of(cur_struct,4,'0') + ".pdb";
		loop_outfile.open(outname.c_str(), std::ios::out | std::ios::binary);
	}

	// set cutpoint variants for correct chainbreak scoring
	Size const nres( pose.total_residue() );
	utility::vector1< bool > is_loop( nres, false );

	for ( Loops::const_iterator it=loops()->begin(), it_end=loops()->end();
			it != it_end; ++it ) {
		for ( Size i= it->start(); i<= it->stop(); ++i ) {
			is_loop[i] = true;
		}
		Size const loop_cut(it->cut());

		if ( loop_cut != nres ) { //c-terminal loop
			if ( ! pose.residue(loop_cut).has_variant_type(chemical::CUTPOINT_LOWER) ) {
				core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, loop_cut );
			}
			if ( ! pose.residue(loop_cut+1).has_variant_type(chemical::CUTPOINT_UPPER) ) {
				core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, loop_cut+1 );
			}
		}
	}

	// scheduler
	int const fast = option[OptionKeys::loops::fast];
	int outer_cycles = option[ OptionKeys::loops::refine_outer_cycles ]();
	if ( option[ OptionKeys::run::test_cycles ]() ) {
		outer_cycles = 3;
	}
	int max_inner_cycles( 200 );
	if ( option[ OptionKeys::loops::max_inner_cycles ].user() ) {
		max_inner_cycles = option[ OptionKeys::loops::max_inner_cycles ]();
	}
	if ( option[ OptionKeys::run::test_cycles ]() ) {
		max_inner_cycles = 3;
	}

	int inner_cycles = std::min( Size(max_inner_cycles), fast ? (int)loops()->loop_size() : 10 * (Size)( loops()->loop_size() ) );
	int repack_period = 20; // should be an option
	if ( option[ OptionKeys::loops::repack_period ].user() ) {
		repack_period = option[ OptionKeys::loops::repack_period ]();
	}

	// scorefxn
	scoring::ScoreFunctionOP local_scorefxn;
	scoring::ScoreFunctionOP min_scorefxn;
	if ( scorefxn() != 0 ) local_scorefxn = scorefxn()->clone();
	else {
		local_scorefxn = get_fa_scorefxn();
	}

	// For testing JK's new solvation term. The exact form isn't yet differentiable
	min_scorefxn = local_scorefxn->clone();
	if ( min_scorefxn->get_weight( core::scoring::occ_sol_exact ) > 0.0 ) {
		min_scorefxn->set_weight( core::scoring::fa_sol, 0.65 );
		min_scorefxn->set_weight( core::scoring::occ_sol_exact, 0.0 );
	}

	//scoring::ScoreFunctionOP scorefxn( get_score_function() );
	//scorefxn->set_weight( chainbreak, 1.0 ); // confirm that chainbreak weight is set

	// scorefxn->set_weight( chainbreak, 1.0*10.0/3.0 ); // confirm that chainbreak weight is set
	// min_scorefxn->set_weight( chainbreak, 1.0*10.0/3.0 ); // confirm that chainbreak weight is set
	loop_mover::loops_set_chainbreak_weight( scorefxn(), 1 );
	loop_mover::loops_set_chainbreak_weight( min_scorefxn, 1 );
	loop_mover::loops_set_chainbreak_weight( local_scorefxn, 1 ); // AS April 25, 2013 -- the local_scorefxn didn't have a chainbreak weight, restoring that

	//scorefxn->set_weight(core::scoring::mm_bend, 1.0);

	// AS: add rama2b weight in case we're sampling with rama2b
	if ( option[ OptionKeys::loops::kic_rama2b ]() ) {
		min_scorefxn->set_weight( rama2b, min_scorefxn->get_weight( rama ) );
		min_scorefxn->set_weight( rama, 0);

		local_scorefxn->set_weight( rama2b, local_scorefxn->get_weight( rama ) );
		local_scorefxn->set_weight( rama, 0);
	}


	// monte carlo
	float const init_temp( option[ OptionKeys::loops::refine_init_temp ]() );
	float const final_temp( option[ OptionKeys::loops::refine_final_temp ]() );
	float const gamma = std::pow( (final_temp/init_temp), 1.0f/(outer_cycles*inner_cycles) );
	float temperature = init_temp;
	protocols::moves::MonteCarlo mc( pose, *local_scorefxn, temperature );

	// AS Feb 6 2013: rewriting the minimizer section to use the MinMover, which allows seamless integration of cartesian minimization
	protocols::simple_moves::MinMoverOP min_mover;

	const std::string min_type = "dfpmin";
	core::Real dummy_tol( 0.001 );
	bool use_nblist( true ), deriv_check( false ), use_cartmin ( option[ OptionKeys::loops::kic_with_cartmin ]() ); // true ); // false );
	if ( use_cartmin ) runtime_assert( min_scorefxn->get_weight( core::scoring::cart_bonded ) > 1e-3 );
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		min_mover = protocols::simple_moves::MinMoverOP( new simple_moves::symmetry::SymMinMover() );
	} else {
		min_mover = protocols::simple_moves::MinMoverOP( new protocols::simple_moves::MinMover() );
	}
	min_mover->score_function( min_scorefxn );
	min_mover->min_type( min_type );
	min_mover->tolerance( dummy_tol );
	min_mover->nb_list( use_nblist );
	min_mover->deriv_check( deriv_check );
	min_mover->cartesian( use_cartmin );


	// show temps
	tr() << "refine init temp: " << init_temp << std::endl;
	tr() << "refine final temp: " << final_temp << std::endl;

	// Set up the packer tasks: one for rotamer trials, one for repacking (with design if resfile supplied)
	using namespace pack::task;
	if ( task_factory == 0 ) {
		task_factory = core::pack::task::TaskFactoryOP( new TaskFactory );
		// TaskOperations replace the following kind of code:
		// base_packer_task->initialize_from_command_line().or_include_current( true );
		task_factory->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
		task_factory->push_back( operation::TaskOperationCOP( new operation::IncludeCurrent ) );
		if ( option[ OptionKeys::packing::resfile ].user() ) {
			// Note - resfile is obeyed, so use NATAA as default to maintain protocol behavior
			task_factory->push_back( operation::TaskOperationCOP( new core::pack::task::operation::ReadResfile ) );
			redesign_loop_ = true;
		}
	}

	PackerTaskOP repack_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
	repack_packer_task->set_bump_check( true );
	pack::task::PackerTaskOP rottrials_packer_task;
	if ( !option[ OptionKeys::packing::resfile ].user() && !redesign_loop_ ) {
		// Not designing -- just repack
		repack_packer_task->restrict_to_repacking();
		tr() << "Not designing" << std::endl;
	} else tr() << "Activating design" << std::endl;

	// setting redes loop, but no resfile specified. all non-loop positions only repack. loop positions can design.
	if ( redesign_loop_ && !option[ OptionKeys::packing::resfile ].user() ) {
		tr() << "Auto-setting loop design for residues:";
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( !is_loop[i] ) repack_packer_task->nonconst_residue_task( i ).restrict_to_repacking();
			else tr() << " " << i;
		}
		tr() << std::endl;
	}

	// setup kinematic mover
	//protocols::loops::kinematic_closure::KinematicMover myKinematicMover( init_temp );
	KinematicMoverOP myKinematicMover_op( new loop_closure::kinematic_closure::KinematicMover );
	KinematicMoverCAP myKinematicMover_cap( myKinematicMover_op );
	KinematicMover & myKinematicMover = *myKinematicMover_op;

	// AS Oct 3, 2012: create appropriate perturber here, depending on input flags
	if ( option[ OptionKeys::loops::vicinity_sampling ]() ) {
		loop_closure::kinematic_closure::VicinitySamplingKinematicPerturberOP perturber( new loop_closure::kinematic_closure::VicinitySamplingKinematicPerturber( myKinematicMover_cap ) );
		perturber->set_vary_ca_bond_angles(  ! option[OptionKeys::loops::fix_ca_bond_angles ]()  );
		myKinematicMover.set_perturber( perturber );
		perturber->set_degree_vicinity( option[ OptionKeys::loops::vicinity_degree ]() );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::restrict_kic_sampling_to_torsion_string ].user() ||
			basic::options::option[ basic::options::OptionKeys::loops::derive_torsion_string_from_native_pose ]() ) { // torsion-restricted sampling

		core::conformation::torsion_bin_string torsion_bins;
		try {
			std::string user_torbin_string = basic::options::option[ basic::options::OptionKeys::loops::restrict_kic_sampling_to_torsion_string ]();
			if ( ! user_torbin_string.empty() ) {
				core::conformation::map_string_to_torsion_bin_string( user_torbin_string );
			}
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			throw utility::excn::EXCN_Msg_Exception( "Error in converting the string given in the option loops::restrict_kic_sampling_to_torsion_string\n"
				"into a valid torsion set of torsion bins.  Check your inputs\n" + e.msg() );
		}

// derive torsion string from native/input pose, if requested -- warning: this overwrites the externally provided one
		if ( basic::options::option[ basic::options::OptionKeys::loops::derive_torsion_string_from_native_pose ]() ) {
			torsion_bins = torsion_features_string( native_pose );
		}
		// loop boundaries are not available here -- maybe the perturber should check the length?
		//runtime_assert(torsion_bins.size() != loop_end - loop_begin + 1);

		loop_closure::kinematic_closure::TorsionRestrictedKinematicPerturberOP perturber( new loop_closure::kinematic_closure::TorsionRestrictedKinematicPerturber ( myKinematicMover_cap, torsion_bins ) );
		// to make the code cleaner this should really be a generic function, even though some classes may not use it
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_rama2b ]() &&
			basic::options::option[ basic::options::OptionKeys::loops::taboo_in_fa ]() ) {
		// TabooSampling with rama2b (neighbor-dependent phi/psi lookup) -- note that Taboo Sampling will
		// only be active during the first half of the full-atom stage, after that we follow the energy
		// landscape / Monte Carlo
		loop_closure::kinematic_closure::NeighborDependentTabooSamplingKinematicPerturberOP perturber( new loop_closure::kinematic_closure::NeighborDependentTabooSamplingKinematicPerturber( myKinematicMover_cap ) );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::taboo_in_fa ] ) { // TabooSampling (std rama) -- note that Taboo Sampling will only be active during the first half of the full-atom stage, after that we follow the energy landscape / Monte Carlo
		loop_closure::kinematic_closure::TabooSamplingKinematicPerturberOP perturber( new loop_closure::kinematic_closure::TabooSamplingKinematicPerturber ( myKinematicMover_cap ) );

		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_rama2b ]() ) { // neighbor-dependent phi/psi selection
		loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturberOP
			perturber( new loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturber( myKinematicMover_cap ) );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );
	} else { // standard KIC
		loop_closure::kinematic_closure::TorsionSamplingKinematicPerturberOP perturber( new loop_closure::kinematic_closure::TorsionSamplingKinematicPerturber( myKinematicMover_cap ) );
		perturber->set_vary_ca_bond_angles(  ! option[OptionKeys::loops::fix_ca_bond_angles ]()  );
		myKinematicMover.set_perturber( perturber );
	}

	myKinematicMover.set_vary_bondangles( true ); // why is this hard-coded?
	myKinematicMover.set_sample_nonpivot_torsions( option[ OptionKeys::loops::nonpivot_torsion_sampling ]());
	myKinematicMover.set_rama_check( true );
	Size kic_start, kic_middle, kic_end; // three pivot residues for kinematic loop closure

	// AS: adding option to change the number of rotamer trials -- just one (as in the current implementation) may be way too little
	core::Size num_rot_trials = Size( option[ OptionKeys::loops::kic_num_rotamer_trials ]() );

	// AS: setting up weights for ramping rama[2b] and/or fa_rep
	core::Real orig_local_fa_rep_weight = local_scorefxn->get_weight( fa_rep );
	core::Real orig_min_fa_rep_weight = min_scorefxn->get_weight( fa_rep );

	core::Real orig_local_rama_weight = local_scorefxn->get_weight( rama );
	core::Real orig_min_rama_weight = min_scorefxn->get_weight( rama );

	core::Real orig_local_rama2b_weight = local_scorefxn->get_weight( rama2b );
	core::Real orig_min_rama2b_weight = min_scorefxn->get_weight( rama2b );

	if ( basic::options::option [ basic::options::OptionKeys::loops::ramp_fa_rep ]() ) {
		local_scorefxn->set_weight( fa_rep, orig_local_fa_rep_weight/(outer_cycles + 1) );
		min_scorefxn->set_weight( fa_rep, orig_min_fa_rep_weight/(outer_cycles + 1) );
	}
	//AS: ramp rama -- currently only in refine mode
	if ( basic::options::option [ basic::options::OptionKeys::loops::ramp_rama ]() ) {
		local_scorefxn->set_weight( rama, orig_local_rama_weight/(outer_cycles + 1) );
		min_scorefxn->set_weight( rama, orig_min_rama_weight/(outer_cycles + 1) );

		// ramp rama2b instead if we're using that for sampling
		if ( option[ OptionKeys::loops::kic_rama2b ]() ) {
			//TR << "ramp repulsive: rama2b set to " << orig_local_rama2b_weight/(outer_cycles + 1) << std::endl;
			local_scorefxn->set_weight( rama2b, orig_local_rama2b_weight/(outer_cycles + 1) );
			min_scorefxn->set_weight( rama2b, orig_min_rama2b_weight/(outer_cycles + 1) );
		}
	}


	// perform initial repack trial
	utility::vector1<bool> allow_sc_move_all_loops( nres, false );
	(*local_scorefxn)(pose); // update 10A nbr graph, silly way to do this
	// here we'll optimize all side-chains within neighbor_dist_ of any loop (unless fix_natsc_ is true)
	select_loop_residues( pose, *loops(), !fix_natsc_, allow_sc_move_all_loops, neighbor_dist_);
	core::pose::symmetry::make_residue_mask_symmetric( pose, allow_sc_move_all_loops );
	repack_packer_task->restrict_to_residues( allow_sc_move_all_loops );
	core::pack::pack_rotamers( pose, *local_scorefxn, repack_packer_task );

	// setup rottrials packer task here to ensure we're using current sequence if design was active
	rottrials_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
	rottrials_packer_task->restrict_to_repacking();
	rottrials_packer_task->set_bump_check( true );
	rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );

	// minimize after initial repack
	std::string move_type = "repack";
	pose.update_residue_neighbors(); // to update 10A nbr graph
	kinematics::MoveMap mm_all_loops; // DJM tmp
	loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_);
	if ( flank_residue_min_ ) { add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
	MoveMapOP mm_all_loops_OP( new kinematics::MoveMap( mm_all_loops ) ); // AS, required for using MinMover
	min_mover->movemap( mm_all_loops_OP );
	min_mover->score_function( min_scorefxn ); // AS: needs to be adapted in case we ramp any weights (or does this happen automatically through pointers?)
	// minimizer->run( pose, mm_all_loops, *min_scorefxn, options ); // DJM tmp
	min_mover->apply( pose );
	mc.boltzmann( pose, move_type );
	mc.show_scores();
	if ( redesign_loop_ ) {
		tr() << "Sequence after design step: "
			<< pose.sequence() << std::endl;
	}

	// Get a vector of move_maps and allow_sc_move vectors, one for each loop. This way if we have multiple loops
	// we can optimize only the side-chains around the loop selected in the inner_cycle
	utility::vector1< kinematics::MoveMap > move_maps ( loops()->size() );
	utility::vector1< utility::vector1< bool > > allow_sc_vectors ( loops()->size() );
	update_movemap_vectors( pose, move_maps );
	update_allow_sc_vectors( pose, allow_sc_vectors );

	if ( local_movie ) {
		// this assumes there is only one loops_.
		Loops::const_iterator one_loop( loops()->begin() );
		Size begin_loop=one_loop->start();
		Size end_loop=one_loop->stop();
		loop_outfile << "MODEL" << std::endl;
		utility::vector1<Size> indices(end_loop - begin_loop + 3);
		for ( Size i=begin_loop-1, j=1; i<=end_loop+1; i++, j++ ) {
			indices[j]=i;
		}
		pose.dump_pdb(loop_outfile, indices, "initial_repack");
		loop_outfile << "ENDMDL" << std::endl;
	}

	//outer_cycles = 1;
	//inner_cycles = 20;
	core::Size num_kinematic_trials = 2;

	for ( int i=1; i<=outer_cycles; ++i ) {
		// increase CHAINBREAK weight and update monte carlo
		//scorefxn->set_weight( chainbreak, float(i)*10.0/3.0 );
		//min_scorefxn->set_weight( chainbreak, float(i)*10.0/3.0 );
		loop_mover::loops_set_chainbreak_weight( scorefxn(), i );
		loop_mover::loops_set_chainbreak_weight( min_scorefxn, i );
		loop_mover::loops_set_chainbreak_weight( local_scorefxn, i ); // AS April 25, 2013
		//scorefxn->set_weight( chainbreak, float(i) );

		// AS: try ramping the fa_rep over the outer cycles
		if ( option [ OptionKeys::loops::ramp_fa_rep ]() ) {
			//TR << "ramp repulsive: fa_rep set to " << orig_local_fa_rep_weight/(outer_cycles - i + 1) << std::endl;
			local_scorefxn->set_weight( fa_rep, orig_local_fa_rep_weight/(outer_cycles - i + 1) );
			min_scorefxn->set_weight( fa_rep, orig_min_fa_rep_weight/(outer_cycles - i + 1) );
		}

		if ( option [ OptionKeys::loops::ramp_rama ]() ) {
			//TR << "ramp repulsive: rama set to " << orig_local_rama_weight/(outer_cycles - i + 1) << std::endl;
			local_scorefxn->set_weight( rama, orig_local_rama_weight/(outer_cycles - i + 1) );
			min_scorefxn->set_weight( rama, orig_min_rama_weight/(outer_cycles - i + 1) );

			if ( option[ OptionKeys::loops::kic_rama2b ]() ) {
				//TR << "ramp repulsive: rama2b set to " << orig_local_rama2b_weight/(outer_cycles - i + 1) << std::endl;
				local_scorefxn->set_weight( rama2b, orig_local_rama2b_weight/(outer_cycles - i + 1) );
				min_scorefxn->set_weight( rama2b, orig_min_rama2b_weight/(outer_cycles - i + 1) );

			}
		}


		// more than 50% done? -> switch off Taboo Sampling (if applicable)
		if ( i > outer_cycles/2 && option[ OptionKeys::loops::taboo_in_fa ]() ) {
			// if rama2b is active, generate a new rama2b perturber, otherwise use a standard torsion sampling perturber
			if ( basic::options::option[ basic::options::OptionKeys::loops::kic_rama2b ]() ) { // neighbor-dependent phi/psi selection
				loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturberOP
					perturber( new loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturber( myKinematicMover_cap ) );
				perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
				myKinematicMover.set_perturber( perturber );
			} else { // standard KIC
				loop_closure::kinematic_closure::TorsionSamplingKinematicPerturberOP perturber( new loop_closure::kinematic_closure::TorsionSamplingKinematicPerturber( myKinematicMover_cap ) );
				perturber->set_vary_ca_bond_angles(  ! option[OptionKeys::loops::fix_ca_bond_angles ]()  );
				myKinematicMover.set_perturber( perturber );
			}
		}


		mc.score_function( *local_scorefxn );
		// recover low
		if ( recover_low_ ) {
			mc.recover_low( pose );
		}
		// score info
		if ( verbose ) tr() << "cycle: " << i << "  " << (*local_scorefxn)(pose) << std::endl;
		//pose.energies().show( tr );
		for ( int j=1; j<=inner_cycles; ++j ) {
			temperature *= gamma;
			mc.set_temperature( temperature );
			if ( verbose ) {
				tr() << "refinement cycle (outer/inner): "
					<< i << "/" << outer_cycles << " "
					<< j << "/" << inner_cycles << " "
					<< std::endl;
			}

			// choose a random loop for both rounds
			//Loops::const_iterator it( loops_.one_random_loop() );
			Size loop_ind = numeric::random::rg().random_range(1, loops()->size());
			Loops one_loop;
			//one_loop.add_loop( it );
			one_loop.add_loop( ( *loops() )[ loop_ind ] );
			// get loop endpoints
			Size begin_loop=one_loop.begin()->start();
			Size end_loop=one_loop.begin()->stop();

			myKinematicMover.set_loop_begin_and_end( begin_loop, end_loop ); // AS -- for restricted torsion bin sampling, the mover needs to know about the start of the defined loop, not just the segment that is sampled in a given move

			// for Taboo Sampling we need to update the sequence here every time, in case it changes (e.g. when modeling multiple loops)
			utility::vector1< core::chemical::AA > loop_sequence;
			loop_sequence.resize(0); // make sure it is empty
			for ( core::Size cur_res = begin_loop; cur_res <= end_loop; cur_res++ ) {
				loop_sequence.push_back(pose.aa(cur_res));
			}
			myKinematicMover.update_sequence( loop_sequence );

			// get the movemap for minimization and allowed side-chains for rottrials for this loop
			kinematics::MoveMap cur_mm = move_maps[ loop_ind ];
			MoveMapOP cur_mm_OP( new kinematics::MoveMap( cur_mm ) );
			utility::vector1<bool> cur_allow_sc_move = allow_sc_vectors[ loop_ind ];
			rottrials_packer_task->restrict_to_residues( cur_allow_sc_move );

			// AS Thu Oct 25 19:41:14 PDT 2012
			// rewriting the kinematic trials so that the 1st and 2nd round are handled by a loop instead of code duplication
			// then incorporating the possibility to repack after each loop move (and before calling mc.boltzmann)

			for ( core::Size kinematic_trial = 1; kinematic_trial <= num_kinematic_trials; kinematic_trial++ ) {

				// AS: the previous implementation had a "history bias" towards the N-terminus of the loop, as the start pivot can be anywhere between begin_loop and end_loop-2, while the choice of the end pivot depends on the start pivot
				if ( option[ OptionKeys::loops::legacy_kic ]() || j % 2 == 0 ) {
					kic_start = numeric::random::rg().random_range(begin_loop,end_loop-2);
					// choose a random end residue so the length is >= 3, <= min(loop_end, start+maxlen)
					kic_end = numeric::random::rg().random_range(kic_start+2, std::min((kic_start+max_seglen_ - 1), end_loop));
					Size middle_offset = (kic_end - kic_start) / 2;
					kic_middle = kic_start + middle_offset;
				} else {
					kic_end = numeric::random::rg().random_range(begin_loop+2,end_loop);
					kic_start = numeric::random::rg().random_range(std::max((kic_end - std::min(max_seglen_, kic_end) + 1), begin_loop), kic_end-2);
					Size middle_offset = (kic_end - kic_start) / 2;
					kic_middle = kic_start + middle_offset;
				}

				myKinematicMover.set_pivots(kic_start, kic_middle, kic_end);
				myKinematicMover.set_temperature(temperature);

				core::pose::Pose last_accepted_pose = pose; // backup in case of undetected chain breaks
				myKinematicMover.apply( pose );

				// AS April 25, 2013
				// There seems to be a bug in KIC that makes it occasionally report open conformations but stating they are closed. We're looking into why exactly this happens -- for the moment this is a workaround that discards structures with a chainbreak score > 0, as this indicates that the chain probably isn't closed after all

				// note that overall scores here may be very high, due to side chain clashes (rotamer trials and/or repacking below)

				/* // moved to KinematicMover.cc
				(*local_scorefxn)(pose);
				core::Real current_chainbreak = pose.energies().total_energies()[ core::scoring::chainbreak ] * local_scorefxn->get_weight( chainbreak );

				//getPoseExtraScore(pose, "chainbreak", current_chainbreak);

				tr() << "AS_DEBUG -- chainbreak score " << current_chainbreak << std::endl;

				if ( current_chainbreak > option[ OptionKeys::loops::kic_chainbreak_threshold ]() ) { // todo: make this option-dependent -- AS_DEBUG
				tr() << "AS_DEBUG -- this structure probably has a chain break, skipping it " << kic_start << "-" << kic_end << std::endl;
				local_scorefxn->show(pose); // AS_DEBUG


				pose.dump_pdb("high_chainbreak_score_" + utility::to_string(i) + "_" + utility::to_string(j) + "_" + utility::to_string(kic_start) + "_" + utility::to_string(kic_end) + ".pdb");

				pose = last_accepted_pose;

				} else
				*/

				if ( myKinematicMover.last_move_succeeded() ) {
					// get the movemap and allowed side-chains for the current loop


					//fpd symmetrize mm_all_loops
					if ( core::pose::symmetry::is_symmetric( pose ) )  {
						core::pose::symmetry::make_residue_mask_symmetric( pose, cur_allow_sc_move );
					}

					// do optimizations
					if ( optimize_only_kic_region_sidechains_after_move_ ) {
						set_rottrials_from_kic_segment( pose, rottrials_packer_task, kic_start, kic_end ); // for rottrials
						set_movemap_from_kic_segment( pose, cur_mm, kic_start, kic_end ); // for minimization
					}

					// AS Oct 25, 2012: considering the order of magnitude that full-scale KIC modeling can introduce, it may be very valuable to repack before testing for acceptance
					// * note that this probably is much less relevant in vicinity sampling -- Rotamer Trials may be enough
					// * also note that this allows for design, which we're not testing at the moment in the standard loop modeling benchmark -- it should be checked whether this leads to collapsing structures / strong preference for small residues, and in the latter case, perhaps this part of the code should only repack, but not design
					if ( !option[ OptionKeys::loops::legacy_kic ]() && (j%repack_period)==0 ) {
						// implementation copied from "main_repack_trial" below

						//tr() << " AS -- repacking after loop closure but before mc.boltzmann " << std::endl; // AS_DEBUG

						update_movemap_vectors( pose, move_maps );
						update_allow_sc_vectors( pose, allow_sc_vectors );

						// the repack/design and subsequent minimization within this main_repack_trial block apply
						// to all loops (and their neighbors, if requested)
						loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_);
						if ( flank_residue_min_ ) { add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
						select_loop_residues( pose, *loops(), !fix_natsc_, allow_sc_move_all_loops, neighbor_dist_);

						core::pose::symmetry::make_residue_mask_symmetric( pose, allow_sc_move_all_loops );  //fpd symmetrize res mask -- does nothing if pose is not symm
						repack_packer_task->restrict_to_residues( allow_sc_move_all_loops );
						rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );
						pack::pack_rotamers( pose, *local_scorefxn, repack_packer_task ); // design here if resfile supplied
						if ( verbose ) tr() << "energy after design: " << (*local_scorefxn)(pose) << std::endl; // DJM remove
						if ( redesign_loop_ ) { // need to make new rottrials packer task with redesigned sequence
							rottrials_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
							rottrials_packer_task->restrict_to_repacking();
							rottrials_packer_task->set_bump_check( true );
							rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );
							if ( verbose ) tr() << "energy after design repack: " << (*local_scorefxn)(pose) << std::endl; // DJM remove
						}

						// minimize after repack if requested
						// AS: note that this may not be necessary here, as we'll minimize after the rotamer trials step anyway
						if ( min_after_repack_ ) {
							if ( core::pose::symmetry::is_symmetric( pose ) )  {
								//fpd  minimizing with the reduced movemap seems to cause occasional problems
								//     in the symmetric case ... am looking into this
								loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_ );
								if ( flank_residue_min_ ) { add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
							}
							min_mover->movemap( mm_all_loops_OP ); // will symmetry changes automatically be transferred through the pointer?
							min_mover->score_function( min_scorefxn );
							min_mover->apply( pose );
						}

						std::string move_type;
						if ( redesign_loop_ ) move_type = "repack+design";
						else move_type = "repack";
						mc.boltzmann( pose, move_type );
						mc.show_scores();
						if ( redesign_loop_ ) {
							tr() << "Sequence after design step: "
								<< pose.sequence() << std::endl;
						}
						if ( verbose ) tr() << "energy after repack: " << (*local_scorefxn)(pose) << std::endl;
					}

					// AS: first do repacking (if applicable), then rotamer trials -- we've shown that this increases recovery of native side chain conformations
					for ( Size i = 1; i <= num_rot_trials; i++ ) {
						pack::rotamer_trials( pose, *local_scorefxn, rottrials_packer_task );
						pose.update_residue_neighbors(); // to update 10A nbr graph
					}
					//pack::pack_rotamers_loop( pose, *local_scorefxn, rottrials_packer_task, num_rot_trials ); // can we use this instead? is this really RT, or actual packing?


					// AS: for non-legacy-KIC one might consider putting this into an ELSE branch -- however, min_after_repack_ defaults to false, in which case we wouldn't have minimization here, which could strongly affect the likelihood of acceptance
					min_mover->movemap( cur_mm_OP );
					min_mover->score_function( min_scorefxn );
					min_mover->apply( pose );

					// test for acceptance
					std::stringstream k_trial;
					k_trial << kinematic_trial;
					std::string move_type = "kic_refine_r" + k_trial.str();
					bool accepted = mc.boltzmann( pose, move_type );

					if ( accepted ) {
						tr() << "RMS to native after accepted kinematic round " << kinematic_trial << " move on loop "
							<< loops()->size() + 1 - loop_ind << ": " // reverse the order so it corresponds with the loop file
							<< loop_rmsd(pose, native_pose, *loops()) << std::endl;
						//tr << "after accepted move res " << kic_start + 2 << " omega is " << pose.omega(kic_start+2) << std::endl;
						if ( local_movie ) {
							loop_outfile << "MODEL" << std::endl;
							utility::vector1<Size> indices(end_loop - begin_loop + 3);
							for ( Size i=begin_loop-1, j=1; i<=end_loop+1; i++, j++ ) {
								indices[j]=i;
							}
							pose.dump_pdb(loop_outfile, indices, "refine_r" + k_trial.str());
							loop_outfile << "ENDMDL" << std::endl;
						}
						//tr << "temperature: " << temperature << std::endl;
						//tr << "chainbreak: " << pose.energies().total_energies()[ core::scoring::chainbreak ] << std::endl;
						if ( verbose ) tr() << "energy after accepted move: " << (*local_scorefxn)(pose) << std::endl;
					}

					//mc.show_scores();
				}
			}

			{ //main_repack_trial
				if ( ( option[ OptionKeys::loops::legacy_kic ]() && (j%repack_period)==0 ) || j==inner_cycles ) { // AS: always repack at the end of an outer cycle, to make sure the loop environment is repacked before returning from this function

					//tr() << " AS -- repacking after loop closure and mc.boltzmann (legacy KIC behavior) " << std::endl; // AS_DEBUG

					// repack trial
					// DJM: we have found that updating the move maps and rottrials/repack sets once per repack period
					//      gives better performance than updating them after every move, so we do so here for
					//  subsequent cycles until the next repack period
					update_movemap_vectors( pose, move_maps );
					update_allow_sc_vectors( pose, allow_sc_vectors );

					// the repack/design and subsequent minimization within this main_repack_trial block apply
					// to all loops (and their neighbors, if requested)
					loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_);
					if ( flank_residue_min_ ) { add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
					select_loop_residues( pose, *loops(), !fix_natsc_, allow_sc_move_all_loops, neighbor_dist_);

					core::pose::symmetry::make_residue_mask_symmetric( pose, allow_sc_move_all_loops );  //fpd symmetrize res mask -- does nothing if pose is not symm
					repack_packer_task->restrict_to_residues( allow_sc_move_all_loops );
					rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );

					pack::pack_rotamers( pose, *local_scorefxn, repack_packer_task ); // design here if resfile supplied
					if ( verbose ) tr() << "energy after design: " << (*local_scorefxn)(pose) << std::endl; // DJM remove

					if ( redesign_loop_ ) { // need to make new rottrials packer task with redesigned sequence
						rottrials_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
						rottrials_packer_task->restrict_to_repacking();
						rottrials_packer_task->set_bump_check( true );
						rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );
						if ( verbose ) tr() << "energy after design repack: " << (*local_scorefxn)(pose) << std::endl; // DJM remove
					}

					// minimize after repack if requested
					if ( min_after_repack_ ) {
						if ( core::pose::symmetry::is_symmetric( pose ) )  {
							//fpd  minimizing with the reduced movemap seems to cause occasional problems
							//     in the symmetric case ... am looking into this
							loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_ );
							if ( flank_residue_min_ ) { add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
						}
						min_mover->movemap( mm_all_loops_OP );
						min_mover->score_function( min_scorefxn );
						min_mover->apply( pose );
					}

					std::string move_type;
					if ( redesign_loop_ ) move_type = "repack+design";
					else move_type = "repack";
					mc.boltzmann( pose, move_type );
					mc.show_scores();
					if ( redesign_loop_ ) {
						tr() << "Sequence after design step: "
							<< pose.sequence() << std::endl;
					}
					if ( verbose ) tr() << "energy after repack: " << (*local_scorefxn)(pose) << std::endl;
				}
			}
			if ( verbose || local_debug ) tr() << std::flush;
		} //inner_cycle
	} //outer_cycle

	mc.show_counters();
	if ( recover_low_ ) {
		pose = mc.lowest_score_pose();
	} else {
		pose = mc.last_accepted_pose();
	}
	if ( local_movie ) {
		// this assumes there is only one loops_.
		Loops::const_iterator one_loop( loops()->begin() );
		Size begin_loop=one_loop->start();
		Size end_loop=one_loop->stop();
		loop_outfile << "MODEL" << std::endl;
		utility::vector1<Size> indices(end_loop - begin_loop + 3);
		for ( Size i=begin_loop-1, j=1; i<=end_loop+1; i++, j++ ) {
			indices[j]=i;
		}
		pose.dump_pdb(loop_outfile, indices, "final_refine");
		loop_outfile << "ENDMDL" << std::endl;
	}
}

std::string
LoopMover_Refine_KIC::get_name() const {
	return "LoopMover_Refine_KIC";
}

void
LoopMover_Refine_KIC::show(std::ostream & output) const
{
	Mover::show(output);
	output <<   "Neighbor distance of the KIC segment:   " << neighbor_dist_ <<
		"\nMaximum KIC segment length:             " << max_seglen_ <<
		"\nRecover the lowest energy conformation: " << (recover_low_ ? "True" : "False") <<
		"\nMinimize after repacking:               " << (min_after_repack_ ? "True" : "False") <<
		"\nFix sidechains neighboring the loop:    " << (fix_natsc_ ? "True" : "False") <<
		"\nOptimize only within the neighbor distance after every KIC move: " <<
		(optimize_only_kic_region_sidechains_after_move_ ? "True" : "False");
}

/// detailed
/// Update the MoveMaps that define, for each loop in loops_, which residues may be minimized.
/// Allows for minimization to apply only to the neighborhood around the segment that KIC moved, in cases
/// where multiple segments (loops) are defined.
void
LoopMover_Refine_KIC::update_movemap_vectors(
	core::pose::Pose & pose,
	utility::vector1<core::kinematics::MoveMap> & move_maps ) {

	for ( Size i = 1; i <= loops()->size(); i++ ) {
		protocols::loops::Loops cur_loop;
		cur_loop.add_loop( ( *loops() )[ i ] );
		core::kinematics::MoveMap cur_mm;
		loops_set_move_map( pose, cur_loop, fix_natsc_, cur_mm, neighbor_dist_ );
		move_maps[ i ] = cur_mm;
	}
}

/// detailed
/// Update the allow_sc_vectors that define, for each loop in loops_, which residue side-chains may be subject to
/// rotamer trials. Allows for rotamer trials to apply only to the neighborhood around the segment that KIC moved,
/// in cases where multiple segments (loops) are defined.
void
LoopMover_Refine_KIC::update_allow_sc_vectors(
	core::pose::Pose & pose,
	utility::vector1< utility::vector1< bool > > & allow_sc_vectors ) {

	for ( Size i = 1; i <= loops()->size(); i++ ) {
		protocols::loops::Loops cur_loop;
		cur_loop.add_loop( ( *loops() )[ i ] );
		utility::vector1<bool> cur_allow_sc_move( pose.total_residue(), false );
		select_loop_residues( pose, cur_loop, !fix_natsc_, cur_allow_sc_move, neighbor_dist_);
		allow_sc_vectors[ i ] = cur_allow_sc_move;
	}
}

/// detailed
/// Sets the rotamer trials packer task used in LoopMover_Refine_KIC::apply() to only rotamer trial the
/// neighborhood around KIC segment that moved. May slightly hurt performance in de novo loop reconstruction,
/// but may yield a speedup of 2X or more in cases of modeling very large or multiple segments (loop definitions),
/// in conjunction with set_movemap_from_kic_segment(), which are both used when
/// optimize_only_kic_region_sidechains_after_move_ is set to true.
void
LoopMover_Refine_KIC::set_rottrials_from_kic_segment(
	core::pose::Pose & pose,
	pack::task::PackerTaskOP & rottrials_packer_task,
	Size kic_start,
	Size kic_end ) {

	// we'll exploit the loops classes and the select_loop_residues() function to setup the PackerTask
	protocols::loops::Loop kic_seg( kic_start, kic_end, kic_start );
	utility::vector1<bool> to_trials( pose.total_residue(), false );
	select_loop_residues( pose, kic_seg, true, to_trials, neighbor_dist_ );
	rottrials_packer_task->restrict_to_residues( to_trials );
}

/// detailed
/// Sets the MoveMap for minimization used in LoopMover_Refine_KIC::apply() to only minimize the
/// neighborhood around KIC segment that moved. May slightly hurt performance in de novo loop reconstruction,
/// but may yield a speedup of 2X or more in cases of modeling very large or multiple segments (loop definitions),
/// in conjunction with set_movemap_from_kic_segment(), which are both used when
/// optimize_only_kic_region_sidechains_after_move_ is set to true.
void
LoopMover_Refine_KIC::set_movemap_from_kic_segment(
	core::pose::Pose & pose,
	kinematics::MoveMap & cur_mm,
	Size kic_start,
	Size kic_end ) {

	protocols::loops::Loops kic_seg;
	kic_seg.add_loop( kic_start, kic_end, kic_start );
	loops_set_move_map( pose, kic_seg, fix_natsc_, cur_mm, neighbor_dist_);
}

std::ostream &operator<< ( std::ostream &os, LoopMover_Refine_KIC const &mover )
{
	mover.show(os);
	return os;
}

basic::Tracer & LoopMover_Refine_KIC::tr() const
{
	return TR;
}

LoopMover_Refine_KICCreator::~LoopMover_Refine_KICCreator() {}

moves::MoverOP LoopMover_Refine_KICCreator::create_mover() const {
	return moves::MoverOP( new LoopMover_Refine_KIC() );
}

std::string LoopMover_Refine_KICCreator::keyname() const {
	return "LoopMover_Refine_KIC";
}

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
