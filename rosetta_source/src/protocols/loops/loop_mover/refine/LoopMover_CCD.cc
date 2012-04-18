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
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCDCreator.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/conformation/Residue.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// Rosetta Headers
#include <core/chemical/VariantType.hh>

// AUTO-REMOVED #include <core/conformation/util.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <protocols/simple_moves/FragmentMover.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

#include <protocols/loops/util.hh>
#include <basic/Tracer.hh> // tracer output

//Utility and numeric Headers
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>

// option key includes

#include <basic/options/keys/MonteCarlo.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>

//Auto Headers


//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

///////////////////////////////////////////////////////////////////////////////
using namespace core;

static numeric::random::RandomGenerator RG(84846);
static basic::Tracer TR("protocols.loops.loop_mover.refine.LoopMover_Refine_CCD");

//constructors
LoopMover_Refine_CCD::LoopMover_Refine_CCD()
	: LoopMover(),
		outer_cycles_(3),
		max_inner_cycles_(200),
		repack_period_(20),
		temp_initial_(1.5),
		temp_final_(0.5),
		set_fold_tree_from_loops_(false)
{
	read_options();
	set_scorefxn( get_fa_scorefxn() );
	protocols::moves::Mover::type("LoopMover_Refine_CCD");
	set_default_settings();
}

LoopMover_Refine_CCD::LoopMover_Refine_CCD(
	protocols::loops::LoopsOP  loops_in
) : LoopMover( loops_in ),
		outer_cycles_(3),
		max_inner_cycles_(200),
		repack_period_(20),
		temp_initial_(1.5),
		temp_final_(0.5),
		set_fold_tree_from_loops_(false)
{
	read_options();
	set_scorefxn( get_fa_scorefxn() );
	protocols::moves::Mover::type("LoopMover_Refine_CCD");
	set_default_settings();
}


LoopMover_Refine_CCD::LoopMover_Refine_CCD(
	protocols::loops::LoopsOP  loops_in,
	core::scoring::ScoreFunctionOP  scorefxn
) : LoopMover( loops_in ),
		outer_cycles_(3),
		max_inner_cycles_(200),
		repack_period_(20),
		temp_initial_(1.5),
		temp_final_(0.5),
		set_fold_tree_from_loops_(false)
{
	read_options();
	set_scorefxn( scorefxn );
	protocols::moves::Mover::type("LoopMover_Refine_CCD");
	set_default_settings();
}

//destructor
LoopMover_Refine_CCD::~LoopMover_Refine_CCD() {}

std::string
LoopMover_Refine_CCD::get_name() const {
	return "LoopMover_Refine_CCD";
}

//clone
protocols::moves::MoverOP LoopMover_Refine_CCD::clone() const {
		return new LoopMover_Refine_CCD(*this);
}


void LoopMover_Refine_CCD::set_default_settings()
{
		redesign_loop_ = false;
		packing_isolated_to_active_loops_ = false;
}
loop_mover::LoopMover::MoveMapOP LoopMover_Refine_CCD::move_map() const
{
    return move_map_;
}
void LoopMover_Refine_CCD::move_map( LoopMover::MoveMapOP mm )
{
    move_map_ = mm;
}

void
LoopMover_Refine_CCD::read_options()
{
	using namespace basic::options;
	if ( option[ OptionKeys::loops::outer_cycles ].user() )
		outer_cycles_ = option[ OptionKeys::loops::outer_cycles ]();
	if ( option[ OptionKeys::loops::max_inner_cycles ].user() )
		max_inner_cycles_ = option[ OptionKeys::loops::max_inner_cycles ]();
	if ( option[ OptionKeys::loops::repack_period ].user() )
		repack_period_ = option[ OptionKeys::loops::repack_period ]();
	if ( option[ OptionKeys::MonteCarlo::temp_initial ].user() )
		temp_initial_ = option[ OptionKeys::MonteCarlo::temp_initial ]();
	if ( option[ OptionKeys::MonteCarlo::temp_final ].user() )
		temp_final_ = option[ OptionKeys::MonteCarlo::temp_final ]();
}

void LoopMover_Refine_CCD::set_task_factory(
	core::pack::task::TaskFactoryCOP task_factory_in
)
{
	// make local, non-const copy from const input
	runtime_assert( task_factory_in );
	task_factory_ = new core::pack::task::TaskFactory( *task_factory_in );
}

core::pack::task::TaskFactoryCOP LoopMover_Refine_CCD::get_task_factory() const { return task_factory_; }


void
LoopMover_Refine_CCD::parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose ){
  packing_isolated_to_active_loops_ = false;
	//using parser implies that the fold tree probably isn't set correctly
	set_fold_tree_from_loops( tag->getOption< bool >( "set_fold_tree_from_loops", true ) );
	utility::vector1< utility::tag::TagPtr > const branch_tags( tag->getTags() );
	bool specified_movemap( false );
	foreach( utility::tag::TagPtr const tag, branch_tags ){
		if( tag->getName() == "MoveMap" ) specified_movemap = true;
		break;
	}
	if( specified_movemap ){
		move_map_ = new core::kinematics::MoveMap;
		move_map_->set_bb( false );
		move_map_->set_chi( false );
		move_map_->set_jump( false );
		protocols::rosetta_scripts::parse_movemap( tag, pose, move_map_, data, false/*don't reset movemap, keep falses, unless stated otherwise*/ );
	}
	if( tag->hasOption( "loops" ) ){
		std::string const loops_str( tag->getOption< std::string >( "loops" ) );
		loops( loops_from_string( loops_str, pose ) );
	}
	if( tag->hasOption( "scorefxn" ) ) this->set_scorefxn( new core::scoring::ScoreFunction( *data.get< core::scoring::ScoreFunction * >( "scorefxns", tag->getOption<std::string>( "scorefxn" ) ) ) );

	if( tag->hasOption("task_operations") ){
		core::pack::task::TaskFactoryOP task_factory = protocols::rosetta_scripts::parse_task_operations( tag, data );
		this->set_task_factory( task_factory );
	}
	else task_factory_ = NULL;

	if( tag->hasOption( "loops_from_cache" ) ) set_use_loops_from_observer_cache( tag->getOption<bool>( "loops_from_cache", 1 ) );

	if( tag->hasOption( "outer_cycles" ) ) outer_cycles_ = tag->getOption<core::Size>( "outer_cycles", 3 );
	if( tag->hasOption( "max_inner_cycles" ) ) max_inner_cycles_ = tag->getOption<core::Size>( "max_inner_cycles", 250 );
	temp_initial( tag->getOption< core::Real >( "temp_initial", 1.5 ) );
	temp_final( tag->getOption< core::Real >( "temp_final", 0.5 ) );
}

void LoopMover_Refine_CCD::apply(
	core::pose::Pose & pose
) {
	using namespace optimization;
	using namespace scoring;
	using namespace basic::options;

	core::pose::Pose native_pose;
	core::kinematics::FoldTree f_orig;
	if ( get_native_pose() ) native_pose = *get_native_pose();
	else native_pose = pose;

	if( use_loops_from_observer_cache() ) this->set_loops_from_pose_observer_cache( pose );

	if( set_fold_tree_from_loops_ ){
		core::kinematics::FoldTree f_new;
		f_orig = pose.fold_tree();
		loops::fold_tree_from_loops( pose, *( this->loops() ), f_new);
		pose.fold_tree( f_new );
		loops::add_cutpoint_variants( pose );
	}

	// 'verbose' should be an option, or better yet Tracer's output filtering capabilities should be used instead
	bool const verbose( true );
	bool const local_debug( option[ OptionKeys::loops::debug ].user() );

	// set cutpoint variants for correct chainbreak scoring
	Size const nres( pose.total_residue() );   //fpd fix for symmetry
	utility::vector1< bool > is_loop( nres, false );

	for( Loops::const_iterator it=loops()->begin(), it_end=loops()->end();
			 it != it_end; ++it ) {
		for ( Size i= it->start(); i<= it->stop(); ++i ) {
			is_loop[i] = true;
		}
		Size const loop_cut(it->cut());
		//if ( loop_cut == nres ) { //c-terminal loop
		//	utility_exit_with_message("This doesnt make sense; ask Chu about it, Phil.");
		//	if ( ! pose.residue(loop_cut - 1).has_variant_type(chemical::CUTPOINT_LOWER) )
		//		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, loop_cut - 1);
		//	if ( ! pose.residue(loop_cut).has_variant_type(chemical::CUTPOINT_UPPER) )
		//		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, loop_cut );
		//} else {

	//if ( loop_cut != nres ) { //c-terminal loop  //fpd fix for symmetry
	if ( !pose.residue(loop_cut).is_upper_terminus() ) {
			if ( ! pose.residue(loop_cut).has_variant_type(chemical::CUTPOINT_LOWER) )
				core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, loop_cut );
			if ( ! pose.residue(loop_cut+1).has_variant_type(chemical::CUTPOINT_UPPER) )
				core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, loop_cut+1 );
		}
	}

	// scheduler
	int const fast( option[ OptionKeys::loops::fast ] ); // why is this an int?
	Size const inner_cycles(
		std::min( max_inner_cycles_, fast ? loops()->loop_size() : Size(10)*loops()->loop_size() ) );

	// scorefxn
	scoring::ScoreFunctionOP scorefxn;
	if ( scorefxn() != 0 ) scorefxn = scorefxn()->clone();
	else scorefxn = get_fa_scorefxn();

	// confirm that chainbreak weight is set
	scorefxn->set_weight( chainbreak, 1. * 10. / 3. );

	// monte carlo
	Real const gamma(
		std::pow( ( temp_final_ / temp_initial_ ), Real( 1.0f / ( outer_cycles_ * inner_cycles ) ) ) );

	Real temperature( temp_initial_ );
	// need to make sure we have scored before we ask for a tenA graph later...
	(*scorefxn)(pose);

	protocols::moves::MonteCarlo mc( pose, *scorefxn, temperature );
	// minimizer
	AtomTreeMinimizerOP minimizer;
	MinimizerOptions options( "dfpmin", 0.001, true /*use_nblist*/, false /*deriv_check*/ );
	bool const repack_neighbors( (! option[ OptionKeys::loops::fix_natsc ])  || task_factory_ );
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		minimizer = dynamic_cast<AtomTreeMinimizer*>( new core::optimization::symmetry::SymAtomTreeMinimizer );
	} else {
		minimizer = new core::optimization::AtomTreeMinimizer;
	}

	pack::task::PackerTaskOP base_packer_task;

	// create default Packer behavior if none has been set via TaskFactory
	if ( task_factory_ == 0 ) {
		// the default Packer behavior is defined here
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		task_factory_ = new TaskFactory;
		task_factory_->push_back( new InitializeFromCommandline );
		task_factory_->push_back( new IncludeCurrent );
		task_factory_->push_back( new NoRepackDisulfides );
		base_packer_task = task_factory_->create_task_and_apply_taskoperations( pose );
		if ( redesign_loop_ ) {
			// allow design at loop positions
			for ( Size i=1; i<= nres; ++i ) {
				if ( !is_loop[i] ) base_packer_task->nonconst_residue_task( i ).restrict_to_repacking();
			}
		} else {
			// restrict to repacking at all positions
			base_packer_task->restrict_to_repacking();
		}
		// additional default behavior: packing restricted to active loops
		packing_isolated_to_active_loops_ = true;
	} else {
		base_packer_task = task_factory_->create_task_and_apply_taskoperations( pose );
	}
	base_packer_task->set_bump_check( true );

	// repack trial
	pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
	utility::vector1<bool> allow_repacked( nres, false );
	if ( packing_isolated_to_active_loops_ ) {
		select_loop_residues( pose, *loops(), repack_neighbors, allow_repacked, 10.0 /* neighbor_cutoff */ );
		core::pose::symmetry::make_residue_mask_symmetric( pose, allow_repacked );
		             // does nothing if pose is not symm
		this_packer_task->restrict_to_residues( allow_repacked );
	}
	pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	std::string move_type = "repack";
	mc.boltzmann( pose, move_type );
	mc.show_scores();

	if (local_debug) {
		pose.dump_pdb("tmp_fa_repack.pdb");
		std::ofstream out("score.tmp_repack_fa");
		out << "scoring for repack_fa " << (*scorefxn)(pose) << std::endl;
		scorefxn->show( out );
		out << pose.energies();
	}

	// remember all packable/repacked residues
	for ( Size index(1); index <= nres; ++index ) {
		if ( !this_packer_task->residue_task(index).being_packed() ) allow_repacked[index] = false;
		else allow_repacked[index] = true;
	}
	// set minimization degrees of freedom for all loops
	kinematics::MoveMapOP mm_all_loops = new core::kinematics::MoveMap;
	setup_movemap( pose, *loops(), allow_repacked, mm_all_loops );

	// small/shear move parameters
	// should be options
	Size const nmoves = { 1 };
	std::map< char, Real > angle_max;
	angle_max.insert( std::make_pair('H', 0.0) );
	angle_max.insert( std::make_pair('E', 5.0) );
	angle_max.insert( std::make_pair('L', 6.0) );

	for (Size i=1; i<=outer_cycles_; ++i) {
		// increase CHAINBREAK weight and update monte carlo
		scorefxn->set_weight( chainbreak, Real(i)*10.0/3.0 );
		mc.score_function( *scorefxn );
		// recover low
		mc.recover_low( pose );
		// score info
		if ( verbose ) tr() << "cycle: " << i << "  " << (*scorefxn)(pose) << std::endl;
		for ( Size j=1; j<=inner_cycles; ++j ) {
			temperature *= gamma;
			mc.set_temperature( temperature );
			if ( verbose ) tr() << "refinement cycle (outer/inner): "
								<< i << "/" << outer_cycles_ << " "
								<< j << "/" << inner_cycles << " "
								<< std::endl;
			{// small_CCD_min_trial
				Loops::const_iterator it( loops()->one_random_loop() );
				Loops one_loop;
				one_loop.add_loop( it );
				// set up movemap properly
				kinematics::MoveMapOP mm_one_loop( new kinematics::MoveMap() );
				setup_movemap( pose, one_loop, allow_repacked, mm_one_loop );

				if (local_debug) {
					tr() << "chutmp-debug small_move-0: " << "  " << (*scorefxn)(pose) << std::endl;
					tr() << "small_move-0: " << pose.energies().total_energies().weighted_string_of( scorefxn->weights() )
										<< " rmsd: " << F(9,3,loop_rmsd( pose, native_pose, *loops() )) << std::endl;
					pose.dump_pdb("small_move-0.pdb");
				}

				protocols::simple_moves::SmallMover small_moves( mm_one_loop, temperature, nmoves );
				small_moves.apply( pose );

				if (local_debug) {
					tr() << "chutmp-debug small_move-1: " << "  " << (*scorefxn)(pose) << std::endl;
					tr() << "small_move-1: " << pose.energies().total_energies().weighted_string_of( scorefxn->weights() )
										<< " rmsd: " << F(9,3,loop_rmsd( pose, native_pose, *loops() )) << std::endl;
					pose.dump_pdb("small_move-1.pdb");
					std::ofstream out("score.small_move_1");
					out << "scoring of input_pose " << (*scorefxn)(pose) << std::endl;
					scorefxn->show( out );
					out << pose.energies();
				}

				if (! it->is_terminal( pose ) ) ccd_close_loops( pose, one_loop, *mm_one_loop);

				if (local_debug) {
					tr() << "chutmp-debug small_move-2: " << "  " << (*scorefxn)(pose) << std::endl;
					tr() << "small_move-2: " << pose.energies().total_energies().weighted_string_of( scorefxn->weights() )
										<< " rmsd: " << F(9,3,loop_rmsd( pose, native_pose, *loops() )) << std::endl;
					pose.dump_pdb("small_move-2.pdb");
					std::ofstream out("score.small_move_2");
					out << "scoring of input_pose " << (*scorefxn)(pose) << std::endl;
					scorefxn->show( out );
					out << pose.energies();
				}

				pack::rotamer_trials( pose, *scorefxn, this_packer_task );

				if (local_debug) {
					tr() << "chutmp-debug small_move-3: " << "  " << (*scorefxn)(pose) << std::endl;
					tr() << "small_move-3: " << pose.energies().total_energies().weighted_string_of( scorefxn->weights() )
										<< " rmsd: " << F(9,3,loop_rmsd( pose, native_pose, *loops() )) << std::endl;
					pose.dump_pdb("small_move-3.pdb");
				}

				setup_movemap( pose, *loops(), allow_repacked, mm_all_loops );
				minimizer->run( pose, *mm_all_loops, *scorefxn, options );

				if (local_debug) {
					tr() << "chutmp-debug small_move-4: " << "  " << (*scorefxn)(pose) << std::endl;
					tr() << "small_move-4: " << pose.energies().total_energies().weighted_string_of( scorefxn->weights() )
										<< " rmsd: " << F(9,3,loop_rmsd( pose, native_pose, *loops() )) << std::endl;
					pose.dump_pdb("small_move-4.pdb");
				}

				std::string move_type = "small_ccd_min";
				mc.boltzmann( pose, move_type );

				if (local_debug) {
					tr() << "chutmp-debug small_move-5: " << "  " << (*scorefxn)(pose) << std::endl;
					tr() << "small_move-5: " << pose.energies().total_energies().weighted_string_of( scorefxn->weights() )
										<< " rmsd: " << F(9,3,loop_rmsd( pose, native_pose, *loops() )) << std::endl;
					pose.dump_pdb("small_move-5.pdb");
				}

				mc.show_scores();
			}
			{// shear_CCD_min_trial
				Loops::const_iterator it( loops()->one_random_loop() );
				Loops one_loop;
				one_loop.add_loop( it );
				// set up movemap properly
				kinematics::MoveMapOP mm_one_loop( new kinematics::MoveMap() );
				setup_movemap( pose, one_loop, allow_repacked, mm_one_loop );
				protocols::simple_moves::ShearMover shear_moves( mm_one_loop, temperature, nmoves );
				shear_moves.apply( pose );
				if (! it->is_terminal( pose ) ) ccd_close_loops( pose, one_loop, *mm_one_loop);
				pack::rotamer_trials( pose, *scorefxn, this_packer_task );
				(*scorefxn)(pose); // update 10A nbr graph, silly way to do this
				setup_movemap( pose, *loops(), allow_repacked, mm_all_loops );
				minimizer->run( pose, *mm_all_loops, *scorefxn, options );
				std::string move_type = "shear_ccd_min";
				mc.boltzmann( pose, move_type );
				mc.show_scores();
			}
			{ //main_repack_trial
				if ( (j%repack_period_)==0 || j==inner_cycles ) {
					// repack trial

					if ( packing_isolated_to_active_loops_ ) {
						select_loop_residues( pose, *loops(), repack_neighbors, allow_repacked, 10.0 /* neighbor_cutoff */ );
					}
					core::pose::symmetry::make_residue_mask_symmetric( pose, allow_repacked );  //fpd symmetrize res mask -- does nothing if pose is not symm
					this_packer_task->restrict_to_residues( allow_repacked );
					pack::pack_rotamers( pose, *scorefxn, this_packer_task );
					std::string move_type = "repack";
					mc.boltzmann( pose, move_type );
					mc.show_scores();
				}
			}
			if ( verbose || local_debug ) tr() << std::flush;
		} //inner_cycle
	} //outer_cycle

	mc.show_counters();
	pose = mc.lowest_score_pose();

	if( set_fold_tree_from_loops_ ){ //if requested, put back old foldtree
		loops::remove_cutpoint_variants( pose );
		pose.fold_tree( f_orig );
		(*scorefxn)(pose);
	}
}


/// @brief setup an appropriate movemap for the given loops
/// @param[in] loops The loops to model.
/// @param[in] allow_repack Indicates whether or not to allow a position to
///  repack.
/// @param[out] movemap Output movemap, all settings added here.
/// @remarks will enforce the false movemap
void LoopMover_Refine_CCD::setup_movemap(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops,
	utility::vector1< bool > const & allow_repack,
	core::kinematics::MoveMapOP & movemap
)
{
	if( move_map_ ){
		movemap = move_map_;
		return;
	}
	loops_set_move_map( loops, allow_repack, *movemap );
	enforce_false_movemap( movemap );
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	}
}

basic::Tracer & LoopMover_Refine_CCD::tr() const
{
    return TR;
}

LoopMover_Refine_CCDCreator::~LoopMover_Refine_CCDCreator() {}

moves::MoverOP LoopMover_Refine_CCDCreator::create_mover() const {
  return new LoopMover_Refine_CCD();
}

std::string LoopMover_Refine_CCDCreator::keyname() const {
  return "LoopMover_Refine_CCD";
}

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
