// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: LoopModeler
///
///
/// @author Vatsan Raman

//devel headers
#include <devel/integrated_loop/LoopModeler.hh>
#include <devel/rbsegment_Moves/RBSegmentRelax.hh>
#include <devel/rbsegment_Moves/RBSegment.fwd.hh>
#include <devel/integrated_loop/CcdCloseMover.hh>

//core headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>


//protocols
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/loops/looprelax_protocols.hh>

//numeric headers
#include <numeric/random/random.hh>

//Auto Headers
#include <core/pose/util.hh>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;

namespace protocols {
namespace moves {

static THREAD_LOCAL basic::Tracer TR( "devel.IntegratedLoop.LoopModeler" );


//v typedef utility::vector1< protocols::Loop > Loops;
//v typedef utility::vector1< protocols::Loop >::iterator LoopsIt;

////////////////////////////LoopMover/////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
/// @brief LoopModeler set_default
//////////////////////////////////////////////////////////////////////////
void LoopMover::set_default()
{
	nmoves_ = 1;
	mc_temperature_ = 2.0;
	//	movemap_ = new core::kinematics::MoveMap();
	//	movemap_->set_bb( false );
	//	movemap_->set_chi( false );

}
//////////////////////////////////////////////////////////////////////////
/// @brief set_default_mc
//////////////////////////////////////////////////////////////////////////
void LoopMover::set_default_mc( core::pose::Pose & pose )
{

	mc_ = new moves::MonteCarlo( pose, *scorefxn_, mc_temperature_ );

}
//////////////////////////////////////////////////////////////////////////
/// @brief get_mc for viewer
//////////////////////////////////////////////////////////////////////////
protocols::moves::MonteCarloOP LoopMover::get_mc( core::pose::Pose & pose )
{
	set_default_mc( pose );
	mc_created = true;
	return mc_;

}
//////////////////////////////////////////////////////////////////////////
/// @brief set movemap
//////////////////////////////////////////////////////////////////////////
void LoopMover::set_movemap(
	Loops const & LoopsToPerturb,
	core::kinematics::MoveMapOP movemap_
)
{

	using namespace core::id;
	//	movemap_ = new core::kinematics::MoveMap();
	movemap_->set_bb( false );
	movemap_->set_chi( false );
	movemap_->set_jump( false );

	for ( LoopsConsIt it = LoopsToPerturb.begin(), it_end = LoopsToPerturb.end();
				it != it_end; ++it ) {
		for ( Size i = it->loop_begin(); i <= it->loop_end(); ++i ) {
			movemap_->set_bb(i, true);
			movemap_->set(TorsionID(i,BB,3), false); // omega is fixed
			movemap_->set_chi(i, true); // chi of loop residues
		}

	}
}

//////////////////////////////////////////////////////////////////////////
/// @brief set movemap overloaded
//////////////////////////////////////////////////////////////////////////
void LoopMover::set_movemap(
	protocols::Loop const & ThisLoop,
	core::kinematics::MoveMapOP movemap_one_loop
)

{
	using namespace core::id;
	movemap_one_loop->set_bb( false );
	movemap_one_loop->set_chi( false );
	movemap_one_loop->set_jump( false );

	for ( Size i = ThisLoop.loop_begin(); i <= ThisLoop.loop_end(); ++i ) {
		movemap_one_loop->set_bb(i, true);
		movemap_one_loop->set(TorsionID(i,BB,3), false); // omega is fixed
		movemap_one_loop->set_chi(i, true); // chi of loop residues
	}

}

////////////////////////////LoopModeler//////////////////////////////////
void LoopMover::set_one_loop_fold_tree(
	core::pose::Pose & pose,
	protocols::Loop const & ThisLoop
)
{
	using namespace kinematics;

	FoldTree f;
	Size nres( pose.size() );
	Size loop_begin( ThisLoop.loop_begin() );
	Size loop_end( ThisLoop.loop_end() );
	Size cutpoint( ThisLoop.cutpoint() );
	f.clear();
	if( loop_begin == 1 || loop_end == nres ) {
		f.add_edge( 1, nres, Edge::PEPTIDE ); //simple fold tree
	} else {
		f.add_edge( 1, loop_begin - 1, Edge::PEPTIDE ); //one jump fold tree
		f.add_edge( loop_begin - 1, cutpoint, Edge::PEPTIDE );
		f.add_edge( cutpoint + 1, loop_end + 1, Edge::PEPTIDE );
		f.add_edge( loop_end + 1, nres, Edge::PEPTIDE );
		f.add_edge( loop_begin - 1, loop_end + 1, 1 );
	}
	f.reorder( 1 );
	TR << "Pose fold tree " << f << "\n";

	pose.fold_tree( f );


}
//////////////////////////////////////////////////////////////////////////
/// @brief LoopModeler apply
//////////////////////////////////////////////////////////////////////////
void LoopMover::set_loops_fold_tree(
	core::pose::Pose & pose,
	Loops const & LoopsToPerturb
)
{

	using namespace kinematics;
	TR << "Do nothing yet !" << "\n";
	FoldTree f;
	Size const total_residue( pose.size() );
	Loops tmp_loops = LoopsToPerturb;
	Size jump_num = 0;


	for ( LoopsConsIt it = LoopsToPerturb.begin(), it_end = LoopsToPerturb.end(), it_next; it != it_end; ++it )
		{
		it_next = it;
		it_next++;
		Size const jump_start =
			( it->loop_begin() == 1 ) ? it->loop_begin() : it->loop_begin() - 1;
		Size const jump_stop  =
			( it->loop_end() == total_residue ) ? it->loop_end() : it->loop_end() + 1;
		Size const jump_cut   = it->cutpoint();
		Size const jump_next_start =
			( it_next == it_end ) ? total_residue : it_next->loop_begin()-1;
		if(  it->loop_begin() == 1 ){
			f.add_edge( jump_start, jump_stop, Edge::PEPTIDE );
			f.add_edge( jump_stop, jump_next_start, Edge::PEPTIDE );
			continue;
		}else if( it->loop_end() == total_residue ){
			f.add_edge( jump_start, jump_stop, Edge::PEPTIDE );
			continue;
		}
		jump_num++;
		f.add_edge( jump_start, jump_stop, jump_num );
		f.add_edge( jump_start, jump_cut,  Edge::PEPTIDE );
		f.add_edge( jump_cut+1, jump_stop, Edge::PEPTIDE );
		f.add_edge( jump_stop, jump_next_start, Edge::PEPTIDE );
		}

	Size const first_start =
		( tmp_loops.begin()->loop_begin() == 1 ) ? tmp_loops.begin()->loop_begin() : tmp_loops.begin()->loop_begin() - 1;
	//	if ( first_start != 1 )
	f.add_edge( 1, first_start, Edge::PEPTIDE );

	// reorder
	Size root;
	if( tmp_loops.begin()->loop_begin() == 1 &&
			tmp_loops.begin()->loop_end() != total_residue )
		root = tmp_loops.begin()->loop_end()+1;
	else root = 1;
	f.reorder(root);

	TR << "fold tree " << f << "\n";

	pose.fold_tree( f );

}

//////////////////////////////////////////////////////////////////////////
/// @brief LoopModeler apply
//////////////////////////////////////////////////////////////////////////
void LoopModeler::apply( core::pose::Pose & pose )
{

	T("devel.LoopModeler") << pose.size() << "\n";

}


//////////////////////////////////////////////////////////////////////////
/// @brief LoopModeler apply
//////////////////////////////////////////////////////////////////////////
void LoopModeler::apply_mod( core::pose::Pose & pose )
{

	T("devel.LoopModeler") << pose.size() << "\n";

}


////////////////////////////LoopRefiner//////////////////////////////////


//////////////////////////////////////////////////////////////////////////
/// @brief LoopRefiner apply ( DONT USE THIS )
//////////////////////////////////////////////////////////////////////////
void LoopRefiner::apply_mod( core::pose::Pose & pose )
{


	TR << "Inside LoopRefiner apply " << "\n";
	for ( LoopsIt it = LoopList_.begin(), it_end = LoopList_.end();
				it != it_end; ++it ) {
		TR << "Loop test " << it->loop_begin() << "\n";
	}


	scorefxn_->show( std::cout );


	Real init_temp( 2.0 );
	Real final_temp( 1.0 );
	Size cycles( 10 );


	movemap_ = new core::kinematics::MoveMap();
	mc_ = new protocols::moves::MonteCarlo( pose, *scorefxn_, init_temp );

	std::string min_type( "linmin" );
	Real const tolerance( 0.001 );
	bool const nblist( true );

	//	small_move_rot_trial_mover();
	//	shear_move_rot_trial_mover();
	LoopManager LM( pose, LoopList_ );
	Loops tmp = LM.LoopsToPerturb();


	for( LoopsIt it = tmp.begin(), it_end = tmp.end();
			 it != it_end; ++it ) {
		TR << "LoopRefiner apply " << it->loop_begin() << " " << it->loop_end() << " " << it->cutpoint() << " " << it->skip_rate() << "\n";
		set_one_loop_fold_tree( pose, *it );
		//v		set_movemap( *it );
		min_mover_ = new protocols::simple_moves::MinMover( movemap_, scorefxn_, min_type, tolerance, nblist );
		protocols::moves::SequenceMoverOP BBPerturb( new protocols::moves::SequenceMover() );
		//		BBPerturb->add_mover( SmallMovesRotTrial );
		//		BBPerturb->add_mover( ShearMovesRotTrial );
		BBPerturb->add_mover( min_mover_ );

		protocols::moves::TrialMoverOP BBPerturbMC( new protocols::moves::TrialMover( BBPerturb, mc_ ) );
		devel::RBSegment::SimAnnealMoverOP BBPerturbCycle( new devel::RBSegment::SimAnnealMover( BBPerturbMC, mc_, init_temp, final_temp, cycles ) );

		BBPerturbCycle->apply( pose );
		BBPerturb->apply( pose );

		pose.dump_pdb("before_recover_low.pdb");
	}
	T("devel.LoopRefiner") << pose.size() << "\n";

}

//////////////////////////////////////////////////////////////////////////
/// @brief LoopRefiner apply
//////////////////////////////////////////////////////////////////////////
void LoopRefiner::apply( core::pose::Pose & pose )
{


	//parameters - move to private data later
	Real init_temp( 2.0 );
	Real final_temp( 1.0 );
	Size cycles( 10 );
	std::string min_type( "linmin" );
	Real const tolerance( 0.001 );
	bool const nblist( true );

	LoopManager AllLoopsLM( pose, LoopList_ );

	//ROTAMER TRIALS has not been added !!


	movemap_ = new core::kinematics::MoveMap();
	set_movemap( LoopList_, movemap_ );//all loops minimized
	min_mover_ = new protocols::simple_moves::MinMover( movemap_, scorefxn_, min_type, tolerance, nblist );

	if( !mc_created )
		set_default_mc( pose );

	Loops LoopsToPerturb = AllLoopsLM.LoopsToPerturb();
	set_loops_fold_tree( pose, LoopsToPerturb );

	//	LoopManager LoopsToPerturbLM( pose, LoopsToPerturb );
	for( LoopsIt it = LoopsToPerturb.begin(), it_end = LoopsToPerturb.end();
			 it != it_end; ++it ) {

		std::cout << "Loop_begin, loop_end and cutpoint " << it->loop_begin() << " " << it->loop_end() << " " << it->cutpoint() << std::endl;

		core::kinematics::MoveMapOP movemap_one_loop( new core::kinematics::MoveMap() );
		set_movemap( *it, movemap_one_loop );

		if ( it->cutpoint() == pose.size()) {
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, it->cutpoint() - 1 );
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, it->cutpoint() );
	} else {
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, it->cutpoint() );
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, it->cutpoint() + 1 );
	}

		protocols::moves::SequenceMoverOP BBPerturb( new protocols::moves::SequenceMover() );
		protocols::moves::CcdCloseMoverOP ccd_close_mover( new protocols::moves::CcdCloseMover( pose, *movemap_one_loop, *it ) );
		BBPerturb->add_mover( small_move_rot_trial_mover( movemap_one_loop ) );
		BBPerturb->add_mover( shear_move_rot_trial_mover( movemap_one_loop ) );
		BBPerturb->add_mover( ccd_close_mover );
		BBPerturb->add_mover( min_mover_ );

		protocols::moves::TrialMoverOP BBPerturbMC( new protocols::moves::TrialMover( BBPerturb, mc_ ) );
		devel::RBSegment::SimAnnealMoverOP BBPerturbCycle( new devel::RBSegment::SimAnnealMover( BBPerturbMC, mc_, init_temp, final_temp, cycles ) );
		//		BBPerturbMC->apply( pose );
 	 	BBPerturbCycle->apply( pose );
		//		BBPerturb->apply( pose );

		if ( it->cutpoint() == pose.size() ) {
			core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, it->cutpoint() - 1 );
			core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, it->cutpoint() );
	} else {
			core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, it->cutpoint() );
			core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, it->cutpoint() + 1 );
	}

	}

	//	utility::exit( EXIT_FAILURE, __FILE__, __LINE__);

}
//////////////////////////////////////////////////////////////////////////
/// @brief LoopRefiner small_move_perturbation
//////////////////////////////////////////////////////////////////////////
protocols::moves::SequenceMoverOP LoopRefiner::small_move_rot_trial_mover(
	core::kinematics::MoveMapOP movemap_one_loop
)
{

	using namespace moves;
	//setup move objects
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( movemap_one_loop, mc_temperature_, nmoves_ ) );
	small_mover->angle_max( 'H', 0.0 );
	small_mover->angle_max( 'E', 5.0 );
	small_mover->angle_max( 'L', 6.0 );

	moves::SequenceMoverOP SmallMovesRotTrial( new moves::SequenceMover() );
	SmallMovesRotTrial->add_mover( small_mover );
	//v	SmallMovesRotTrial->add_mover( RotTrial );

	return SmallMovesRotTrial;
}

//////////////////////////////////////////////////////////////////////////
/// @brief LoopRefiner shear_move_perturbation
//////////////////////////////////////////////////////////////////////////
protocols::moves::SequenceMoverOP LoopRefiner::shear_move_rot_trial_mover(
	core::kinematics::MoveMapOP movemap_one_loop
)
{

	//setup move objects
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover( movemap_one_loop, mc_temperature_, nmoves_ ) );
	shear_mover->angle_max( 'H', 0.0 );
	shear_mover->angle_max( 'E', 5.0 );
	shear_mover->angle_max( 'L', 6.0 );

	protocols::moves::SequenceMoverOP ShearMovesRotTrial( new moves::SequenceMover() );
	ShearMovesRotTrial->add_mover( shear_mover );
	//v	ShearMovesRotTrial->add_mover( RotTrial );

	return ShearMovesRotTrial;

}

//////////////////////////////////////////////////////////////////////////
/// @brief LoopRefiner Rotamer_trials
//////////////////////////////////////////////////////////////////////////
protocols::simple_moves::RotamerTrialsMoverOP LoopRefiner::rotamer_trial_mover(
	core::pose::Pose & pose
)
{

	core::pack::task::PackerTaskOP task_;
	task_ = pack::task::TaskFactory::create_packer_task( pose );
	task_->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	task_->set_bump_check( true );
	protocols::simple_moves::RotamerTrialsMoverOP RotTrial = new protocols::simple_moves::RotamerTrialsMover( scorefxn_, *task_ );
	return RotTrial;

}


}//moves
}//protocols
