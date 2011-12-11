// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody2/LoopRlxMover.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)


#include <protocols/antibody2/LoopRlxMover.hh>

#include <basic/Tracer.hh>

#include <protocols/simple_moves/BackboneMover.fwd.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/vector1.hh>

#include <core/conformation/Conformation.hh>
#include <protocols/antibody2/CDRH3Modeler2.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/moves/RepeatMover.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <protocols/loops/CcdLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <core/pose/util.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/ResFilters.hh>


namespace protocols {
namespace antibody2 {
using namespace core;

static basic::Tracer TRL("protocols.antibody2.LoopRlxMover");

LoopRlxMover::LoopRlxMover( Size query_start, Size query_end ) : Mover( "LoopRlxMover" ) {
	loop_start_ = query_start;
	loop_end_ = query_end;
	set_default();
} // LoopRlxMover default constructor

void LoopRlxMover::set_default() {
	highres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( "standard", "score12" );
	highres_scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	highres_scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );

	movemap_= new kinematics::MoveMap();
	movemap_->set_chi( false );
	movemap_->set_bb( false );

	mc_ = new moves::MonteCarlo( *highres_scorefxn_, 0.8 );
	tf_ = new core::pack::task::TaskFactory;

	Real const init_temp( 2.0 );
	Real const last_temp( 0.5 );
	inner_cycles_ = (loop_end_ - loop_start_) + 1;
	outer_cycles_ = 2;
	gamma_ = std::pow( (last_temp/init_temp), (1.0/inner_cycles_));
	temperature_ = init_temp;
} // LoopRlxMover::set_default

void LoopRlxMover::setup_objects( pose::Pose & pose ) {
	using namespace protocols::moves;
	using namespace protocols::loops;

	//setting MoveMap
	utility::vector1< bool> allow_bb_move( pose.total_residue(), false );
	for( Size ii = loop_start_; ii <= loop_end_; ii++ )
		allow_bb_move[ ii ] = true;
	movemap_->set_bb( allow_bb_move );
	movemap_->set_jump( 1, false );

	Size loop_size = ( loop_end_ - loop_start_ ) + 1;
	Size cutpoint = loop_start_ + int(loop_size/2);

	one_loop_ = new Loop( loop_start_, loop_end_, cutpoint, 0, false );
	simple_one_loop_fold_tree( pose, *one_loop_ );

	// set cutpoint variants for correct chainbreak scoring
	if( !pose.residue( cutpoint ).is_upper_terminus() ) {
		if( !pose.residue( cutpoint ).has_variant_type( chemical::CUTPOINT_LOWER ) ) core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpoint );
		if( !pose.residue( cutpoint + 1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpoint + 1 );
	}

	Real min_tolerance = 0.001;
	if( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = "dfpmin_armijo_nonmonotone";
	bool nb_list = true;
	loop_min_mover_ = new simple_moves::MinMover( movemap_, highres_scorefxn_, min_type, min_tolerance, nb_list );

	// more params
	Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
	if( benchmark_ ) {
		n_small_moves = 1;
		inner_cycles_ = 1;
		outer_cycles_ = 1;
	}

	Real high_move_temp = 2.00;
	// minimize amplitude of moves if correct parameter is set
    simple_moves::BackboneMoverOP small_mover = new simple_moves::SmallMover( movemap_, high_move_temp, n_small_moves );
    simple_moves::BackboneMoverOP shear_mover = new simple_moves::ShearMover( movemap_, high_move_temp, n_small_moves );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 5.0 );
	small_mover->angle_max( 'L', 6.0 );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 5.0 );
	shear_mover->angle_max( 'L', 6.0 );

	CcdMoverOP ccd_moves = new CcdMover( *one_loop_, movemap_ );
	RepeatMoverOP ccd_cycle = new RepeatMover(ccd_moves, n_small_moves);

	wiggle_loop_ = new protocols::moves::SequenceMover();
	wiggle_loop_->add_mover( small_mover );
	wiggle_loop_->add_mover( shear_mover );
	wiggle_loop_->add_mover( ccd_cycle );
	wiggle_loop_->add_mover( loop_min_mover_ );
}


/// @brief relaxes the region specified
/// @detailed This is all done in high resolution.Hence there are no rigid
///           body moves relative to the docking partners. Only small moves
///           are carried out here to see if there are better fits.
///           Repacking is carried out extensively after each move.
///
/// @param[in] pose, loop begin position, loop end position
/// @global_read none
/// @global_write none

void LoopRlxMover::apply( pose::Pose & pose_in ) {
	using namespace protocols;
	using namespace protocols::loops;
	using namespace protocols::moves;
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;

	TRL<<"I am here 7.4.3"<<std::endl;
	TRL << "LoopRlxMover: Apply" << std::endl;

	if ( !pose_in.is_fullatom() )
		utility_exit_with_message("Fullatom poses only");

	setup_objects( pose_in );

	// storing starting fold tree
	kinematics::FoldTree tree_in( pose_in.fold_tree() );

	utility::vector1< bool> allow_repack( pose_in.total_residue(), false );
	select_loop_residues( pose_in, *one_loop_, true /*include_neighbors*/, allow_repack);
	movemap_->set_chi( allow_repack );

    simple_moves::PackRotamersMoverOP loop_repack = new simple_moves::PackRotamersMover(highres_scorefxn_);
	setup_packer_task( pose_in );
	( *highres_scorefxn_ )( pose_in );
	tf_->push_back( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) );
	loop_repack->task_factory(tf_);
	loop_repack->apply( pose_in );

	// rotamer trials - select loop with new neighbors after packing occurs
	select_loop_residues( pose_in, *one_loop_, true /*include_neighbors*/, allow_repack);
	movemap_->set_chi( allow_repack );
	setup_packer_task( pose_in );
	( *highres_scorefxn_ )( pose_in );
	tf_->push_back( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) );
    simple_moves::RotamerTrialsMoverOP pack_rottrial = new simple_moves::RotamerTrialsMover( highres_scorefxn_, tf_ );

	pack_rottrial->apply( pose_in );

	mc_->reset( pose_in ); // monte carlo reset

	// outer cycle
	for(Size i = 1; i <= outer_cycles_; i++) {
		mc_->recover_low( pose_in );

		// inner cycle
		for ( Size j = 1; j <= inner_cycles_; j++ ) {
			temperature_ *= gamma_;
			mc_->set_temperature( temperature_ );
			wiggle_loop_->apply( pose_in );

			// rotamer trials- select loop with new neighbors
			select_loop_residues( pose_in, *one_loop_, true /*include_neighbors*/, allow_repack);
			movemap_->set_chi( allow_repack );
			setup_packer_task( pose_in );
			( *highres_scorefxn_ )( pose_in );
			tf_->push_back( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) );
            simple_moves::RotamerTrialsMoverOP pack_rottrial = new simple_moves::RotamerTrialsMover( highres_scorefxn_, tf_ );
			pack_rottrial->apply( pose_in );
			mc_->boltzmann( pose_in );

			if ( numeric::mod(j,Size(20))==0 || j==inner_cycles_ ) {
				// repack trial
				loop_repack = new simple_moves::PackRotamersMover( highres_scorefxn_ );
				setup_packer_task( pose_in );
				( *highres_scorefxn_ )( pose_in );
				tf_->push_back( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) );
				loop_repack->task_factory( tf_ );
				loop_repack->apply( pose_in );
				mc_->boltzmann( pose_in );
			}
		} // inner cycles
	} // outer cycles
	mc_->recover_low( pose_in );

	// minimize
	if( !benchmark_ )
		loop_min_mover_->apply( pose_in );

	// Restoring pose stuff
	pose_in.fold_tree( tree_in ); // Tree

	TRL << "LoopRlxMover: Finished Apply" << std::endl;
} // LoopRlxMover::apply

std::string
LoopRlxMover::get_name() const { return "LoopRlxMover"; }



void
LoopRlxMover::setup_packer_task( pose::Pose & pose_in ) {
	using namespace pack::task;
	using namespace pack::task::operation;
	tf_ = new TaskFactory;

	TRL << "LoopRlxMover Setting Up Packer Task" << std::endl;

	tf_->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
	tf_->push_back( new InitializeFromCommandline );
	tf_->push_back( new IncludeCurrent );
	tf_->push_back( new RestrictToRepacking );
	tf_->push_back( new NoRepackDisulfides );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new pack::rotamer_set::UnboundRotamersOperation();
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation = new operation::AppendRotamerSet( unboundrot );
	tf_->push_back( unboundrot_operation );
	// adds scoring bonuses for the "unbound" rotamers, if any
	core::pack::dunbrack::load_unboundrot( pose_in );

	TRL << "LoopRlxMover Done: Setting Up Packer Task" << std::endl;

} // LoopRlxMover::setup_packer_task

} // antibody2
} // protocols
