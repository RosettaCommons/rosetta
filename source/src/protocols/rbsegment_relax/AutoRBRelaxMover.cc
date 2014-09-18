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
/// @author Srivatsan Raman
/// @author Frank DiMaio

#include <core/types.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>

#include <core/scoring/Energies.hh>

#include <protocols/rbsegment_relax/AutoRBRelaxMover.hh>
#include <protocols/rbsegment_relax/RBSegmentMover.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <protocols/rbsegment_relax/util.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/relax/FastRelax.hh>

#include <core/scoring/electron_density/util.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>


#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <core/fragment/FragSet.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <numeric/random/random.hh>

// C++ headers
#include <iostream>
#include <string>


//options
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>

#include <core/kinematics/MoveMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace rbsegment_relax {

static thread_local basic::Tracer tr( "protocols.rbsegment_relax.RBSegment.AutoRBRelaxMover" );

////////////////
// stupid mover
class CCDMoveWrapper : public protocols::moves::Mover {
public:
	CCDMoveWrapper( core::kinematics::MoveMapOP movemap, core::Size start, core::Size stop, core::Size cut) :
		movemap_(movemap),
		start_(start),
		stop_(stop),
		cut_(cut) { }

	void apply( Pose & pose ) {
		protocols::loops::loop_closure::ccd::CCDLoopClosureMover ccd_mover(
				protocols::loops::Loop( start_, stop_, cut_ ), movemap_ );
		ccd_mover.max_cycles( 125 );
		ccd_mover.apply( pose );
	}

	virtual std::string get_name() const {
		return ("CCDMoveWrapper");
	}

private:
	core::kinematics::MoveMapOP movemap_;
	core::Size start_, stop_, cut_;
};

typedef utility::pointer::owning_ptr< CCDMoveWrapper >  CCDMoveWrapperOP;



///////////////////
///
AutoRBMover::AutoRBMover() {
	using namespace basic::options;

	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
		option[ OptionKeys::RBSegmentRelax::rb_scorefxn ]() );

	movemap_ = new core::kinematics::MoveMap();

	nouter_cycles_ = option[ OptionKeys::RBSegmentRelax::nrboutercycles ]();
	ninner_cycles_ = option[ OptionKeys::RBSegmentRelax::nrbmoves ]();

	loop_melt_ = 3;

	// load frags
	loops::read_loop_fragments( frag_libs_ );

	allowSeqShiftMoves_ = !(option[ OptionKeys::RBSegmentRelax::skip_seqshift_moves]());

	allowSSFragInserts_ = false;

	// fa stuff
	tf_ = new core::pack::task::TaskFactory();
	tf_->push_back( new core::pack::task::operation::RestrictToRepacking );
	tf_->push_back( new core::pack::task::operation::InitializeFromCommandline );
	tf_->push_back( new core::pack::task::operation::IncludeCurrent );
	tf_->push_back( new core::pack::task::operation::NoRepackDisulfides );

	fa_scorefxn_ = core::scoring::get_score_function();
	fa_scorefxn_->set_weight( core::scoring::chainbreak, 10.0/3.0);

	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn_ );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *fa_scorefxn_ );
	}
}

void
AutoRBMover::apply( core::pose::Pose & pose ) {
	// ensure pose is in centroid mode
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	to_centroid.apply( pose );

	// ensure pose is rooted on VRT
	core::pose::addVirtualResAsRoot( pose );

	// Get DSSP parse; use it to:
	//    build fold tree (star topology)
	//    set movemap
	//    set variant types
	setup_topology( pose );

	// load constraints
	core::scoring::constraints::add_constraints_from_cmdline( pose, *scorefxn_ );

	core::Size nres =  pose.total_residue()-1;
	bool loops_closed = false;
	grow_flexible( loop_melt_, nres, 1 );

	core::pose::Pose start_pose = pose;

	while( !loops_closed ) {
		pose = start_pose;

		// setup fragment movers
		utility::vector1< protocols::simple_moves::FragmentMoverOP > fragmover;
		for ( utility::vector1< core::fragment::FragSetOP >::const_iterator
					it = frag_libs_.begin(), it_end = frag_libs_.end();
					it != it_end; it++ ) {
			protocols::simple_moves::ClassicFragmentMoverOP cfm = new protocols::simple_moves::ClassicFragmentMover( *it, movemap_ );
			cfm->set_check_ss( false );
			cfm->enable_end_bias_check( false );
			fragmover.push_back( cfm );
		}

		// mc object
		core::Real init_temp = 2.0;
		core::Real temperature = init_temp;
		core::Real final_temp = 0.6;
		core::Real const gamma = std::pow(final_temp/init_temp, (1.0/(nouter_cycles_*ninner_cycles_)) );
		moves::MonteCarloOP mc_ = new moves::MonteCarlo( pose, *scorefxn_, init_temp );

		// movement loop
		float final_chain_break_weight = 1.0;
		float delta_weight( final_chain_break_weight/nouter_cycles_ );

		// random mover
		protocols::moves::RandomMover random_move;

		// loop fragment insertion
		for ( std::vector< protocols::simple_moves::FragmentMoverOP >::const_iterator
				it = fragmover.begin(),it_end = fragmover.end(); it != it_end; it++ )
			random_move.add_mover(*it, rb_chunks_.size());

		// rigid-body move
		for (int i=1; i<=(int)rb_chunks_.size(); ++i)
			random_move.add_mover(new rigid::RigidBodyPerturbMover( i , 3.0 , 1.0 ));

		//TODO rigid-chunk fragment insertion
		//if (allowSSFragInserts_) ;

		// sequence shift
		if (allowSeqShiftMoves_) {
			for (int i=1; i<=(int)rb_chunks_.size(); ++i)
			for (int j=1; j<=(int)rb_chunks_[i].nContinuousSegments(); ++j) {
				protocols::moves::SequenceMoverOP seq_shift_move = new protocols::moves::SequenceMover;
				seq_shift_move->add_mover( new SequenceShiftMover(rb_chunks_[i][j]) );

				// find adjacent loops
				for (core::Size k=1; k<=loops_.size(); ++k) {
					bool adjLoopN = (loops_[k].stop() >= rb_chunks_[i][j].start()-1) && (loops_[k].stop() <= rb_chunks_[i][j].end()+1);
					bool adjLoopC = (loops_[k].start() >= rb_chunks_[i][j].start()-1) && (loops_[k].start() <= rb_chunks_[i][j].end()+1);
					if ( adjLoopN || adjLoopC ) {
						seq_shift_move->add_mover( new CCDMoveWrapper(movemap_, loops_[k].start(), loops_[k].stop(), loops_[k].cut() ) );
					}
				}
				random_move.add_mover(seq_shift_move, 0.5);
			}
		}

		scorefxn_->set_weight( core::scoring::chainbreak, 0.0 );
		for( int n1 = 1; n1 <= (int)nouter_cycles_; ++n1 ) {
			mc_->recover_low( pose );
			scorefxn_->set_weight( core::scoring::chainbreak, n1*delta_weight );

			(*scorefxn_)(pose);
			if ( tr.visible() ) { scorefxn_->show_line( tr.Info , pose ); }
			tr.Info << std::endl;
			mc_->score_function( *scorefxn_ );

			for( int n2 = 1; n2 <= (int)ninner_cycles_; ++n2 ) {
				// cool
				temperature *= gamma;
				mc_->set_temperature( temperature );

				// randomly do something
				if( numeric::random::uniform()*nouter_cycles_ > n1 ) {
					random_move.apply(pose);
				} else {
					protocols::loops::Loops::const_iterator it( loops_.one_random_loop() );
					protocols::loops::loop_closure::ccd::CCDLoopClosureMover ccd_mover( *it, movemap_ );
					ccd_mover.max_cycles( 25 );
					ccd_mover.apply( pose );
				}

				mc_->boltzmann( pose );
			}
		}
		mc_->recover_low( pose );
		mc_->show_counters();

		scorefxn_->set_weight( core::scoring::chainbreak, 1.0 );
		(*scorefxn_)(pose);

		loops_closed = ( pose.energies().total_energies()[ core::scoring::chainbreak ] ) <= loops_.size()*0.5;
		if (!loops_closed) {
			tr << "Loops not closed! ("
			   << pose.energies().total_energies()[ core::scoring::chainbreak ]
			   << " > " << loops_.size()*0.5 << ")" << std::endl;
			grow_flexible( loop_melt_, nres );
		} else {
			tr << "Loops closed! ("
			   << pose.energies().total_energies()[ core::scoring::chainbreak ]
			   << " <= " << loops_.size()*0.5 << ")" << std::endl;
		}
	}

	//////////////
	// fastrelax -- keep foldtree
	relax::FastRelax fast_relax( fa_scorefxn_ );
	fast_relax.apply( pose );
}

////
//// grow loops
//// dont allow jump residues to be flexible
void
AutoRBMover::grow_flexible( core::Size maxlen , core::Size nres , core::Size minlen ) {
	if (maxlen == 0) return;

	tr << "EXTENDING LOOPS:" << std::endl;
	for ( core::Size i=1; i <= loops_.size(); i++ ) {
		core::Size extend_start = (core::Size) numeric::random::random_range(minlen, maxlen-minlen);
		core::Size extend_stop  = (core::Size) numeric::random::random_range(minlen, maxlen-minlen);
		if ( ( extend_start == 0 ) && ( extend_stop == 0 ) ) {
			if ( numeric::random::uniform() > 0.5) extend_start = 1;
			else extend_stop  = 1;
		}

		// dont go past termini
		if (loops_[i].start()  < 1 + extend_start) extend_start = loops_[i].start()-1;
		if (loops_[i].stop() + extend_stop > nres ) extend_stop = nres - loops_[i].stop();

		// make sure we dont go over a jump point
		for (core::Size j=1; j<=jumps_.size(); ++j) {
			if (loops_[i].start()-extend_start <= jumps_[j] && loops_[i].start() > jumps_[j])
				extend_start = loops_[i].start() - jumps_[j] + 1;
			if (loops_[i].stop()+extend_stop >= jumps_[j] && loops_[i].stop() < jumps_[j])
				extend_stop = jumps_[j] - loops_[i].stop() +1;
		}

		loops_[i].set_start( loops_[i].start()-extend_start );
		loops_[i].set_stop( loops_[i].stop()+extend_stop );
	}
	tr << loops_ << std::endl;
}

////
//// set up foldtree, variants, movemap, etc.
void
AutoRBMover::setup_topology( core::pose::Pose & pose ) {
	//core::Size nres = pose.total_residue()-1; // terminal VRT

	rigid_segs_.clear();
	rb_chunks_.clear();
	jumps_.clear();
	loops_.clear();

	guess_rbsegs_from_pose( pose, rigid_segs_, rb_chunks_, loops_ );
	jumps_ = setup_pose_rbsegs_keep_loops( pose,  rigid_segs_ , loops_,  movemap_ );
}

}
}

