// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief RBSegmentRelax protocol
/// @details
///
///
///
/// @author Srivatsan Raman
/// @author Frank DiMaio

#include <core/pose/util.hh>
#include <protocols/rbsegment_relax/FragInsertAndAlignMover.hh>
#include <protocols/rbsegment_relax/RBSegmentRelax.hh>
#include <protocols/rbsegment_relax/RBSegment.fwd.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <protocols/rbsegment_relax/RBSegmentMover.hh>
#include <protocols/rbsegment_relax/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>

#include <core/fragment/FragSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

//Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

//Map
#include <map>

#include <core/id/SequenceMapping.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rbsegment_relax {

static thread_local basic::Tracer TS( "protocols.moves.RBSegmentRelax" );
using namespace core;
using basic::Error;

RBSegmentRelax::RBSegmentRelax() {}
RBSegmentRelax::~RBSegmentRelax() {}

RBSegmentRelax::RBSegmentRelax(
	core::scoring::ScoreFunctionOP scorefxn,
	utility::vector1< RBSegment > const & rbsegs_input,
	protocols::loops::Loops      const & loops_input
) :
	Mover("RBSegmentRelax"),
	scorefxn_( scorefxn ),
	rbsegs_input_( rbsegs_input ),
	loops_input_( loops_input )
{
	using namespace basic::options;

	// Default parameter settings
	init_temp = 2.0;
	final_temp = 1.0;

	// take defaults from the command line
	utility::vector1< core::Real > helix_def   = option[ OptionKeys::RBSegmentRelax::helical_movement_params ]();
	utility::vector1< core::Real > strand_def  = option[ OptionKeys::RBSegmentRelax::strand_movement_params ]();
	utility::vector1< core::Real > generic_def = option[ OptionKeys::RBSegmentRelax::default_movement_params ]();
	cycles_ = option[ OptionKeys::RBSegmentRelax::nrbmoves ]();

	helical_sigR = helix_def[1]; helical_sigT = helix_def[2];
	helical_sigOffAxisR = helix_def[3]; helical_sigOffAxisT = helix_def[4];

	strand_sigR = strand_def[1]; strand_sigT = strand_def[2];
	strand_sigOffAxisR = strand_def[3]; strand_sigOffAxisT = strand_def[4];

	genericRB_sigR = generic_def[1]; genericRB_sigT = generic_def[2];

	cst_weight_ = option[ OptionKeys::constraints::cst_weight ]();
	cst_width_  = option[ OptionKeys::relax::coord_cst_width ]();
	cst_stdev_  = option[ OptionKeys::relax::coord_cst_stdev ]();
	cst_seqwidth_ = option[ OptionKeys::RBSegmentRelax::cst_seqwidth ]();

	rand_ = 0;
	bootstrap_ = no_lr_ = false;
}


//////////////////////////////////////////////////////////////////////////
/// @brief setup_RBSegmentRelax; initialize movesets to defaults
//////////////////////////////////////////////////////////////////////////
void RBSegmentRelax::initialize( utility::vector1< core::fragment::FragSetOP > const &frag_libs , core::Real rnd)
{
	using namespace basic::options;

	// set up default movesets ( start and end positions are set in apply() )
	if ( ! option[ OptionKeys::RBSegmentRelax::skip_seqshift_moves ]() ) {
		HelixMoveSet_.push_back( RBSegmentMoverOP( new SequenceShiftMover(  ) ) );
		StrandMoveSet_.push_back( RBSegmentMoverOP( new SequenceShiftMover(  ) ) );
	}

	if ( ! option[ OptionKeys::RBSegmentRelax::skip_rb_moves ]() ) {
		HelixMoveSet_.push_back( RBSegmentMoverOP( new HelicalGaussianMover( helical_sigR, helical_sigT, helical_sigOffAxisR, helical_sigOffAxisT ) ) );
		StrandMoveSet_.push_back( RBSegmentMoverOP( new StrandTwistingMover( strand_sigR, strand_sigT, strand_sigOffAxisR, strand_sigOffAxisT ) ) );
		GenericRBMoveSet_.push_back( RBSegmentMoverOP( new GaussianRBSegmentMover( genericRB_sigR , genericRB_sigT ) ) );
		CompositeSegmentMoveSet_.push_back( RBSegmentMoverOP( new GaussianRBSegmentMover( genericRB_sigR , genericRB_sigT ) ) );
	}

	// frags
	frag_libs_ = frag_libs;
	randomness_ = rnd;
}


////////////////////////////
void RBSegmentRelax::set_temperature( Real start, Real final ) {
	init_temp = start;
	final_temp = final;
}

////////////////////////////
void RBSegmentRelax::set_helicalMoveStepsize( Real onAxisTrans, Real onAxisRot, Real offAxisTrans, Real offAxisRot ) {
	helical_sigT = onAxisTrans;
	helical_sigR = onAxisRot;
	helical_sigOffAxisT = offAxisTrans;
	helical_sigOffAxisR = offAxisRot;
}

////////////////////////////
void RBSegmentRelax::set_genericRBMoveStepsize( Real trans, Real rot ) {
	genericRB_sigT = trans;
	genericRB_sigR = rot;
}

////////////////////////////
void RBSegmentRelax::set_ncycles( int ncycles ) {
	cycles_ = ncycles;      // cycles PER rb-segment
}


////////////////////////////
void RBSegmentRelax::set_cst_weight( core::Real wt ) {
	cst_weight_ = wt;
}

////////////////////////////
void RBSegmentRelax::set_cst_width ( core::Real width ) {
	cst_width_ = width;
}


//////////////////////////////////////////////////////////////////////////
/// @brief apply method
//////////////////////////////////////////////////////////////////////////
void RBSegmentRelax::apply( core::pose::Pose & pose ) {
	using namespace basic::options;

	protocols::moves::RandomMoverOP SegmentRandomizeMover( new protocols::moves::RandomMover() );

	// Remove loops from pose, connect with jumps
	core::pose::Pose pose_noloops, pose_input = pose;
	core::id::SequenceMapping resmap;
	core::kinematics::MoveMap mm;
	bool fullatom = pose.is_fullatom();

	// make loopless pose
	setup_pose_from_rbsegs( rbsegs_input_ , pose , pose_noloops , resmap, mm, fix_ligands_ );
	TS << "LOOPS-REMOVED fold tree " << pose_noloops.fold_tree() << std::endl;
	TS << "LOOPS-REMOVED secstruct " << pose_noloops.secstruct() << std::endl;
	protocols::viewer::add_conformation_viewer( pose_noloops.conformation() );

	// remap rbsegs
	remap_rb_segments( rbsegs_input_, rbsegs_remap_, resmap);

	// set up constraints
	set_rb_constraints( pose_noloops, pose, rbsegs_input_ , resmap, cst_width_, cst_stdev_, cst_seqwidth_ );  // does the mapping internally
	scorefxn_->set_weight( core::scoring::coordinate_constraint, cst_weight_ );  // add constraint term

	bool doFragInserts=false;
	FragInsertAndAlignMoverOP frag_ins;
	if ( option[ OptionKeys::loops::vall_file ].user() && !option[ OptionKeys::RBSegmentRelax::skip_fragment_moves ]() ) {
		frag_ins = FragInsertAndAlignMoverOP( new
			FragInsertAndAlignMover(rbsegs_remap_, pose_noloops, randomness_ ) );
		doFragInserts = true;

		if ( bootstrap_ ) {
			frag_ins->bootstrapCATrace( pose_noloops ); // Bootstrap a model from the Ca trace (idealizing in the process)

			// if fullatom do a repack now
			core::pack::task::PackerTaskOP task_fast = core::pack::task::TaskFactory::create_packer_task( pose_noloops );
			task_fast->restrict_to_repacking(); task_fast->or_include_current(false);
			protocols::simple_moves::PackRotamersMover pack_fast( scorefxn_, task_fast );
			if ( fullatom ) {
				pack_fast.apply( pose_noloops );
			}
		}
	}

	// loop over segments, add a mover for each
	int nmovers=0;
	for ( RBIt it_seg = rbsegs_remap_.begin(), it_seg_end = rbsegs_remap_.end();
			it_seg != it_seg_end; ++it_seg ) {
		// hardcoded since it needs a few extra things
		// initialize the frag-insert mover
		if ( doFragInserts ) {
			SegmentRandomizeMover->add_mover( frag_ins );
			nmovers++;
		}

		// add the whole structure movers
		// add a copy for each RB segment
		for ( std::vector< protocols::moves::MoverOP >::iterator it_mover = WholeStructureMoveSet_.begin(),
				it_mover_end = WholeStructureMoveSet_.end();
				it_mover != it_mover_end; ++it_mover ) {
			SegmentRandomizeMover->add_mover( *it_mover );
			nmovers++;
		}

		//////////////////////
		// HELIX SEGMENTS
		if ( it_seg->isHelix() ) {
			// add each mover from the helix moveset
			// resids refer to those in the loop-removed pose
			for ( std::vector< RBSegmentMoverOP >::iterator it_mover = HelixMoveSet_.begin(),
					it_mover_end = HelixMoveSet_.end(); it_mover != it_mover_end; ++it_mover ) {
				RBSegmentMoverOP moverToAdd = utility::pointer::dynamic_pointer_cast< RBSegmentMover >((*it_mover)->clone());  // make a deep copy of the mover
				moverToAdd->setResidueRange( *it_seg );
				if ( it_seg->initialized() ) moverToAdd->set_movement(*it_seg);  // override default movement params
				SegmentRandomizeMover->add_mover( moverToAdd );
				nmovers++;
			}

			//////////////////////
			// SHEET SEGMENTS
		} else if ( it_seg->isSheet() ) {
			// add each mover from the strand moveset
			// resids refer to those in the loop-removed pose
			for ( std::vector< RBSegmentMoverOP >::iterator it_mover = StrandMoveSet_.begin(),
					it_mover_end = StrandMoveSet_.end(); it_mover != it_mover_end; ++it_mover ) {
				RBSegmentMoverOP moverToAdd = utility::pointer::dynamic_pointer_cast< RBSegmentMover >((*it_mover)->clone()); // make a deep copy of the mover
				if ( it_seg->initialized() ) moverToAdd->set_movement(*it_seg);  // override default movement params
				moverToAdd->setResidueRange( *it_seg );
				SegmentRandomizeMover->add_mover( moverToAdd );
				nmovers++;
			}

			//////////////////////
			// GENERIC RB SEGMENTS
		} else if ( it_seg->isGenericRB() ) {
			// Neither helix nor strand, we can include more logic here to be more specific
			// add each genericRB mover
			// resids refer to those in the loop-removed pose
			for ( std::vector< RBSegmentMoverOP >::iterator it_mover = GenericRBMoveSet_.begin(),
					it_mover_end = GenericRBMoveSet_.end(); it_mover != it_mover_end; ++it_mover ) {
				RBSegmentMoverOP moverToAdd = utility::pointer::dynamic_pointer_cast< RBSegmentMover >((*it_mover)->clone());  // make a deep copy of the mover
				moverToAdd->setResidueRange( *it_seg );
				if ( it_seg->initialized() ) moverToAdd->set_movement(*it_seg);  // override default movement params
				SegmentRandomizeMover->add_mover( moverToAdd );
				nmovers++;
			}

			//////////////////////
			// COMPOUND SEGMENTS
		} else if ( it_seg->isCompound() ) {
			if ( ! it_seg->initialized() ) {
				it_seg->set_movement( genericRB_sigT, genericRB_sigR );
			}
			// add each genericRB mover
			// resids refer to those in the loop-removed pose
			for ( std::vector< RBSegmentMoverOP >::iterator it_mover = CompositeSegmentMoveSet_.begin(),
					it_mover_end = CompositeSegmentMoveSet_.end(); it_mover != it_mover_end; ++it_mover ) {
				RBSegmentMoverOP moverToAdd = utility::pointer::dynamic_pointer_cast< RBSegmentMover >((*it_mover)->clone()); // make a deep copy of the mover
				moverToAdd->setResidueRange( *it_seg );
				if ( it_seg->initialized() ) moverToAdd->set_movement(*it_seg);  // override default movement params
				SegmentRandomizeMover->add_mover( moverToAdd );
				nmovers++;
			}

			/////////////////////
			//  error?
		} else {
			TS << "[ ERROR ] Unknown segment type\n";
			exit(1);
		}
	}

	if ( nmovers == 0 ) {
		TS << "[ ERROR ] RBSegmentRelax::apply() caled with no segments or movers defined!\n";
		return;
	}

	// repacking task
	// to do: use factory
	core::pack::task::PackerTaskOP taskstd = core::pack::task::TaskFactory::create_packer_task( pose_noloops );
	taskstd->initialize_from_command_line();
	taskstd->restrict_to_repacking();
	taskstd->or_include_current(true);
	protocols::simple_moves::PackRotamersMover pack( scorefxn_, taskstd );

	// rand
	if ( rand_>0 ) {
		for ( int i=1; i<=(int)rand_; ++i ) {
			SegmentRandomizeMover->apply( pose_noloops );
		}

		// repack
		pack.apply( pose_noloops );
	}

	// set up move-minimization cycle

	// create MC object
	//  -  Currently, this is not exposed in the interface, since the MC object needs
	//     the pose object _without_ loops, which is not created until apply is called
	//  -  In the future, removal of loops should be separated from apply, which would
	//     allow the MC object to be exposed in the interface
	//  -  A side effect of this is that the viewer -- that needs the MC object --
	//     can only be called here.  Again, this should be fixed in the future
	mc_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo( pose_noloops , *scorefxn_ , init_temp ) );

	// wrap in a MC search for 'cycles' iteration
	int total_cycles = cycles_;
	core::Real temperature = init_temp;
	core::Real gamma = std::pow( final_temp/init_temp, 1.0/total_cycles );
	mc_->reset( pose_noloops );
	mc_->set_temperature( temperature );

	//Minimizer object
	optimization::AtomTreeMinimizer mzr;
	optimization::MinimizerOptions options1( "linmin", 0.1, true, false );
	core::optimization::MinimizerOptions options2( "dfpmin_armijo_nonmonotone", 1e-3, true, false );

	// set movemaps
	if ( fullatom ) {
		mm.set_bb( true );
		mm.set_chi( true );

		// fix ligands
		if ( fix_ligands_ ) {
			for ( int i=1; i<=(int)pose.total_residue(); ++i ) {
				if ( pose.residue(i).is_ligand() ) {
					mm.set_bb( false );
					mm.set_chi( false );
				}
			}
		}
	} else {
		//mm.set_bb( true );
		mm.set_bb( false );  // ? think about this
		mm.set_chi( false );
	}

	scorefxn_->show( TS.Debug , pose_noloops );
	for ( int i=0; i<total_cycles; ++i ) {
		temperature *= gamma;
		mc_->set_temperature( temperature );

		// 1: random move
		SegmentRandomizeMover->apply( pose_noloops );

		// 2: minimization
		if ( fullatom ) {
			// if fullatom mode
			//    repack + SC + RB minimize (as in docking)
			// TO DO restrict packing to interface reses only
			pack.apply( pose_noloops );
			mzr.run( pose_noloops, mm, *scorefxn_, options2 );
		} else {
			// quick BB + RB linmin
			mzr.run( pose_noloops, mm, *scorefxn_, options1 );
		}

		// 3: MC
		(*scorefxn_)(pose_noloops);
		scorefxn_->show_line( TS.Debug , pose_noloops );
		//scorefxn_->show( std::cout , pose_noloops );
		if (  mc_->boltzmann( pose_noloops )  ) TS.Debug << " *";
		TS.Debug << std::endl;

		mc_->boltzmann( pose_noloops );
	}

	// recover low-energy pose
	mc_->recover_low( pose_noloops );
	mc_->show_counters();

	if ( no_lr_ || loops_input_.size() == 0 ) {
		pose = pose_noloops;
		// rescore
		(*scorefxn_)(pose);
	} else {
		//pose_noloops.dump_pdb("rb_noloops.pdb");
		restore_pose_from_rbsegs( rbsegs_input_ , pose_noloops, pose );

		// cut loops & restore original loop conformation
		// dilate loops
		protocols::loops::LoopsOP loops( new protocols::loops::Loops( loops_input_ ) );
		for ( loops::Loops::iterator it = loops->v_begin(), it_end = loops->v_end(); it != it_end; ++it ) {
			if ( !pose.fold_tree().is_cutpoint( it->start() - 1 ) ) {
				it->set_start( it->start() - 1 );
			}
			if ( !pose.fold_tree().is_cutpoint( it->stop() ) ) {
				it->set_stop( it->stop() + 1 );
			}
		}

		core::kinematics::FoldTree f, f_in = pose.fold_tree();
		loops->auto_choose_cutpoints( pose );
		protocols::loops::fold_tree_from_loops( pose, *loops, f);

		//pose.dump_pdb("precopy.pdb");
		//protocols::viewer::add_conformation_viewer( pose.conformation() );   // <<<< added in looprelax mover
		pose.fold_tree( f );
		for ( loops::Loops::const_iterator it = loops->begin(), it_end = loops->end(); it != it_end; ++it ) {
			core::Size lstart = it->start(), lstop = it->stop();
			idealize_loop( pose, *it );
			for ( core::Size k=lstart; k<lstop; ++k ) {
				pose.set_phi( k, pose_input.phi(k) );
				pose.set_psi( k, pose_input.psi(k) );
				pose.set_omega( k, pose_input.omega(k) );
			}
		}
		pose.fold_tree( f_in );

		// get LR params
		std::string remodel, intermedrelax, refine, relax;

		// defaults
		remodel = "quick_ccd";
		refine = "no";
		intermedrelax = "no";
		relax = "fastrelax";

		// override defaults is specified on the command line
		if ( option[ OptionKeys::loops::remodel ].user() ) {
			remodel = option[ OptionKeys::loops::remodel ]();
		}
		if ( option[ OptionKeys::loops::intermedrelax ].user() ) {
			intermedrelax = option[ OptionKeys::loops::intermedrelax ]();
		}
		if ( option[ OptionKeys::loops::refine ].user() ) {
			refine = option[ OptionKeys::loops::refine ]();
		}
		if ( option[ OptionKeys::loops::relax ].user() ) {
			relax = option[ OptionKeys::loops::relax ]();
		}

		// set up LR constraints (if provided in input pose)
		set_constraints( pose, pose_input, cst_width_, cst_stdev_, cst_seqwidth_ );

		// setup scorefunction for loopbuilding
		core::scoring::ScoreFunctionOP lr_fa_scorefxn = core::scoring::get_score_function();
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *lr_fa_scorefxn );

		// call looprelax
		protocols::comparative_modeling::LoopRelaxMover mover;
		mover.fa_scorefxn( lr_fa_scorefxn );
		mover.frag_libs( frag_libs_ );
		mover.loops( loops );
		mover.relax( relax );
		mover.refine( refine );
		mover.remodel( remodel );
		mover.intermedrelax( intermedrelax );
		mover.cmd_line_csts( false );   // csts come from this protocol, not the cmd line
		// however, use wts from cmd line
		mover.apply( pose );
		//TS << "FINAL fold tree " << pose.fold_tree() << std::endl;
	}
}


std::string
RBSegmentRelax::get_name() const {
	return "RBSegmentRelax";
}

//////////////////////////////////////////////////////////////////////////
/// @brief
//////////////////////////////////////////////////////////////////////////
void RBSegmentRelax::add_helixMover( RBSegmentMoverOP newMover ) {
	HelixMoveSet_.push_back( newMover );
}


//////////////////////////////////////////////////////////////////////////
/// @brief
//////////////////////////////////////////////////////////////////////////
void RBSegmentRelax::add_strandMover( RBSegmentMoverOP newMover ) {
	StrandMoveSet_.push_back( newMover );
}


//////////////////////////////////////////////////////////////////////////
/// @brief
//////////////////////////////////////////////////////////////////////////
void RBSegmentRelax::add_genericRBMover( RBSegmentMoverOP newMover ) {
	GenericRBMoveSet_.push_back( newMover );
}


//////////////////////////////////////////////////////////////////////////
/// @brief
//////////////////////////////////////////////////////////////////////////
void RBSegmentRelax::add_compositeSegmentMover( RBSegmentMoverOP newMover ) {
	CompositeSegmentMoveSet_.push_back( newMover );
}

//////////////////////////////////////////////////////////////////////////
/// @brief
//////////////////////////////////////////////////////////////////////////
void RBSegmentRelax::add_wholeStructureMover( protocols::moves::MoverOP newMover ) {
	WholeStructureMoveSet_.push_back( newMover );
}


//////////////////////////////////////////////////////////////////////////
/// @brief
//////////////////////////////////////////////////////////////////////////
void RBSegmentRelax::clear_movesets() {
	HelixMoveSet_.clear();
	StrandMoveSet_.clear();
	GenericRBMoveSet_.clear();
	CompositeSegmentMoveSet_.clear();
	WholeStructureMoveSet_.clear();
}

}// moves
}//protocols
