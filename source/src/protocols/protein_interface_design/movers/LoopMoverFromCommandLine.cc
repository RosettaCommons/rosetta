// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/Loops/LoopMoverFromCommandLine.cc
/// @brief Parseable class to do full loop remodeling with input fragment files from command line
/// @author Jordan Willis (jordan.r.willis@vanderbilt.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/LoopMoverFromCommandLine.hh>
#include <protocols/protein_interface_design/movers/LoopMoverFromCommandLineCreator.hh>

// Package headers

// Project headers
// AUTO-REMOVED #include <protocols/loops/kinematic_closure/KinematicWrapper.hh>
// AUTO-REMOVED #include <protocols/loops/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_QuickCCD.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loops_main.hh> // for various loop utility fxns
#include <protocols/loops/Loops.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/forge/methods/util.hh>
// AUTO-REMOVED #include <protocols/forge/remodel/RemodelLoopMover.hh>


#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh>
#include <protocols/rosetta_scripts/util.hh>

// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>

// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/dssp/Dssp.hh> //getting SS from frag files

#include <core/scoring/ScoreFunction.hh>


#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
// AUTO-REMOVED #include <core/pack/task/operation/OperateOnCertainResidues.hh>
// AUTO-REMOVED #include <core/pack/task/operation/ResLvlTaskOperations.hh>

// AUTO-REMOVED #include <protocols/simple_moves/PackRotamersMover.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/packing.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/loops.OptionKeys.gen.hh>
//create option keys for loop movers

// AUTO-REMOVED #include <core/fragment/ConstantLengthFragSet.hh>
// AUTO-REMOVED #include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragSet.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

// AUTO-REMOVED #include <core/fragment/Frame.hh>
// AUTO-REMOVED #include <core/fragment/picking_old/vall/util.hh>
// AUTO-REMOVED #include <core/fragment/picking_old/FragmentLibraryManager.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {

static thread_local basic::Tracer TR( "protocols.moves.LoopRemodelFromCommandLine" );
static thread_local basic::Tracer TR_report( "protocols.moves.LoopRemodelFromCommandLine.REPORT" );

std::string
LoopMoverFromCommandLineCreator::keyname() const
{
	return LoopMoverFromCommandLineCreator::mover_name();
}

protocols::moves::MoverOP
LoopMoverFromCommandLineCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopMoverFromCommandLine );
}

std::string	LoopMoverFromCommandLineCreator::mover_name()
{
	return "LoopMoverFromCommandLine";
}

LoopMoverFromCommandLine::~LoopMoverFromCommandLine() {}


protocols::moves::MoverOP
LoopMoverFromCommandLine::clone() const
{
	return( protocols::moves::MoverOP( new LoopMoverFromCommandLine( *this ) ) );
}


//call on empty constructor
LoopMoverFromCommandLine::LoopMoverFromCommandLine() :
	simple_moves::DesignRepackMover( LoopMoverFromCommandLineCreator::mover_name() ),
	intermedrelax_( "no" ),
	remodel_( "no" ),
	relax_( "no" ),
	string_refine_( "no" )
{
	design(false);
}

//full member variables defined in constructor

LoopMoverFromCommandLine::LoopMoverFromCommandLine(
		std::string const protocol,
		bool const perturb,
		bool const refine,
		core::scoring::ScoreFunctionOP & hires_score,
		core::scoring::ScoreFunctionOP & lores_score,
		std::string const loop_file_name,
		protocols::loops::LoopsCOP loops
		) :
		simple_moves::DesignRepackMover ( LoopMoverFromCommandLineCreator::mover_name()),
		protocol_ ( protocol ),
		perturb_( perturb),
		refine_(refine),
		intermedrelax_( "no" ),
		remodel_( "no" ),
		relax_( "no" ),
		string_refine_( "no" )
{
		hires_score_ = hires_score;
		lores_score = lores_score->clone();
		loop_file_name_= loop_file_name;
		loops_ = protocols::loops::LoopsOP( new protocols::loops::Loops( *loops ) );
		design(false);
}



//apply to pose
void
LoopMoverFromCommandLine::apply ( core::pose::Pose & pose)
{
	using namespace protocols::loops;
	using core::pack::task::operation::TaskOperationCOP;
	pose.conformation().detect_disulfides(); // I don't think that this is important but just in case
	core::pose::Pose native_pose = pose;
	loops::set_secstruct_from_psipred_ss2( pose );
	LoopsOP loops( new protocols::loops::Loops( loop_file_name_ ) );
	loops->verify_against(pose);
	loops->auto_choose_cutpoints(pose);
	if( loops->size() == 0)  {
		TR << "No loops found!" << std::endl;
		return; // bounce out if we didn't define any loops
	}
	else
	{
		TR << *loops << std::endl;
	}
	if( loops->size() > 0 ) {
		core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::IncludeCurrent ) );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::NoRepackDisulfides ) );
		// set up temporary fold tree for loop closure
		TR.Debug << "Original FoldTree " << pose.fold_tree() << std::endl;
		core::kinematics::FoldTree old_ft( pose.fold_tree() );
		for( Loops::iterator it = loops->v_begin(); it != loops->v_end(); ++it ) {
						it->set_extended( true ); // set all loops to extended (needed for CCD mover to really perturb)
						protocols::loops::LoopsOP single_loop( new protocols::loops::Loops() );
						single_loop->add_loop(*it);
		core::kinematics::FoldTree new_ft = protocols::forge::methods::fold_tree_from_loops( pose, *single_loop );
		pose.fold_tree( new_ft );
		add_cutpoint_variants( pose );
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		//pose will always start full atom
		//protocols::moves::MonteCarlo mc( pose, *scorefxn_repack_, mc_kt );
		protocols::protein_interface_design::movers::SaveAndRetrieveSidechains retrieve_sc( pose );
		retrieve_sc.allsc( true );

		if( protocol_ == "automatic" ){
			utility::vector1< core::fragment::FragSetOP > frag_libs;
			protocols::loops::read_loop_fragments( frag_libs );

			protocols::comparative_modeling::LoopRelaxMover lrm;
			lrm.frag_libs( frag_libs );
			lrm.loops( single_loop );
			lrm.relax( relax() );
			lrm.refine( string_refine() );
			lrm.remodel( remodel() );
			lrm.intermedrelax( intermedrelax() );
			lrm.scorefxns( lores_score_, hires_score_ );
			lrm.apply( pose );
			return;
		}
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID);
		if( protocol_ == "kinematic" ) {
						if( perturb_ ) {
							protocols::loops::loop_mover::perturb::LoopMover_Perturb_KIC perturb(single_loop, lores_score_ );
							perturb.set_native_pose( PoseCOP( new core::pose::Pose ( native_pose ) ) );
							perturb.apply( pose );
						}
						core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
						retrieve_sc.apply( pose ); // recover sidechains from pre-centroid pose
						if( refine_ ) {
							protocols::loops::loop_mover::refine::LoopMover_Refine_KIC refine( single_loop, hires_score_ );
							refine.set_redesign_loop(false); // design?
							refine.set_native_pose( PoseCOP( new core::pose::Pose ( native_pose ) ) );
							pose.update_residue_neighbors();
							refine.apply( pose );
						}
		} // protocol == kinematic
		else if( protocol_ == "ccd" ) {
			TR << "Task Factory =" << task_factory;
			TR << "ccd protocol" << std::endl;
			pose.update_residue_neighbors();
			core::scoring::dssp::Dssp dssp( pose );
			dssp.insert_ss_into_pose( pose );
			std::string const full_ss = pose.secstruct();
			std::string const full_sequence = pose.sequence();
			utility::vector1< core::fragment::FragSetOP > frag_libs;
			protocols::loops::read_loop_fragments( frag_libs );
				if( perturb_ ) {
					protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCD perturb(single_loop, lores_score_ );
					for ( core::Size i = 1; i <= frag_libs.size(); ++i ) {
						perturb.add_fragments( frag_libs[i] );
					}
					perturb.set_strict_loops( true );
					perturb.set_native_pose( PoseCOP( new core::pose::Pose ( native_pose ) ) );
					perturb.apply( pose );
				}
			core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
			retrieve_sc.apply( pose ); // recover sidechains from pre-centroid pose
			if( refine_ ) {
				protocols::loops::loop_mover::refine::LoopMover_Refine_CCD refine(single_loop, hires_score_ );
				for ( core::Size i = 1; i <= frag_libs.size(); ++i ) {
						refine.add_fragments( frag_libs[i] );
				}
				//core::pack::task::PackerTaskOP task = task_factory->create_task_and_apply_taskoperations( pose );
				refine.set_redesign_loop( false );
				refine.set_native_pose( PoseCOP( new core::pose::Pose ( native_pose ) ) );
				refine.apply( pose );
			}//refine
		}//ccd
	  }//end single loop
	}//loops>0
}
std::string
LoopMoverFromCommandLine::get_name() const {
	return LoopMoverFromCommandLineCreator::mover_name();
}
void
LoopMoverFromCommandLine::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	protocol_ = tag->getOption<std::string>( "protocol", "ccd" );
	perturb_ = tag->getOption<bool>( "perturb", 1 );
	if( protocol_ == "automatic" ) // ugly, but LoopRemodelMover accepts string whereas the other movers accept bool
		refine( tag->getOption< std::string >( "refine", "no" ) );
	else
		refine( tag->getOption<bool>( "refine", 1 ) );
	intermedrelax( tag->getOption< std::string >( "intermedrelax", "no" ) );
	remodel( tag->getOption< std::string >( "remodel", "no" ) );
	relax( tag->getOption< std::string > ("relax", "no" ) );

	hires_score_ = protocols::rosetta_scripts::parse_score_function( tag, "refine_score", data )->clone();
	lores_score_ = protocols::rosetta_scripts::parse_score_function( tag, "perturb_score", data, "score4L" )->clone();

	loop_file_name_ = tag->getOption<std::string>("loop_file", "loops.loops");
//	task_factory(protocols::rosetta_scripts::parse_task_operations( tag, data ));

}//parsemytags
}//movers
}//protein_interface_design
}//protocols
