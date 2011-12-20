// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/Splice.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/Splice.hh>
#include <protocols/protein_interface_design/movers/SpliceCreator.hh>
#include <utility/string_util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/ScoreFunction.hh>
// Package headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>
#include <protocols/loops/LoopMover_CCD.hh>
#include <numeric/xyzVector.hh>
#include <protocols/loops/FoldTreeFromLoopsWrapper.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

static basic::Tracer TR( "protocols.protein_interface_design.movers.Splice" );

std::string
SpliceCreator::keyname() const
{
	return SpliceCreator::mover_name();
}

protocols::moves::MoverOP
SpliceCreator::create_mover() const {
	return new Splice;
}

std::string
SpliceCreator::mover_name()
{
	return "Splice";
}

Splice::Splice() :
	Mover( SpliceCreator::mover_name() ),
	from_res_( 0 ),
	to_res_( 0 ),
	source_pdb_( "" ),
	ccd_( true ),
	rms_cutoff_( 999999 ),
	res_move_( 4 )
{
}


Splice::~Splice() {}


/// @brief Return the number of the residue on source that is nearest to res on target. If the distance
/// is greater than 2.0 returns 0 to indicate error
core::Size
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res ){
	core::Real min_dist( 100000 ); core::Size nearest_res( 0 );
	for( core::Size i = 1; i <= source.total_residue(); ++i ){
	  core::Real const dist( target.residue( res ).xyz( "CA" ).distance( source.residue( i ).xyz( "CA" ) ) );
		if( dist <= min_dist ){
			min_dist = dist;
			nearest_res = i;
		}
	}
	if( min_dist <= 2.0 ) return nearest_res;
	else return 0;
}

void
Splice::apply( core::pose::Pose & pose )
{
	core::pose::Pose source_pose;
	core::import_pose::pose_from_pdb( source_pose, source_pdb_ );

	core::Size const nearest_to_from( find_nearest_res( source_pose, pose, from_res() ) );
	core::Size const nearest_to_to( find_nearest_res( source_pose, pose, to_res() ) );

	if( nearest_to_from == 0 || nearest_to_to == 0 ){
		TR<<"nearest_to_from: "<<nearest_to_from<<" nearest_to_to: "<<nearest_to_to<<". Failing"<<std::endl;
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
		return;
	}
/// Don't do anything with loops that contain disulfides. At one point, might be a good idea to copy these disulfides onto new pose.
/// Keep Gly/Pro residues as those have torsion angle strangeness
	std::string threaded_seq( "" );
	for( core::Size i = nearest_to_from; i <= nearest_to_to; ++i ){
		if( source_pose.residue( i ).name3() == "CYD" ){
			TR<<"Residue "<<i<<" in "<<source_pdb()<<" is a disulfide. Failing"<<std::endl;
			set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
			return;
		}
		if( source_pose.residue( i ).name3() == "GLY" )
			threaded_seq += "G";
		else if( source_pose.residue( i ).name3() == "PRO" )
			threaded_seq += "P";
		else{
			core::Size const nearest_on_target( find_nearest_res( pose, source_pose, i ) );
			if( nearest_on_target > 0 && source_pose.residue( i ).name3() == pose.residue( nearest_on_target ).name3() )
				threaded_seq += source_pose.residue(i).name1();
			else
				threaded_seq += "A";
		}
	}

/// make fold tree compatible with the loop (starts and ends 6 residue away from the start points, cuts at loop terminus
	protocols::loops::FoldTreeFromLoops ffl;
	using namespace utility;
	protocols::loops::Loop loop( from_res() - 6/*start*/, to_res() + 6/*stop*/, to_res()/*cut*/ );
	protocols::loops::LoopsOP loops = new protocols::loops::Loops();
	loops->push_back( loop );
	ffl.loops( loops );
	ffl.apply( pose );
	core::Size const residue_diff( nearest_to_to - nearest_to_from - ( to_res() - from_res()) );
/// change the loop length
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	llc.loop_start( from_res() );
	llc.loop_end( to_res());
	llc.delta( residue_diff );
	llc.apply( pose );
/// set torsions
	core::Size const total_residue_new( nearest_to_to - nearest_to_from + 1 );
	for( core::Size i = 0; i < total_residue_new; ++i ){
		pose.set_phi( from_res() + i, source_pose.phi( nearest_to_from + i ) );
		pose.set_psi( from_res() + i, source_pose.psi( nearest_to_from + i ) );
		pose.set_omega( from_res() + i, source_pose.omega( nearest_to_from + i ) );
	}
	if( ccd() ){
/// Set ccd to minimize 4 residues at each loop terminus including the first residue of the loop. This way,
/// the torsion in the loop are maintained. Allow repacking around the loop.
/// If disulfide occurs in the range that is allowed to minimize, adjust that region to not include disulf
		core::scoring::ScoreFunctionOP scorefxn_local( scorefxn()->clone() );
//		scorefxn_local->set_weight( core::scoring::sheet, 5.0 );
		protocols::loops::LoopMover_Refine_CCD ccd_mover( loops, scorefxn_local );
		ccd_mover.temp_initial( 1.5 );
		ccd_mover.temp_final( 0.5 );
		core::kinematics::MoveMapOP mm;
		using namespace protocols::toolbox::task_operations;
		DesignAroundOperationOP dao = new DesignAroundOperation;
		ThreadSequenceOperationOP tso = new ThreadSequenceOperation;
		dao->design_shell( 5 ); // threaded sequence operation needs to design, and will restrict design to the loop only
		dao->repack_shell( 8.0 );
		tso->target_sequence( threaded_seq );
		tso->start_res( from_res() );
		tso->allow_design_around( false );
		TR<<"Threading sequence: "<<threaded_seq<<" starting from "<<from_res()<<std::endl;

		mm = new core::kinematics::MoveMap;
		mm->set_chi( false ); mm->set_bb( false ); mm->set_jump( false );
	/// First look for disulfides. Those should never be moved.
		core::Size disulfn( 0 ), disulfc( 0 );
		for( core::Size i = from_res() - 3; i <= from_res(); ++i ){
			if( pose.residue( i ).name3() == "CYD" ){
				disulfn = i;
			}
		}
		for( core::Size i = from_res() + total_residue_new - 1; i <= from_res() + total_residue_new + 2; ++i ){
			if( pose.residue( i ).name3() == "CYD" ){
				disulfc = i;
				break;
			}
		}
		core::Size const startn( disulfn > 0 ? disulfn + 1 : from_res() - 3 );
		core::Size const startc( disulfc > 0 ? disulfc - 6 : from_res() + total_residue_new - ( res_move() - 3 ) );
		for( core::Size i = startn; i <= startn + res_move() - 1; ++i ){
			mm->set_chi( i, true );
			mm->set_bb( i, true );
		}
		for( core::Size i = startc; i <= startc + res_move() - 1; ++i ){
			mm->set_chi( i, true );
			mm->set_bb( i, true );
		}
		for( core::Size i = from_res() - 3; i <= from_res() + total_residue_new + 2; ++i ){
			if( pose.residue( i ).name3() != "CYD" )
				dao->include_residue( i );
		}
	  core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory;
		tf->push_back( dao );
		tf->push_back( new core::pack::task::operation::InitializeFromCommandline );
		tf->push_back( tso );
		ccd_mover.set_task_factory( tf );
		ccd_mover.move_map( mm );
		ccd_mover.apply( pose );
	}
	core::Real rms( 0 );
	for( core::Size i = 0; i <= total_residue_new - 1; ++i ){
		core::Real const dist( pose.residue( from_res() + i ).xyz( "CA" ).distance( source_pose.residue( nearest_to_from+ i ).xyz("CA" ) ) );
		rms += dist;
	}
	core::Real const average_rms( rms / total_residue_new );
	TR<<"Average distance of spliced segment to original: "<< average_rms<<std::endl;
	if( average_rms >= rms_cutoff() ){
		TR<<"Failing."<<std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}

	protocols::protein_interface_design::movers::AddChainBreak acb;
	acb.resnum( utility::to_string( from_res() + total_residue_new - 1 ));
	acb.find_automatically( false );
	acb.change_foldtree( false );
	acb.apply( pose );
	TR<<"Adding chainbreak at: "<<from_res() + total_residue_new - 1<<std::endl;
}

std::string
Splice::get_name() const {
	return SpliceCreator::mover_name();
}

void
Splice::parse_my_tag( TagPtr const tag, protocols::moves::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{
	from_res( protocols::rosetta_scripts::parse_resnum( tag->getOption< std::string >( "from_res" ), pose ) );
	to_res( protocols::rosetta_scripts::parse_resnum( tag->getOption< std::string >( "to_res" ), pose ) );
	source_pdb( tag->getOption< std::string >( "source_pdb" ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	ccd( tag->getOption< bool >( "ccd", 1 ) );
	rms_cutoff( tag->getOption< core::Real >( "rms_cutoff", 999999 ) );
	res_move( tag->getOption< core::Size >( "res_move", 4 ) );
	TR<<"from_res: "<<from_res()<<" to_res: "<<to_res()<<" source_pdb: "<<source_pdb()<<" ccd: "<<ccd()<<" rms_cutoff: "<<rms_cutoff()<<" res_move: "<<res_move()<<std::endl;
}

protocols::moves::MoverOP
Splice::clone() const {
    return( protocols::moves::MoverOP( new Splice( *this ) ));
}

void
Splice::scorefxn( core::scoring::ScoreFunctionOP sf ){
	scorefxn_ = sf;
}

core::scoring::ScoreFunctionOP
Splice::scorefxn() const{
	return scorefxn_;
}

} //movers
} //protein_interface_design
} //protocols
