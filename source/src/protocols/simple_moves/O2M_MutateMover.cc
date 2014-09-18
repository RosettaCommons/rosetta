// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/O2M_MutateMover.cc
/// @brief Mover that takes a single starting structure and some task ops, and produces every single possible point mutation
/// DOES NOT WORK WITH ROSETTASCRIPTS
/// @author Ken Jung

// Unit headers
#include <protocols/simple_moves/O2M_MutateMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <protocols/elscripts/util.hh>

#include <utility/string_util.hh>

//Auto Headers
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>

namespace protocols{
namespace simple_moves{

static thread_local basic::Tracer TR( "protocols.simple_moves.O2M_MutateMover" );

void O2M_MutateMover::apply( core::io::serialization::PipeMap & pmap)
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;
	using namespace core::conformation;

	if( ! task_factory_) {
		TR << "Task factory not initialized" << std::endl;
		return;
	}
	if( ! scorefxn_) {
		TR << "Score function not initialized" << std::endl;
		return;
	}

	PoseSP starting_pose = (*pmap["input"])[0];
	core::pose::add_comment( *starting_pose, "mut_pos", "AA0" );
	PackerTaskCOP starting_task = task_factory_->create_task_and_apply_taskoperations( *starting_pose );

	for( core::Size resi = 1; resi <= starting_pose->total_residue(); ++resi ){
		if( starting_task->residue_task( resi ).being_designed() && starting_pose->residue(resi).is_protein() ) {
			std::list<ResidueTypeCOP> const & allowed( starting_task->residue_task( resi ).allowed_residue_types() );
			for( std::list<ResidueTypeCOP>::const_iterator itr=allowed.begin(); itr != allowed.end(); itr++ ){
				if( (*itr)->aa() != starting_pose->residue( resi ).aa() ) {
					PoseSP working_pose( new Pose( *starting_pose) );

					PackerTaskOP mutate_task( starting_task->clone() );

					// Create the new residue and replace it
					ResidueOP new_res = ResidueFactory::create_residue(
						**itr, working_pose->residue(resi),
						working_pose->conformation());
					// Make sure we retain as much info from the previous res as possible
					copy_residue_coordinates_and_rebuild_missing_atoms( working_pose->residue(resi),
						*new_res, working_pose->conformation() );
					working_pose->replace_residue(resi, *new_res, false );


					TR << "Mutated pos " << resi << " from " << starting_pose->residue( resi ).name3() << " to " << working_pose->residue( resi ).name3() << std::endl;
					// although normal mutation notation is A20Q, by doing AQ20 we can use index to get residue type easily
					std::string mut_pos;
					mut_pos += starting_pose->residue( resi ).name1();
					mut_pos += working_pose->residue(resi).name1();
					mut_pos += utility::to_string(resi);
					core::pose::add_comment( *working_pose, "mut_pos", mut_pos );
					pmap["input"]->push_back( working_pose );	
				}
			}
		}
	}
}

void O2M_MutateMover::parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks,
		protocols::moves::MoverCacheSP ) {

	if( def["scorefxn"] ) {
		scorefxn_ = protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns );
	} else {
		scorefxn_ = score_fxns["score12"].to<core::scoring::ScoreFunctionSP>()->clone();
	}
	if( def["tasks"] ) {
		task_factory_ = protocols::elscripts::sp_parse_taskdef( def["tasks"], tasks );
	} else {
		task_factory_ = core::pack::task::TaskFactorySP( new core::pack::task::TaskFactory );
	}
}

}//simple_moves
}//protocols
