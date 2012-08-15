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

#include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>

#include <protocols/elscripts/util.hh>


//Auto Headers
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>

namespace protocols{
namespace simple_moves{

static basic::Tracer TR( "protocols.simple_moves.O2M_MutateMover" );

void O2M_MutateMover::apply( core::io::serialization::PipeMap & pmap)
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	if( ! task_factory_) {
		TR << "Task factory not initialized" << std::endl;
		return;
	}
	if( ! scorefxn_) {
		TR << "Score function not initialized" << std::endl;
		return;
	}

	PoseSP starting_pose = (*pmap["input"])[0];
	PackerTaskCOP starting_task = task_factory_->create_task_and_apply_taskoperations( *starting_pose );

  utility::vector1< bool > allowed_aas;
  allowed_aas.assign( num_canonical_aas, false );
	for( core::Size resi = 1; resi <= starting_pose->total_residue(); ++resi ){
		if( starting_task->residue_task( resi ).being_designed() && starting_pose->residue(resi).is_protein() ) {
			std::list<ResidueTypeCAP> const & allowed( starting_task->residue_task( resi ).allowed_residue_types() );
			for( std::list<ResidueTypeCAP>::const_iterator itr=allowed.begin(); itr != allowed.end(); itr++ ){
				if( (*itr)->aa() != starting_pose->residue( resi ).aa() ) {
					allowed_aas[ (*itr)->aa() ] = true;
					PackerTaskOP mutate_task( starting_task->clone() );
					for( core::Size resj = 1; resj <= starting_pose->total_residue(); ++resj ){
						if( resj != resi )
							mutate_task->nonconst_residue_task( resj ).restrict_to_repacking();
						else
							mutate_task->nonconst_residue_task( resj ).restrict_absent_canonical_aas( allowed_aas );
					}
					protocols::simple_moves::PackRotamersMoverOP pack;
					if( core::pose::symmetry::is_symmetric( *starting_pose ) )
						pack =  new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn_, mutate_task );
					else
						pack = new protocols::simple_moves::PackRotamersMover( scorefxn_, mutate_task );

					PoseSP working_pose( new Pose( *starting_pose) );
					pack->apply( *working_pose );
					TR << "Mutated pos " << resi << " from " << starting_pose->residue( resi ).name3() << " to " << working_pose->residue( resi ).name3() << std::endl;
					// although normal mutation notation is A20Q, by doing AQ20 we can use index to get residue type easily
					std::string mut_pos;
					mut_pos += starting_pose->residue( resi ).name1();
					mut_pos += working_pose->residue(resi).name1();
					mut_pos += utility::to_string(resi);
					core::pose::add_comment( *working_pose, "mut_pos", mut_pos );
					pmap["input"]->push_back( working_pose );	
					allowed_aas[ (*itr)->aa() ] = false;
				}
			}
		}
	}
}

void O2M_MutateMover::parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks,
		protocols::moves::MoverCacheSP cache ) {

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
