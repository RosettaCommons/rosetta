// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief 
/// @author Neil King ( neilking@uw.edu )

// Unit headers
#include <protocols/simple_moves/symmetry/TaskAwareSymMinMover.hh>
#include <protocols/simple_moves/symmetry/TaskAwareSymMinMoverCreator.hh>

// project headers
#include <core/kinematics/MoveMap.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>

#include <core/chemical/ResidueConnection.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/ScoreFunction.hh>

#include <basic/datacache/DataMap.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

using basic::Warning;
static thread_local basic::Tracer TR( "protocols.simple_moves.symmetry.TaskAwareSymMinMover" );

namespace protocols {
namespace simple_moves {
namespace symmetry {

using namespace core;
using namespace utility;

// -------------  Mover Creator -------------
std::string
TaskAwareSymMinMoverCreator::keyname() const
{
	return TaskAwareSymMinMoverCreator::mover_name();
}

protocols::moves::MoverOP
TaskAwareSymMinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TaskAwareSymMinMover );
}

std::string
TaskAwareSymMinMoverCreator::mover_name()
{
	return "TaskAwareSymMinMover";
}
// -------------  Mover Creator -------------

TaskAwareSymMinMover::TaskAwareSymMinMover() :
  min_chi_(false),
  min_bb_(false),
  min_rb_(false),
  min_type_("dfpmin_armijo_nonmonotone"),
  tolerance_(1e-5),
	minmover_(/* 0 */),
	factory_(/* 0 */),
	designable_only_(true)
{ }
	
// constructor with TaskFactory
TaskAwareSymMinMover::TaskAwareSymMinMover(
	protocols::simple_moves::MinMoverOP minmover_in,
	core::pack::task::TaskFactoryCOP factory_in
) :
	min_chi_(false),
	min_bb_(false),
	min_rb_(false),
	min_type_("dfpmin_armijo_nonmonotone"),
	tolerance_(1e-5),
	minmover_(minmover_in),
	factory_(factory_in),
	designable_only_(true)
	{}
	//protocols::simple_moves::TaskAwareMinMover(
	//	minmover_in, factory_in) {}

TaskAwareSymMinMover::TaskAwareSymMinMover(const TaskAwareSymMinMover& rval) :
//	protocols::moves::Mover(),
	protocols::simple_moves::TaskAwareMinMover(), //Jeliazkov experimenat
  scorefxn_( rval.scorefxn_ ),
  min_chi_( rval.min_chi_),
  min_bb_( rval.min_bb_),
  min_rb_( rval.min_rb_),
  min_type_( rval.min_type_),
  tolerance_( rval.tolerance_),
	minmover_( rval.minmover_ ),
  factory_( rval.factory_),
	designable_only_( rval.designable_only_ )
{ }

protocols::moves::MoverOP 
TaskAwareSymMinMover::clone() const {
	return protocols::moves::MoverOP( new TaskAwareSymMinMover( *this ) );
}

protocols::moves::MoverOP 
TaskAwareSymMinMover::fresh_instance() const {
	return protocols::moves::MoverOP( new TaskAwareSymMinMover() );
}


void
TaskAwareSymMinMover::apply(Pose & pose) {

	runtime_assert( factory_ != 0 );

	// Initialize a MoveMap, set all moves to false
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_jump(false);
	movemap->set_bb(false);
	movemap->set_chi(false);

	// Modify MoveMap to set allowable move types at positions defined by TaskOperations
	core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	core::Size nres_monomer = sym_info->num_independent_residues();
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
  if ( factory_ != 0 ) {
    task = factory_->create_task_and_apply_taskoperations( pose );
  } else {
		TR.Warning << "Warning: You have not provided any TaskOperations. A default will be used." << std::endl;
	}
	for (core::Size i = 1; i <= nres_monomer; i++) {
		if ( designable_only_ ) {
			if ( task->design_residue( i ) ) {
				movemap->set_bb (i, min_bb_);
				movemap->set_chi(i, min_chi_);
			}
		} else {
			if ( task->pack_residue( i ) ) {
				movemap->set_bb (i, min_bb_);
				movemap->set_chi(i, min_chi_);
			}
		}
	}

	// Make MoveMap symmetric, apply it to minimize the pose
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	protocols::simple_moves::symmetry::SymMinMover m1( movemap, scorefxn_, min_type_, tolerance_, true, false, false );
	m1.apply(pose);

	// If rigid body minimization is desired, minimize again with that
  if (min_rb_) {
    movemap->set_jump(true);
    core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
    protocols::simple_moves::symmetry::SymMinMover m2( movemap, scorefxn_, min_type_, tolerance_, true, false, false );
    m2.apply(pose);
  }

}

void 
TaskAwareSymMinMover::parse_my_tag( utility::tag::TagCOP tag,
										 basic::datacache::DataMap & data,
										 protocols::filters::Filters_map const &,
										 protocols::moves::Movers_map const &,
										 core::pose::Pose const & ) {

	scorefxn_name_ = tag->getOption< std::string >( "scorefxn", "score12_symm" );
	min_chi_ = tag->getOption< bool >( "chi", true );
	min_bb_ = tag->getOption< bool >( "bb", false );
	min_rb_ = tag->getOption< bool >( "rb", false );
  min_type_ = tag->getOption< std::string >( "type", "dfpmin_armijo_nonmonotone" );
  tolerance_ = tag->getOption< core::Real >( "tolerance", 1e-5 );
	designable_only_ = tag->getOption< bool >( "designable_only", true );
	// Get the ScoreFunction and TaskOperations from the basic::datacache::DataMap
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data, scorefxn_name_ );
	factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );

}

void TaskAwareSymMinMover::parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks,
		protocols::moves::MoverCacheSP ){
	scorefxn_name_ = def["scorefxn"] ? def["scorefxn"].to<std::string>() : "score12_symm";
	min_chi_ = def["chi"] ? def["chi"].to<bool>() : true;
	min_bb_ = def["bb"] ? def["bb"].to<bool>() : false;
	min_bb_ = def["rb"] ? def["rb"].to<bool>() : false;
	min_type_ = def["type"] ? def["type"].to<std::string>() : "dfpmin_armijo_nonmonotone";
	tolerance_ = def["tolerance"] ? def["tolerance"].to<core::Real>() : 1e-5;
	// Get the ScoreFunction and TaskOperations from the basic::datacache::DataMap
	scorefxn_ = protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns );
	factory_ = protocols::elscripts::parse_taskdef( def["tasks"], tasks );
}


} // symmetry
} // simple_moves
} // protocols
