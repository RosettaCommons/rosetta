// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/TaskAwareMinMover.cc
/// @brief TaskAwareMinMover methods implemented
/// @author Steven Lewis


// Unit Headers
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMoverCreator.hh>

// Package Headers

// Project Headers
#include <core/kinematics/MoveMap.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


// C++ Headers

using basic::Error;
using basic::Warning;


static basic::Tracer TR( "protocols.simple_moves.TaskAwareMinMover" );

namespace protocols {
namespace simple_moves {

// XRW TEMP std::string
// XRW TEMP TaskAwareMinMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return TaskAwareMinMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP TaskAwareMinMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new TaskAwareMinMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP TaskAwareMinMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "TaskAwareMinMover";
// XRW TEMP }

TaskAwareMinMover::TaskAwareMinMover()
: protocols::moves::Mover("TaskAwareMinMover"),
	minmover_(/* nullptr */),
	base_movemap_(/* nullptr */),
	factory_(/* nullptr */),
	chi_(true),
	bb_(false),
	jump_(false)
{}

/// @brief constructor with TaskFactory
TaskAwareMinMover::TaskAwareMinMover(
	protocols::simple_moves::MinMoverOP minmover_in,
	core::pack::task::TaskFactoryCOP factory_in
) : protocols::moves::Mover("TaskAwareMinMover"),
	minmover_(std::move(minmover_in)),
	factory_(std::move(factory_in)),
	chi_(true),
	bb_(false),
	jump_(false)
{
	protocols::moves::Mover::type( "TaskAwareMinMover" );
	if ( minmover_ ) {
		base_movemap_ = minmover_->explicitly_set_movemap(); // Not ideal, as MinMover::movemap() should really be called with a Pose
	}
}

TaskAwareMinMover::~TaskAwareMinMover()= default;

/// @details apply will extract the movemap from your minmover, modify it to include sidechain DOFs that are packable according to some TaskFactory, run the minmover with this movemap, and revert the minmover to its original movemap.
void TaskAwareMinMover::apply( core::pose::Pose & pose ){
	runtime_assert( minmover_ != nullptr );
	runtime_assert( factory_ != nullptr );

	using core::kinematics::MoveMapOP;
	using core::kinematics::MoveMap;

	// Make a default MoveMap
	MoveMapOP mm;
	if ( base_movemap_ ) {
		mm = base_movemap_->clone();
	} else {
		mm = MoveMapOP( new MoveMap );
	}

	//generate task
	using core::pack::task::PackerTaskOP;
	PackerTaskOP task( factory_->create_task_and_apply_taskoperations( pose ) );

	//modify movemap by task
	//  core::kinematics::modify_movemap_from_packertask( *mm, *task );
	Size const nres( task->total_residue() );

	mm->set_jump( jump_ );
	for ( Size i(1); i <= nres; ++i ) {
		if ( task->design_residue( i ) ) {
			// the MoveMap initializes to false for all degrees of freedom
			// this class only turns on minimization at packable dofs, it does not turn them off
			if ( chi_ ) mm->set_chi( i, chi_ );
			if ( bb_ ) mm->set_bb( i, bb_ );
		} else if ( task->pack_residue( i ) ) {
			if ( chi_ ) mm->set_chi( i, chi_ );
		}
	}

	//pass the modified map into the MinMover
	minmover_->movemap( mm );

	//now run MinMover
	minmover_->apply( pose );

	// We don't bother resetting the MinMover's original movemap, as we'll be resetting it later from the stored base

}//apply

// XRW TEMP std::string
// XRW TEMP TaskAwareMinMover::get_name() const {
// XRW TEMP  return TaskAwareMinMover::mover_name();
// XRW TEMP }

protocols::moves::MoverOP TaskAwareMinMover::fresh_instance() const { return protocols::moves::MoverOP( new TaskAwareMinMover ); }
protocols::moves::MoverOP TaskAwareMinMover::clone() const { return protocols::moves::MoverOP( new protocols::simple_moves::TaskAwareMinMover( *this ) ); }

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
TaskAwareMinMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->hasOption("chi") ) chi_ = tag->getOption<bool>("chi");
	if ( tag->hasOption("bb") ) bb_ = tag->getOption<bool>("bb");
	if ( tag->hasOption( "jump" ) ) {
		std::string const jump( tag->getOption< std::string >( "jump" ) );
		if ( jump != "0" && jump != "1" ) {
			utility_exit_with_message( "TaskAwareMinMover only knows how to interpret jump=1(all jumps true) or jump=0 (false). I got jump = "+jump );
		}
	}
	jump_ = tag->getOption< bool >( "jump", false );
	minmover_ = protocols::simple_moves::MinMoverOP( new MinMover );
	// call to MinMover::parse_my_tag avoided here to prevent collision of movemap tag options
	minmover_->parse_opts( tag, datamap, filters, movers );
	parse_task_operations( tag, datamap, filters, movers, pose );
	minmover_->score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap) );
}

/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
TaskAwareMinMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	using namespace core::pack::task;
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == nullptr ) return;
	factory_ = new_task_factory;
}

std::string TaskAwareMinMover::get_name() const {
	return mover_name();
}

std::string TaskAwareMinMover::mover_name() {
	return "TaskAwareMinMover";
}

void TaskAwareMinMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// AMW: Some of below relates to options parsed by MinMover::parse_opts
	// which is called in this parse_my_tag in an effort to avoid reparsing
	// chi and bb.

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "chi", xsct_rosetta_bool, "Allow chi degrees of freedom to minimize" )
		+ XMLSchemaAttribute( "bb", xsct_rosetta_bool, "Allow bb degrees of freedom to minimize" )
		+ XMLSchemaAttribute( "jump", xsct_rosetta_bool, "Allow jump degrees of freedom to minimize" );

	attlist + XMLSchemaAttribute::attribute_w_default( "max_iter", xsct_positive_integer, "Number of iterations", "200" )
		+ XMLSchemaAttribute::attribute_w_default( "type", xsct_minimizer_type, "Minimizer type, chosen from a long list of algorithms", "lbfgs_armijo_nonmonotone" )
		+ XMLSchemaAttribute::attribute_w_default( "tolerance", xsct_real, "Minimization tolerance (absolute)", "0.01" )
		+ XMLSchemaAttribute::attribute_w_default( "cartesian", xsct_rosetta_bool, "Use cartesian minimization (not internal coordinate)", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "bondangle", xsct_rosetta_bool, "Minimize bond angle degrees of freedom", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "bondlength", xsct_rosetta_bool, "Minimize bond length degrees of freedom", "0" );

	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string TaskAwareMinMoverCreator::keyname() const {
	return TaskAwareMinMover::mover_name();
}

protocols::moves::MoverOP
TaskAwareMinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TaskAwareMinMover );
}

void TaskAwareMinMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TaskAwareMinMover::provide_xml_schema( xsd );
}


}//moves
}//protocols

