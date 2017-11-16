// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_mover/refine/LoopRefineInnerCycleContainer.cc
/// @brief This class is a LoopRefineInnerCycle that contains one or more other LoopRefineInnerCycles to allow a developer to
/// quickly string together existing LoopRefineInnerCycles in new ways to create new loop refinement protocols.
/// @details
///
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )

// Unit headers
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleContainer.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleContainerCreator.hh>

// Package headers

// Project headers
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.loops.loop_mover.refine.LoopRefineInnerCycleContainer" );

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

using utility::excn::EXCN_Msg_Exception;

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
LoopRefineInnerCycleContainer::LoopRefineInnerCycleContainer() : LoopRefineInnerCycle()
{
	init();
}

/// @brief copy constructor
LoopRefineInnerCycleContainer::LoopRefineInnerCycleContainer( LoopRefineInnerCycleContainer const & rhs ) : LoopRefineInnerCycle(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

/// @brief assignment operator
LoopRefineInnerCycleContainer & LoopRefineInnerCycleContainer::operator=( LoopRefineInnerCycleContainer const & rhs ){
	//abort self-assignment
	if ( this == &rhs ) return *this;
	LoopRefineInnerCycle::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

//destructor
LoopRefineInnerCycleContainer::~LoopRefineInnerCycleContainer() {}

/// @brief Each derived class must specify its name.
// XRW TEMP std::string LoopRefineInnerCycleContainer::get_name() const
// XRW TEMP {
// XRW TEMP  return type();
// XRW TEMP }

//@brief clone operator, calls the copy constructor
protocols::moves::MoverOP
LoopRefineInnerCycleContainer::clone() const
{
	return protocols::moves::MoverOP( new LoopRefineInnerCycleContainer( *this ) );
}

/// @brief fresh_instance returns a default-constructed object for JD2
protocols::moves::MoverOP
LoopRefineInnerCycleContainer::fresh_instance() const
{
	return protocols::moves::MoverOP( new LoopRefineInnerCycleContainer() );
}

/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
bool LoopRefineInnerCycleContainer::reinitialize_for_new_input() const
{
	return true;
}

void LoopRefineInnerCycleContainer::register_options()
{
	///  PUT THE LIST OF OPTIONS THAT ARE USED HERE  ///

	///  RECURSIVELY CALL REGISTER OPTIONS ON ALL MOVERS THAT THIS CLASS HAS AN OWNING_PTR TO  ///
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

void LoopRefineInnerCycleContainer::apply( Pose & pose )
{
	// TR << "Beginning apply function of " + get_name() + "." << std::endl;

	setup_objects( pose );

	// show( TR );
	// TR << "Applying each refinement step in this LoopInnerRefineCycleContainer..." << std::endl;

	inner_cycle_steps_->apply( pose );
}

void LoopRefineInnerCycleContainer::setup_objects( Pose const & pose )
{
	// TR << "Setting up data for " + get_name() + "." << std::endl;

	LoopRefineInnerCycle::setup_objects( pose );
}

void LoopRefineInnerCycleContainer::init()
{
	type( "LoopRefineInnerCycleContainer" );
	inner_cycle_list_  = InnerCycleList();
	inner_cycle_steps_ = moves::SequenceMoverOP( new moves::SequenceMover );

	init_options();
}

void LoopRefineInnerCycleContainer::init_for_equal_operator_and_copy_constructor(
	LoopRefineInnerCycleContainer & lhs,
	LoopRefineInnerCycleContainer const & rhs
)
{
	// copy all data members from rhs to lhs
	lhs.inner_cycle_list_ = rhs.inner_cycle_list_;
	lhs.inner_cycle_steps_ = rhs.inner_cycle_steps_;
}

void LoopRefineInnerCycleContainer::init_options()
{
	/* UNCOMMENT WHEN THERE ARE ACTUALLY OPTIONS TO PROCESS
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	*/
	// Set options here.
}

void LoopRefineInnerCycleContainer::add_inner_cycle_step( LoopRefineInnerCycleOP inner_cycle_step )
{
	inner_cycle_list_.push_back( inner_cycle_step );
	inner_cycle_steps_->add_mover( inner_cycle_step );
}

// This could probably be simplified if I passed around function pointers, but who wants to deal with all that bullshit.
void LoopRefineInnerCycleContainer::set_mc( moves::MonteCarloOP mc)
{
	LoopRefineInnerCycle::set_mc( mc );

	for ( InnerCycleList::iterator it = inner_cycle_list_.begin(); it != inner_cycle_list_.end(); ++it ) {
		(*it)->set_mc( mc );
	}
}

void LoopRefineInnerCycleContainer::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn)
{
	LoopRefineInnerCycle::set_scorefxn( scorefxn );

	for ( InnerCycleList::iterator it = inner_cycle_list_.begin(); it != inner_cycle_list_.end(); ++it ) {
		(*it)->set_scorefxn( scorefxn );
	}
}

void LoopRefineInnerCycleContainer::set_task_factory( core::pack::task::TaskFactoryOP tf )
{
	LoopRefineInnerCycle::set_task_factory( tf );

	for ( InnerCycleList::iterator it = inner_cycle_list_.begin(); it != inner_cycle_list_.end(); ++it ) {
		(*it)->set_task_factory( tf );
	}
}

void LoopRefineInnerCycleContainer::set_loop_mover( LoopMover_Refine_CCDAP new_owner_in_town )
{
	LoopRefineInnerCycle::set_loop_mover( new_owner_in_town );

	for ( InnerCycleList::iterator it = inner_cycle_list_.begin(); it != inner_cycle_list_.end(); ++it ) {
		(*it)->set_loop_mover( new_owner_in_town );
	}
}

void LoopRefineInnerCycleContainer::set_native_pose( PoseCOP pose )
{
	moves::Mover::set_native_pose( pose );

	for ( InnerCycleList::iterator it = inner_cycle_list_.begin(); it != inner_cycle_list_.end(); ++it ) {
		(*it)->set_native_pose( pose );
	}
}
void
LoopRefineInnerCycleContainer::show( std::ostream & out ) const
{
	out << *this;
}

std::ostream & operator<<(std::ostream& out, LoopRefineInnerCycleContainer const & loop_refine_inner_cycle_container )
{
	out << loop_refine_inner_cycle_container.get_name() << " with the following LoopRefineInnerCycles: " << std::endl;
	// Iterate over contained LoopRefineInnerCycles and print their names
	LoopRefineInnerCycleContainer::InnerCycleList inner_cycle_list = loop_refine_inner_cycle_container.inner_cycle_list_;

	for ( LoopRefineInnerCycleContainer::InnerCycleList::const_iterator it = inner_cycle_list.begin();
			it != inner_cycle_list.end(); ++it ) {
		out << "    "  << **it << std::endl;
	}

	return out;
}

// XRW TEMP LoopRefineInnerCycleContainerCreator::~LoopRefineInnerCycleContainerCreator() {}

// XRW TEMP moves::MoverOP LoopRefineInnerCycleContainerCreator::create_mover() const {
// XRW TEMP  return moves::MoverOP( new LoopRefineInnerCycleContainer() );
// XRW TEMP }

// XRW TEMP std::string LoopRefineInnerCycleContainerCreator::keyname() const {
// XRW TEMP  return "LoopRefineInnerCycleContainer";
// XRW TEMP }

std::string LoopRefineInnerCycleContainer::get_name() const {
	return mover_name();
}

std::string LoopRefineInnerCycleContainer::mover_name() {
	return "LoopRefineInnerCycleContainer";
}

void LoopRefineInnerCycleContainer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"This class is a LoopRefineInnerCycle that contains one or more other LoopRefineInnerCycles"
		" to allow a developer to quickly string together existing LoopRefineInnerCycles"
		" in new ways to create new loop refinement protocols."
		, attlist );
}

std::string LoopRefineInnerCycleContainerCreator::keyname() const {
	return LoopRefineInnerCycleContainer::mover_name();
}

protocols::moves::MoverOP
LoopRefineInnerCycleContainerCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopRefineInnerCycleContainer );
}

void LoopRefineInnerCycleContainerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopRefineInnerCycleContainer::provide_xml_schema( xsd );
}


} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
