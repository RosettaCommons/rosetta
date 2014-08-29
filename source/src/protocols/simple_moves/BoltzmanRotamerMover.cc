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
/// @author Noah Ollikainen

// Unit headers
#include <protocols/simple_moves/BoltzmannRotamerMover.hh>
#include <protocols/simple_moves/BoltzmannRotamerMoverCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/rotamer_trials.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>


#include <protocols/filters/Filter.fwd.hh>


// Parser headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

// Utility Headers
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <utility/string_util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>



using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.simple_moves.BoltzmannRotamerMover");

namespace protocols {
namespace simple_moves {

std::string
BoltzmannRotamerMoverCreator::keyname() const {
	return BoltzmannRotamerMoverCreator::mover_name();
}

protocols::moves::MoverOP
BoltzmannRotamerMoverCreator::create_mover() const {
	return new BoltzmannRotamerMover;
}

std::string
BoltzmannRotamerMoverCreator::mover_name() {
	return "BoltzmannRotamerMover";
}

// default constructor
BoltzmannRotamerMover::BoltzmannRotamerMover() : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "BoltzmannRotamer" );
}

// constructor with arguments
BoltzmannRotamerMover::BoltzmannRotamerMover(
	ScoreFunctionCOP scorefxn_in,
	PackerTask & task_in
) : protocols::moves::Mover(), scorefxn_( scorefxn_in ), factory_( NULL ), show_packer_task_( false )
{
	protocols::moves::Mover::type( "BoltzmannRotamer" );
	task_ = task_in.clone();
}

// constructor with arguments
BoltzmannRotamerMover::BoltzmannRotamerMover(
	ScoreFunctionCOP scorefxn_in,
	TaskFactoryCOP factory_in
) : protocols::moves::Mover(), scorefxn_( scorefxn_in ), task_( NULL ), factory_( factory_in ), show_packer_task_( false )
{
	protocols::moves::Mover::type( "BoltzmannRotamer" );
}

// copy constructor
BoltzmannRotamerMover::BoltzmannRotamerMover( BoltzmannRotamerMover const & rval ):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( rval ),
	scorefxn_( rval.scorefxn_ ),
	task_( rval.task_ ),
	factory_( rval.factory_ ),
	show_packer_task_( rval.show_packer_task_ )
{}

// destructor
BoltzmannRotamerMover::~BoltzmannRotamerMover(){}

// clone this object
BoltzmannRotamerMover::MoverOP
BoltzmannRotamerMover::clone() const
{
	return new protocols::simple_moves::BoltzmannRotamerMover( *this );
}

// create this type of object
BoltzmannRotamerMover::MoverOP
BoltzmannRotamerMover::fresh_instance() const
{
	return new protocols::simple_moves::BoltzmannRotamerMover();
}

// setters
void BoltzmannRotamerMover::score_function( core::scoring::ScoreFunctionCOP sf ) { scorefxn_ = sf; }
void BoltzmannRotamerMover::task_factory( core::pack::task::TaskFactoryCOP tf ) { factory_ = tf; }

void
BoltzmannRotamerMover::apply( core::pose::Pose & pose )
{
	//task() contains the call to the TaskFactory
	//TR << *(task(pose)) << std::flush;
	PackerTaskCOP ptask = task( pose );
	if( show_packer_task_ ) {
		TR << *ptask;
	}
	core::pack::rotamer_trials( pose, *scorefxn_, ptask );
}

std::string
BoltzmannRotamerMover::get_name() const {
	return "BoltzmannRotamerMover";
}

void
BoltzmannRotamerMover::show(std::ostream & output) const
{
	Mover::show(output);
	if ( scorefxn() != 0 ) {
		output << "Score function: " << scorefxn()->get_name() << std::endl;
	}
	else { output << "Score function: none" << std::endl; }
}

/// @brief read access for derived classes
BoltzmannRotamerMover::ScoreFunctionCOP
BoltzmannRotamerMover::scorefxn() const
{
	return scorefxn_;
}

/// @brief read access for derived classes
BoltzmannRotamerMover::PackerTaskCOP
BoltzmannRotamerMover::task( core::pose::Pose const & pose ) const
{
	//if we have a factory, generate and return a new task
	if(factory_) return factory_->create_task_and_apply_taskoperations( pose );
	//else runtime_assert( task_is_valid( pose ) );

	//else return the unsafe one
	return task_;
}

/// @brief parse xml
void
BoltzmannRotamerMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{
	using core::scoring::ScoreFunction;
	using core::pack::task::operation::TaskOperation;
	using core::pack::task::TaskFactoryOP;
	using core::pack::task::TaskFactory;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	show_packer_task_ = tag->getOption<bool>( "show_packer_task", 0 );
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
}

std::ostream &operator<< (std::ostream &os, BoltzmannRotamerMover const &mover)
{
	mover.show(os);
	return os;
}

} // moves
} // protocols
