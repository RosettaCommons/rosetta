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
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov

// Unit headers
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMoverCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_trials.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>


#include <protocols/filters/Filter.fwd.hh>


// Parser headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.simple_moves.RotamerTrialsMover" );

namespace protocols {
namespace simple_moves {

std::string
RotamerTrialsMoverCreator::keyname() const {
	return RotamerTrialsMoverCreator::mover_name();
}

protocols::moves::MoverOP
RotamerTrialsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RotamerTrialsMover );
}

std::string
RotamerTrialsMoverCreator::mover_name() {
	return "RotamerTrialsMover";
}

// default constructor
RotamerTrialsMover::RotamerTrialsMover() : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "RotamerTrials" );
}

// constructor with arguments
RotamerTrialsMover::RotamerTrialsMover(
	ScoreFunctionCOP scorefxn_in,
	PackerTask & task_in
) : protocols::moves::Mover(), scorefxn_( scorefxn_in ), factory_( /* NULL */ ), show_packer_task_( false )
{
	protocols::moves::Mover::type( "RotamerTrials" );
	task_ = task_in.clone();
}

// constructor with arguments
RotamerTrialsMover::RotamerTrialsMover(
	ScoreFunctionCOP scorefxn_in,
	TaskFactoryCOP factory_in
) : protocols::moves::Mover(), scorefxn_( scorefxn_in ), task_( /* NULL */ ), factory_( factory_in ), show_packer_task_( false )
{
	protocols::moves::Mover::type( "RotamerTrials" );
}

// copy constructor
RotamerTrialsMover::RotamerTrialsMover( RotamerTrialsMover const & rval ):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( rval ),
	scorefxn_( rval.scorefxn_ ),
	task_( rval.task_ ),
	factory_( rval.factory_ ),
	show_packer_task_( rval.show_packer_task_ )
{}

// destructor
RotamerTrialsMover::~RotamerTrialsMover(){}

// clone this object
RotamerTrialsMover::MoverOP
RotamerTrialsMover::clone() const
{
	return RotamerTrialsMover::MoverOP( new protocols::simple_moves::RotamerTrialsMover( *this ) );
}

// create this type of object
RotamerTrialsMover::MoverOP
RotamerTrialsMover::fresh_instance() const
{
	return RotamerTrialsMover::MoverOP( new protocols::simple_moves::RotamerTrialsMover() );
}

// setters
void RotamerTrialsMover::score_function( core::scoring::ScoreFunctionCOP sf ) { scorefxn_ = sf; }
void RotamerTrialsMover::task_factory( core::pack::task::TaskFactoryCOP tf ) { factory_ = tf; }

void
RotamerTrialsMover::apply( core::pose::Pose & pose )
{
	//task() contains the call to the TaskFactory
	//TR << *(task(pose)) << std::flush;
	PackerTaskCOP ptask = task( pose );
	if ( show_packer_task_ ) {
		TR << *ptask;
	}
	core::pack::rotamer_trials( pose, *scorefxn_, ptask );
}

std::string
RotamerTrialsMover::get_name() const {
	return "RotamerTrialsMover";
}

void
RotamerTrialsMover::show(std::ostream & output) const
{
	Mover::show(output);
	if ( scorefxn() != 0 ) {
		output << "Score function: " << scorefxn()->get_name() << std::endl;
	} else { output << "Score function: none" << std::endl; }
}

/// @brief read access for derived classes
RotamerTrialsMover::ScoreFunctionCOP
RotamerTrialsMover::scorefxn() const
{
	return scorefxn_;
}

/// @brief read access for derived classes
RotamerTrialsMover::PackerTaskCOP
RotamerTrialsMover::task( core::pose::Pose const & pose ) const
{
	//if we have a factory, generate and return a new task
	if ( factory_ ) return factory_->create_task_and_apply_taskoperations( pose );
	//else runtime_assert( task_is_valid( pose ) );

	//else return the unsafe one
	return task_;
}

/// @brief parse xml
void
RotamerTrialsMover::parse_my_tag(
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

std::ostream &operator<< (std::ostream &os, RotamerTrialsMover const &mover)
{
	mover.show(os);
	return os;
}


// default constructor
EnergyCutRotamerTrialsMover::EnergyCutRotamerTrialsMover() :
	protocols::simple_moves::RotamerTrialsMover()
{
	protocols::moves::Mover::type( "EnergyCutRotamerTrials" );
}

// constructor with arguments
EnergyCutRotamerTrialsMover::EnergyCutRotamerTrialsMover(
	ScoreFunctionCOP scorefxn_in,
	PackerTask & task_in,
	protocols::moves::MonteCarloOP mc_in,
	core::Real energycut_in
) : protocols::simple_moves::RotamerTrialsMover(scorefxn_in, task_in), mc_( mc_in ), energycut_( energycut_in )
{
	protocols::moves::Mover::type( "EnergyCutRotamerTrials" );
}

// constructor with arguments
EnergyCutRotamerTrialsMover::EnergyCutRotamerTrialsMover(
	ScoreFunctionCOP scorefxn_in,
	TaskFactoryCOP factory_in,
	protocols::moves::MonteCarloOP mc_in,
	core::Real energycut_in
) : protocols::simple_moves::RotamerTrialsMover(scorefxn_in, factory_in), mc_( mc_in ), energycut_( energycut_in )
{
	protocols::moves::Mover::type( "EnergyCutRotamerTrials" );
}

EnergyCutRotamerTrialsMover::~EnergyCutRotamerTrialsMover() {}

void
EnergyCutRotamerTrialsMover::apply( core::pose::Pose & pose )
{
	PackerTaskOP rottrial_task( task(pose)->clone() );
	( *scorefxn() )(pose);
	/// Now handled automatically.  scorefxn()->accumulate_residue_total_energies( pose );
	setup_energycut_task( pose, *mc_, *rottrial_task );
	/// This call is dangerous.  If sequence or length has changed since task was created, it will crash.
	/// Not a problem if you used a TaskFactory
	core::pack::rotamer_trials( pose, *scorefxn(), rottrial_task );
}

std::string
EnergyCutRotamerTrialsMover::get_name() const {
	return "EnergyCutRotamerTrialsMover";
}

/// @details starting from a fresh task, it reduces the number of residues to be repacked to only
/// those whose energy has increased by energycut_ since the application of the last move.
void
EnergyCutRotamerTrialsMover::setup_energycut_task(
	core::pose::Pose const & pose,
	protocols::moves::MonteCarlo const & mc,
	core::pack::task::PackerTask & task_in
) const
{
	using namespace core;
	using core::scoring::total_score;

	//Size count_fixed( 0 ), count_repacked( 0 );

	task_in.restrict_to_repacking();

	for ( int i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
		core::Real const resE ( pose.energies().residue_total_energy(i) );
		core::Real const lowest_resE( mc.lowest_score_pose().energies().residue_total_energy(i) );
		core::Real const deltaE ( resE - lowest_resE );
		if ( deltaE < energycut_ ) {
			task_in.nonconst_residue_task(i).prevent_repacking();
			//++count_fixed;
		} else {
			// let this residue be repacked
			//++count_repacked;
		}
	}
}
// protected accessor function for derived mover
protocols::moves::MonteCarloOP
EnergyCutRotamerTrialsMover::mc()
{
	return mc_;
}

} // moves
} // protocols
