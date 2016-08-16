// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RotamerTrialsMinMover.cc
/// @brief protocols::moves::Mover for Rotamer-Trials with Minimization (based on RotamerTrialsMover)
/// @author Barak Raveh

// Unit headers
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMoverCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pack/rtmin.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.RotamerTrialsMinMover" );

namespace protocols {
namespace simple_moves {

// default constructor
RotamerTrialsMinMover::RotamerTrialsMinMover() : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "RotamerTrialsMin" );
	init();
}

// constructor with arguments
RotamerTrialsMinMover::RotamerTrialsMinMover(
	ScoreFunctionCOP scorefxn_in,
	PackerTask & task_in
) : protocols::moves::Mover(), scorefxn_( scorefxn_in ), factory_( /* NULL */ )
{
	protocols::moves::Mover::type( "RotamerTrialsMin" );
	task_ = task_in.clone();
	init();
}

// constructor with arguments
RotamerTrialsMinMover::RotamerTrialsMinMover(
	ScoreFunctionCOP scorefxn_in,
	TaskFactoryCOP factory_in
) : protocols::moves::Mover(), scorefxn_( scorefxn_in ), task_( /* NULL */ ), factory_( factory_in )
{
	protocols::moves::Mover::type( "RotamerTrialsMin" );
	init();
}

void
RotamerTrialsMinMover::init()
{
	nonideal_ = basic::options::option[ basic::options::OptionKeys::optimization::scmin_nonideal ]();
	cartesian_ = basic::options::option[ basic::options::OptionKeys::optimization::scmin_cartesian ]();
}


RotamerTrialsMinMover::~RotamerTrialsMinMover() {}

// setters
void RotamerTrialsMinMover::score_function( core::scoring::ScoreFunctionCOP sf ) { scorefxn_ = sf; }
void RotamerTrialsMinMover::task_factory( core::pack::task::TaskFactoryCOP tf ) { factory_ = tf; }

void
RotamerTrialsMinMover::apply( core::pose::Pose & pose )
{
	//task() contains the call to the TaskFactory
	//TR << *(task(pose)) << std::flush;
	( *scorefxn_ )(pose); // Ensure scorefunction data is appropriately initialized
	core::pack::RTMin RTMin;
	RTMin.set_nonideal(nonideal_);
	RTMin.set_cartesian(cartesian_);
	RTMin.rtmin( pose, *scorefxn_, task(pose) );
}

std::string
RotamerTrialsMinMover::get_name() const {
	return "RotamerTrialsMinMover";
}

void
RotamerTrialsMinMover::show(std::ostream & output) const
{
	Mover::show(output);
	if ( scorefxn() != 0 ) {
		output << "Score function: " << scorefxn()->get_name() << std::endl;
	} else { output << "Score function: none" << std::endl; }
}

/// @brief read access for derived classes
RotamerTrialsMinMover::ScoreFunctionCOP
RotamerTrialsMinMover::scorefxn() const
{
	return scorefxn_;
}

/// @brief read access for derived classes
RotamerTrialsMinMover::PackerTaskOP
RotamerTrialsMinMover::task( core::pose::Pose const & pose ) const
{
	//if we have a factory, generate and return a new task
	if ( factory_ ) return factory_->create_task_and_apply_taskoperations( pose );
	//else assert( task_is_valid( pose ) );

	//else return the unsafe one
	return task_->clone();
}

/// @brief Parse XML for RosettaScripts
void
RotamerTrialsMinMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	assert( tag->getName() == "RotamerTrialsMinMover" );

	core::scoring::ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == 0 ) {
		TR << "Using default score function for RotamerTrialsMinMover." << std::endl;
		new_score_function = core::scoring::get_score_function();
	}
	score_function( new_score_function );

	core::pack::task::TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0 ) {
		TR << "Using default Task Operations for RotamerTrialsMinMover." << std::endl;
		new_task_factory = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );
	}
	task_factory( new_task_factory );

	if ( tag->hasOption( "nonideal" ) ) {
		nonideal_ = tag->getOption<bool>( "nonideal" );
	}
	if ( tag->hasOption( "cartesian" ) ) {
		cartesian_ = tag->getOption<bool>( "cartesian" );
	}
}

/// @brief Return a new mover instance (for RosettaScripts)
protocols::moves::MoverOP
RotamerTrialsMinMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RotamerTrialsMinMover );
}

/// @brief Return a copy of this mover instance (for RosettaScripts)
protocols::moves::MoverOP
RotamerTrialsMinMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::RotamerTrialsMinMover( *this ) );
}

std::string
RotamerTrialsMinMoverCreator::keyname() const
{
	return RotamerTrialsMinMoverCreator::mover_name();
}

protocols::moves::MoverOP
RotamerTrialsMinMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new RotamerTrialsMinMover );
}

std::string
RotamerTrialsMinMoverCreator::mover_name()
{
	return "RotamerTrialsMinMover";
}

std::ostream &operator<< (std::ostream &os, RotamerTrialsMinMover const &mover)
{
	mover.show(os);
	return os;
}


// default constructor
EnergyCutRotamerTrialsMinMover::EnergyCutRotamerTrialsMinMover() :
	protocols::simple_moves::RotamerTrialsMinMover()
{
	protocols::moves::Mover::type( "EnergyCutRotamerTrialsMin" );
	init();
}

// constructor with arguments
EnergyCutRotamerTrialsMinMover::EnergyCutRotamerTrialsMinMover(
	ScoreFunctionCOP scorefxn_in,
	PackerTask & task_in,
	protocols::moves::MonteCarloOP mc_in,
	core::Real energycut_in
) : protocols::simple_moves::RotamerTrialsMinMover(scorefxn_in, task_in), mc_( mc_in ), energycut_( energycut_in )
{
	protocols::moves::Mover::type( "EnergyCutRotamerTrialsMin" );
	init();
}

// constructor with arguments
EnergyCutRotamerTrialsMinMover::EnergyCutRotamerTrialsMinMover(
	ScoreFunctionCOP scorefxn_in,
	TaskFactoryCOP factory_in,
	protocols::moves::MonteCarloOP mc_in,
	core::Real energycut_in
) : protocols::simple_moves::RotamerTrialsMinMover(scorefxn_in, factory_in), mc_( mc_in ), energycut_( energycut_in )
{
	protocols::moves::Mover::type( "EnergyCutRotamerTrialsMin" );
	init();
}

EnergyCutRotamerTrialsMinMover::~EnergyCutRotamerTrialsMinMover() {}

void
EnergyCutRotamerTrialsMinMover::apply( core::pose::Pose & pose )
{
	PackerTaskOP rottrial_task( task(pose) );
	( *scorefxn() )(pose);
	setup_energycut_task( pose, *mc_, *rottrial_task );
	/// This call is dangerous.  If sequence or length has changed since task was created, it will crash.
	/// Not a problem if you used a TaskFactory
	core::pack::RTMin RTMin;
	RTMin.rtmin( pose, *scorefxn(), rottrial_task );
}


std::string
EnergyCutRotamerTrialsMinMover::get_name() const {
	return "EnergyCutRotamerTrialsMinMover";
}

/// @details starting from a fresh task, it reduces the number of residues to be repacked to only
/// those whose energy has increased by energycut_ since the application of the last move.
void
EnergyCutRotamerTrialsMinMover::setup_energycut_task(
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
		core::Real const resE ( pose.energies().residue_total_energies(i)[ total_score ] );
		core::Real const lowest_resE( mc.lowest_score_pose().energies().residue_total_energies(i)[ total_score ] );
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

} // moves
} // protocols
