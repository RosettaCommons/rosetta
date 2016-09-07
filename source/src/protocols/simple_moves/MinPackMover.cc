// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/simple_moves/MinPackMover.cc
/// @brief  Implementation of the MinPackMover class; a wrapper class for invoking core::pack::min_pack
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/simple_moves/MinPackMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/pack/min_pack.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// option key includes

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>


namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;

using basic::Warning;
using basic::t_warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.MinPackMover" );


std::string
MinPackMoverCreator::keyname() const
{
	return MinPackMoverCreator::mover_name();
}

protocols::moves::MoverOP
MinPackMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MinPackMover );
}

std::string
MinPackMoverCreator::mover_name()
{
	return "MinPackMover";
}

MinPackMover::MinPackMover() :
	protocols::moves::Mover("MinPackMover"),
	scorefxn_(/* 0 */),
	task_(/* 0 */),
	task_factory_(/* 0 */),
	off_rotamer_pack_( false )
{
	init();
}

MinPackMover::MinPackMover( std::string const & type_name ) :
	protocols::moves::Mover( type_name ),
	scorefxn_(/* 0 */),
	task_(/* 0 */),
	task_factory_(/* 0 */),
	off_rotamer_pack_( false )
{
	init();
}

// constructors with arguments
MinPackMover::MinPackMover(
	ScoreFunctionCOP scorefxn
) :
	protocols::moves::Mover("MinPackMover"),
	scorefxn_(std::move( scorefxn )),
	task_( /* 0 */ ),
	task_factory_(/* 0 */),
	off_rotamer_pack_( false )
{
	init();
}


MinPackMover::MinPackMover(
	ScoreFunctionCOP scorefxn,
	PackerTaskCOP task
) :
	protocols::moves::Mover("MinPackMover"),
	scorefxn_(std::move( scorefxn )),
	task_(std::move( task )),
	task_factory_(/* 0 */),
	off_rotamer_pack_( false )
{
	init();
}

MinPackMover::~MinPackMover()= default;

void
MinPackMover::init()
{
	nonideal_ = basic::options::option[ basic::options::OptionKeys::optimization::scmin_nonideal ]();
	cartesian_ = basic::options::option[ basic::options::OptionKeys::optimization::scmin_cartesian ]();
}

MinPackMover::MinPackMover( MinPackMover const & other ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( other ),
	off_rotamer_pack_( other.off_rotamer_pack_ )
{
	scorefxn_ = other.score_function();
	task_ = other.task();
	task_factory_ = other.task_factory();
	nonideal_ = other.nonideal_;
	cartesian_ = other.cartesian_;
}

void
MinPackMover::apply( Pose & pose )
{
	if ( scorefxn_ == nullptr ) {
		Warning() << "undefined ScoreFunction -- creating a default one" << std::endl;
		scorefxn_ = get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	}

	core::pack::task::PackerTaskCOP task;
	if ( task_factory_ ) {
		task = task_factory_->create_task_and_apply_taskoperations( pose );
	} else {
		runtime_assert( task_ != nullptr );
		runtime_assert( task_is_valid( pose ) );
		task = task_;
	}

	if ( off_rotamer_pack_ ) {
		core::pack::off_rotamer_pack( pose, *scorefxn_, task );
	} else {
		core::pack::min_pack( pose, *scorefxn_, task, cartesian_, nonideal_ );
	}

}

std::string
MinPackMover::get_name() const {
	return MinPackMoverCreator::mover_name();
}

/// @brief when the PackerTask was not generated locally, verify compatibility with pose
/// @details the pose residue types must be equivalent to the ones used to generate the ResidueLevelTasks, because of the way that prevent_repacking and its associated flags work
bool
MinPackMover::task_is_valid( Pose const & pose ) const
{
	if ( task_->total_residue() != pose.total_residue() ) return false;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		chemical::ResidueTypeCOP r( pose.residue_type(i).get_self_ptr() );
		if ( ! task_->residue_task(i).is_original_type( r ) ) return false;
	}
	return true;
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
MinPackMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->getName() != "MinPackMover" ) {
		TR(t_warning) << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );

	if ( tag->hasOption( "nonideal" ) ) {
		nonideal_ = tag->getOption<bool>( "nonideal" );
	}
	if ( tag->hasOption( "cartesian" ) ) {
		cartesian_ = tag->getOption<bool>( "cartesian" );
	}
	if ( tag->hasOption( "off_rotamer_pack" ) ) {
		off_rotamer_pack_ = tag->getOption<bool>( "off_rotamer_pack" );
	}
}

/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
MinPackMover::parse_score_function(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == nullptr ) return;
	score_function( new_score_function );
}

/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
MinPackMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == nullptr ) return;
	task_factory( new_task_factory );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MinPackMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new MinPackMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MinPackMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::MinPackMover( *this ) );
}


// setters
void MinPackMover::score_function( ScoreFunctionCOP sf )
{
	runtime_assert( sf != nullptr );
	scorefxn_ = sf;
}

void MinPackMover::task( task::PackerTaskCOP t ) { task_ = t; }

void MinPackMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf != nullptr );
	task_factory_ = tf;
}

// accessors
ScoreFunctionCOP MinPackMover::score_function() const { return scorefxn_; }
PackerTaskCOP MinPackMover::task() const { return task_; }
TaskFactoryCOP MinPackMover::task_factory() const { return task_factory_; }

void MinPackMover::off_rotamer_pack( bool setting ) { off_rotamer_pack_ = setting; }
bool MinPackMover::off_rotamer_pack() const { return off_rotamer_pack_; }

} // moves
} // protocols

