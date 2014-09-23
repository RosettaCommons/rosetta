// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Assigns a ConstraintSet to a pose. Reads and creats ConstraintSet from file via command line option -constraints::cst_file, unless a ConstraintSet is supplied via the constructor or the constraint_set() method.
/// @author ashworth

#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/ConstraintSetMoverCreator.hh>

// AUTO-REMOVED #include <basic/datacache/DataMap.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.simple_moves.ConstraintSetMover" );

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace scoring;
using namespace constraints;
using namespace utility::tag;

std::string
ConstraintSetMoverCreator::keyname() const
{
	return ConstraintSetMoverCreator::mover_name();
}

protocols::moves::MoverOP
ConstraintSetMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ConstraintSetMover );
}

std::string
ConstraintSetMoverCreator::mover_name()
{
	return "ConstraintSetMover";
}

ConstraintSetMover::ConstraintSetMover()
	: protocols::moves::Mover( ConstraintSetMoverCreator::mover_name() )
{
	read_options();
}

ConstraintSetMover::~ConstraintSetMover(){}

ConstraintSetMover::ConstraintSetMover( std::string const & type )
	: protocols::moves::Mover(type)
{
	read_options();
}

void
ConstraintSetMover::read_options()
{
	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		cst_file_ = option[ OptionKeys::constraints::cst_file ]().front();
	}

	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		cst_fa_file_ = option[ OptionKeys::constraints::cst_fa_file ]().front();
	} else {
		cst_fa_file_ = cst_file_;
	}

	add_constraints_ = false; // This is never used, lol.
}

void
ConstraintSetMover::constraint_set( ConstraintSetCOP cst_set )
{
	constraint_set_low_res_ = ConstraintSetOP( new ConstraintSet( *cst_set ) );
	constraint_set_high_res_ = ConstraintSetOP( new ConstraintSet( *cst_set ) );
}

void
ConstraintSetMover::constraint_file( std::string const & cst_file )
{
	cst_file_ = cst_file;
	cst_fa_file_ = cst_file;
}


ConstraintSetOP ConstraintSetMover::constraint_set() { return constraint_set_low_res_; }
ConstraintSetCOP ConstraintSetMover::constraint_set() const { return constraint_set_low_res_; }

void
ConstraintSetMover::apply( Pose & pose )
{
	if ( !constraint_set_low_res_ && !pose.is_fullatom() ) {
		// uninitialized filename not tolerated, in order to avoid potential confusion
		if ( cst_file_.empty() ) utility_exit_with_message("Can\'t read constraints from empty file!");
		// special case: set cst_file_ to "none" to effectively remove constraints from Pose
		else if ( cst_file_ == "none" ) constraint_set_low_res_ = ConstraintSetOP( new ConstraintSet );
		else {
			constraint_set_low_res_ =
				ConstraintIO::get_instance()->read_constraints( cst_file_, ConstraintSetOP( new ConstraintSet ), pose );
			//ConstraintIO::get_instance()->read_constraints_new( cst_file_, new ConstraintSet, pose );
		}
	}

	if ( !constraint_set_high_res_ && pose.is_fullatom() ) {
		// uninitialized filename not tolerated, in order to avoid potential confusion
		if ( cst_fa_file_.empty() ) utility_exit_with_message("Can\'t read constraints from empty file!");
		// special case: set cst_file_ to "none" to effectively remove constraints from Pose
		else if ( cst_fa_file_ == "none" ) constraint_set_high_res_ = ConstraintSetOP( new ConstraintSet );
		else {
			constraint_set_high_res_ =
				ConstraintIO::get_instance()->read_constraints( cst_fa_file_, ConstraintSetOP( new ConstraintSet ), pose );
			//ConstraintIO::get_instance()->read_constraints_new( cst_fa_file_, new ConstraintSet, pose );
		}
	}

	if ( pose.is_fullatom() ) {
		if( add_constraints() )
			pose.add_constraints( constraint_set_high_res_->get_all_constraints());
		else
			pose.constraint_set( constraint_set_high_res_ );
		if (TR.Debug.visible()) {
			TR.Debug << "High-res constraints:" << std::endl;
			constraint_set_high_res_->show_definition(TR.Debug, pose);
			TR.Debug.flush(); // Make sure the tracer is output if the definitions don't use a std::endl
		}
	} else {
		if( add_constraints() )
			pose.add_constraints( constraint_set_low_res_->get_all_constraints() );
		else
			pose.constraint_set( constraint_set_low_res_ );
		if (TR.Debug.visible()) {
			TR.Debug << "Low-res constraints:" << std::endl;
			constraint_set_low_res_->show_definition(TR.Debug, pose);
			TR.Debug.flush(); // Make sure the tracer is output if the definitions don't use a std::endl
		}
	}
}

std::string
ConstraintSetMover::get_name() const {
	return ConstraintSetMoverCreator::mover_name();
}

protocols::moves::MoverOP ConstraintSetMover::clone() const { return protocols::moves::MoverOP( new protocols::simple_moves::ConstraintSetMover( *this ) ); }
protocols::moves::MoverOP ConstraintSetMover::fresh_instance() const { return protocols::moves::MoverOP( new ConstraintSetMover ); }

void
ConstraintSetMover::register_options()
{
	option.add_relevant( OptionKeys::constraints::cst_file );
	option.add_relevant( OptionKeys::constraints::cst_fa_file );
}

void
ConstraintSetMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	if ( tag->hasOption("cst_file") ) cst_file_ = tag->getOption<std::string>("cst_file");
	if ( tag->hasOption("cst_fa_file") ) cst_fa_file_ = tag->getOption<std::string>("cst_fa_file");
	else cst_fa_file_=cst_file_;
	add_constraints( tag->getOption< bool >( "add_constraints", false ) );
	TR << "of type ConstraintSetMover with constraint file: "<< cst_file_ <<std::endl;
	if ( cst_fa_file_ != cst_file_ ) {
		TR << "of type ConstraintSetMover with fullatom constraint file: "<< cst_fa_file_ <<std::endl;
	}
}
void ConstraintSetMover::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & /*score_fxns*/,
				utility::lua::LuaObject const & /*tasks*/,
				protocols::moves::MoverCacheSP /*cache*/ ) {
	if ( def["cst_file"] ) cst_file_ = def["cst_file"].to<std::string>();
	if ( def["cst_fa_file"] ) cst_fa_file_ = def["cst_fa_file"].to<std::string>();
	else cst_fa_file_=cst_file_;
	TR << "of type ConstraintSetMover with constraint file: "<< cst_file_ <<std::endl;
	if ( cst_fa_file_ != cst_file_ ) {
		TR << "of type ConstraintSetMover with fullatom constraint file: "<< cst_fa_file_ <<std::endl;
	}
}
} // moves
} // protocols
