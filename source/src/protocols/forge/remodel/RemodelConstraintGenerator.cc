// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelConstraintGenerator.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009
/// @modified Tom Linsky, tlinsky@uw.edu

// AUTO-REMOVED #include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/id/SequenceMapping.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.forge.remodel.remodelconstraintgenerator" );
namespace protocols {
namespace forge {
namespace remodel {

std::map<std::string,core::scoring::constraints::ConstraintCOPs> RemodelConstraintGenerator::cst_map_;

/// @details Auto-generated virtual destructor
RemodelConstraintGenerator::~RemodelConstraintGenerator()
{}


RemodelConstraintGenerator::RemodelConstraintGenerator()
	: Mover("RemodelConstraintGenerator"),
		id_( "" ),
		seqmap_(NULL),
		vlb_(NULL)
{}

RemodelConstraintGenerator::RemodelConstraintGenerator( RemodelConstraintGenerator const & rval )
	: Mover( rval ),
		id_( rval.id_ ),
		seqmap_( rval.seqmap_ ),
		vlb_( rval.vlb_ )
{}

/// @details When called, generates constraints for the pose, and adds them to the pose
void
RemodelConstraintGenerator::apply( core::pose::Pose & pose )
{
	init( pose );
	add_remodel_constraints_to_pose( pose );
}

/// @brief This is called if this mover is instantiated from XML
void
RemodelConstraintGenerator::parse_my_tag( utility::tag::TagCOP const tag,
																					basic::datacache::DataMap &,
																					protocols::filters::Filters_map const &,
																					protocols::moves::Movers_map const &,
																					core::pose::Pose const & )
{
	// if there are any options that we might want, they will go here...
	id_ = tag->getOption<std::string>( "name", id_ );
	TR << "Setting id for type " << this->get_name() << " = " << id_ << std::endl;
}

void
RemodelConstraintGenerator::add_remodel_constraints_to_pose(
	core::pose::Pose & pose )
{
	//TR << "Clearing old constraints in RCG" << std::endl;
	if ( csts_.size() > 0 ){
		clear_constraints();
	}

	//TR << "Generating remodel constraints in RCG" << std::endl;
	generate_remodel_constraints(	pose );
	//TR << "Done generating remodel csts in RCG" << std::endl;

	//safeguard against an RCG not generating anything
	if( csts_.size() == 0 ) return;

	//TR << this->get_name() << " generated " << csts_.size() << " constraints." << std::endl;

	csts_ = pose.add_constraints( csts_ );
	store_constraints();
}

void
RemodelConstraintGenerator::remove_remodel_constraints_from_pose(
	core::pose::Pose & pose
) const
{
	core::scoring::constraints::ConstraintCOPs remodel_csts( csts_ );

	//safeguard against an RCG not generating anything
	if( remodel_csts.size() == 0 ) return;

	//TR << this->get_name() << " is about to try to remove " << remodel_csts.size() << " constraints." << std::endl;

	if( ! pose.remove_constraints( remodel_csts, true ) ){
		utility_exit_with_message("Remodel constraints somehow got lost among the way");
	}
}

void
RemodelConstraintGenerator::add_constraint(	core::scoring::constraints::ConstraintCOP cst )
{
	csts_.push_back( cst );
}

void
RemodelConstraintGenerator::add_constraints( core::scoring::constraints::ConstraintCOPs csts )
{
	for( core::scoring::constraints::ConstraintCOPs::const_iterator cst_it = csts.begin();
			 cst_it != csts.end(); ++cst_it ){
		add_constraint( (*cst_it) );
	}
}

void
RemodelConstraintGenerator::clear_constraints()
{
	csts_.clear();
	clear_stored_constraints();
}

RemodelConstraintGenerator::VarLengthBuildAP
RemodelConstraintGenerator::vlb() const {
	return vlb_;
}

void
RemodelConstraintGenerator::set_vlb(
	VarLengthBuildAP vlb )
{
	vlb_ = vlb;
}

std::string
RemodelConstraintGenerator::id() const
{
	return id_;
}

void
RemodelConstraintGenerator::set_id( std::string const & id )
{
	id_ = id;
}

void
RemodelConstraintGenerator::set_seqmap(
	core::id::SequenceMappingCOP seqmap )
{
	seqmap_ = seqmap;
}

core::id::SequenceMappingCOP
RemodelConstraintGenerator:: seqmap() const
{
	return seqmap_;
}

void
RemodelConstraintGenerator::clear_stored_constraints()
{
	if ( id_ == "" ) return;
	// find id for this class in the map
	std::map< std::string, core::scoring::constraints::ConstraintCOPs >::iterator cst_it( cst_map_.find( id_ ) );
	// if the constraints are found, warn the user and erase the old data
	if ( cst_it != cst_map_.end() ) {
		TR << "Overwriting constraints for " << this->get_name() << " named " << id_ << std::endl;
		cst_map_.erase( cst_it );
	}
}

void
RemodelConstraintGenerator::store_constraints()
{
	if ( id_ == "" ) {
		//TR << "ID is not set for this " << this->get_name() << " object. Constraints will not be removable by XML." << std::endl;
		return;
	}
	// clear any stored csts
	clear_stored_constraints();

	// store the csts
	cst_map_.insert( std::pair< std::string, core::scoring::constraints::ConstraintCOPs >( id_, csts_ ) );
}

core::scoring::constraints::ConstraintCOPs const
RemodelConstraintGenerator::lookup_stored_constraints( std::string const & id )
{
	if ( id == "" ) {
		utility_exit_with_message( "ID is not set! The constraint set returned will be empty. Something is being mis-used." );
	}
	// find the constraint set in map
	std::map< std::string, core::scoring::constraints::ConstraintCOPs >::const_iterator cst_it( cst_map_.find( id ) );
	if ( cst_it == cst_map_.end() ) {
		utility_exit_with_message( "Tried to remove constraints that aren't stored for ID=" + id );
	}
	return cst_it->second;
}

} //namespace remodel
} //namespace forge
} //namespace protocols
