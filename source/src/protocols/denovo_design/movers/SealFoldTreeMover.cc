// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/movers/SealFoldTreeMover.cc
/// @brief Creates a sealed foldtree, and removes all cutpoints from the pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/movers/SealFoldTreeMover.hh>
#include <protocols/denovo_design/movers/SealFoldTreeMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.SealFoldTreeMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

SealFoldTreeMover::SealFoldTreeMover():
	protocols::moves::Mover( SealFoldTreeMover::class_name() ),
	cutpoints_(),
	fg_(),
	roots_()
{
}

SealFoldTreeMover::SealFoldTreeMover( components::StructureData const & sd, protocols::loops::Loops const & loops ):
	protocols::moves::Mover( SealFoldTreeMover::class_name() ),
	cutpoints_(),
	fg_( new components::FoldGraph( sd ) ),
	roots_()
{
	for ( protocols::loops::Loops::const_iterator l=loops.begin(); l!=loops.end(); ++l ) {
		if ( l->cut() > 0 ) cutpoints_.push_back( l->cut() );
	}
	if ( sd.segments_begin() != sd.segments_end() ) {
		roots_.push_back( *sd.segments_begin() );
	}
}

SealFoldTreeMover::~SealFoldTreeMover(){}

void
SealFoldTreeMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
}

protocols::moves::MoverOP
SealFoldTreeMover::clone() const
{
	return protocols::moves::MoverOP( new SealFoldTreeMover( *this ) );
}

protocols::moves::MoverOP
SealFoldTreeMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new SealFoldTreeMover );
}

std::string
SealFoldTreeMover::get_name() const
{
	return SealFoldTreeMover::class_name();
}

std::string
SealFoldTreeMover::class_name()
{
	return "SealFoldTreeMover";
}

void
SealFoldTreeMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, SealFoldTreeMover const & mover )
{
	mover.show(os);
	return os;
}

void
SealFoldTreeMover::apply( core::pose::Pose & pose )
{
	remove_cutpoints( pose );

	if ( roots_.empty() ) {
		TR.Warning << "No new root segments found in roots_... not doing anything" << std::endl;
		return;
	}

	// create new fold tree
	pose.fold_tree( fg_->fold_tree( roots_ ) );
	TR << "Set fold tree to " << pose.fold_tree() << std::endl;
}

void
SealFoldTreeMover::remove_cutpoints( core::pose::Pose & pose ) const
{
	// remove cutpoints
	for ( Cutpoints::const_iterator cut=cutpoints_.begin(); cut!=cutpoints_.end(); ++cut ) {
		core::Size const resid = *cut;

		TR.Debug << "Removing cutpoint at " << pose.residue( resid ).name() << " " << resid << std::endl;
		// remove cutpoint variants from this residue if they are present
		if ( pose.residue( resid ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, resid );
			if ( ( resid + 1 <= pose.total_residue() ) && ( pose.residue( resid ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) ) {
				core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, resid + 1 );
			}
		}

		/*
		// this will catch rogue upper cutpoint variants
		if ( pose.residue( resid ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, resid );
		}
		*/
	}

}

/////////////// Creator ///////////////

protocols::moves::MoverOP
SealFoldTreeMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new SealFoldTreeMover );
}

std::string
SealFoldTreeMoverCreator::keyname() const
{
	return SealFoldTreeMover::class_name();
}

} //protocols
} //denovo_design
} //movers

