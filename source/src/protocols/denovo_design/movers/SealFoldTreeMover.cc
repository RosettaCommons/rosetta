// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/SealFoldTreeMover.cc
/// @brief Creates a sealed foldtree, and removes all cutpoints from the pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/movers/SealFoldTreeMover.hh>
#include <protocols/denovo_design/movers/SealFoldTreeMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/kinematics/FoldTree.hh>

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
	components::StructureData sd = components::StructureDataFactory::get_instance()->get_from_pose( pose );
	remove_cutpoints( pose );
	remove_cutpoints( sd );
	components::StructureDataFactory::get_instance()->save_into_pose( pose, sd );

	if ( roots_.empty() ) {
		TR.Warning << "No new root segments found in roots_... not doing anything" << std::endl;
		return;
	}

	// create new fold tree
	core::kinematics::FoldTree const ft = fg_->fold_tree( roots_ );

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::kinematics::FoldTree const symm_ft = symmetric_fold_tree( pose, ft );
		pose.fold_tree( symm_ft );
	} else {
		pose.fold_tree( ft );
	}
	TR << "Set fold tree to " << pose.fold_tree() << std::endl;
}

void
SealFoldTreeMover::remove_cutpoints( components::StructureData & sd ) const
{
	for ( Cutpoints::const_iterator cut=cutpoints_.begin(); cut!=cutpoints_.end(); ++cut ) {
		core::Size const resid = *cut;
		SegmentName const seg_name = sd.segment_name( resid );
		if ( sd.segment( seg_name ).cutpoint() == resid ) sd.set_cutpoint( seg_name, 0 );
		TR << "Removed cutpoint at " << resid << " (segment " << seg_name << ") from StructureData." << std::endl;
	}
}

void
SealFoldTreeMover::remove_cutpoints( core::pose::Pose & pose ) const
{
	components::StructureData sd = components::StructureDataFactory::get_instance()->get_from_pose( pose );
	// remove cutpoints
	for ( Cutpoints::const_iterator cut=cutpoints_.begin(); cut!=cutpoints_.end(); ++cut ) {
		core::Size const resid = *cut;

		TR.Debug << "Removing cutpoint at " << pose.residue( resid ).name() << " " << resid << std::endl;
		// remove cutpoint variants from this residue if they are present
		if ( pose.residue( resid ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, resid );
			if ( ( resid + 1 <= pose.size() ) && ( pose.residue( resid ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) ) {
				core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, resid + 1 );
			}
		}

		std::string const segment = sd.segment_name( resid );
		if ( sd.segment( segment ).cutpoint() == resid ) {
			components::Segment newseg = sd.segment( segment );
			newseg.set_cutpoint( 0 );
			sd.replace_segment( segment, newseg );
		}
		/*
		// this will catch rogue upper cutpoint variants
		if ( pose.residue( resid ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, resid );
		}
		*/
	}
	components::StructureDataFactory::get_instance()->save_into_pose( pose, sd );

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

