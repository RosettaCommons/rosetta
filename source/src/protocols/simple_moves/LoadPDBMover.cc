// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LoadPDBMover.cc
/// @brief simple mover that loads a pdb from file and replaces the current pdb with it. Useful for checkpointing

// Unit headers
#include <protocols/simple_moves/LoadPDBMover.hh>
#include <protocols/simple_moves/LoadPDBMoverCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.LoadPDBMover" );
#include <utility/tag/Tag.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

namespace protocols {
namespace simple_moves {

std::string
LoadPDBMoverCreator::keyname() const
{
	return LoadPDBMoverCreator::mover_name();
}

protocols::moves::MoverOP
LoadPDBMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoadPDBMover );
}

std::string
LoadPDBMoverCreator::mover_name()
{
	return "LoadPDB";
}

LoadPDBMover::LoadPDBMover()
: moves::Mover("LoadPDB"),
	filename_( "" ),
	append_( false )
{
}

void
LoadPDBMover::apply( Pose & pose )
{
	TR<<"Loading pdb file "<<filename_<<std::endl;
	core::pose::Pose loaded_pose = *core::import_pose::pose_from_file( filename_, false/*read foldtree*/ , core::import_pose::PDB_file);
	core::pose::read_comment_pdb( filename_, loaded_pose ); //read comments from pdb file

	if ( append() ) {
		pose.append_pose_by_jump( loaded_pose, pose.size() );
		TR.Debug << "Fold tree after insertion " << pose.fold_tree() << std::endl;
	} else {
		pose = loaded_pose;
	}
}

std::string
LoadPDBMover::get_name() const {
	return LoadPDBMoverCreator::mover_name();
}

moves::MoverOP
LoadPDBMover::clone() const
{
	return moves::MoverOP( new LoadPDBMover( *this ) );
}

moves::MoverOP
LoadPDBMover::fresh_instance() const
{
	return moves::MoverOP( new LoadPDBMover );
}

void
LoadPDBMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	filename( tag->getOption< std::string >( "filename" ) );
	TR<<"filename: "<<filename_<<std::endl;

	append( tag->getOption< bool >( "append", false ) );
}

void
LoadPDBMover::filename( std::string const & s ){
	filename_ = s;
}

std::string
LoadPDBMover::filename() const{
	return filename_;
}

} // simple_moves
} // protocols
