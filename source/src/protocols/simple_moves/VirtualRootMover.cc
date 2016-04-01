// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/VirtualRootMover.cc
/// @brief  Manipulate virtual roots on poses. This is a seperate mover mainly to be RosettaScriptable
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// unit headers
#include <protocols/simple_moves/VirtualRootMover.hh>
#include <protocols/simple_moves/VirtualRootMoverCreator.hh>

// type headers
#include <core/types.hh>

// project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.VirtualRootMover" );

std::string
VirtualRootMoverCreator::keyname() const
{
	return VirtualRootMoverCreator::mover_name();
}

protocols::moves::MoverOP
VirtualRootMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new VirtualRootMover );
}

std::string
VirtualRootMoverCreator::mover_name()
{
	return "VirtualRoot";
}

VirtualRootMover::VirtualRootMover() :
	moves::Mover( VirtualRootMoverCreator::mover_name() ),
	remove_(false),
	removable_(false)
{}

VirtualRootMover::~VirtualRootMover() {}

moves::MoverOP VirtualRootMover::clone() const {
	return moves::MoverOP( new VirtualRootMover( *this ) );
}
moves::MoverOP VirtualRootMover::fresh_instance() const {
	return moves::MoverOP( new VirtualRootMover );
}

void
VirtualRootMover::apply( core::pose::Pose & pose ) {
	if ( remove_ ) {
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			TR << "Can't remove virtual root - pose is not rooted on virtual residue." << std::endl;
			return;
		}
		core::Real virtroot;
		if ( ! core::pose::getPoseExtraScore( pose, "VirtualRootMover_root", virtroot ) ) {
			TR << "Virtual root not added as removable by VirtualRootMover - not removing root." << std::endl;
			return;
		}
		if ( pose.residue( virtroot ).aa() != core::chemical::aa_vrt || virtroot != pose.fold_tree().root() ) {
			TR.Warning << "Residue added by VirtualRootMover no longer virtual root - not removing root." << std::endl;
		} else {
			pose.conformation().delete_residue_slow(virtroot);
		}
		core::pose::clearPoseExtraScore( pose, "VirtualRootMover_root");
	} else {
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::pose::addVirtualResAsRoot( pose );
			if ( removable_ ) {
				core::Size virtroot = pose.fold_tree().root();
				assert( virtroot == pose.n_residue() && pose.residue(virtroot).aa() == core::chemical::aa_vrt );
				core::pose::setPoseExtraScore( pose, "VirtualRootMover_root", virtroot );
			}
		} else {
			TR << "Pose already virtually rooted. Not changing root." << std::endl;
		}
	}
}

std::string
VirtualRootMover::get_name() const {
	return VirtualRootMoverCreator::mover_name();
}

void VirtualRootMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	moves::Movers_map const &,
	Pose const &
) {
	set_remove( tag->getOption<bool>( "remove",  remove_ ) );
	set_removable( tag->getOption<bool>( "removable", removable_ ) );
}

void
VirtualRootMover::set_removable( bool const removable )
{
	removable_ = removable;
}

void
VirtualRootMover::set_remove( bool const remove )
{
	remove_ = remove;
}

} // simple_moves
} // protocols
