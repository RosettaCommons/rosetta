// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.VirtualRootMover" );

VirtualRootMover::VirtualRootMover() :
	moves::Mover( VirtualRootMover::mover_name() ),
	remove_(false),
	removable_(false)
{}

VirtualRootMover::~VirtualRootMover() = default;

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
				debug_assert( virtroot == pose.size() && pose.residue(virtroot).aa() == core::chemical::aa_vrt );
				core::pose::setPoseExtraScore( pose, "VirtualRootMover_root", virtroot );
			}
		} else {
			TR << "Pose already virtually rooted. Not changing root." << std::endl;
		}
	}
}

// XRW TEMP std::string
// XRW TEMP VirtualRootMover::get_name() const {
// XRW TEMP  return VirtualRootMover::mover_name();
// XRW TEMP }

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

std::string VirtualRootMover::get_name() const {
	return mover_name();
}

std::string VirtualRootMover::mover_name() {
	return "VirtualRoot";
}

void VirtualRootMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"remove", xsct_rosetta_bool,
		"Removes the virtual root from the pose. Useful for subsequent use to a previous VirtualRoot call")
		+ XMLSchemaAttribute(
		"removable", xsct_rosetta_bool,
		"Set this to true of you want the virtual root to be removable. See remove option.");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Reroot the pose foldtree on a (new) virtual residue. Useful for minimization "
		"in the context of absolute frames (coordinate constraints, electron density information, etc.) "
		"By default, the mover will add a virtual root residue to the pose if one "
		"does not already exist. If you wish to later remove the virtual root, "
		"add the root with a mover with removable set to true, and then later "
		"use a separate VirtualRoot mover with remove set to true to remove it.",
		attlist );
}

std::string VirtualRootMoverCreator::keyname() const {
	return VirtualRootMover::mover_name();
}

protocols::moves::MoverOP
VirtualRootMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new VirtualRootMover );
}

void VirtualRootMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	VirtualRootMover::provide_xml_schema( xsd );
}


} // simple_moves
} // protocols
