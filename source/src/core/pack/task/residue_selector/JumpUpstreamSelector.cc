// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/JumpUpstreamSelector.hh
/// @brief  The JumpUpstreamSelector selects residues downstream of a given jump in a FoldTree
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <core/pack/task/residue_selector/JumpUpstreamSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/kinematics/FoldTree.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <utility/assert.hh>


namespace core {
namespace pack {
namespace task {
namespace residue_selector {

JumpUpstreamSelector::JumpUpstreamSelector():
jump_(0) {}

JumpUpstreamSelector::JumpUpstreamSelector( int jump )
{
	jump_ = jump;
}


JumpUpstreamSelector::~JumpUpstreamSelector() {}

ResidueSubset
JumpUpstreamSelector::apply( core::pose::Pose const & pose ) const
{
debug_assert( jump_ > 0 );
	ResidueSubset subset( pose.total_residue(), false );

	ObjexxFCL::FArray1D_bool upstream( pose.total_residue() );
	pose.fold_tree().partition_by_jump( jump_, upstream );

	for( core::Size ii = 1; ii < upstream.size(); ++ii ) {
		subset[ ii ] = upstream( ii );
	}
	return subset;
}

void
JumpUpstreamSelector::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &)
{
	try {
		set_jump( tag->getOption< int >( "jump" ) );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream err_msg;
		err_msg << "Failed to access required option 'jump' from JumpUpstreamSelector::parse_my_tag.\n";
		err_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}
}

void
JumpUpstreamSelector::set_jump( int jump )
{
	jump_ = jump;
}

std::string JumpUpstreamSelector::get_name() const {
	return JumpUpstreamSelector::class_name();
}

std::string JumpUpstreamSelector::class_name() {
		  return "JumpUpstream";
}

ResidueSelectorOP
JumpUpstreamSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new JumpUpstreamSelector );
}

std::string
JumpUpstreamSelectorCreator::keyname() const {
	return JumpUpstreamSelector::class_name();
}

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core
