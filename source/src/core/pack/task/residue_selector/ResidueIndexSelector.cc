// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueIndexSelector.hh
/// @brief  The ResidueIndexSelector selects residues using a string containing pose indices
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <core/pack/task/residue_selector/ResidueIndexSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// C++ headers
#include <utility/assert.hh>


namespace core {
namespace pack {
namespace task {
namespace residue_selector {

ResidueIndexSelector::ResidueIndexSelector():
	index_str_() {}

ResidueIndexSelector::ResidueIndexSelector( std::string const & index_str )
{
	index_str_ = index_str;
}


ResidueIndexSelector::~ResidueIndexSelector() {}

ResidueSubset
ResidueIndexSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( !index_str_.empty() );

	ResidueSubset subset( pose.total_residue(), false );
	std::set< Size > const res_set( get_resnum_list( index_str_, pose ) );

	for ( std::set< Size >::const_iterator it = res_set.begin();
			it != res_set.end(); ++it ) {
		if ( *it == 0 || *it > subset.size() ) {
			std::stringstream err_msg;
			err_msg << "Residue " << *it << " not found in pose!\n";
			throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
		}
		subset[ *it ] = true;
	}
	return subset;
}

void
ResidueIndexSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	try {
		set_index( tag->getOption< std::string >( "resnums" ) );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream err_msg;
		err_msg << "Failed to access required option 'resnums' from ResidueIndexSelector::parse_my_tag.\n";
		err_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}
}

void
ResidueIndexSelector::set_index( std::string const &index_str )
{
	index_str_ = index_str;
}

std::string ResidueIndexSelector::get_name() const {
	return ResidueIndexSelector::class_name();
}

std::string ResidueIndexSelector::class_name() {
	return "Index";
}

ResidueSelectorOP
ResidueIndexSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ResidueIndexSelector );
}

std::string
ResidueIndexSelectorCreator::keyname() const {
	return ResidueIndexSelector::class_name();
}

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core
