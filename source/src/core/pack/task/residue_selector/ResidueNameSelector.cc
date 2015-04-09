// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueNameSelector.hh
/// @brief  The ResidueNameSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

// Unit headers
#include <core/pack/task/residue_selector/ResidueNameSelector.hh>
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

ResidueNameSelector::ResidueNameSelector():
res_name_str_() {}

ResidueNameSelector::ResidueNameSelector( std::string const & res_name_str )
{
	res_name_str_ = res_name_str;
}


ResidueNameSelector::~ResidueNameSelector() {}

ResidueSubset
ResidueNameSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( !res_name_str_.empty() );

	ResidueSubset subset( pose.total_residue(), false );
	utility::vector1< std::string > const res_name_vec = utility::string_split( res_name_str_, ',' );
	std::set< std::string > res_set;
	for (	utility::vector1< std::string >::const_iterator n=res_name_vec.begin(), endn=res_name_vec.end(); n!=endn; ++n ) {
		res_set.insert( *n );
	}

	for ( core::Size i=1, endi=pose.total_residue(); i<=endi; ++i ) {
		if ( res_set.find( pose.residue(i).name() ) != res_set.end() ) {
			subset[i] = true;
		}
	}

	return subset;
}

void
ResidueNameSelector::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &)
{
	try {
		set_residue_names( tag->getOption< std::string >( "residue_names" ) );
	} catch ( utility::excn::EXCN_Msg_Exception e ) {
		std::stringstream err_msg;
		err_msg << "Failed to access required option 'residue_names' from ResidueNameSelector::parse_my_tag.\n";
		err_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}
}

void
ResidueNameSelector::set_residue_names( std::string const &res_name_str )
{
	res_name_str_ = res_name_str;
}

std::string ResidueNameSelector::get_name() const {
	return ResidueNameSelector::class_name();
}

std::string ResidueNameSelector::class_name() {
	return "ResidueName";
}

ResidueSelectorOP
ResidueNameSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ResidueNameSelector );
}

std::string
ResidueNameSelectorCreator::keyname() const {
	return ResidueNameSelector::class_name();
}

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core

