// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/OrResidueSelector.cc
/// @brief  The OrResidueSelector combines logic from multiple ResidueSelectors
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)
/// @see    core::pack::task::residue_selector::AndResidueSelector

// Unit headers
#include <core/pack/task/residue_selector/OrResidueSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/pack/task/residue_selector/ResidueSelectorFactory.hh>
#include <core/pose/Pose.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cassert>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {


OrResidueSelector::OrResidueSelector() {}
OrResidueSelector::~OrResidueSelector() {}

OrResidueSelector::OrResidueSelector( ResidueSelectorCOP selector1, ResidueSelectorCOP selector2 )
{
	add_residue_selector( selector1 );
	add_residue_selector( selector2 );
}

void OrResidueSelector::apply( core::pose::Pose const & pose, ResidueSubset & subset ) const
{
	assert( subset.size() == pose.total_residue() );
	assert( num_selectors() > 0 );

	// make subset neutral for OR operations
	subset = ResidueSubset( pose.total_residue(), false );
	for(std::list< ResidueSelectorCOP >::const_iterator rs = selectors_.begin();
			rs != selectors_.end();
			++rs) {
		ResidueSubset tmp(subset.size());
		(*rs)->apply(pose, tmp);
		apply_or_to_subset(tmp, subset);
	}
}

void OrResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	// grab the comma-separated list of residue selectors that should be ORed together
	// from the tag, and then grab each of the indicated residue selectors from the datamap.
	std::list< ResidueSelectorCOP > local_selectors;
	if( tag->hasOption( "selectors" ) ) {
		std::string selectors_str;
		try {
			selectors_str = tag->getOption< std::string >( "selectors" );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selectors' from OrResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
		utility::vector1< std::string > selector_names = utility::string_split( selectors_str, ',' );
	
	
		for ( core::Size ii = 1; ii <= selector_names.size(); ++ii ) {
			try {
				ResidueSelectorCOP selector = datamap.get< ResidueSelector const * >( "ResidueSelector", selector_names[ ii ] );
				local_selectors.push_back( selector );
			} catch ( utility::excn::EXCN_Msg_Exception e ) {
				std::stringstream error_msg;
				error_msg << "Failed to find ResidueSelector named '" << selector_names[ ii ] << "' from the Datamap from OrResidueSelector::parse_my_tag.\n";
				error_msg << e.msg();
				throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
			}
		}
	}
	// add selectors from tags
	for(utility::vector0< utility::tag::TagCOP >::const_iterator itag = tag->getTags().begin();
			itag != tag->getTags().end(); ++itag) {
		ResidueSelectorCOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
				(*itag)->getName(),
				(*itag),
				datamap
			);
		local_selectors.push_back( rs );	
	}


	if ( local_selectors.size() == 0 ) {
		std::stringstream error_msg;
		error_msg << "No ResidueSelectors given to the OrResidueSelector; OrResidueSelector requires at least one ResidueSelector as input\n";
		throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
	}

	for ( std::list< ResidueSelectorCOP >::const_iterator
			iter = local_selectors.begin(), iter_end = local_selectors.end();
			iter != iter_end; ++iter ) {
		add_residue_selector( *iter );
	}
}

void OrResidueSelector::add_residue_selector( ResidueSelectorCOP selector )
{
	selectors_.push_back(selector);
}

Size OrResidueSelector::num_selectors() const
{
	return selectors_.size();
}

void
OrResidueSelector::apply_or_to_subset(ResidueSubset const & newSubset, ResidueSubset & existingSubset) const
{
	assert( existingSubset.size() == newSubset.size() );
	for( Size ii = 1; ii <= existingSubset.size(); ++ii) {
		existingSubset[ ii ] = existingSubset[ ii ] || newSubset[ ii ];
	}
}

std::string OrResidueSelector::get_name() const {
	return OrResidueSelector::class_name();
}

std::string OrResidueSelector::class_name() {
	return "Or";
}

ResidueSelectorOP
OrResidueSelectorCreator::create_residue_selector() const {
	return new OrResidueSelector;
}

std::string
OrResidueSelectorCreator::keyname() const {
	return OrResidueSelector::class_name();
}


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core

