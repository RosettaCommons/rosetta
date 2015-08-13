// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ChainSelector.cc
/// @brief  The ChainSelector combines logic from multiple ResidueSelectors
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/task/residue_selector/ChainSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

ChainSelector::ChainSelector() {}
ChainSelector::ChainSelector( std::string chains ) : chain_strings_( utility::string_split( chains, ',' ) ) {}
ChainSelector::~ChainSelector() {}

ResidueSubset
ChainSelector::apply(
	core::pose::Pose const & pose
) const
{
	ResidueSubset subset( pose.total_residue(), false );
	for ( core::Size ii = 1; ii <= chain_strings_.size(); ++ii ) {
		std::string const & iichain_string = chain_strings_[ ii ];
		core::Size ii_num = 0;
		try {
			ii_num = boost::lexical_cast< core::Size > ( iichain_string );
			select_chain_by_index( pose, subset, ii_num );
		} catch ( boost::bad_lexical_cast & ) {
			if ( iichain_string.size() != 1 ) {
				std::stringstream error_message;
				error_message << "ChainSelector was given a string identifier '" << iichain_string << "' that cannot be parsed as an unsigned integer and is also longer than a single character\n";
				throw utility::excn::EXCN_Msg_Exception( error_message.str() );
			}
			select_chain_by_pdb_chain_char( pose, subset, iichain_string[0] );
		}
	}
	return subset;
}

/// @details Chains are given either as single characters (matched against the chain in the input PDB) or as integers
/// (matched against the pose-assigned index for the chain) in a comma-separated list.
void ChainSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	std::string chain_string;
	try {
		chain_string = tag->getOption< std::string >( "chains" );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream error_message;
		error_message << "ChainSelector::parse_my_tag was not able to find the required option 'chains' in the input Tag\n";
		error_message << e.msg();
		throw utility::excn::EXCN_Msg_Exception( error_message.str() );
	}
	chain_strings_ = utility::string_split( chain_string, ',' );
}

std::string ChainSelector::get_name() const {
	return ChainSelector::class_name();
}

std::string ChainSelector::class_name() {
	return "Chain";
}

utility::vector1< std::string > const &
ChainSelector::chain_strings() const {
	return chain_strings_;
}

void
ChainSelector::set_chain_strings(
	utility::vector1< std::string > const & setting
)
{
	chain_strings_ = setting;
}

void ChainSelector::select_chain_by_index(
	core::pose::Pose const & pose,
	ResidueSubset & subset,
	core::Size chain_index
) const
{
	runtime_assert( chain_index <= pose.conformation().num_chains() );
	for ( core::Size ii = pose.conformation().chain_begin( chain_index ),
			ii_end = pose.conformation().chain_end( chain_index ); ii <= ii_end; ++ii ) {
		subset[ ii ] = true;
	}
}

void ChainSelector::select_chain_by_pdb_chain_char(
	core::pose::Pose const & pose,
	ResidueSubset & subset,
	char ch
) const
{
  if ( !pose.pdb_info() ){
    std::ostringstream err;
    err << get_name() << "Selector recieved a pose without a valid PDBInfo--chains cannot be selected.";
    throw utility::excn::EXCN_NullPointer( err.str() );
  }

	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( pose.pdb_info()->chain( ii ) == ch ) subset[ ii ] = true;
	}
}


ResidueSelectorOP
ChainSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ChainSelector );
}

std::string
ChainSelectorCreator::keyname() const {
	return ChainSelector::class_name();
}

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core
