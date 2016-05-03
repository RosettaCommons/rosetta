// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ChainSelector.cc
/// @brief  The ChainSelector combines logic from multiple ResidueSelectors
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

ChainSelector::ChainSelector():
	chain_strings_()
{}

/// @brief Copy constructor
///
ChainSelector::ChainSelector( ChainSelector const &src) :
	chain_strings_( src.chain_strings_ )
{}

ChainSelector::ChainSelector( std::string chains ) : chain_strings_( utility::string_split( chains, ',' ) ) {}
ChainSelector::~ChainSelector() {}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP ChainSelector::clone() const { return ResidueSelectorOP( new ChainSelector(*this) ); }

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

void
ChainSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	// first define the chain_cslist restriction -- this could be moved elsewhere
	XMLSchemaRestriction chain_cslist;
	chain_cslist.name( "chain_cslist" );
	chain_cslist.base_type( xs_string );
	chain_cslist.add_restriction( xsr_pattern, "[A-Z](,[A-Z])*" );
	xsd.add_top_level_element( chain_cslist);

	utility::tag::AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute( "chains", "chain_cslist" ) );
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
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
	if ( !pose.pdb_info() ) {
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

void
ChainSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ChainSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::ChainSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( chain_strings_ ) ); // utility::vector1<std::string>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::ChainSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( chain_strings_ ); // utility::vector1<std::string>
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::ChainSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::ChainSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_ChainSelector )
#endif // SERIALIZATION
