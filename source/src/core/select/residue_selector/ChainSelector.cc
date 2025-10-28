// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers

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

static basic::Tracer TR("core.select.residue_selector.ChainSelector");

ChainSelector::ChainSelector():
	chain_strings_()
{}

/// @brief Copy constructor
///
ChainSelector::ChainSelector( ChainSelector const &src) :
	chain_strings_( src.chain_strings_ )
{}

ChainSelector::ChainSelector( utility::vector1< std::string > const & chains ) : chain_strings_( chains ) {}
ChainSelector::ChainSelector( std::string const & chains ) : chain_strings_( utility::string_split( chains, ',' ) ) {}
ChainSelector::ChainSelector( char chain ) : chain_strings_( 1, std::string( 1, chain ) ) {}
ChainSelector::ChainSelector( int chainid ) : chain_strings_( 1, std::to_string( chainid ) ) {}
ChainSelector::ChainSelector( core::Size chainid ) : chain_strings_( 1, std::to_string( chainid ) ) {}
ChainSelector::~ChainSelector() = default;

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP ChainSelector::clone() const { return utility::pointer::make_shared< ChainSelector >(*this); }

ResidueSubset
ChainSelector::apply(
	core::pose::Pose const & pose
) const
{
	ResidueSubset subset( pose.size(), false );
	for ( std::string const & iichain_string : chain_strings_ ) {
		core::Size ii_num = 0;
		try {
			ii_num = boost::lexical_cast< core::Size > ( iichain_string );
			select_chain_by_index( pose, subset, ii_num );
		} catch ( boost::bad_lexical_cast & ) {
			select_chain_by_pdb_chain_char( pose, subset, iichain_string );
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
	} catch ( utility::excn::Exception & e ) {
		std::stringstream error_message;
		error_message << "ChainSelector::parse_my_tag was not able to find the required option 'chains' in the input Tag\n";
		error_message << e.msg();
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
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

	utility::tag::AttributeList attributes;
	attributes + XMLSchemaAttribute(
		"chains", xsct_chain_cslist,
		"The string given for the \"chains\" option should be a "
		"comma-separated list of chain identifiers. "
		"Each chain identifier should be either an integer, so that the Pose "
		"chain index will be used, or a single character, so that the PDB chain ID can be used.");
	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The ChainSelector sets the positions corresponding to all the residues "
		"in the given set of chains to true, and all the other positions to false.",
		attributes );
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
	std::string const & ch
) const
{
	if ( !pose.pdb_info() ) {
		std::ostringstream err;
		err << get_name() << "Selector received a pose without a valid PDBInfo--chains cannot be selected.";
		throw CREATE_EXCEPTION(utility::excn::NullPointerError,  err.str() );
	}

	bool found = false;

	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( pose.pdb_info()->chain( ii ) == ch ) {
			subset[ ii ] = true;
			found = true;
		}
	}

	if ( !found && ch.size() != 1 ) {
		TR.Warning << "ChainSelector was unable to interpret `" << ch << "` as a valid chain letter specifier. " << std::endl;
	}
}


ResidueSelectorOP
ChainSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared< ChainSelector >();
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
