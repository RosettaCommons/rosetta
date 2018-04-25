// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/ResidueInSequenceMotifSelector.hh
/// @brief  Select residues by motif search
/// @author raemisch (raemisch@scripps.edu)

// Unit headers
#include <core/select/residue_selector/ResidueInSequenceMotifSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>
#include <regex>
#include <string>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.select.residue_selector.ResidueInSequenceMotifSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
ResidueInSequenceMotifSelector::ResidueInSequenceMotifSelector():
	core::select::residue_selector::ResidueSelector(),
	regex_(""),
	position_in_motif_(1)
{
}

/// @brief Destructor.
ResidueInSequenceMotifSelector::~ResidueInSequenceMotifSelector() {}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
ResidueInSequenceMotifSelector::clone() const {
	return ResidueSelectorOP( utility::pointer::dynamic_pointer_cast<ResidueSelector>( ResidueSelectorOP( new ResidueInSequenceMotifSelector(*this) ) ) );
}

ResidueSelectorOP
ResidueInSequenceMotifSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
		ResidueInSequenceMotifSelectorOP( new ResidueInSequenceMotifSelector )
		)
	);
}


/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResidueInSequenceMotifSelector::ResidueSubset
ResidueInSequenceMotifSelector::apply(  core::pose::Pose const & pose ) const {

	TR << "Use ResidueInSequenceMotifSelector with regular expression " << std::endl;
	// test if regex can work
	if ( !regex_usable() ) {
		utility_exit_with_message("Your build does not support regular expression. Try to use e.g. gcc4.9. ResidueInSequenceMotifSelector will not work.");
	}
	std::string sequence = pose.sequence();
	const std::string full_sequence = sequence;
	std::regex re(regex_);
	ResidueSubset subset( sequence.length(), false );
	std::smatch string_match;
	Size pos = 0;
	int nr_positions = 0;
	while ( regex_search(sequence,string_match,re) ) {
		TR << "searching in: " << sequence << std::endl;
		pos += ( string_match.position( 0 ) + position_in_motif_ );
		std::string match = string_match[0];
		TR << "Found motif " << " (" << match << ") ";
		TR << "at sequence position " << pos << std::endl;
		subset[ pos ] = true;
		sequence = full_sequence.substr(pos);
		// increment pos by the length of match minus how far inside the motif we aleady are
		//pos += match.size() - position_in_motif_;
		++nr_positions;
		//sequence = string_match.suffix().str();
	}
	TR << "Found " << nr_positions << " motifs in sequence." << std::endl;

	return subset;
}

/// @brief XML tarse.
/// @details Parse RosettaScripts tags and set up this mover.
void
ResidueInSequenceMotifSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	if ( tag->hasOption("motif") ) {
		utility::vector1< std::string > option = utility::string_split_multi_delim( tag->getOption< std::string >("motif"), ",");
		regex_ = option[1];
		position_in_motif_ = utility::string2Size( option[2] );
	} else {
		std::string error_message = "ResidueInSequenceMotifSelector needs a string that describes the motif and the position within that motif e.g. (N[TSAG]RT),1\n";
		throw CREATE_EXCEPTION( utility::excn::Exception, error_message );
	}

}

std::string ResidueInSequenceMotifSelector::get_name() const
{
	return ResidueInSequenceMotifSelector::class_name();
}

std::string ResidueInSequenceMotifSelector::class_name()
{
	return "ResidueInSequenceMotifSelector";
}

void ResidueInSequenceMotifSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "motif",  xsct_string_cslist, "Comma-separated string: Regular a exression that describes the sequence motif, position within the motif (e.g motif=\"N.[ST],1)\"" );
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Select a residue in a all occurances of a user-specified sequence motif.", attributes );
}

std::string
ResidueInSequenceMotifSelectorCreator::keyname() const {
	return ResidueInSequenceMotifSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
ResidueInSequenceMotifSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	ResidueInSequenceMotifSelector::provide_xml_schema( xsd );
}

} //core
} //select
} //residue_selector

#ifdef SERIALIZATION

template< class Archive >
void
core::select::residue_selector::ResidueInSequenceMotifSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( regex_ ) );
	arc( CEREAL_NVP( position_in_motif_ ) );
}

template< class Archive >
void
core::select::residue_selector::ResidueInSequenceMotifSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( regex_ );
	arc( position_in_motif_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::ResidueInSequenceMotifSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::ResidueInSequenceMotifSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_ResidueInSequenceMotifSelector )
#endif // SERIALIZATION
