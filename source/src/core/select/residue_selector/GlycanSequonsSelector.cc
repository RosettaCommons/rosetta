// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/GlycanSequonsSelector.hh
/// @brief  Find glycosylation sequons in pose
/// @detail Select amino acids based on a sequence motif that is known to trigger glycosylation. The best known one is Asn-notPro-Ser/Thr.
/// @author raemisch (raemisch@scripps.edu)


// Unit headers
#include <core/select/residue_selector/GlycanSequonsSelector.hh>
#include <core/select/residue_selector/ResidueInSequenceMotifSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>
// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.select.residue_selector.GlycanSequonsSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
/// Select N-notP-S/T by default
GlycanSequonsSelector::GlycanSequonsSelector():
	ResidueInSequenceMotifSelector("N[^P][ST]",1)
{}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
ResidueSelectorOP
GlycanSequonsSelector::clone() const {
	return ResidueSelectorOP( utility::pointer::dynamic_pointer_cast<ResidueSelector>( GlycanSequonsSelectorOP( new GlycanSequonsSelector(*this) ) ) );
}
ResidueSelectorOP
GlycanSequonsSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
		GlycanSequonsSelectorOP( new GlycanSequonsSelector )
		)
	);
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
GlycanSequonsSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	NxST_ = tag->getOption< bool >("NxST", NxST_);
	NxC_ = tag->getOption< bool >("NxC", NxC_);
	WxxW_ = tag->getOption< bool >("WxxW", WxxW_);
	WSTxC_ = tag->getOption< bool >("WSTxC", WSTxC_);

	int nr_trues = 0;
	if ( NxST_ ) { ++nr_trues; }
	if ( NxC_ ) {
		++nr_trues;
		regex_ = "N[A_Z]C";
	}
	if ( WxxW_ ) {
		++nr_trues;
		regex_ = "W[A-Z][A-Z]W";
	}
	if ( WSTxC_ ) {
		++nr_trues;
		regex_ = "W[ST][A-Z]C";
	}

	if ( nr_trues == 0 ) {
		throw CREATE_EXCEPTION( utility::excn::Exception, "No glycosylation motif specified.");
	}
	if ( nr_trues > 1 ) {
		throw CREATE_EXCEPTION( utility::excn::Exception, "Too many glycosylation motifs specified. To select any other sequon than NxST, you need to specify NxST=\"0\". Please specify only one motif and combine selections using an OrResidueSelector");
	}
}

std::string GlycanSequonsSelector::get_name() const
{
	return GlycanSequonsSelector::class_name();
}

std::string GlycanSequonsSelector::class_name()
{
	return "GlycanSequonsSelector";
}

void GlycanSequonsSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute::attribute_w_default("NxST", xsct_rosetta_bool, "Use N-notP-S/T glycosylation motif", "true")
		+ XMLSchemaAttribute::attribute_w_default("NxC", xsct_rosetta_bool, "Use N-X-C glycosylation motif", "false")
		+ XMLSchemaAttribute::attribute_w_default("WxxW", xsct_rosetta_bool, "Use W-X-X-W mannosylation motif", "false")
		+ XMLSchemaAttribute::attribute_w_default("WSTxC", xsct_rosetta_bool, "Use W-S/T-X-C mannosylation motif", "false");
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Select residues that can be glycosylated, based on a known sequence motif. Default: N-(not P)-S/T", attributes );
}

std::string
GlycanSequonsSelectorCreator::keyname() const {
	return GlycanSequonsSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
GlycanSequonsSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	GlycanSequonsSelector::provide_xml_schema( xsd );
}

} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION

template< class Archive >
void
core::select::residue_selector::GlycanSequonsSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( NxST_) );
	arc( CEREAL_NVP( NxC_) );
	arc( CEREAL_NVP( WxxW_) );
	arc( CEREAL_NVP( WSTxC_) );
}

template< class Archive >
void
core::select::residue_selector::GlycanSequonsSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( NxST_);
	arc( NxC_);
	arc( WxxW_);
	arc( WSTxC_);
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::GlycanSequonsSelector );

CEREAL_REGISTER_TYPE( core::select::residue_selector::GlycanSequonsSelector )
CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_GlycanSequonsSelector )


#endif // SERIALIZATION
