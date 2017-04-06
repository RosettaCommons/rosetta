// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/GlycanLayerSelector.hh
/// @brief  A selector for choosing glycan residues based on their layer - as measured by the residue distance to the start of the glycan tree.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/select/residue_selector/GlycanLayerSelector.hh>
#include <core/select/residue_selector/GlycanLayerSelectorCreator.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/chemical/ResidueType.hh>

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


static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.GlycanLayerSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
GlycanLayerSelector::GlycanLayerSelector():
	ResidueSelector()
{
}

/// @brief Destructor.
///
GlycanLayerSelector::~GlycanLayerSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//GlycanLayerSelector::GlycanLayerSelector(GlycanLayerSelector const & src):
// ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
GlycanLayerSelector::ResidueSelectorOP
GlycanLayerSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast<core::select::residue_selector::ResidueSelector>(
		GlycanLayerSelectorOP( new GlycanLayerSelector(*this) )
		)
	);
}


/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
GlycanLayerSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )

{
	if ( tag->hasOption("start") && tag->hasOption("end") ) {
		range_set_ = true;
		start_ = tag->getOption< core::Size >("start", start_);
		end_ = tag->getOption< core::Size >("end", end_);
	} else if ( tag->hasOption( "start ") || tag->hasOption( "end" ) ) {
		utility_exit_with_message("For range selection, start and end must be set!");
	} else {
		start_from_as_layer_ = tag->getOption< core::Size >("layer_as_greater_than_or_equal_to", start_from_as_layer_);
		end_for_layer_ = tag->getOption< core::Size >("layer_as_less_than_or_equal_to", end_for_layer_);
		range_set_ = false;
	}
}

std::string GlycanLayerSelector::get_name() const
{
	return GlycanLayerSelector::class_name();
}

std::string GlycanLayerSelector::class_name()
{
	return "GlycanLayerSelector";
}

void GlycanLayerSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes

		+ XMLSchemaAttribute(  "start",      xsct_non_negative_integer ,
		"Start of the glycan layer we are selecting from."
		"0-N" )

		+ XMLSchemaAttribute(  "end",    xsct_non_negative_integer ,
		"End of the glycan layer we are selecting from."
		"1-N" )

		+ XMLSchemaAttribute(  "layer_as_greater_than_or_equal_to", xsct_non_negative_integer ,
		"Set the layer as all residues greater than or equal to this number (such as the end of the tree)" )

		+ XMLSchemaAttribute(  "layer_as_less_than_or_equal_to", xsct_non_negative_integer ,
		"Set the layer as all residue less or equal to this number (the beginning of the tree)" );

	xsd_type_definition_w_attributes( xsd, class_name(),
		"A selector for choosing glycan residues based on their layer - as measured by the residue distance to the start of the glycan tree."
		"If no layer is set, will select all glycan residues.", attributes );

}

void
GlycanLayerSelector::set_layer(core::Size start, core::Size end){
	start_ = start;
	end_ = end;
	start_from_as_layer_ = 0;
	end_for_layer_ = 0;
	range_set_ = true;
}

void
GlycanLayerSelector::set_layer_as_greater_than_or_equal_to(core::Size start ){
	start_from_as_layer_ = start;
	start_ = 0;
	end_ = 0;
	end_for_layer_ = 0;
	range_set_ = false;
}

void
GlycanLayerSelector::set_layer_as_less_than_or_equal_to( core::Size end ){
	end_for_layer_ = end;
	start_ = 0;
	end_ = 0;
	end_for_layer_ = 0;
	range_set_ = false;
}


/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
GlycanLayerSelector::ResidueSubset
GlycanLayerSelector::apply(
	core::pose::Pose const & pose
) const {

	debug_assert( start_ <= end_);
	utility::vector1< bool > layer_subsets(pose.total_residue(), false );

	for ( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {
		if ( pose.residue_type( resnum ).is_carbohydrate() ) {
			core::Size layer = pose.glycan_tree_set()->get_distance_to_start( resnum );

			if ( range_set_ ) {

				if ( layer >= start_ && layer <= end_ ) {
					TR.Debug << "Adding " << resnum << " at layer " << layer << std::endl;
					layer_subsets[ resnum ] = true;
				}
			} else if ( end_for_layer_ != 0 ) {
				if ( layer <= end_for_layer_ ) {
					TR.Debug << "2-Adding " << resnum << " at layer " << layer << std::endl;
					layer_subsets[ resnum ] = true;
				}
			} else {
				if ( layer >= start_from_as_layer_ ) {
					TR.Debug << "3-Adding " << resnum << " at layer " << layer << std::endl;
					layer_subsets[ resnum ] = true;
				}
			}
		}
	}
	TR.Trace << "Layer Subset: " << layer_subsets << std::endl;
	return layer_subsets;
}





GlycanLayerSelector::ResidueSelectorOP
GlycanLayerSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
		GlycanLayerSelectorOP( new GlycanLayerSelector )
		)
	);
}

std::string
GlycanLayerSelectorCreator::keyname() const {
	return GlycanLayerSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
GlycanLayerSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	GlycanLayerSelector::provide_xml_schema( xsd );
}









} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION
template< class Archive >
void
core::select::residue_selector::GlycanLayerSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( start_) );
	arc( CEREAL_NVP( end_) );
	arc( CEREAL_NVP( end_for_layer_ ) );
	arc( CEREAL_NVP( start_from_as_layer_ ) );
	arc( CEREAL_NVP( range_set_ ));
}

template< class Archive >
void
core::select::residue_selector::GlycanLayerSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( start_ );
	arc( end_ );
	arc( end_for_layer_ );
	arc( start_from_as_layer_ );
	arc( range_set_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::GlycanLayerSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::GlycanLayerSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_GlycanLayerSelector )
#endif // SERIALIZATION
