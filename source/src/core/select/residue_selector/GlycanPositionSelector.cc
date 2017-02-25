// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/GlycanPositionSelector.hh
/// @brief  A Residue Selector for selecting specific parts of arbitrary glycans.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/select/residue_selector/GlycanPositionSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/carbohydrates/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>
#include <utility/string_util.hh>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.GlycanPositionSelector" );


namespace core {
namespace select {
namespace residue_selector {

using namespace core::select::residue_selector;

/// @brief Constructor.
///
GlycanPositionSelector::GlycanPositionSelector():
	ResidueSelector(),
	from_residue_(0),
	to_residue_(0)
{



}

/// @brief Destructor.
///
GlycanPositionSelector::~GlycanPositionSelector() {}


GlycanPositionSelector::GlycanPositionSelector( GlycanPositionSelector const & src ):
	ResidueSelector(src),
	ranges_(src.ranges_),
	positions_(src.positions_),
	from_residue_(src.from_residue_),
	to_residue_(src.to_residue_)

{

}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
GlycanPositionSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast<core::select::residue_selector::ResidueSelector>(
		GlycanPositionSelectorOP( new GlycanPositionSelector(*this) )
		)
	);
}

ResidueSelectorOP
GlycanPositionSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
		GlycanPositionSelectorOP( new GlycanPositionSelector )
		)
	);
}


std::string GlycanPositionSelector::get_name() const
{
	return GlycanPositionSelector::class_name();
}

std::string GlycanPositionSelector::class_name()
{
	return "GlycanPositionSelector";
}

std::string
GlycanPositionSelectorCreator::keyname() const {
	return GlycanPositionSelector::class_name();
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
GlycanPositionSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	from_residue_ = tag->getOption< Size >("from", from_residue_);
	to_residue_ = tag->getOption< Size >("to", to_residue_);

	if ( tag->hasOption("positions") ) {
		utility::vector1< std::string > position_strings = utility::string_split_multi_delim( tag->getOption< std::string >("positions"), ",'`~+*&|;. ");
		for ( Size index = 1; index <= position_strings.size(); ++index ) {
			Size position = utility::string2Size( position_strings[ index ]);
			positions_.push_back( position );
		}
	}
	if ( tag->hasOption("position") ) {
		Size position = tag->getOption< Size >("position");
		positions_.push_back(position);
	}
	if ( tag->hasOption("range") | tag->hasOption("ranges") ) {
		utility::vector1< std::string > range_strings;
		if ( tag->hasOption("range") ) {
			range_strings.push_back( tag->getOption< std::string>( "range"));
		} else if ( tag->hasOption("ranges") ) {
			range_strings = utility::string_split_multi_delim( tag->getOption< std::string >("ranges"), ", ");
		}
		for ( Size index = 1; index <= range_strings.size(); ++index ) {
			utility::vector1< std::string> const str_limits( utility::string_split( range_strings[ index ] , '-' ) );
			if ( str_limits.size() == 2 ) {
				Size const start ( utility::string2Size(str_limits[1]) );
				Size const end ( utility::string2Size(str_limits[2] ) );

				ResidueRange range = ResidueRange(start, end);
				ranges_.push_back(range);
			}
		}
	}
}







void
GlycanPositionSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;

	std::string description = " "
		" @brief A Residue Selector for selecting specific parts of arbitrary glycans by 'position'."
		" @details."

		" Background"

		"  Lets assume that the ASN that a particular N-linked Glycan is attached to, starts from Residue 0."
		"  The residues off this continue with 1 being the next residue, 2, and so on.  Each branch corresponds to a number."

		"  This allows you to choose parts of the glycan, without knowing the actual glycan residue numbers.  "
		"   For example, maybe you want to select the outer part of all glycans or between specific positions."

		" Tips For use"

		"\t This Selector works on all glycans of the pose at once."

		"    Settings are:"

		"      range (start, end)"
		"      positions (specific resnums, such as the carbohydrate at position 4)"

		"      from_residue (all glycan foliage from this and including this residue.)"
		"      to_residue (all glycan foliage up to and including this residue.)"
		" "
		"    Use the 'glycan_info' application to determine the glycan position numbers from a pose."
		"    Combine with the GlycanResidueSelector to get unions of specific glycans."
		"     (such as the leaf of all Man5 residues or the stem of the glycan that starts at ASN85.)";


	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "from",      xsct_positive_integer,
		"all glycan foliage from this and including this residue.", "0" )

		+ XMLSchemaAttribute::attribute_w_default(  "to", xsct_positive_integer,
		"all glycan foliage up to and including this residue.", "0" )

		+ XMLSchemaAttribute(  "position",      xs_string,
		"a specific position to select on.")

		+ XMLSchemaAttribute(  "positions",      xs_string,
		"a specific set of positions to select on.")

		+ XMLSchemaAttribute(  "range",      xs_string,
		"Select a the glycans using a range of positions" )

		+ XMLSchemaAttribute(  "ranges",      xs_string,
		"Select glycans using a range of positions" );

	xsd_type_definition_w_attributes( xsd, class_name(), description, attributes );

}

void
GlycanPositionSelector::add_range(ResidueRange const & range){
	ranges_.push_back(range);
}

void
GlycanPositionSelector::set_range(utility::vector1<ResidueRange> const & ranges){
	ranges_ = ranges;
}

void
GlycanPositionSelector::set_positions(utility::vector1< Size > const & positions){
	positions_ = positions;
}

void
GlycanPositionSelector::set_select_from_residue_position( Size const select_from_residue_position){
	from_residue_ = select_from_residue_position;
}

void
GlycanPositionSelector::set_select_to_residue_position( Size const select_to_residue_position){
	to_residue_ = select_to_residue_position;
}





/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
GlycanPositionSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	GlycanPositionSelector::provide_xml_schema( xsd );
}



/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResidueSubset
GlycanPositionSelector::apply(
	pose::Pose const & pose
) const
{

	utility::vector1< Size > glycan_positions;

	//Range Case
	if ( ranges_.size() > 0 ) {
		for ( Size index = 1; index <= ranges_.size(); ++index ) {
			for ( Size resnum = ranges_[ index ].start(); resnum <= ranges_[ index ].stop(); ++resnum ) {
				if ( ! glycan_positions.has_value( resnum ) ) glycan_positions.push_back( resnum );

			}
		}
	}

	//Positions Case
	if ( positions_.size() > 0 ) {
		for ( Size index = 1; index <= positions_.size(); ++index ) {
			Size resnum = positions_[ index ];
			if ( ! glycan_positions.has_value( resnum ) ) glycan_positions.push_back( resnum );
		}
	}

	//
	if ( to_residue_ == from_residue_ ) {
		GlycanResidueSelector all_selector = GlycanResidueSelector();
		return all_selector.apply(pose);
	} else if ( to_residue_ > from_residue_ ) {
		utility_exit_with_message("GlycanPositionSelector: to_residue > from_residue, resulting in an overlap.  This is probably a mistake!  Please correct!");
	} else {

		Size end_residue;
		if ( to_residue_ <= 0 ) {
			end_residue = to_residue_;
		}
		if ( from_residue_ <= 0 ) {
			end_residue = pose.glycan_tree_set()->get_largest_glycan_tree_length();
		}

		for ( Size i = 1; i <= end_residue; ++i ) {
			if ( ! glycan_positions.has_value( i ) ) glycan_positions.push_back( i );
		}

	}
	return core::pose::carbohydrates::get_resnums_from_glycan_positions( pose, glycan_positions );
}


} //core
} //select
} //residue_selector
