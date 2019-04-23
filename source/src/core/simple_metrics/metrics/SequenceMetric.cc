// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SequenceMetric.cc
/// @brief A SimpleMetric to output the single-letter OR three-letter sequence of a protein or subset of positions/regions using a ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added support for writing full residue type names or basenames.

// Unit headers
#include <core/simple_metrics/metrics/SequenceMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/util.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ Includes
#include <map>

static basic::Tracer TR( "core.simple_metrics.metrics.SequenceMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/// @brief The map of enum to name.
static const std::map< SequenceMetricMode, std::string > mode_to_name_map_(
	{
	{ SMM_ONELETTER_CODE, "oneletter" },
	{ SMM_THREELETTER_CODE, "threeletter" },
	{ SMM_BASE_NAME, "basename" },
	{ SMM_FULL_NAME, "fullname" },
	{ SMM_INVALID_MODE, "INVALID" }
	}
);

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SequenceMetric::SequenceMetric():
	core::simple_metrics::StringMetric(),
	selector_(nullptr),
	output_mode_(SMM_ONELETTER_CODE)
{}

SequenceMetric::SequenceMetric( select::residue_selector::ResidueSelectorCOP selector ):
	core::simple_metrics::StringMetric(),
	selector_(nullptr), //Set below
	output_mode_(SMM_ONELETTER_CODE)
{
	set_residue_selector( selector );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SequenceMetric::~SequenceMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SequenceMetric::SequenceMetric( SequenceMetric const & src ):
	core::simple_metrics::StringMetric( src ),
	output_mode_(src.output_mode_)
{
	selector_ = src.selector_;
}

core::simple_metrics::SimpleMetricOP
SequenceMetric::clone() const {
	return utility::pointer::make_shared< SequenceMetric >( *this );

}

std::string
SequenceMetric::name() const {
	return name_static();
}

std::string
SequenceMetric::name_static() {
	return "SequenceMetric";

}
std::string
SequenceMetric::metric() const {
	return "sequence";
}

/// @brief Returns all allowed output modes.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
std::string
SequenceMetric::allowed_output_modes() {
	std::stringstream outstring;
	for ( core::Size i(1); i < static_cast<core::Size>( SMM_END_OF_LIST ); ++i ) {
		if ( i > 1 ) {
			outstring << ", ";
			if ( i == static_cast< core::Size >( SMM_END_OF_LIST )-1 ) {
				outstring << "or ";
			}
		}
		outstring << mode_name_from_enum( static_cast< SequenceMetricMode >(i) );
	}
	return outstring.str();
}

/// @brief Returns all allowed output modes as a vector of strings.
/// @details Note that this is a bit inefficient.  It generates the vector each time, and returns it by copy.
/// Not intended for repeated calls.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< std::string >
SequenceMetric::allowed_output_modes_as_vector() {
	utility::vector1< std::string > outvec;
	outvec.reserve( static_cast< core::Size >( SMM_END_OF_LIST - 1 ) );
	for ( core::Size i(1); i < static_cast<core::Size>( SMM_END_OF_LIST ); ++i ) {
		outvec.push_back( mode_name_from_enum( static_cast<SequenceMetricMode>(i) ) );
	}
	return outvec;
}


/// @brief Given an output mode enum, get its string representation.
/// @details Returns "INVALID" if invalid.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
std::string const &
SequenceMetric::mode_name_from_enum(
	SequenceMetricMode const mode_enum
) {
	if ( static_cast<core::Size>(mode_enum) == 0 || mode_enum >= SMM_END_OF_LIST ) {
		return mode_to_name_map_.at( SMM_INVALID_MODE );
	}
	return mode_to_name_map_.at( mode_enum );
}


/// @brief Given an output mode string, get its enum.
/// @details Returns SMM_INVALID_MODE if invalid.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
SequenceMetricMode
SequenceMetric::mode_enum_from_name(
	std::string const & mode_string
) {
	for ( core::Size i(1); i < static_cast<core::Size>(SMM_END_OF_LIST); ++i ) {
		if ( !mode_string.compare( mode_name_from_enum( static_cast< SequenceMetricMode >( i ) ) ) ) {
			return static_cast< SequenceMetricMode >(i);
		}
	}
	return SMM_INVALID_MODE;
}

void
SequenceMetric::set_residue_selector(select::residue_selector::ResidueSelectorCOP selector ){
	selector_ = selector;
}

/// @brief Set the output mode -- one-letter code (e.g. Y), three-letter code (e.g. DTY), residue base name (e.g. DTYR), or
/// full residue name (e.g. DTYR:CtermProteinFull).
/// @details Throws an error if invalid enum provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
void
SequenceMetric::set_output_mode(
	SequenceMetricMode const mode_in
) {
	runtime_assert_string_msg( mode_in > 0 && mode_in < SMM_END_OF_LIST, "Error in SequenceMetric::set_output_mode(): An unknown output mode was provided." );
	output_mode_ = mode_in;
}

/// @brief Set the output mode using the string corresponding to the output mode.
/// @details Throws an error if string is invalid.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
void
SequenceMetric::set_output_mode(
	std::string const & mode_in
) {
	set_output_mode( mode_enum_from_name(mode_in) ); //Will throw an error if string is invalid.
}

void
SequenceMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap)
{
	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(core::select::residue_selector::parse_residue_selector( tag, datamap ));
	}
	runtime_assert_string_msg( !tag->hasOption("three_letter"), "Error in SequenceMetric::parse_my_tag(): The \"three_letter\" option has been deprecated.  Please use the \"output_mode\" option to specify the output mode." );

	set_output_mode( tag->getOption<std::string>( "output_mode", mode_name_from_enum(SMM_ONELETTER_CODE) ) );
}

void
SequenceMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default ("three_letter", xsct_rosetta_bool, "DEPRECATED.  Use of this option will trigger a runtime error.  Use the output_mode option instead.", "false");

	utility::vector1< std::string > const allowed_modes_vec( allowed_output_modes_as_vector() );
	utility::tag::add_schema_restrictions_for_strings( xsd, "SequenceMetric_output_modes", allowed_modes_vec );
	attlist + XMLSchemaAttribute::attribute_w_default( "output_mode", "SequenceMetric_output_modes", "The format for the sequence.  Allowed output formats are: " + allowed_output_modes() + ".", mode_name_from_enum( SMM_ONELETTER_CODE ) );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Output the sequence of only the selected residues." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"A metric for measuring ... and adding it to the resulting score file.", attlist);
}

std::string
SequenceMetric::calculate(const pose::Pose & pose ) const {

	utility::vector1< bool > subset( pose.size(), true);
	if ( selector_ ) {
		subset = selector_->apply( pose );
	}

	std::string sequence = "";
	utility::vector1< core::Size > selection = select::get_residues_from_subset( subset );
	//std::cout << "Selection"  << utility::to_string(selection);

	for ( core::Size resnum : selection ) {
		switch( output_mode_ ) {
		case SMM_ONELETTER_CODE :
			sequence += pose.residue_type( resnum ).name1();
			break;
		case SMM_THREELETTER_CODE :
			sequence += pose.residue_type( resnum ).name3();
			if ( resnum != selection.size() ) {
				sequence += ",";
			}
			break;
		case SMM_BASE_NAME :
			sequence += pose.residue_type( resnum ).base_name();
			if ( resnum != selection.size() ) {
				sequence += ",";
			}
			break;
		case SMM_FULL_NAME :
			sequence += pose.residue_type( resnum ).name();
			if ( resnum != selection.size() ) {
				sequence += ",";
			}
			break;
		case SMM_INVALID_MODE :
			utility_exit_with_message( "Invalid configuration for SequenceMetric!" );
			break;
		}
	}
	return sequence;
}

void
SequenceMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SequenceMetric::provide_xml_schema( xsd );
}

std::string
SequenceMetricCreator::keyname() const {
	return SequenceMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SequenceMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< SequenceMetric >();

}

} //core
} //simple_metrics
} //metrics






