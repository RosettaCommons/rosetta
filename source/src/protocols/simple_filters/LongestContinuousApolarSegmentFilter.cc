// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/LongestContinuousApolarSegmentFilter.cc
/// @brief This filter computes the longest continuous stretch of apolar residues within a pose or selection.
/// @author Yang Hsia (yhsia@uw.edu)

#include <protocols/simple_filters/LongestContinuousApolarSegmentFilter.hh>
#include <protocols/simple_filters/LongestContinuousApolarSegmentFilterCreator.hh>

//Core includes
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>

//Protocols includes
#include <protocols/rosetta_scripts/util.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR( "protocols.simple_filters.LongestContinuousApolarSegmentFilter" );

namespace protocols {
namespace simple_filters {

/// @brief Constructor of parent
LongestContinuousApolarSegmentFilter::LongestContinuousApolarSegmentFilter():
	LongestContinuousPolarSegmentFilter( "LongestContinuousApolarSegmentFilter" )
{}

/// @brief Destructor.
LongestContinuousApolarSegmentFilter::~LongestContinuousApolarSegmentFilter() = default;

protocols::filters::FilterOP
LongestContinuousApolarSegmentFilter::clone() const
{
	return protocols::filters::FilterOP( new LongestContinuousApolarSegmentFilter( *this ) );
}

protocols::filters::FilterOP
LongestContinuousApolarSegmentFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new LongestContinuousApolarSegmentFilter );
}

std::string LongestContinuousApolarSegmentFilter::name() const {
	return class_name();
}

std::string LongestContinuousApolarSegmentFilter::class_name() {
	return "LongestContinuousApolarSegment";
}

void LongestContinuousApolarSegmentFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default("exclude_chain_termini", xsct_rosetta_bool, "If true, apolar stretches at the ends of chains are not counted.  If false, they are.  True by default.  Note that if this is set to true, an entirely apolar chain will not be counted.", "true")
		+ XMLSchemaAttribute::attribute_w_default("filter_out_high", xsct_rosetta_bool, "If true, poses with more than the cutoff number of residues in the longest apolar stretch will be rejected.  If false, poses with fewer than the cutoff number of residues in the longest apolar stretch will be rejected.  True by default.", "true")
		+ XMLSchemaAttribute::attribute_w_default( "cutoff", xsct_non_negative_integer, "The maximum (or minimum, if \"filter_out_high\" is set to \"false\") number of residues in the longest apolar stretch that will still allow the pose to pass this filter.  Default 5.", "5" )
		;

	protocols::rosetta_scripts::attributes_for_parse_residue_selector( attlist, "An optional, previously-defined residue selector.  If provided, the filter will only consider stretches of apolar residues that have at least one residue in the selection.  Not used if not specified." );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"The LongestContinuousApolarSegment filter counts the number of amino acid residues in the longest continuous stretch of apolar amino acids in any chain of your pose (or, optionally, in a residue selection).  Optionally, it can ignore the polar stretches at the N- and C-termini of chains.",
		attlist );
}

/// @brief Given a residue type, determine whether it's one of the types that this filter
/// should count.
/// @details Based on whether the residue type has the APOLAR property. APOLAR is defined as "not POLAR".
bool
LongestContinuousApolarSegmentFilter::is_counted(
	core::chemical::ResidueType const &restype
) const {
	if ( ! restype.is_polar() ) return true;
	return false;
}

std::string
LongestContinuousApolarSegmentFilter::counted_residue_description() const {
	return "apolar";
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
LongestContinuousApolarSegmentFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new LongestContinuousApolarSegmentFilter );
}

std::string
LongestContinuousApolarSegmentFilterCreator::keyname() const
{
	return LongestContinuousApolarSegmentFilter::class_name();
}

void LongestContinuousApolarSegmentFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LongestContinuousApolarSegmentFilter::provide_xml_schema( xsd );
}

} //protocols
} //simple_filters
