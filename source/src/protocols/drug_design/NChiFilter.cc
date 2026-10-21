// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/NChiFilter.cc
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)


//Unit Headers
#include <protocols/drug_design/NChiFilter.hh>
#include <protocols/drug_design/NChiFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace drug_design {

static basic::Tracer TR( "protocols.drug_design.NChiFilter" );

protocols::filters::FilterOP
NChiFilterCreator::create_filter() const { return protocols::filters::FilterOP( new NChiFilter ); }

std::string
NChiFilterCreator::keyname() const { return NChiFilter::class_name(); }

void NChiFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NChiFilter::provide_xml_schema( xsd );
}

std::string NChiFilter::class_name() {
	return "NChi";
}

void NChiFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "residue", xsct_refpose_enabled_residue_number, "Residue for which to calculate the number of rotatable bonds" )
		+ XMLSchemaAttribute::attribute_w_default( "exclude_proton", xsct_rosetta_bool, "Don't count proton chis", "true" )
		+ XMLSchemaAttribute( "threshold", xsct_real, "Fail if the number of rotatable bonds is greater than this." );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(),
		"Count the number of rotatable bonds (as determined by Rosetta).",
		attlist );
}

void
NChiFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	residue_ = tag->getOption< std::string >( "residue" );
	threshold_ = tag->getOption<core::Size>( "threshold", threshold_ );
	exclude_proton_ = tag->getOption<bool>( "exclude_proton", exclude_proton_ );
	TR << "NChi filter for residue " << residue_ << " with an upper threshold of  " << threshold_ << (exclude_proton_?" excluding":" including") << " proton chis." << std::endl;
}

bool
NChiFilter::apply( core::pose::Pose const & pose ) const {
	core::Size const value( compute( pose ) );
	if ( value > threshold_ ) {
		TR << "Failing NChi  filter on residue " << residue_ << " with value " << value << " which is above the threshold of " << threshold_ << std::endl;
		return false;
	} else {
		return true;
	}
}

void
NChiFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size const value( compute( pose ) );
	out << "NChi value for residue " << residue_ << " is " << value << std::endl;
}

core::Real
NChiFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Size
NChiFilter::compute( core::pose::Pose const & pose ) const {
	core::Size resnum( core::pose::parse_resnum( residue_, pose ) );
	if ( resnum == 0 || resnum > pose.total_residue() ) {
		TR.Warning << "Attempted to access residue " << residue_ << " in pose with " <<  pose.total_residue() << " residues. Failing filter. " << std::endl;
		utility_exit_with_message("Cannot apply NChi filter on non-existant residue!");
	}
	core::chemical::ResidueType const & restype( pose.residue(resnum).type() );

	core::Real retval( restype.nchi() );

	if ( exclude_proton_ ) {
		retval -= restype.n_proton_chi();
	}

	TR << "NChi for residue " << residue_ << ", of type " << restype.name() << " is " << retval << ( exclude_proton_ ? " excluding":"including") << " proton chis." << std::endl;

	return retval;
}

}
}
