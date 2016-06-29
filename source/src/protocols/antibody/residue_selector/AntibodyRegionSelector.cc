// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/AntibodyRegionSelector.hh
/// @brief  A simple selector to select residues of particular antibody regions.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/antibody/residue_selector/AntibodyRegionSelector.hh>
#include <protocols/antibody/residue_selector/AntibodyRegionSelectorCreator.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/util.hh>

#include <core/pose/Pose.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/util.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.residue_selector.AntibodyRegionSelector" );


namespace protocols {
namespace antibody {
namespace residue_selector {

using namespace core::select::residue_selector;
using namespace basic::options;

/// @brief Constructor.
///
AntibodyRegionSelector::AntibodyRegionSelector():
	ResidueSelector(),
	ab_info_(/* NULL */)
{

	set_defaults();

}

AntibodyRegionSelector::AntibodyRegionSelector( AntibodyInfoCOP ab_info ):
	ResidueSelector(),
	ab_info_(ab_info)
{
	set_defaults();
}

AntibodyRegionSelector::AntibodyRegionSelector( AntibodyInfoCOP ab_info, AntibodyRegionEnum region ):
	ResidueSelector(),
	ab_info_(ab_info)
{
	set_defaults();
	region_ = region;
}


void
AntibodyRegionSelector::set_defaults(){

	region_ = unknown_ab_region;

	AntibodyEnumManager manager = AntibodyEnumManager();
	std::string numbering_scheme = option [OptionKeys::antibody::numbering_scheme]();
	std::string cdr_definition = option [OptionKeys::antibody::cdr_definition]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);
	cdr_definition_ = manager.cdr_definition_string_to_enum(cdr_definition);
}


/// @brief Destructor.
///
AntibodyRegionSelector::~AntibodyRegionSelector() {}

AntibodyRegionSelector::AntibodyRegionSelector( AntibodyRegionSelector const & src):
	ResidueSelector(),
	ab_info_(src.ab_info_),
	region_(src.region_),
	numbering_scheme_(src.numbering_scheme_),
	cdr_definition_(src.cdr_definition_)

{


}



/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
AntibodyRegionSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast<core::select::residue_selector::ResidueSelector>(
		AntibodyRegionSelectorOP( new AntibodyRegionSelector(*this) )
		)
	);
}



/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
AntibodyRegionSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	AntibodyEnumManager manager = AntibodyEnumManager();

	if ( tag->hasOption("region") ) {
		region_ = manager.antibody_region_string_to_enum(tag->getOption< std::string >( "region" ));
	}

	if ( tag->hasOption("cdr_definition") && tag->hasOption("numbering_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();

		cdr_definition_ = manager.cdr_definition_string_to_enum(tag->getOption<std::string>("cdr_definition"));
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("numbering_scheme"));

	} else if ( tag->hasOption("cdr_definition") || tag->hasOption("numbering_scheme") ) {
		TR <<"Please pass both cdr_definition and numbering_scheme.  These can also be set via cmd line options of the same name." << std::endl;

	}

}

std::string AntibodyRegionSelector::get_name() const
{
	return AntibodyRegionSelector::class_name();
}

std::string AntibodyRegionSelector::class_name()
{
	return "AntibodyRegion";
}

void AntibodyRegionSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "region",      xs_boolean, "true" )
		+ XMLSchemaAttribute::attribute_w_default(  "cdr_definition", xs_string, option [OptionKeys::antibody::cdr_definition]() )
		+ XMLSchemaAttribute::attribute_w_default(  "numbering_scheme", xs_string, option [OptionKeys::antibody::numbering_scheme]());
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );

}

ResidueSelectorOP
AntibodyRegionSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
		AntibodyRegionSelectorOP( new AntibodyRegionSelector )
		)
	);
}

std::string
AntibodyRegionSelectorCreator::keyname() const {
	return AntibodyRegionSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
AntibodyRegionSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	AntibodyRegionSelector::provide_xml_schema( xsd );
}

void
AntibodyRegionSelector::set_region(AntibodyRegionEnum region){
	region_ = region;
}

void
AntibodyRegionSelector::set_ab_info(AntibodyInfoCOP ab_info){
	ab_info_ = ab_info;
}


/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResidueSubset
AntibodyRegionSelector::apply(
	core::pose::Pose const & pose
) const {

	if ( region_ == unknown_ab_region ) {
		utility_exit_with_message("AntibodyRegionSelector: Please set the region to select a region.");
	}

	utility::vector1< bool > subset( pose.total_residue(), false);

	//This is due to const apply and no pose in parse_my_tag.
	AntibodyInfoOP local_ab_info;
	if ( ! ab_info_ ) {
		local_ab_info = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, cdr_definition_));
	} else {
		local_ab_info = ab_info_->clone();
	}


	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( local_ab_info->get_region_of_residue(pose, i) == region_ ) {
			subset[ i ] = true;
		}
	}


	return subset;


}

} //protocols
} //antibody
} //residue_selector
