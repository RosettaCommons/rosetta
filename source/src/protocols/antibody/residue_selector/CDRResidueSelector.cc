// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/CDRResidueSelector.hh
/// @brief  Select CDR residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/antibody/residue_selector/CDRResidueSelector.hh>
#include <protocols/antibody/residue_selector/CDRResidueSelectorCreator.hh>

#include <core/select/residue_selector/util.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/util.hh>

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

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>



static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.residue_selector.CDRResidueSelector" );


namespace protocols {
namespace antibody {
namespace residue_selector {

using namespace core::pack::task::operation;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility::tag;

using namespace core::select::residue_selector;

/// @brief Constructor.
///
CDRResidueSelector::CDRResidueSelector():
	ResidueSelector(),
	ab_info_(/* NULL */)
{
	set_defaults();
}

CDRResidueSelector::CDRResidueSelector( AntibodyInfoCOP ab_info ):
	ResidueSelector(),
	ab_info_(ab_info)
{
	set_defaults();
}

CDRResidueSelector::CDRResidueSelector( AntibodyInfoCOP ab_info, utility::vector1< CDRNameEnum > cdrs ):
	ResidueSelector(),
	ab_info_(ab_info)
{
	set_defaults();
	set_cdrs( cdrs );
}

CDRResidueSelector::CDRResidueSelector( AntibodyInfoCOP ab_info, utility::vector1< bool > cdrs ):
	ResidueSelector(),
	ab_info_(ab_info)
{
	set_defaults();
	set_cdrs( cdrs );
}





/// @brief Destructor.
///
CDRResidueSelector::~CDRResidueSelector() {}

CDRResidueSelector::CDRResidueSelector(CDRResidueSelector const & src):
	ResidueSelector(),
	ab_info_(src.ab_info_),
	cdrs_(src.cdrs_),
	numbering_scheme_(src.numbering_scheme_),
	cdr_definition_(src.cdr_definition_)
{

}


/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
CDRResidueSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast<core::select::residue_selector::ResidueSelector>(
		CDRResidueSelectorOP( new CDRResidueSelector(*this) )
		)
	);
}



void
CDRResidueSelector::set_cdrs( utility::vector1< bool > cdrs ){
	cdrs_ = cdrs;
}

void
CDRResidueSelector::set_cdrs( utility::vector1< CDRNameEnum > cdrs ){
	cdrs_.clear();
	cdrs_.resize(8, false);
	for ( core::Size i = 1; i < cdrs.size(); ++i ) {
		cdrs_[ cdrs[ i ] ] = true;
	}
}

void
CDRResidueSelector::set_ab_info(AntibodyInfoCOP ab_info){
	ab_info_ = ab_info;
}

void
CDRResidueSelector::set_defaults() {
	cdrs_.clear();
	cdrs_.resize(6, true);

	AntibodyEnumManager manager = AntibodyEnumManager();
	std::string numbering_scheme = option [OptionKeys::antibody::numbering_scheme]();
	std::string cdr_definition = option [OptionKeys::antibody::cdr_definition]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);
	cdr_definition_ = manager.cdr_definition_string_to_enum(cdr_definition);

}




/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
CDRResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	AntibodyEnumManager manager = AntibodyEnumManager();

	if ( tag->hasOption("cdrs") ) {
		TR << "Setting CDRs from settings" << std::endl;
		cdrs_ = get_cdr_bool_from_tag(tag, "cdrs", true /* include_cdr4*/);
	}

	if ( tag->hasOption("cdr_definition") && tag->hasOption("numbering_scheme") ) {
		AntibodyEnumManager manager = AntibodyEnumManager();
		cdr_definition_ = manager.cdr_definition_string_to_enum(tag->getOption<std::string>("cdr_definition"));
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("numbering_scheme"));
	} else if ( tag->hasOption("cdr_definition") || tag->hasOption("numbering_scheme") ) {
		TR <<"Please pass both cdr_definition and numbering_scheme.  These can also be set via cmd line options of the same name." << std::endl;
	}

}

std::string CDRResidueSelector::get_name() const
{
	return CDRResidueSelector::class_name();
}

std::string CDRResidueSelector::class_name()
{
	return "CDRResidueSelector";
}

void CDRResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "cdrs",      xs_string, "L1,L2,L3,H1,H2,H3" )
		+ XMLSchemaAttribute::attribute_w_default(  "cdr_definition", xs_string, option [OptionKeys::antibody::cdr_definition]() )
		+ XMLSchemaAttribute::attribute_w_default(  "numbering_scheme", xs_string, option [OptionKeys::antibody::numbering_scheme]());
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );

}

ResidueSelectorOP
CDRResidueSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
		CDRResidueSelectorOP( new CDRResidueSelector )
		)
	);
}

std::string
CDRResidueSelectorCreator::keyname() const {
	return CDRResidueSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
CDRResidueSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	CDRResidueSelector::provide_xml_schema( xsd );
}



/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResidueSubset
CDRResidueSelector::apply(
	core::pose::Pose const & pose
) const {

	utility::vector1< bool > subset(pose.total_residue(), false);

	//This is due to const apply and no pose in parse_my_tag.
	AntibodyInfoOP local_ab_info;
	if ( ! ab_info_ ) {
		local_ab_info = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, cdr_definition_));
	} else {
		local_ab_info = ab_info_->clone();
	}

	for ( core::Size i = 1; i <= cdrs_.size(); ++i ) {

		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );


		if ( ! cdrs_[ i ] ) continue;
		if ( local_ab_info->is_camelid() && local_ab_info->get_CDR_chain( cdr ) == 'L' ) continue;



		core::Size start = local_ab_info->get_CDR_start( cdr, pose );
		core::Size end = local_ab_info->get_CDR_end( cdr, pose );

		for ( core::Size resnum = start; resnum <= end; ++resnum ) {
			subset[ resnum ] = true;
		}
	}


	return subset;


}






} //protocols
} //antibody
} //residue_selector
