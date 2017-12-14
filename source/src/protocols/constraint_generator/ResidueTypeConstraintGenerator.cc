// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/ResidueTypeConstraintGenerator.cc
/// @brief Generates residue type constraints for a set of residues from the current or reference pose
/// @author Sharon Guffy (guffy@email.unc.edu)



// Unit headers
#include <protocols/constraint_generator/ResidueTypeConstraintGenerator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>
#include <protocols/constraint_generator/ResidueTypeConstraintGeneratorCreator.hh>
#include <protocols/constraint_generator/util.hh> //I'll need this later for sequence mapping to/from reference poses
// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/types.hh>
#include <core/id/SequenceMapping.hh>
//Basic headers
#include <basic/Tracer.hh>
// Utility headers
#include <utility/fixedsizearray1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>


static basic::Tracer TR( "protocols.constraint_generator.ResidueTypeConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

///@brief Generates atom pair constraints for a set of residues from the current or reference pose
ResidueTypeConstraintGenerator::ResidueTypeConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( ResidueTypeConstraintGenerator::class_name() ),
	selector_( core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::TrueResidueSelector ) ),
	rsd_type_name3_( "" )
{}

ResidueTypeConstraintGenerator::~ResidueTypeConstraintGenerator() = default;

ConstraintGeneratorOP
ResidueTypeConstraintGenerator::clone() const{
	return ConstraintGeneratorOP( new ResidueTypeConstraintGenerator );
}

core::scoring::constraints::ConstraintCOPs
ResidueTypeConstraintGenerator::apply( core::pose::Pose const & pose ) const{
	runtime_assert( selector_ != nullptr );
	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );
	core::scoring::constraints::ConstraintCOPs csts;
	if ( ref_pose_ == nullptr ) {
		TR << "No reference pose detected!" << std::endl;
		csts = generate_residue_type_constraints( pose, subset );
	} else {
		core::id::SequenceMapping seqmap = generate_seqmap_from_poses( pose, *ref_pose_ );
		TR << "Generating constraints from native reference" << std::endl;
		csts = generate_residue_type_constraints( pose, *ref_pose_, subset, seqmap );
	}
	return csts;
}

void
ResidueTypeConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ){
	bool const use_native = tag->getOption< bool >( "native", false );
	if ( use_native ) {
		set_reference_pose( get_native_pose() );
		if ( ! ref_pose_ ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "'native' option for ResidueTypeConstraintGenerator specified, but no native pose is availible." );
		}
		//Any other source of reference poses?
	}
	favor_native_bonus_ = tag->getOption< core::Real >( "favor_native_bonus", 1.0 );
	rsd_type_name3_ = tag->getOption< std::string >( "rsd_type_name3", "" );
	//Option name for selector is "residue_selector"
	core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) {
		selector_ = selector;
	}
	//Otherwise we'll constrain all residues
	//TEMP for debugging purposes
	if ( !tag->hasOption( "residue_selector" ) ) {
		runtime_assert( selector_->get_name() == "True" );
	}
	//END TEMP

}

//Setters
void
ResidueTypeConstraintGenerator::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector ){
	selector_ = selector->clone();
}

void
ResidueTypeConstraintGenerator::set_favor_native_bonus( core::Real bonus ){
	favor_native_bonus_ = bonus;
}

void
ResidueTypeConstraintGenerator::set_rsd_type_name3( std::string name3 ){
	rsd_type_name3_ = name3;
}

void
ResidueTypeConstraintGenerator::set_reference_pose( core::pose::PoseCOP ref ){
	ref_pose_ = ref;
}

//Getters
core::select::residue_selector::ResidueSelectorCOP
ResidueTypeConstraintGenerator::get_residue_selector() const{
	return selector_;
}


core::Real
ResidueTypeConstraintGenerator::get_favor_native_bonus() const{
	return favor_native_bonus_;
}


std::string
ResidueTypeConstraintGenerator::get_rsd_type_name3() const{
	return rsd_type_name3_;
}

core::pose::PoseCOP
ResidueTypeConstraintGenerator::get_reference_pose() const{
	return ref_pose_;
}



core::scoring::constraints::ConstraintCOPs
ResidueTypeConstraintGenerator::generate_residue_type_constraints(
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose,
	core::select::residue_selector::ResidueSubset const & subset,
	core::id::SequenceMapping const & seqmap ) const
{
	core::scoring::constraints::ConstraintCOPs csts;
	core::Size nres = compute_nres_in_asymmetric_unit( pose );
	for ( core::Size i = 1; i <= nres; ++i ) {
		//Skip residues that are not being constrained
		if ( !subset[ i ] ) {
			continue;
		}
		//Skip virtual residues
		if ( pose.residue( i ).aa() == core::chemical::aa_vrt ) {
			continue;
		}
		//Otherwise add constraint to this residue
		core::scoring::constraints::ConstraintCOP new_cst;
		if ( rsd_type_name3_ != "" ) { //In this case the native is actually not even used
			TR << "WARNING: You specified both use_native and a specific residue identity to use for your constraint. Only the specified residue identity will be used; ignoring native." << std::endl;
			new_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::ResidueTypeConstraint( pose, i, rsd_type_name3_, favor_native_bonus_ ) );
		} else {
			//TR << "Getting desired identity from reference pose" << std::endl;
			core::Size j = seqmap.get_corresponding_residue_in_current( i );
			if ( j == 0 ) {
				TR << "Residue " << i << " does not have a  corresponding native residue! Skipping." << std::endl;
				continue;
			}
			//TR << "Pose residue type: " << pose.residue_type( i ).name3() << " Refpose residue type: " << ref_pose.residue_type( j ).name3() << std::endl;
			new_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::ResidueTypeConstraint( pose, i, ref_pose.residue_type( j ).name3(), favor_native_bonus_ ) );
		}
		csts.push_back( new_cst );
	}//end for i <= nres
	return csts;
}



void
ResidueTypeConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "favor_native_bonus", xsct_real, "Unweighted score bonus the pose will receive for having the specified residue type at the specified position", "1.0" )
		+ XMLSchemaAttribute( "rsd_type_name3", xs_string, "Three-letter code for the amino acid to which you want to constrain this residue. If unspecified, this defaults to the native amino acid at this position." )
		+ XMLSchemaAttribute::attribute_w_default( "use_native", xsct_rosetta_bool, "Use native structure (provided with in:file:native) as reference pose for defining desired residue identities", "false" );
	core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name( attributes, "Selector specifying residues to be constrained. When not provided, all residues are selected" );
	ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		"Generates residue type constraints (either to native or to a specific residue type) for the specified residues",
		attributes );
}
//Call this function if we do not have a reference pose
core::scoring::constraints::ConstraintCOPs
ResidueTypeConstraintGenerator::generate_residue_type_constraints(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & subset ) const
{
	core::scoring::constraints::ConstraintCOPs csts;
	core::Size nres = compute_nres_in_asymmetric_unit( pose );
	for ( core::Size i = 1; i <= nres; ++i ) {
		//Skip residues that are not being constrained
		if ( !subset[ i ] ) {
			continue;
		}
		//Skip virtual residues
		if ( pose.residue( i ).aa() == core::chemical::aa_vrt ) {
			continue;
		}
		//Otherwise add constraint to this residue
		core::scoring::constraints::ConstraintCOP new_cst;
		if ( rsd_type_name3_ != "" ) {
			new_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::ResidueTypeConstraint( pose, i, rsd_type_name3_, favor_native_bonus_ ) );
		} else {
			new_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::ResidueTypeConstraint( pose, i, favor_native_bonus_ ) );
		}
		csts.push_back( new_cst );
	}//end for i <= nres
	return csts;
}


//************Creator Functions***************

protocols::constraint_generator::ConstraintGeneratorOP
ResidueTypeConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new ResidueTypeConstraintGenerator );
}

std::string
ResidueTypeConstraintGeneratorCreator::keyname() const
{
	return ResidueTypeConstraintGenerator::class_name();
}

void
ResidueTypeConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	ResidueTypeConstraintGenerator::provide_xml_schema( xsd );
}

//**********End Creator Functions*************




} //protocols
} //constraint_generator
