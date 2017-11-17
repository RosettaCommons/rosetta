// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/SetupMetalsMover.hh
/// @brief  class definition for SetupMetalsMover
/// @author Sharon Guffy (guffy@email.unc.edu)


// Unit headers
#include <protocols/simple_moves/SetupMetalsMover.hh>
#include <protocols/simple_moves/SetupMetalsMoverCreator.hh>

// Protocols headers
#include <protocols/moves/mover_schemas.hh>
//Core headers
#include <core/util/metalloproteins_util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
//Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// Utility Headers
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace simple_moves {

SetupMetalsMover::SetupMetalsMover():
	protocols::moves::Mover()
{
	//prevent_setup_metal_bb_variants_ = false; //default false
	metals_detection_LJ_multiplier_ = 1.0; //default 1.0
	metals_distance_constraint_multiplier_ = 1.0; //default 1.0
	metals_angle_constraint_multiplier_ = 1.0; //default 1.0
	metal_selector_ = nullptr;
	metal_resnums_string_ = "";
	remove_hydrogens_ = true;
	constraints_only_ = false;
	add_constraints_ = true;
	this->set_defaults_from_command_line();
}
/// @brief Copy constructor
SetupMetalsMover::SetupMetalsMover( SetupMetalsMover const & src ):
	protocols::moves::Mover( src )
{
	//prevent_setup_metal_bb_variants_ = src.get_prevent_setup_metal_bb_variants();
	metals_detection_LJ_multiplier_ = src.get_metals_detection_LJ_multiplier();
	metals_distance_constraint_multiplier_ = src.get_metals_distance_constraint_multiplier();
	metals_angle_constraint_multiplier_ = src.get_metals_angle_constraint_multiplier();
	metal_selector_ = src.get_metal_selector();
	metal_resnums_string_ = src.get_metal_resnums_string();
	remove_hydrogens_ = src.get_remove_hydrogens();
	constraints_only_ = src.get_constraints_only();
	add_constraints_ = src.add_constraints_;
	contact_selector_ = src.contact_selector_;
	contact_resnums_string_ = src.contact_resnums_string_;
}

SetupMetalsMover::~SetupMetalsMover(){}


//General mover methods
void
SetupMetalsMover::apply( core::pose::Pose & pose ){
	//If no residue selector or resnum string were provided, just call the standard methods to set up all the metals
	if ( metal_resnums_string_ == "" && !metal_selector_ ) {
		if ( !constraints_only_ ) {
			core::util::auto_setup_all_metal_bonds( pose, metals_detection_LJ_multiplier_, remove_hydrogens_ );
		}
		core::util::auto_setup_all_metal_constraints( pose, metals_distance_constraint_multiplier_, metals_angle_constraint_multiplier_ );
	} else {
		//This takes the provided selector/resnums into account
		utility::vector1< core::Size > metals = find_metal_resnums( pose );
		for ( core::Size resn: metals ) {
			//Get contacts
			utility::vector1< core::id::AtomID > binders;
			utility::vector1< core::id::AtomID > binders_pre = core::util::find_metalbinding_atoms( pose, resn, metals_detection_LJ_multiplier_ );
			std::set< core::Size > contact_resnums = find_contact_resnums( pose );
			if ( contact_resnums.size() == 0 ) {
				binders = binders_pre;
			} else {
				for ( core::id::AtomID id: binders_pre ) {
					if ( contact_resnums.count( id.rsd() ) != 0 ) {
						binders.push_back( id );
					}
				}
			}
			if ( !constraints_only_ ) {
				//Make chemical bonds
				core::util::add_covalent_linkages_to_metal( pose, resn, binders, remove_hydrogens_ );
			}
			//Add constraints
			if ( add_constraints_ ) {
				//TODO!!
				core::util::add_constraints_to_metal( pose, resn, metals_distance_constraint_multiplier_, metals_angle_constraint_multiplier_, binders );
			}
		}
	}
}

void
SetupMetalsMover::show(std::ostream & output) const{
	Mover::show( output );
	output << "Metals detection LJ multiplier: " << metals_detection_LJ_multiplier_ << std::endl;
	output << "Metals distance constraint multiplier: " << metals_distance_constraint_multiplier_ << std::endl;
	output << "Metals angle constraint multiplier: " << metals_angle_constraint_multiplier_ << std::endl;
	//output << "Prevent setup metal backbone variants: " << prevent_setup_metal_bb_variants_ << std::endl;
	output<< "Remove hydrogens: " << remove_hydrogens_ << std::endl;
	output << "Constraints only: " << constraints_only_ << std::endl;
	if ( metal_resnums_string_ != "" ) {
		output << "Metal resnums: " << metal_resnums_string_ << std::endl;
	}
	if ( metal_selector_ ) { //There is no virtual show method for all residue selectors, but this will at least show that a selector is being used
		output << "Metal selector: " << metal_selector_->get_name() << std::endl;
	}
}

protocols::moves::MoverOP
SetupMetalsMover::clone() const{
	return protocols::moves::MoverOP( new SetupMetalsMover( *this ) );
}
protocols::moves::MoverOP
SetupMetalsMover::fresh_instance() const{
	return protocols::moves::MoverOP( new SetupMetalsMover );
}

void
SetupMetalsMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{
	//Tag should be able to take either a named selector or resnum string for metal or contact
	if ( tag->hasOption("metal_residue_selector" ) ) {
		if ( tag->hasOption("metal_resnums") ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "SetupMetalsMover takes EITHER a residue selector or resnum list, not both!\n" );
		}
		metal_selector_ = core::select::residue_selector::parse_residue_selector( tag, data, "metal_residue_selector" );
	} else if ( tag->hasOption( "metal_resnums" ) ) {
		metal_resnums_string_ = tag->getOption< std::string >( "metal_resnums", "" );
	}


	if ( tag->hasOption("contact_residue_selector" ) ) {
		if ( tag->hasOption("contact_resnums") ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "SetupMetalsMover takes EITHER a residue selector or resnum list, not both!\n" );
		}
		contact_selector_ = core::select::residue_selector::parse_residue_selector( tag, data, "contact_residue_selector" );
	} else if ( tag->hasOption( "contact_resnums" ) ) {
		contact_resnums_string_ = tag->getOption< std::string >( "contact_resnums", "" );
	}

	//Now set other options from tag
	remove_hydrogens_ = tag->getOption< bool >( "remove_hydrogens", true );
	constraints_only_ = tag->getOption< bool >( "constraints_only", false );
	add_constraints_ = tag->getOption< bool >( "add_constraints", true );
	//Defaults for these options are taken from the command line
	//this->set_defaults_from_command_line(); Moved to constructor
	if ( tag->hasOption( "metals_detection_LJ_multiplier" ) ) {
		metals_detection_LJ_multiplier_ = tag->getOption< core::Real >( "metals_detection_LJ_multiplier", 1.0 );
	}
	if ( tag->hasOption( "metals_distance_constraint_multiplier" ) ) {
		metals_distance_constraint_multiplier_ = tag->getOption< core::Real >( "metals_distance_constraint_multiplier", 1.0 );
	}
	if ( tag->hasOption( "metals_angle_constraint_multiplier" ) ) {
		metals_angle_constraint_multiplier_ = tag->getOption< core::Real >( "metals_angle_constraint_multiplier", 1.0 );
	}

}

//Static methods
void
SetupMetalsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "remove_hydrogens", xsct_rosetta_bool, "Should hydrogens on metal binding atoms be removed?", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "constraints_only", xsct_rosetta_bool, "Should we ONLY add the constraints to the metal and not set up covalent bonds? Useful for repeated applications of the mover if constraints are deleted.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "add_constraints", xsct_rosetta_bool, "Should we add constraints to the pose? May be used if, for example, the user wants to add their own custom constraints elsewhere.", "true" )
		//+ XMLSchemaAttribute::attribute_w_default( "prevent_setup_metal_bb_variants", xsct_rosetta_bool, "Should we prevent the setup of variant types for all possible residue types in cases where backbone atoms coordinate the metal? False indicates that the residue types WILL be set up.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "metals_detection_LJ_multiplier", xsct_real, "Multiplier for the distance cutoff for contact detection, which is based on the Lennard-Jones radii of the two atoms. This option overrides default value set through the command line.", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "metals_distance_constraint_multiplier", xsct_real, "Multiplier for distance constraints between metal and metal-binding atoms. This option overrides default value set through the command line.", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "metals_angle_constraint_multiplier", xsct_real, "Multiplier for angle constraints between metal and metal-binding atoms. This option overrides default value set through the command line.", "1.0" )
		+ utility::tag::optional_name_attribute();
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "metal_residue_selector", "Selector used to indicate which metal residues should be set up. If no selector or resnum string is provided, all metals will be set up." );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "contact_residue_selector", "Selector used to indicate which metal contacts should be set up. If no selector or resnum string is provided, all contacts for the specified metals will be set up." );
	core::pose::attributes_for_get_resnum_selector( attlist, xsd, "metal_resnums" );
	core::pose::attributes_for_get_resnum_selector( attlist, xsd, "contact_resnums" );
	//We can provide 0 or 1 subselector
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.element_name( SetupMetalsMover::mover_name() )
		.description( "Mover that adds chemical bonds and distance/angle constraints between metal ions and their coordinating atoms. If a residue selector or resnum list are provided, only sets up the specified metals. Otherwise, sets up all metals."  )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );
}


std::string
SetupMetalsMover::mover_name(){
	return "SetupMetalsMover";
}

//Getters
/*
bool
SetupMetalsMover::get_prevent_setup_metal_bb_variants() const{
return prevent_setup_metal_bb_variants_;
}
*/

bool
SetupMetalsMover::get_constraints_only() const{
	return constraints_only_;
}

core::Real
SetupMetalsMover::get_metals_detection_LJ_multiplier() const{
	return metals_detection_LJ_multiplier_;
}

core::Real
SetupMetalsMover::get_metals_distance_constraint_multiplier() const{
	return metals_distance_constraint_multiplier_;
}

core::Real
SetupMetalsMover::get_metals_angle_constraint_multiplier() const{
	return metals_angle_constraint_multiplier_;
}

core::select::residue_selector::ResidueSelectorCOP
SetupMetalsMover::get_metal_selector() const{
	return metal_selector_;
}

std::string
SetupMetalsMover::get_metal_resnums_string() const{
	return metal_resnums_string_;
}

bool
SetupMetalsMover::get_remove_hydrogens() const{
	return remove_hydrogens_;
}

bool
SetupMetalsMover::get_add_constraints() const{
	return add_constraints_;
}

core::select::residue_selector::ResidueSelectorCOP
SetupMetalsMover::get_contact_selector() const{
	return contact_selector_;
}

std::string
SetupMetalsMover::get_contact_resnums_string() const{
	return contact_resnums_string_;
}



//Setters
/*
void
SetupMetalsMover::set_prevent_setup_metal_bb_variants( bool input ){
prevent_setup_metal_bb_variants_ = input;
}
*/

void
SetupMetalsMover::set_constraints_only( bool input ){
	constraints_only_ = input;
}

void
SetupMetalsMover::set_metals_detection_LJ_multiplier( core::Real input ){
	metals_detection_LJ_multiplier_ = input;
}

void
SetupMetalsMover::set_metals_distance_constraint_multiplier( core::Real input ){
	metals_distance_constraint_multiplier_ = input;
}

void
SetupMetalsMover::set_metals_angle_constraint_multiplier( core::Real input){
	metals_angle_constraint_multiplier_ = input;
}

void
SetupMetalsMover::set_metal_selector( core::select::residue_selector::ResidueSelectorCOP input ){
	metal_selector_ = input;
}

void
SetupMetalsMover::set_metal_resnums_string( std::string input ){
	metal_resnums_string_ = input;
}

void
SetupMetalsMover::set_remove_hydrogens( bool input ){
	remove_hydrogens_ = input;
}

void
SetupMetalsMover::set_add_constraints( bool input){
	add_constraints_ = input;
}

void
SetupMetalsMover::set_contact_selector( core::select::residue_selector::ResidueSelectorCOP input){
	contact_selector_ = input;
}

void
SetupMetalsMover::set_contact_resnums_string(std::string input){
	contact_resnums_string_ = input;
}


//Begin protected methods

///@brief If a residue selector or resnum string is provided, returns any metal ions within the selection.
///Otherwise returns an empty vector (it shouldn't be called in this situation, though)
utility::vector1< core::Size >
SetupMetalsMover::find_metal_resnums( core::pose::Pose const & pose ){
	utility::vector1< core::Size > metals;
	core::select::residue_selector::ResidueSubset residues( pose.total_residue(), false );
	if ( metal_selector_ ) {
		residues = metal_selector_->apply( pose );
	}
	std::set< core::Size > str_residues;
	if ( metal_resnums_string_ != "" ) {
		str_residues = core::pose::get_resnum_list( metal_resnums_string_, pose );
	}
	for ( core::Size resn=1; resn<=pose.size(); ++resn ) { //Loop through all residues.
		if ( pose.residue(resn).is_metal() && (residues[ resn ] || ( str_residues.count( resn ) != 0 ) ) ) {
			metals.push_back( resn );
		}
	}
	return metals;
}

std::set< core::Size >
SetupMetalsMover::find_contact_resnums( core::pose::Pose const & pose ){
	std::set< core::Size > contacts;
	core::select::residue_selector::ResidueSubset residues( pose.total_residue(), false );
	if ( contact_selector_ ) {
		residues = contact_selector_->apply( pose );
	}
	std::set< core::Size > str_residues;
	if ( contact_resnums_string_ != "" ) {
		str_residues = core::pose::get_resnum_list( contact_resnums_string_, pose );
	}

	for ( core::Size resn=1; resn<=pose.size(); ++resn ) { //Loop through all residues.
		if ( residues[ resn ] || ( str_residues.count( resn ) != 0 ) ) {
			contacts.insert( resn );
		}
	}
	return contacts;
}


void
SetupMetalsMover::set_defaults_from_command_line(){
	using namespace basic::options::OptionKeys;
	metals_detection_LJ_multiplier_ = basic::options::option[ in::metals_detection_LJ_multiplier ].value();
	metals_distance_constraint_multiplier_ = basic::options::option[ in::metals_distance_constraint_multiplier ].value();
	metals_angle_constraint_multiplier_ = basic::options::option[ in::metals_angle_constraint_multiplier ].value();
	//prevent_setup_metal_bb_variants_ = basic::options::option[ in::prevent_auto_setup_metal_bb_variants ].value();
}

std::ostream &operator<< (std::ostream &os, SetupMetalsMover const &mover)
{
	mover.show( os );
	return os;
}

//CREATOR METHODS
protocols::moves::MoverOP
SetupMetalsMoverCreator::create_mover() const{
	return protocols::moves::MoverOP( new SetupMetalsMover );
}
std::string
SetupMetalsMoverCreator::keyname() const{
	return SetupMetalsMover::mover_name();
}
void
SetupMetalsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	SetupMetalsMover::provide_xml_schema( xsd );
}
std::string
SetupMetalsMoverCreator::mover_name(){
	return SetupMetalsMover::mover_name();
}

//END CREATOR METHODS

}  // simple_moves
}  // protocols
