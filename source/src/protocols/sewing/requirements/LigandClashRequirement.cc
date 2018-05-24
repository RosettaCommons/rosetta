// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/LigandClashRequirement.cc
/// @brief Assembly requirement that clash checks the assembly backbone with its bound ligands
/// @author Minnie Langlois (minnie@email.unc.edu)

#include <protocols/sewing/requirements/LigandClashRequirement.hh>
#include <protocols/sewing/requirements/LigandClashRequirementCreator.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/sewing/requirements/ClashRequirement.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/conformation/Atom.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.sewing.requirements.LigandClashRequirement" );


namespace protocols {
namespace sewing {
namespace requirements {

LigandClashRequirement::LigandClashRequirement():
	LigandAssemblyRequirement(),
	maximum_clashes_allowed_( 0 ),
	clash_radius_( 4.0 ),
	active_resnum_( 0 )
{
	atom_types_ = core::chemical::ChemicalManager::get_instance()->atom_type_set( "fa_standard" );
}

LigandClashRequirement::~LigandClashRequirement()=default;

LigandClashRequirement::LigandClashRequirement( LigandClashRequirement const & src):
	LigandAssemblyRequirement( src ),
	maximum_clashes_allowed_( src.maximum_clashes_allowed_ ),
	clash_radius_( src.clash_radius_ ),
	atom_types_( src.atom_types_ )
{
}



LigandClashRequirementOP
LigandClashRequirement::clone() const {
	return LigandClashRequirementOP( new LigandClashRequirement( *this ) );
}

std::pair<bool,bool>
LigandClashRequirement::test(data_storage::SmartAssemblyOP assembly) {
	TR << "Beginning ligand clash check" << std::endl;
	current_clashes_ = 0;

	//NEW: Check for ligand clashes with partner IFF this ligand just underwent a conformer switch



	active_segment_ = assembly->get_n_terminal_segment();
	while ( active_segment_ != nullptr ) {
		for ( active_resnum_ = 1; active_resnum_ <= active_segment_->get_length(); ++active_resnum_ ) {
			active_residue_ = active_segment_->get_residue( active_resnum_ );
			if ( !active_residue_->get_atom_vector().size() ) {
				TR << "Residue" << active_resnum_ << " in segment " << active_segment_->get_segment_id() << " reports no atoms!" << std::endl;
			}
			for ( std::pair< const core::Size, data_storage::LigandResidueOP > active_ligand : assembly->get_local_ligands() ) {
				active_ligand_residue_ = active_ligand.second;
				for ( core::conformation::Atom active_ligand_atom : active_ligand_residue_->get_atom_vector() ) {
					if ( ( *atom_types_.lock() )[ active_ligand_atom.type() ].element() == "H" ) {
						continue;
					}
					if ( active_residue_->get_atom( 2 ).xyz().distance( active_ligand_atom.xyz()) <= clash_radius_ ) {
						//is the clash actually a contact?
						bool is_clash_ligand_contact = false;
						for ( data_storage::LigandContactOP active_ligand_contact_ : active_ligand_residue_->get_current_contacts() ) {
							if ( active_ligand_contact_->segment_id == active_segment_->get_segment_id()
									&& active_ligand_contact_->residue_number == active_resnum_ ) {
								is_clash_ligand_contact = true; // if yes, don't count it as a clash
								break;
							}
						}
						if ( !is_clash_ligand_contact ) {
							current_clashes_++;
							if ( current_clashes_ > maximum_clashes_allowed_ ) {
								TR << "LigandClashRequirement failed!" << std::endl;
								test_results_.first = false;
								test_results_.second = false;
								return test_results_;
							}
						}
					}
				}
			}
		}
		active_segment_ = active_segment_->get_c_terminal_neighbor();
	}
	if ( current_clashes_ > maximum_clashes_allowed_ ) { //this shouldn't happen
		TR << "LigandClashRequirement failed!" << std::endl;
		test_results_.first = false;
		test_results_.second = false;
	} else {
		TR << "LigandClashRequirement passed!" << std::endl;
		test_results_.first = true;
		test_results_.second = true;
	}
	return test_results_;
}

//Getters and Setters
core::Size
LigandClashRequirement::get_maximum_clashes_allowed() const{
	return maximum_clashes_allowed_;
}

core::Real
LigandClashRequirement::get_clash_radius() const{
	return clash_radius_;
}

void
LigandClashRequirement::set_maximum_clashes_allowed( core::Size setting ){
	maximum_clashes_allowed_ = setting;
}

void
LigandClashRequirement::set_clash_radius( core::Real setting ){
	clash_radius_ = setting;
}


void
LigandClashRequirement::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	/*
	AttributeList attributes;
	attributes
	+ XMLSchemaAttribute::attribute_w_default( "maximum_clashes_allowed", xsct_non_negative_integer, "Maximum number of clashes to allow in the assembly", "0" )
	+ XMLSchemaAttribute::attribute_w_default( "clash_radius", xs_decimal, "Radius in Angstroms within which two residues are considered to be clashing", "5.0" );
	*/
	AttributeList attributes = ClashRequirement::get_xml_attributes();
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyRequirementFactory::assembly_requirement_ct_namer )
		.element_name( LigandClashRequirement::type_name() )
		.add_attributes( attributes )
		.description( "Checks for clashes between the assembly and its ligands" )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}

void
LigandClashRequirement::set_options_from_tag(
	utility::tag::TagCOP requirement_tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Setting up LigandClashRequirement" << std::endl;
	maximum_clashes_allowed_ = requirement_tag->getOption< core::Size >( "maximum_clashes_allowed", 0 );
	TR << "Number of ligand clashes allowed: " << maximum_clashes_allowed_ << std::endl;
	clash_radius_ = requirement_tag->getOption< core::Real >( "clash_radius", 4.0 );
	TR << "Clash radius (Angstroms): " << clash_radius_ << std::endl;
}

////// Creator Methods /////

void
LigandClashRequirementCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	LigandClashRequirement::provide_xml_schema( xsd );
}



AssemblyRequirementOP
LigandClashRequirementCreator::create_requirement() const{
	return AssemblyRequirementOP( new LigandClashRequirement() );
}

std::string
LigandClashRequirement::type_name(){
	return "LigandClashRequirement";
}

std::string
LigandClashRequirementCreator::keyname() const{
	return LigandClashRequirement::type_name();
}




} //protocols
} //sewing
} //requirements






