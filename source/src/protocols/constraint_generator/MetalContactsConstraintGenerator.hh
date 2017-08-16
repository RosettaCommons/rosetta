// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/MetalContactsConstraintGenerator.hh
/// @brief This constraint generator takes residue selectors for a residue containing metal ion(s) and for residue(s) for which to set up contacts. It allows users to specify which base atoms will be used to define angles/dihedrals to constrain; ideal values for angles/dihedrals/distances; and an option to constrain to native values.
/// @author guffysl (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_constraint_generator_MetalContactsConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_MetalContactsConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/MetalContactsConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <set>
namespace protocols {
namespace constraint_generator {

///@brief This constraint generator takes residue selectors for a residue containing metal ion(s) and for residue(s) for which to set up contacts. It allows users to specify which base atoms will be used to define angles/dihedrals to constrain; ideal values for angles/dihedrals/distances; and an option to constrain to native values.
class MetalContactsConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	MetalContactsConstraintGenerator();

	MetalContactsConstraintGenerator( MetalContactsConstraintGenerator const & src );

	virtual ~MetalContactsConstraintGenerator();

	static std::string
	class_name();

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const;

	virtual core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


public:

	//GETTERS
	std::string
	get_ligand_atom_name() const;

	bool
	get_use_ligand_selector() const;

	core::select::residue_selector::ResidueSelectorCOP
	get_ligand_selector() const;

	std::string
	get_ligand_resnum_string() const;

	bool
	get_use_contact_selector() const;

	core::select::residue_selector::ResidueSelectorCOP
	get_contact_selector() const;

	std::string
	get_contact_resnum_string() const;

	std::string
	get_base_atom_name() const;

	std::string
	get_base_base_atom_name() const;

	core::Real
	get_ideal_distance() const;

	utility::vector1< core::Real >
	get_ideal_angle_about_contact() const;

	utility::vector1< core::Real >
	get_ideal_dihedral_about_contact() const;

	utility::vector1< core::Real >
	get_ideal_angle_about_metal() const;

	utility::vector1< core::Real >
	get_ideal_dihedral_about_metal() const;

	utility::vector1< core::Real >
	get_ideal_dihedral_3() const;

	bool
	get_score_against_internal_contacts() const;

	core::Real
	get_dist_cutoff_multiplier() const;

	bool
	get_constrain_to_closest() const;

	//SETTERS

	void
	set_ligand_atom_name( std::string);

	void
	set_use_ligand_selector( bool );

	void
	set_ligand_selector( core::select::residue_selector::ResidueSelectorCOP );

	void
	set_ligand_resnum_string( std::string);

	void
	set_use_contact_selector( bool );

	void
	set_contact_selector(core::select::residue_selector::ResidueSelectorCOP );

	void
	set_contact_resnum_string( std::string );

	void
	set_base_atom_name(std::string);

	void
	set_base_base_atom_name(  std::string );

	void
	set_ideal_distance( core::Real );

	void
	set_ideal_angle_about_contact( utility::vector1< core::Real > );

	void
	set_ideal_dihedral_about_contact( utility::vector1< core::Real >);

	void
	set_ideal_angle_about_metal( utility::vector1< core::Real > );

	void
	set_ideal_dihedral_about_metal( utility::vector1< core::Real >);

	void
	set_ideal_dihedral_3( utility::vector1< core::Real >);

	void
	set_score_against_internal_contacts(bool);

	void
	set_dist_cutoff_multiplier( core::Real );

	void
	set_constrain_to_closest( bool );
private:
	//Private methods

	//Function to get ligand resnum
	///@brief Uses private data to compute ligand resnums based on selector/resnum string and pose
	core::Size
	get_ligand_resnum( core::pose::Pose const & pose ) const;

	//Function to get selector resnums
	///@brief Uses private data to compute contact residue numbers based on selector/resnum string and pose; inserts them into the provided set
	void
	get_contact_resnums( core::pose::Pose const & pose, std::set< core::Size > & output ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

private:
	core::Real dist_cutoff_multiplier_=1.0;
	bool use_ligand_selector_=false;
	std::string ligand_atom_name_=""; //REQUIRED
	core::select::residue_selector::ResidueSelectorCOP ligand_selector_;
	std::string ligand_resnum_string_="";
	bool use_contact_selector_=false;
	core::select::residue_selector::ResidueSelectorCOP contact_selector_;
	std::string contact_resnum_string_="";
	std::string base_atom_name_="";
	std::string base_base_atom_name_="";
	core::Real ideal_distance_=-1.0;//

	//THE FOLLOWING COULD HAVE MULTIPLE OK VALUES depending on contacts:
	utility::vector1< core::Real > ideal_angle_about_contact_;
	utility::vector1< core::Real > ideal_dihedral_about_contact_;
	utility::vector1< core::Real > ideal_angle_about_metal_;
	utility::vector1< core::Real > ideal_dihedral_about_metal_;
	utility::vector1< core::Real > ideal_dihedral_3_;

	//Should we also score dihedrals/angles about metal against internal contacts?
	bool score_against_internal_contacts_=false;
	bool constrain_to_closest_=true;


};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_MetalContactsConstraintGenerator_hh
