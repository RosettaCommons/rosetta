// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/LigandArea.hh>

// Package Headers
#include <protocols/ligand_docking/LigandDockingLoaders.hh>

// Project Headers
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static basic::Tracer ligand_area_tracer( "protocols.ligand_docking.ligand_options.LigandArea", basic::t_debug );

LigandArea::LigandArea():
	ReferenceCount(),
	chain_('?'),
	cutoff_(0),
	Calpha_restraints_(0),
	minimize_ligand_(0),
	tether_ligand_(0),
	high_res_angstroms_(0.1),
	high_res_degrees_(2.8648),// the old code had 0.05 radians
	add_nbr_radius_(false),
	all_atom_mode_(false)
{}

LigandArea::~LigandArea() = default;

void LigandArea::parse_my_tag(
	utility::tag::TagCOP tag
){
	if ( ! tag->hasOption("chain") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'LigandArea' requires 'chain' tag");
	chain_= tag->getOption<char>("chain");

	if ( ! tag->hasOption("cutoff") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'LigandArea' requires 'cutoff' tag");
	cutoff_= tag->getOption<core::Real>("cutoff");

	if ( tag->hasOption("Calpha_restraints") ) {
		Calpha_restraints_= tag->getOption<core::Real>("Calpha_restraints");
	}

	if ( tag->hasOption("minimize_ligand") ) {
		minimize_ligand_= tag->getOption<core::Real>("minimize_ligand");
	}

	if ( tag->hasOption("tether_ligand") ) {
		tether_ligand_= tag->getOption<core::Real>("tether_ligand");
	}

	if ( tag->getOption< bool >("add_nbr_radius") ) {
		add_nbr_radius_= true;
	} else {
		add_nbr_radius_= false;
	}

	if ( tag->getOption< bool >("all_atom_mode") ) {
		all_atom_mode_ = true;
	} else {
		all_atom_mode_ = false;
	}

	if ( tag->hasOption("high_res_angstroms") ) {
		high_res_angstroms_= tag->getOption<float>("high_res_angstroms");
	}

	if ( tag->hasOption("high_res_degrees") ) {
		high_res_degrees_= tag->getOption<float>("high_res_degrees");
	}

}

std::string LigandArea::element_name() { return "LigandArea"; }
void LigandArea::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using Attr = XMLSchemaAttribute;
	AttributeList attributes;
	attributes + optional_name_attribute();

	attributes + Attr::required_attribute( "chain", xsct_char, "The ligand's chain" )
		+ Attr::required_attribute( "cutoff", xsct_real, "The distance cutoff from the ligand" )
		+ Attr( "Calpha_restraints", xsct_real, "Size of one standard deviation (in Angstroms) for restraints on C-alphas" )
		+ Attr( "minimize_ligand", xsct_real, "Size of one standard deviation (in degrees) for ligand torsion"
		" angles used by the ResidueTorsionRestraints class to create CircularHarmonic restraints on the ligand"
		" dihedrals during minimization to keep these dihedrals near their starting conformation" )
		+ Attr( "tether_ligand", xsct_real, "The standard deviation (in Angstroms) used by the HighResDocker to create"
		" coordinate constraints on the ligand's neighbor atom" )
		+ Attr::required_attribute( "add_nbr_radius", xsct_rosetta_bool, "Used by an InterfaceBuilder for deciding how to define"
		" the distance threshold between the ligand and the protein's residues. If this is 'true', then the neighbor radius of"
		" the protein residue is added into the LigandArea's 'cutoff' parameter. If this is 'false', then the 'cutoff' parameter"
		"is used unaltered" )
		+ Attr::required_attribute( "all_atom_mode", xsct_rosetta_bool, "When deciding whether a protein residue is sufficiently"
		" close to the ligand, should the distance between the protein residue's neighbor atom and every atom in the ligand be"
		" measured? If not, then only the ligand's neighbor atom will be used in that decision." )
		+ Attr( "high_res_angstroms", xsct_real, "The euclidean perturbation magnitude, in Angstroms, used by the HighResDocker" )
		+ Attr( "high_res_degrees", xsct_real, "The rotational perturbation magnitude, in degrees, used by the HighResDocker" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( element_name() )
		.complex_type_naming_func( & LigandAreaLoader::ligand_area_ct_namer )
		.description( "LigandAreas are used to define InterfaceBuilders which in turn are used to define MoveMapBuilders." )
		.add_attributes(attributes)
		.write_complex_type_to_schema( xsd );
}


} //namespace ligand_docking
} //namespace protocols
