// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)


// Project Headers
#include <core/types.hh>
#include <basic/Tracer.hh>

// Unit Headers
#include <protocols/ligand_docking/LigandArea.hh>

// Utility headers
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer ligand_area_tracer( "protocols.ligand_docking.ligand_options.LigandArea", basic::t_debug );

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

LigandArea::~LigandArea() {}

void LigandArea::parse_my_tag(
	utility::tag::TagCOP tag
){
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'LigandArea' requires 'chain' tag");
	chain_= tag->getOption<char>("chain");

	if ( ! tag->hasOption("cutoff") ) throw utility::excn::EXCN_RosettaScriptsOption("'LigandArea' requires 'cutoff' tag");
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

	if ( tag->getOption<std::string>("add_nbr_radius") == "true" ) {
		add_nbr_radius_= true;
	} else if ( tag->getOption<std::string>("add_nbr_radius") != "false" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'add_nbr_radius' option is true or false");
	}

	if ( tag->getOption<std::string>("all_atom_mode") == "true" ) {
		all_atom_mode_= true;
	} else if ( tag->getOption<std::string>("all_atom_mode") != "false" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'all_atom_mode' option is true or false");
	}

	if ( tag->hasOption("high_res_angstroms") ) {
		high_res_angstroms_= tag->getOption<float>("high_res_angstroms");
	}

	if ( tag->hasOption("high_res_degrees") ) {
		high_res_degrees_= tag->getOption<float>("high_res_degrees");
	}

}

} //namespace ligand_docking
} //namespace protocols
