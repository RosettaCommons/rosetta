// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   residue_io.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_residue_io_hh
#define INCLUDED_core_chemical_residue_io_hh

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <utility/io/izstream.fwd.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace chemical {

/// @brief  If polymer, determine a list of main chain atoms by shortest path from LOWER to UPPER.
AtomIndices define_mainchain_atoms( ResidueTypeOP rsd );

/// @brief virtual constructor for ResidueType objects
ResidueTypeOP
read_topology_file(
	std::string const & filename,
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
	// chemical::CSDAtomTypeSetCAP csd_atom_types kwk commenting out csd_atom_types until I have a chance to fully implement them.
	chemical::ResidueTypeSetCAP rsd_type_set
);

ResidueTypeOP
read_topology_file(
	utility::io::izstream & istream,
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
	// chemical::CSDAtomTypeSetCAP csd_atom_types kwk commenting out csd_atom_types until I have a chance to fully implement them.
	chemical::ResidueTypeSetCAP rsd_type_set
);

/// @brief writes a .params file from a given ResidueType object
void
write_topology_file(
	ResidueType const & rsd,
	std::string filename = ""
);

void
setup_atom_type_reassignments_from_commandline(
	std::string const & rsd_type_name,
	std::string const & rsd_type_set_name,
	std::map< std::string, std::string > & atom_type_reassignments
);


void
setup_atomic_charge_reassignments_from_commandline(
	std::string const & rsd_type_name,
	std::string const & rsd_type_set_name,
	std::map< std::string, Real > & atomic_charge_reassignments
);


void
setup_icoor_reassignments_from_commandline(
	std::string const & rsd_type_name,
	std::string const & rsd_type_set_name,
	std::map< std::string, utility::vector1< std::string > > & icoor_reassignments
);


} // chemical
} // core

#endif
