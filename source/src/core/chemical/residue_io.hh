// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/io/izstream.fwd.hh>

// C++ headers
#include <iosfwd>
#include <map>

namespace core {
namespace chemical {

/// @brief  If polymer, determine a list of main chain atoms by shortest path from LOWER to UPPER.
AtomIndices define_mainchain_atoms( ResidueTypeOP rsd );

/// @brief function to convert params files into ResidueType objects (repackages string filename into istream, gets needed subsidiary type sets from rsd_type_set
ResidueTypeOP
read_topology_file(
	std::string const & filename,
	chemical::ResidueTypeSetCOP rsd_type_set
);

/// @brief function to convert params files into ResidueType objects, gets needed subsidiary type sets from rsd_type_set
ResidueTypeOP
read_topology_file(
	utility::io::izstream & istream,
	chemical::ResidueTypeSetCOP rsd_type_set
);

/// @brief function to convert params files into ResidueType objects, gets needed subsidiary type sets from rsd_type_set
ResidueTypeOP
read_topology_file(
	std::istream & istream,
	std::string const & filename, //this may be a fake filename
	chemical::ResidueTypeSetCOP rsd_type_set
);

/// @brief function to convert params files into ResidueType objects (repackages string filename into istream)
ResidueTypeOP
read_topology_file(
	std::string const & filename,
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types
);

/// @brief main function to convert params files into ResidueType objects
ResidueTypeOP
read_topology_file(
	std::istream & data,
	std::string const & filename, //MAY be faux filename if stream is not izstream of file
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types
);

/// @brief writes a .params file from a given ResidueType object
void
write_topology_file(
	ResidueType const & rsd,
	std::string filename = ""
);

/// @brief Produces a graphviz dot representation of the ResidueType to the given output stream
/// If header is true (the default) a line with an explanatory message will be printed first.
void
write_graphviz(
	ResidueType const & rsd,
	std::ostream & out,
	bool header = true
);

/// @brief Certain commandline flags override the default RamaPrePro maps used by the 20
/// canonical amino acids.  This function applies those overrides.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
set_up_mapfile_reassignments_from_commandline(
	ResidueTypeOP rsd
);

void
setup_atom_type_reassignments_from_commandline(
	std::string const & rsd_type_name,
	TypeSetMode rsd_type_set_mode,
	std::map< std::string, std::string > & atom_type_reassignments
);


void
setup_atomic_charge_reassignments_from_commandline(
	std::string const & rsd_type_name,
	TypeSetMode rsd_type_set_mode,
	std::map< std::string, Real > & atomic_charge_reassignments
);


void
setup_icoor_reassignments_from_commandline(
	std::string const & rsd_type_name,
	TypeSetMode rsd_type_set_mode,
	std::map< std::string, utility::vector1< std::string > > & icoor_reassignments
);

/// @brief Symmetrize the glycine params file (if the user has used the -symmetric_gly_tables option).
/// @details Ugh.  Special-case logic.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
void
apply_symm_gly_corrections(
	std::string const &child_atom,
	core::Real &phi,
	core::Real &theta,
	core::Real &d,
	std::string &parent_atom,
	std::string &angle_atom,
	std::string &torsion_atom
);

} // chemical
} // core

#endif
