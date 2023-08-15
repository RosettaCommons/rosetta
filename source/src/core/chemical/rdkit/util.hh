// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/util.hh
/// @brief Utilities for interacting with the RDKit library.
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_rdkit_util_hh
#define INCLUDED_core_chemical_rdkit_util_hh

#include <core/chemical/Bond.fwd.hh>
#include <core/chemical/rdkit/RDKit.fwd.hh>

#ifdef PYROSETTA
#include <rdkit/ForceField/ForceField.h>
#endif

#include <core/chemical/AtomRefMapping.hh>
#include <core/types.hh>

#include <map>

#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h> // For MatchVectType

namespace core {
namespace chemical {
namespace rdkit {

/// @brief Initialize the RDKit random number generator.
/// @details Note that seed is an int to match the seed generated in core/init.cc
void initialize_rdkit_random( int seed );

/// @brief Initialize the RDKit output levels with the Rosetta commandline settings
/// @details You can set the global RDKit output by controlling the "RDKit" tracer.
void initialize_rdkit_tracers();

/// @brief Convert a Rosetta BondName enum to an RDKit BondType value
::RDKit::Bond::BondType
convert_to_rdkit_bondtype( core::chemical::BondName bondtype, bool aro2double = false);

/// @brief Convert an RDKit BondType value to a Rosetta BondName enum
core::chemical::BondName
convert_from_rdkit_bondtype( ::RDKit::Bond::BondType bondtype);

/// @brief Get the name of the RDMol
std::string
get_name(::RDKit::ROMol const & mol);

/// @brief Does the molecule have physical hydrogens?
bool
has_physical_Hs(::RDKit::ROMol const & mol);

/// @brief Does the molecule have "explicit" (but not physical) hydrogens?
bool
has_explicit_Hs(::RDKit::ROMol const & mol);

/// @brief Does the molecule have implicit hydrogens?
bool
has_implicit_Hs(::RDKit::ROMol const & mol);

/// @brief Load an RDKit forcefield. Will prefer MMFF, but will fall back to UFF
::RDKit::ForceFieldOP
get_forcefield(::RDKit::ROMol & mol, int conf_num =-1);

/// @brief Non strict sanitization, useful if working with molecules which aren't 100% acceptable by RDKit (e.g. protonation/kekulization issues)
void
softSanitize(::RDKit::RWMol & mol);

/// @brief Remove any excess hydrogens, where "excess" is defined as any which contribute to a positive formal charge
/// (Assumes a graph-hydrogen removed form.)
void
remove_excess_protons(::RDKit::RWMol & rdmol);

typedef std::vector< std::pair< ::RDKit::RWMolOP, int > > ChargeTransformList;

/// @brief Reset the charges on particular groups
/// (Assumes a graph-hydrogen removed form.)
void
apply_charge_transforms( ::RDKit::RWMol & rdmol, ChargeTransformList const & transforms );

/// @brief Change the protonation state on an rdmol to match physiological pH.
/// Will reset the molecule to an graph-explicit hydrogen form
void
reprotonate_rdmol(::RDKit::RWMol & rdmol);

/// @brief Remove any pH-dependent charges on the molecule
/// If addHs is true, then the resultant molecule with have hydrogens added.
/// If addHs is false, then the resultant molecule will be without hydrogens.
void
neutralize_rdmol(::RDKit::RWMol & rdmol, bool addHs=true);

/// @brief Final neutralization process invoked in constructing an RDMol
/// from a piece of a Pose
void
final_neutralize(
	RDKit::RWMOL_SPTR const & rdmol
);

/// @brief Label a molecule with it's index values (for find_mapping, later)
void
label_with_index(  ::RDKit::ROMol & rdmol, std::string const & index_prop = "Orig_Index" );

/// @brief Convert the MatchVectType to an IndexIndex map, going query->molecule
/// Use -1 as the invalid value, as zero is a valid one.
core::chemical::IndexIndexMapping
convert_match_vect_to_index_index_map( ::RDKit::MatchVectType const & match_vect );

/// @brief Find a mapping from one RDMol to another
/// If index_prop is found, it will be used (if present) as the starting index.
core::chemical::IndexIndexMapping
find_mapping( ::RDKit::ROMOL_SPTR from, ::RDKit::ROMOL_SPTR to, std::string const & index_propr = "" );

/// @brief Load the structures from the given file into the provided vector
/// Structures will be appended to the ones already in the vector
void load_sdf( std::string const & filename, utility::vector1< ::RDKit::ROMolOP > & mol_vector, bool removeHs = true );

/// @brief Return a set containing all the valid names for the rdkit_metric() function mapped to short descriptions
std::map< std::string, std::string >
get_metric_names();

/// @brief Return the value of a given RDKit metric for the given mol
core::Real
rdkit_metric(::RDKit::ROMol const & mol, std::string const & metric );

} // namespace rdkit
} // namespace chemical
} // namespace core

#endif // Include guard.
