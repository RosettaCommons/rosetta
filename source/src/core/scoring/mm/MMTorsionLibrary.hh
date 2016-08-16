// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/MMTorsionLibrary.hh
/// @brief  Molecular mechanics torsion library class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_scoring_mm_MMTorsionLibrary_hh
#define INCLUDED_core_scoring_mm_MMTorsionLibrary_hh

// Unit headers
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.fwd.hh>

// Utility header
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.fwd.hh>

#ifdef PYROSETTA
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.hh>
#endif

namespace core {
namespace scoring {
namespace mm {

// all ints for now
typedef utility::keys::Key4Tuple< int, int, int, int > mm_torsion_atom_quad;
typedef utility::keys::Key3Tuple< double, int, double > mm_torsion_param_set;
typedef std::multimap< mm_torsion_atom_quad, mm_torsion_param_set > mm_torsion_library;
typedef std::multimap< mm_torsion_atom_quad, mm_torsion_param_set >::const_iterator mm_torsion_library_citer;
typedef std::pair< mm_torsion_library_citer, mm_torsion_library_citer > mm_torsion_library_citer_pair;

/// @brief A class to maintain a set of MM torsion paramaters.
///
/// @details The MMTorsionLibrary class contains functions and data structures to initialize, store, and provide lookups
/// for molecular mechanics torsional paramaters. The sets of 4 atom types that make up a dihedral angle are stored in an
/// mm_torsion_atom_quad . The sets 3 parameters that act as values for the actual MM torsion potential are stored in an
/// mm_torsion_param_set . The class maintains 2 libraries; one for set of torsion parameter for which all 4 atom types in
/// the dihedral at defined ( fully_assigned_mm_torsion_library_ ), and a library that is only dependant on the two central
/// atom types ( wildcard_mm_torsion_library_ ). Wildcar atom types are designated by 'X'. Lookups are done with MMAtomType
/// indicies or strings that are converted to MMAtomType indicies. Lookup functions return pairs of iterators to the map.
/// The same set of 4 atom types can corespond to multiple paramter sets with different multiplicities to more accuratly
/// define the potential and is the reason multimaps are used instead of just maps.
///
class MMTorsionLibrary  : public utility::pointer::ReferenceCount
{

public:
	/// @brief ctor
	MMTorsionLibrary( std::string filename, core::chemical::MMAtomTypeSetCOP mm_atom_set );
	virtual ~MMTorsionLibrary();

	/// @brief lookup by atom type ints
	mm_torsion_library_citer_pair
	lookup( int atom1, int atom2, int atom3, int atom4 ) const;

	/// @brief lookup by atom type strings
	mm_torsion_library_citer_pair
	lookup( std::string atom1, std::string atom2,
		std::string atom3, std::string atom4 ) const;

private:

	/// @brief library that contains torsion params for sets in which all mm atom types
	/// are specified for all 4 positions
	mm_torsion_library fully_assigned_mm_torsion_library_;

	/// @brief library that contains torsion params for sets in which all mm atom types
	/// are NOT specified for all 4 positions
	mm_torsion_library wildcard_mm_torsion_library_;

	/// @brief the MMAtomTypeSet associated with the library
	core::chemical::MMAtomTypeSetCAP mm_atom_set_;

};


} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_torsion_library_HH
