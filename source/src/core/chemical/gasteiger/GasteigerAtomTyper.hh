// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a Rosetta Commons Member Institution

/// @file   core/chemical/gasteiger/GasteigerAtomTyper.hh
/// @brief  The type assigner for gasteiger type data.
/// @author Rosetta conversion: Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_gasteiger_GasteigerAtomStandardizer_hh
#define INCLUDED_core_chemical_gasteiger_GasteigerAtomStandardizer_hh

#include <core/chemical/gasteiger/GasteigerAtomTypeSet.fwd.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh> // For enums
#include <core/chemical/Element.fwd.hh>

#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/ResidueType.fwd.hh>

#include <core/types.hh>

#include <utility/vector0.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <list>

namespace core {
namespace chemical {
namespace gasteiger {

class PossibleAtomTypesForAtom; // Forward declaration

void assign_gasteiger_atom_types( core::chemical::ResidueType & restype, GasteigerAtomTypeSetCOP gasteiger_atom_type_set, bool keep_existing, bool allow_unknown = false);

/// @breif Version which gets typeset from ResidueType, or just uses default
void assign_gasteiger_atom_types( core::chemical::ResidueType & restype, bool keep_existing = true, bool allow_unknown = false);

PossibleAtomTypesForAtom GetPossibleTypesForAtom(
	core::chemical::RealResidueGraph const & graph,
	RealResidueVD const & atomVD,
	GasteigerAtomTypeSetCOP gasteiger_atom_type_set,
	core::Size connections );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @class PossibleAtomTypesForAtom
//! @brief A helper class by which AtomTypes can return all possible atom types for a given atom in a structure
//!        that is easily accessed by orbital type
//!
//! @see @link example_chemistry_possible_atom_types_for_atom.cpp @endlink
//! @author mendenjl
//! @date Aug 27, 2010
//!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PossibleAtomTypesForAtom :
	public utility::pointer::ReferenceCount
{
	//private:
	//
	////////////
	//// data //
	////////////

	//vector0 because it's indexed with an enum which starts at zero.
	//!< Number of atom types with each hybridization
	GasteigerAtomTypeSetCOP gasteiger_atom_type_set_;
	std::map< std::string, PossibleAtomTypesForAtom> s_atomic_env_outside_arom_ring_to_types_map;
	std::map< std::string, PossibleAtomTypesForAtom> s_element_bonds_in_arom_ring_to_types_map;
	utility::vector0< core::Size >   m_NumberAtomTypesWithHybridization;
	std::list< GasteigerAtomTypeDataCOP >   m_AtomTypesByDecreasingStability;   //!< Most stable types first
	core::Size m_NumberConjugatedTypes;            //!< Number of conjugated types in the list

	//! Function used to resolve the final atom type; only used in a few ambiguous cases
	//! such as trigonal vs. tetrahedral nitrogen with 3 bonds, otherwise NULL
	//! Parameters are
	//! 1. the atom
	//! 2. the size of the smallest ring that this atom is part of
	void ( PossibleAtomTypesForAtom::*m_FinalizeFunction)( core::chemical::RealResidueGraph const &, core::chemical::RealResidueVD const & );

	//! @brief Create the map from atom environment string to possible atom types
	//! @param IN_AROMATIC_RING whether to only include types that could be in an aromatic ring
	static std::map< std::string, PossibleAtomTypesForAtom> CreateAtomicEnvironmentToTypesMap
	(
		const bool IN_AROMATIC_RING, GasteigerAtomTypeSetCOP GASTEIGER_ATOM_TYPE_SET
	);

public:

	//  //! @brief write out the atom typing scheme
	//  //! @param OSTREAM stream to write the atom typing scheme to
	//  static std::ostream &WriteDetailedScheme( std::ostream &OSTREAM);
	//
	////////////
	//// data //
	////////////
	//
	//  //! single instance of that class
	//  static const util::SiPtr< const util::ObjectInterface> s_Instance;
	//
	////////////////////////////////////
	//// construction and destruction //
	////////////////////////////////////

	//! @brief default constructor
	PossibleAtomTypesForAtom();

	//! @brief constructor from the known information about the atom
	//! @param ELEMENT element type,
	//! @param NUMBER_ELECTRONS_IN_BONDS number of electrons in bonds for the atom type
	//! @param NUMBER_BONDS number of bonds for the atom
	//! @param SUSPECTED_CHARGE; expected charge, ignored if no atom type matching the other criteria if found
	//! @param IN_AROMATIC_RING true iff the atom has bonds of the aromatic unspecified type
	static PossibleAtomTypesForAtom
	FindPossibleAtomTypesForAtom
	(
		GasteigerAtomTypeSetCOP GASTEIGER_ATOM_TYPE_SET,
		const Element &ELEMENT,
		const core::Size NUMBER_ELECTRONS_IN_BONDS,
		const core::Size NUMBER_BONDS,
		const int SUSPECTED_CHARGE,
		const bool IN_AROMATIC_RING
	);

	//  //! @brief Clone function
	//  //! @return pointer to new PossibleAtomTypesForAtom
	//  PossibleAtomTypesForAtom *Clone() const;
	//
	///////////////////
	//// data access //
	///////////////////
	//
	//  //! @brief returns class name of the object behind a pointer or the current object
	//  //! @return the class name
	//  const std::string &GetClassIdentifier() const;

	void gasteiger_atom_type_set( GasteigerAtomTypeSetCOP GASTEIGER_ATOM_TYPE_SET );

	//! @brief tell whether a particular hybrid orbital type is possible given what we know about this atom
	//! @param HYBRID the type of hybrid orbital
	//! @return true iff there is a possible atom type for that hybrid orbital
	bool CouldHaveHybridization( const GasteigerAtomTypeData::HybridOrbitalType HYBRID) const;

	//! @brief return the number of types that the atom has the potential to become
	//! @return the number of types that the atom has the potential to become
	core::Size GetNumberPossibleTypes() const;

	//! @brief return the most stable type
	//! @return the most stable type - NULL if no such type exists
	GasteigerAtomTypeDataCOP GetMostStableType() const;

	//  //! @brief get the alternate atom types
	//  //! @return the alternative atom types
	//  storage::Vector< AtomType> GetAlternateTypes() const;
	//
	//  //! @brief get the alternate atom type with the given charge
	//  //! @param CHARGE the charge desired
	//  //! @return an alternative atom type
	//  AtomType GetAlternateTypeWithCharge( const short &CHARGE) const;
	//
	//  //! @brief returns true if any of the candidate types are conjugated
	//  //! @return true if any of the candidate types are conjugated
	//  bool CouldBeConjugated() const;
	//
	//  //! @brief returns true if all of the candidate types are conjugated
	//  //! @return true if all of the candidate types are conjugated
	//  bool MustBeConjugated() const;
	//
	//  //! @brief determine the maximal # of pi-electrons in the pi-electron system
	//  //! @return the maximal # of pi-electrons in the pi-electron system
	//  std::size_t GetMaxElectronsParticipatingInPiSystem() const;
	//
	//////////////////
	//// operations //
	//////////////////

	//! @brief add an atom type to be considered
	//! @param ATOM_TYPE the type of atom to consider
	void AddAtomType( GasteigerAtomTypeDataCOP ATOM_TYPE);

	//! @brief set this object to only consider the given atom type
	//! @param ATOM_TYPE the atom type desired
	void SetToType( GasteigerAtomTypeDataCOP ATOM_TYPE);

	//! @brief set the final type based on the given atom and smallest ring size
	void Finalize( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD);

	/////////////////
	//// operators //
	/////////////////
	//
	////////////////////////
	//// input and output //
	////////////////////////
	//
	//protected:
	//
	//  //! @brief read from std::istream
	//  //! @param ISTREAM input stream
	//  //! @return istream which was read from
	//  std::istream &Read( std::istream &ISTREAM);
	//
	//  //! @brief write to std::ostream
	//  //! @param OSTREAM outputstream to write to
	//  //! @param INDENT number of indentations
	//  //! @return outputstream which was written to
	//  std::ostream &Write( std::ostream &OSTREAM, const std::size_t INDENT) const;
	//
	////////////////////////
	//// helper functions //
	////////////////////////

	core::Size hybridization_rank( GasteigerAtomTypeData::HybridOrbitalType const hybrid );

	//! @brief remove a particular hybrid orbital type from the possible types, unless that would remove all possibilities
	//! @param HYBRID the type of hybrid orbital to remove
	void RemoveHybridization( const GasteigerAtomTypeData::HybridOrbitalType HYBRID);

	//! @brief add an atom type to the search using a set of rules for atom types in aromatic rings
	//! @param ATOM_TYPE the type of atom to consider
	//! @param DESIRED_CHARGE the charge desired
	//! The atom type will be ordered using the distance from the desired charge as the first sort key, second by
	//! the stability.  Unlike AddAtomType, AddAromaticAtomType always adds the type to the considered list
	void AddAromaticAtomType( GasteigerAtomTypeDataCOP ATOM_TYPE, const int DESIRED_CHARGE);

	//! @brief Select the best choice for the atom type wherever possible
	//! @see @link https://structbio.vanderbilt.edu:8443/display/MeilerLab/RethinkingAtomTypeDetection @endlink
	//! @details the link above contains the statistics and models used to select the current set of rules
	void Finalize();

	//! @brief choose the preferred atom type (using VSEPR theory) assuming that the orbitals do not hybridize
	//! @details This is used for elements in group 1, 2, & 7, which do hybridize in the gasteiger scheme
	void FinalizeUnhybridized();

	//! @brief only keep the most stable types for the atom that span the set of desired pi orbital electrons (0-2)
	//! @param DESIRED_CHARGE the preferred charge
	//! used during construction of the maps when there is no part of standardization that
	//! should edit this class
	void FinalizeAromatic( const int DESIRED_CHARGE);

	//! @brief choose the final atom type for Nitrogen with two single bonds
	void FinalizeNitrogenTwoSingle( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD );

	//! @brief choose the final atom type for a nitrogen with a single and a double bond
	void FinalizeNitrogenSingleDouble( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD );

	//! @brief choose the final atom type for Nitrogen with three single bonds
	void FinalizeNitrogenThreeSingle( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD );

	//! @brief choose the final atom type for Oxygen with two single bonds
	void FinalizeOxygenTwoSingle( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD );

	//! @brief choose the final atom type for Oxygen with a single and a double bond
	void FinalizeOxygenSingleDouble( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD );

	//! @brief choose the final atom type for Oxygen with three single bonds
	void FinalizeOxygenThreeSingle( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD );

	//  //! @brief get connected element types
	//  //! @param ATOM the atom of interest
	//  //! @return a set of the connected element types
	//  static storage::Set< ElementType> GetConnectedElementTypes( const AtomConformationalInterface &ATOM);

	//! @brief test whether a particular atom is unsaturated, without relying on atom types having already been set
	//! @return true if atom has no A. unsaturated bonds or B. is part of an aromatic ring or C. has empty orbitals
	bool IsUnsaturated( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) const;

	//! @brief count unsaturated neighbors
	//! @return the number of unsaturated neighbors around ATOM
	core::Size CountUnsaturatedNeighbors( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) const;

	//! @brief test whether atom is bonded to any halogens
	//! @param ATOM the atom of interest
	//! @return true if the atom is bonded to any halogens
	bool IsBondedToAHalogen( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) const;

}; // class PossibleAtomTypesForAtom

} // gasteiger
} // chemical
} // core

#endif

