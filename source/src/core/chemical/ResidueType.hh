// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResidueType.hh
/// @brief Method declarations and simple accessors/getters for ResidueType
/// @author
/// Phil Bradley
/// Rocco Moretti (rmorettiase@gmail.com)
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to properties, rings, lipids, carbohydrates, and other non-AAs)


#ifndef INCLUDED_core_chemical_ResidueType_hh
#define INCLUDED_core_chemical_ResidueType_hh

// Unit headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeBase.hh>

#include <core/chemical/MutableResidueType.fwd.hh> // For conversion operation only
// IMPORTANT - Don't allow ResidueGraphTypes.hh to be included in this header directly or indirectly.
// It slows compilation immensely.

// Package headers
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomProperty.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/Bond.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/Element.fwd.hh>
#include <core/chemical/Elements.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/AtomProperties.fwd.hh>
#include <core/chemical/ResidueProperties.fwd.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#if defined(__INTEL_COMPILER) || defined(WIN32)
#include <core/chemical/RestypeDestructionEvent.hh>
#else
#include <core/chemical/RestypeDestructionEvent.fwd.hh>
#endif

#ifdef WIN32
#include <core/chemical/Orbital.hh>
#include <core/chemical/ResidueConnection.hh>
#else
#include <core/chemical/Orbital.fwd.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#endif
#include <core/chemical/rna/RNA_Info.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.fwd.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.fwd.hh>
#include <core/chemical/rings/RingConformerSet.fwd.hh>
#include <core/chemical/rings/RingSaturationType.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <map>
//#include <list>

#ifdef MULTI_THREADED
#include <mutex>
#endif

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

typedef utility::keys::Key2Tuple< Size, Size > two_atom_set;
typedef utility::keys::Key3Tuple< Size, Size, Size > three_atom_set;
typedef utility::keys::Key3Tuple< Size, Size, Size > bondangle_atom_set;
typedef utility::keys::Key4Tuple< Size, Size, Size, Size > dihedral_atom_set;

/// @brief A class for defining a type of residue
/// @details
/// This class contains the "chemical" information for residues as well as the ideal xyz and internal coordinates for a
/// residue (generated xyz coordinates are found in core/conformation/Residue.hh).  A ResidueType in Rosetta can be a
/// ligand, DNA, amino acid, or basically anything. ResidueTypes are generated through .params files, which are read
/// from the database chemical/residue_types.  For ligands, a parameter has to be provided to rosetta through the
/// -extra_res_fa flag.  Primary data are set through the residue_io.cc class.  The primary data that are set are:
/// atoms, mmatoms, orbitals, and properties of the particular ResidueType.  These properties can be modified through
/// patches, which create new ResidueTypes, and are controlled through PatchOperations.cc.
///
/// When a ResidueType is initially created, it's created as a MutableResidueType, which represents the data as an atom graph,
/// which facilitates atomistic modification. This class (the plain ResidueType class) is a "finalized" version of the
/// ResidueType, which doesn't permit modification, and is the main class used by most Rosetta simulations.
/// It's created from a MutableResidueType, and for efficiency represents the data as a "struct of vectors" representation.
///
/// Correspondingly, the main representation of atoms is different between the classes. Whereas the MutableResidueType uses
/// an opaque "vertex descriptor" to identify atoms, the plain ResidueType class uses an integer index.
/// The order of atoms is well-defined, with heavy atoms coming before hydrogens, and backbone atoms coming before sidechain atoms.
/// (and hydrogens in the same order as the heavyatoms they're connected to).
///
/// Properties: Properties of a residue include things like DNA, PROTEIN, CHARGED, etc.  These properties indicate the
/// type of residue it is and what properties are associated with the residue. They are set when read in.
/// To add new ResidueProperties, add them to core/chemical/residue_properties/general_properties.list.
///
/// Orbitals: Orbitals are indexed separately from atoms.  They function much the same way as atoms, except for some
/// key differences.  To find atoms bonded to orbitals, you must provide the atom index, not the orbital index.  (I
/// haven't figured out how to get the reverse to work because of the separate indices.)  Orbital xyz coordinates are
/// not updated when atom coordinates are.  This is to keep speed consistent with just having atoms.  To output the
/// orbitals, use the flag -output_orbitals.
class ResidueType : public ResidueTypeBase, public utility::pointer::enable_shared_from_this< ResidueType >
{
	// All constructors are private.
	// The only way to create a ResidueType is with the ResidueType::make(MutableResidueType) factory function
	// But we need a way to create a blank type for serialization code, hence the following ifdef
#ifdef SERIALIZATION
public:
#else
private:
#endif

	ResidueType(); // private, not deleted, as serialization needs it.

private:

	/// @brief Make (most of) a ResidueType from a MutableResidueType
	/// Use ResidueType::make() instead.
	ResidueType( MutableResidueType const & residue_type );

	/// @brief make a copy of this ResidueType
	ResidueType( ResidueType const & residue_type ) = delete;

	/// @brief make a copy
	ResidueTypeOP clone() const = delete;

	/// @brief make a copy
	ResidueTypeOP placeholder_clone() const = delete;

	/// @brief Copies  <src>  into the ResidueType
	ResidueType &
	operator=( ResidueType const & src ) = delete;

public:

	/// @brief Make a ResidueType from a Modifiable residue type.
	/// @details In practice, this is the only way to create (or change)
	/// a ResidueType object.
	/// @details By using a static function and not a constructor,
	/// this ensures that all ResidueTypes are on the heap (and can shared_from_this()),
	/// and allows us to do OP-requiring subobject updates in self_pointer_updates().
	static
	ResidueTypeCOP
	make( MutableResidueType const & residue_type  );

	/// @brief destructor
	~ResidueType() override;

public:

	/// self pointers
	inline ResidueTypeCOP get_self_ptr() const      { return shared_from_this(); }
	inline ResidueTypeOP  get_self_ptr()            { return shared_from_this(); }
	inline ResidueTypeCAP get_self_weak_ptr() const { return ResidueTypeCAP( shared_from_this() ); }
	inline ResidueTypeAP  get_self_weak_ptr()       { return ResidueTypeAP( shared_from_this() ); }

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Summary Functions               ////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief number of atoms
	Size
	natoms() const override
	{
		return natoms_;
	}

	/// @brief number of heavy atoms
	Size
	nheavyatoms() const override
	{
		return nheavyatoms_;
	}

	/// @brief number of hbond_acceptors
	Size
	n_hbond_acceptors() const
	{
		return n_hbond_acceptors_;
	}

	/// @brief number of hbond_donors
	Size
	n_hbond_donors() const
	{
		return n_hbond_donors_;
	}

	/// @brief number of bonds
	Size
	nbonds() const override
	{
		return nbonds_;
	}

	/// @brief Counts the number of virtual atoms and returns the count.
	/// @details The virtual count is not stored in the resiude type.  This count is performed on the fly, and
	///can hurt performance if reapeatedly carried out.  Not intended for use in large loops -- instead, call
	///once and store the value.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	Size
	n_virtual_atoms () const override;

	/// @brief index of the last backbone heavy atom
	Size
	last_backbone_atom() const
	{
		return n_backbone_heavyatoms_;
	}

	/// @brief index of the first sidechain atom (heavy or hydrogen)
	Size
	first_sidechain_atom() const
	{
		return ( ( n_backbone_heavyatoms_ == nheavyatoms_ ) ? natoms() + 1 : n_backbone_heavyatoms_ + 1 );
	}

	/// @brief index of the first sidehchain hydrogen
	Size
	first_sidechain_hydrogen() const
	{
		return first_sidechain_hydrogen_;
	}

	/// @brief get nbr_atom used to define residue-level neighbors
	Size
	nbr_atom() const
	{
		return nbr_atom_;
	}

	/// @brief get nbr_radius_ used to define residue-level neighbors
	Real
	nbr_radius() const
	{
		return nbr_radius_;
	}

	/// @brief get the molecular weight of this residue
	core::Real
	mass() const { return mass_; }

	/// @brief Returns true if this residue has shadow atoms, false otherwise.
	///
	bool has_shadow_atoms() const { return !atom_shadowed_indices_.empty(); }

	/// @brief Return the index of the atom that the "atom_shadowing"
	/// atom is shadowing; returns zero if the "atom_shadowing" atom is
	/// not shadowing anyone.
	Size
	atom_being_shadowed( Size atom_shadowing ) const {
		return atom_shadowed_indices_[ atom_shadowing ];
	}

	////////////////////////////////////////////////////////////////////////
	/////////////////// Summary Atom Lists              ////////////////////
	////////////////////////////////////////////////////////////////////////

	/// @brief for all heavy atoms, index numbers of their first attached Hydrogen
	AtomIndices const &
	attached_H_begin() const
	{
		return attached_H_begin_;
	}

	/// @brief for all heavy atoms, index numbers of their last attached Hydrogen
	AtomIndices const &
	attached_H_end() const
	{
		return attached_H_end_;
	}

	using ResidueTypeBase::get_metal_binding_atoms;

	/// @brief Gets indices of all atoms that can form bonds to metals
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void
	get_metal_binding_atoms( AtomIndices &metal_binding_indices ) const {
		metal_binding_indices.clear();

		for ( std::string const & name: get_metal_binding_atoms() ) {
			if ( has(name) ) {
				metal_binding_indices.push_back( atom_index( name )  );
			}
		}

		return;
	}

	/// @brief Indices of all backbone atoms, hydrogens and heavyatoms
	AtomIndices const &
	all_bb_atoms() const {
		return all_bb_atoms_;
	}

	/// @brief Indices of all sidechain atoms, hydrogens and heavyatoms
	AtomIndices const &
	all_sc_atoms() const {
		return all_sc_atoms_;
	}


	/// @brief return indices of aromatic Hydrogens
	AtomIndices const &
	Haro_index() const
	{
		return Haro_index_;
	}

	/// @brief return indices of polar Hydrogens
	AtomIndices const &
	Hpol_index() const
	{
		return Hpol_index_;
	}


	/// @brief indices of polar hydrogens as Hbond donors
	AtomIndices const &
	Hpos_polar() const
	{
		return Hpos_polar_;
	}

	/// @brief indices of non-polar hydrogens as potential carbon Hbond donors
	AtomIndices const &
	Hpos_apolar() const
	{
		return Hpos_apolar_;
	}

	AtomIndices const &
	Hpos_polar_sc() const
	{
		return Hpos_polar_sc_;
	}

	/// @brief indices of atoms as Hbond acceptors
	AtomIndices const &
	accpt_pos() const
	{
		return accpt_pos_;
	}

	/// @brief indices of atoms as Hbond acceptors
	AtomIndices const &
	accpt_pos_sc() const
	{
		return accpt_pos_sc_;
	}

	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	/////////////////// Atom Functions              ////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////

	/// @brief is this atom present in this residue?
	bool
	has( std::string const & atom_name ) const override {
		return atom_name_to_index_.count( atom_name );
	}

#ifdef WIN32  // Fixes incorrect cast on WIN32 where has("string") actually calls has( VD )
	inline bool has( const char *name ) const
	{
		return has( std::string(name) );
	}
#endif

	/// @brief get atom name by index
	std::string const &
	atom_name( Size const atomno ) const {
		PyAssert((atomno > 0) && (atomno <= atom_names_.size()), "ResidueType " + name() + " does not have an atom numbered " + std::to_string( atomno ));
		return atom_names_[ atomno ];
	}

	/// @brief get atom index by name
	Size
	atom_index( std::string const & atom_name ) const {
		if ( atom_name_to_index_.count( atom_name ) == 0 ) {
#if defined BOINC
			// chu temporary graphic fix for boinc
			if ( atom_name == "CA" && !is_protein() ) return 1;
#endif
			if ( atom_name == "CA" && is_membrane() ) return 2; // THIS IS NOT THE WAY TO HANDLE THIS! (look into atom aliases instead)
			PyAssert(atom_name_to_index_.count( atom_name ), "ResidueType " + name() + " does not have an atom " + atom_name);
		}
		return atom_name_to_index_.at( atom_name );
	}

#ifdef WIN32
	// Fixes incorrect cast on WIN32 where atom_index("string") actually calls atom_index( VD )
	inline
	Size
	atom_index( const char *name ) const { return atom_index( std::string(name) ); }
#endif

	/// @brief Get the chemical atom_type for this atom by it index number in this residue
	/// @details If we want the atom_type index (integer), we get this from
	/// the conformation::Atom itself, as seen in the code below
	AtomType const &
	atom_type( Size const atomno ) const
	{
		PyAssert((atomno > 0) && (atomno <= vec_atom_types_.size()), "ResidueType::atom_type( Size const atomno ): atomno is not in this ResidueType!");
		return vec_atom_types_[ atomno ];
	}

	core::Size
	atom_type_index(core::Size atomno) const {
		PyAssert((atomno > 0) && (atomno <= vec_atom_type_indexes_.size()), "ResidueType::atom_type_index(): atomno is not in this ResidueType!");
		return vec_atom_type_indexes_[ atomno ];
	}

	ElementCOP
	element_type(core::Size atomno) const {
		PyAssert((atomno > 0) && (atomno <= elements_.size()), "ResidueType::element_type(): atomno is not in this ResidueType!");
		return elements_[ atomno ];
	};

	/// @brief Convenience function to go directly to the element enum
	core::chemical::element::Elements element(core::Size atomno) const;

	/// @brief path distance (number of bonds separated) between a pair of atoms
	int
	path_distance( Size at1, Size at2 ) const { return path_distance_[ at1][at2]; }

	/// @brief shortest path distance for an atom to all other residue atoms
	utility::vector1< int > const &
	path_distance( Size atom ) const
	{
		return path_distance_[ atom ];
	}

	/// @brief accessor of path_distance_ data for this residue, which is a 2D array
	utility::vector1< utility::vector1< int > > const &
	path_distances() const
	{
		return path_distance_;
	}

	/// @brief index number of the first attached Hydrogen on an atom
	Size
	attached_H_begin( Size const atom ) const
	{
		return attached_H_begin_[ atom ];
	}

	/// @brief index number of the last attached Hydrogen on an atom
	Size
	attached_H_end( Size const atom ) const
	{
		return attached_H_end_[ atom ];
	}

	bool
	heavyatom_has_polar_hydrogens( Size atomno ) const;

	bool
	heavyatom_is_an_acceptor( Size atomno ) const {
		return vec_atom_types_[ atomno ].is_acceptor();
	}

	bool
	atom_is_polar_hydrogen( Size atomno ) const {
		return vec_atom_types_[ atomno ].is_polar_hydrogen();
	}

	bool
	atom_is_aro_hydrogen( Size atomno ) const {
		return vec_atom_types_[ atomno ].is_haro();
	}


	/// @brief indices of all mainchain atoms
	AtomIndices const &
	mainchain_atoms() const
	{
		return mainchain_atoms_indices_;
	}

	/// @brief index of mainchain atom
	Size
	mainchain_atom( Size const mainchain_index ) const
	{
		return mainchain_atoms_indices_[mainchain_index];
	}

	int
	formal_charge( Size const atomno ) const {
		PyAssert((atomno > 0) && (atomno <= formal_charge_.size()), "ResidueType::formal_charge(): atomno is not in this ResidueType!");
		return formal_charge_[ atomno ];
	}

	core::Real
	atom_charge( Size const atomno ) const {
		PyAssert((atomno > 0) && (atomno <= partial_charge_.size()), "ResidueType::atom_charge(): atomno is not in this ResidueType!");
		return partial_charge_[ atomno ];
	}

	core::Vector
	ideal_xyz( Size const atomno ) const {
		PyAssert((atomno > 0) && (atomno <= ideal_xyz_.size()), "ResidueType::ideal_xyz(): atomno is not in this ResidueType!");
		return ideal_xyz_[ atomno ];
	}
	core::Vector
	ideal_xyz( std::string const & atomname ) const {
		return ideal_xyz( atom_index( atomname ) );
	}

	/// @brief is a backbone atom (heavy or hydrogen)?
	bool
	atom_is_backbone( Size const atomno ) const
	{
		debug_assert( atomno <= natoms() );
		return ( ( atomno <= n_backbone_heavyatoms_ ) ||
			( atomno > nheavyatoms_ && atomno < first_sidechain_hydrogen_ ) );
	}

	bool
	atom_is_sidechain( Size const atomno ) const
	{
		debug_assert( atomno <= natoms() );
		return all_sc_atoms_.has_value( atomno );
	}

	/// @brief quick lookup: is the atom with the given index a hydrogen or not?
	/// Atoms are sorted so that heavy atoms come first and hydrogen atoms come last.
	bool
	atom_is_hydrogen( Size const atomno ) const
	{
		debug_assert( atomno <= natoms() );
		return atomno > nheavyatoms_;
	}

	/// @brief  Check if atom is virtual.
	bool is_virtual( Size const atomno ) const;

	/// @brief  Check if atom is repulsive.
	bool is_repulsive( Size const atomno ) const;

	bool
	atom_has_property( Size const atomno, AtomProperty const property ) const;

	/// @brief Get the AtomicProperities object for the atom.
	/// (Note this is only the manually set properties, not all the properties of the atom.)
	AtomProperties const &
	atom_properties( Size const atomno ) const {
		return *(atomic_properties_[atomno]);
	}

	void
	show_all_atom_names( std::ostream & out ) const override;

	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	/////////////////// Bond Functions              ////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////

	/// @brief number of bonds for given atom
	Size
	nbonds( Size atom ) const
	{
		return bonded_neighbor_[ atom ].size();
	}

	/// @brief indicates how many proton bonded neighbors an atom has
	Size
	number_bonded_hydrogens( Size const atomno ) const
	{
		if ( attached_H_end_[ atomno ] == 0 || atomno > attached_H_end_.size() ) return 0;
		else return attached_H_end_[ atomno ] - attached_H_begin_[ atomno ] + 1;
	}

	/// @brief indicates how many heavyatom bonded neighbors an atom has
	Size
	number_bonded_heavyatoms( Size const atomno ) const;

	AtomIndices const &
	bonded_neighbor( Size const atomno ) const;

	/// @brief Returns a list of all pairs of atoms which are bonded.
	/// This returned list should only have only one entry per bond
	/// (So [1,2] and [2,1] will not both be present.)
	utility::vector1< std::pair< core::Size, core::Size > >
	bonds() const;

	/// @brief Indicates whether or not two atom indices have a chemical bond linking them.
	/// @details Note that this assumes that the Rosetta machinery is set up so that if
	/// atom 1 is bonded to atom 2, atom 2 is bonded to atom 1.  This function breaks if
	/// that assumption breaks.
	/// @author Vikram K. Mulligan
	bool atoms_are_bonded( core::Size const atomindex1, core::Size const atomindex2 ) const;

	utility::vector1<BondName> const & bonded_neighbor_types(Size const atomno) const;

	utility::vector1<BondRingness> const & bonded_neighbor_ringnesses(Size const atomno) const;

	/// @brief Returns the type of bond between the atoms
	/// Will return UnknownBond for something which isn't bonded.
	BondName
	bond_type(core::Size const atomindex1, core::Size const atomindex2 ) const;

	/// @brief Returns whether the bond between the atoms is in a ring.
	/// Will return UnknownRingness for something which isn't bonded.
	BondRingness
	bond_ringness(core::Size const atomindex1, core::Size const atomindex2 ) const;

	/// @brief indices of the bonded neighbors for an atom
	AtomIndices const &
	cut_bond_neighbor( Size const atomno ) const
	{
		return cut_bond_neighbor_indices_[atomno];
	}

	/// @brief indices of the bonded neighbors for an atom, shortcut for bonded_neighbor(atomno)
	AtomIndices const &
	nbrs( Size const atomno ) const
	{
		return bonded_neighbor(atomno);
	}

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Chi Functions               ////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief number of chi angles
	Size
	nchi() const
	{
		return chi_atoms_.size();
	}

	/// @brief Return number of nu (internal ring) angles.
	Size
	n_nus() const
	{
		return nu_atoms_.size();
	}

	/// @brief  Return the number of rings in this residue.
	Size
	n_rings() const
	{
		return ring_atoms_.size();
	}

	/// @brief Return the number of intraresidue dihedrals.  This covers all pairs
	/// of atoms that are separated by four bonds, and all pairs of intervening
	/// atoms.
	Size
	ndihe() const
	{
		return dihedral_atom_sets_.size();
	}

	/// @brief number of proton chis
	Size
	n_proton_chi() const
	{
		return proton_chis_.size();
	}

	/// @brief number of proton chis
	bool
	is_proton_chi( Size const chino ) const
	{
		return ( std::find( proton_chis_.begin(), proton_chis_.end(), chino ) != proton_chis_.end() );
	}

	/// @brief translate proton_chi to global chi
	Size
	proton_chi_2_chi( Size proton_chi_id ) const
	{
		return proton_chis_[ proton_chi_id ];
	}
	Size
	chi_2_proton_chi( Size chi_index ) const
	{
		return chi_2_proton_chi_[ chi_index ];
	}

	utility::vector1< Real > const &
	proton_chi_samples( Size proton_chi ) const
	{
		return proton_chi_samples_[ proton_chi ];
	}

	utility::vector1< Real > const &
	proton_chi_extra_samples( Size proton_chi ) const
	{
		return proton_chi_extra_samples_[ proton_chi ];
	}

	/// @brief all rotamers bins (mean, std) for a given chi angle
	utility::vector1< std::pair< Real, Real > > const &
	chi_rotamers( Size const chino ) const
	{
		return chi_rotamers_[ chino ];
	}
	/// @brief indices of the atoms which are used to define a given chi angle (chino)
	AtomIndices const &
	chi_atoms( Size const chino ) const
	{
		debug_assert(chino <= chi_atoms_.size());
		return chi_atoms_[ chino ];

	}

	/// @brief indices of the atoms which are used to define all the chi angles
	utility::vector1< AtomIndices > const &
	chi_atoms() const
	{
		return chi_atoms_;
	}

	/// @brief Return indices of the atoms used to define a given nu (internal ring) angle.
	AtomIndices const &
	nu_atoms( core::uint const nu_index ) const
	{
		debug_assert( nu_index <= nu_atoms_.size() );
		return nu_atoms_[nu_index];
	}

	/// @brief Return list of indices of the atoms used to define all the nu (internal ring) angles.
	utility::vector1< AtomIndices > const &
	nu_atoms() const
	{
		return nu_atoms_;
	}

	/// @brief Return list of indices of the atoms within this residue's nth cycle, not counting virtual atoms.
	AtomIndices const &
	ring_atoms( uint const ring_num ) const
	{
		debug_assert( ring_num <= ring_atoms_.size() );
		return ring_atoms_[ ring_num ];
	}

	/// @brief Return list of indices of the atoms within this residue's cycles, not counting virtual atoms.
	utility::vector1< AtomIndices > const &
	ring_atoms() const
	{
		return ring_atoms_;
	}

	/// @brief Return whether this atom is in a particular ring
	bool
	is_ring_atom( uint const ring_num, uint const atom_id ) const;

	/// @brief  Return the saturation level of this residue's nth cycle.
	core::chemical::rings::RingSaturationType
	ring_saturation_type( uint const ring_num ) const
	{
		return ring_saturation_types_[ ring_num ];
	}

	/// @brief  Return the saturation level for all of this residue's cycles.
	utility::vector1< core::chemical::rings::RingSaturationType > const &
	ring_saturation_types() const { return ring_saturation_types_; }

	/// @brief For each ring, what's the lowest conformer?
	utility::vector1< std::string > const &
	lowest_ring_conformers() const { return lowest_ring_conformer_; }

	/// @brief   Low-energy ring conformers for each ring
	utility::vector1< utility::vector1< std::string > > const &
	low_ring_conformers() const { return low_ring_conformers_; }

	/// @brief Read access to the last_controlling_chi_ array
	utility::vector1< Size > const &
	last_controlling_chi() const {
		return last_controlling_chi_;
	}

	/// @brief The last_controlling_chi for an atom.  0 if an atom is controlled by no chi.
	Size
	last_controlling_chi( Size atomno ) const {
		return last_controlling_chi_[ atomno ];
	}

	/// @brief Read access to the atoms_last_controlled_by_chi_ array
	utility::vector1< AtomIndices > const &
	atoms_last_controlled_by_chi() const {
		return atoms_last_controlled_by_chi_;
	}

	/// @brief Read access to the Atoms last controlled by a particular chi
	AtomIndices const &
	atoms_last_controlled_by_chi( Size chi ) const {
		return atoms_last_controlled_by_chi_[ chi ];
	}

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Icoord Functions               //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief get root_atom used as the base of the icoor tree.
	core::Size
	root_atom() const
	{
		return root_atom_;
	}

	/// @brief AtomICoord of an atom
	AtomICoor const &
	icoor( Size const atm ) const;

	/// @brief get index of an atom's base atom
	Size
	atom_base( Size const atomno ) const;

	/// @brief get index of an atom's second base atom
	Size
	abase2( Size const atomno ) const;

	/// @brief Dump out atomnames and icoor values
	void
	debug_dump_icoor() const;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// MMAtom Functions              //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Get the MM atom_type for this atom by its index number in this residue
	MMAtomType const &
	mm_atom_type( Size const atomno ) const;

	std::string const &
	mm_name(core::Size atomno) const;

	core::Size
	mm_atom_type_index(core::Size atomno) const {
		debug_assert( (atomno > 0) && (atomno <= natoms_) );
		return mm_atom_type_index_[ atomno ];
	}

	/// @brief Get the Gasteiger atom type for this atom
	/// Can be null if the atom type isn't set.
	gasteiger::GasteigerAtomTypeDataCOP
	gasteiger_atom_type(core::Size atomno) const;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////          Orbital Functions     //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief indices of the orbitals bonded to an atom
	utility::vector1<core::Size> const &
	bonded_orbitals(Size const atomno) const {
		return bonded_orbitals_[ atomno ];
	}

	//@brief indices of atoms with orbitals
	AtomIndices const &
	atoms_with_orb_index() const
	{
		return atoms_with_orb_index_;
	}

	orbitals::ICoorOrbitalData const &
	orbital_icoor_data(Size const orbital_index) const;

	orbitals::ICoorOrbitalData const &
	new_orbital_icoor_data(Size const orbital_index) const;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////  Ring Conformer Set Functions  //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief    Return a pointer to the object containing the set of ring
	/// conformers possible for this residue's nth cycle.
	core::chemical::rings::RingConformerSetCOP ring_conformer_set( core::uint ring_num ) const;

	core::Size
	n_ring_conformer_sets() const;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////      Connection Functions      //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	// Connections
	// Lower
	ResidueConnection const & lower_connect() const;

	Size
	lower_connect_id() const
	{
		return lower_connect_id_;
	}

	/// @brief index number of the atom which connects to the lower connection
	Size lower_connect_atom() const;

	// Upper
	ResidueConnection const & upper_connect() const;

	Size
	upper_connect_id() const
	{
		return upper_connect_id_;
	}

	/// @brief index number of the atom which connects to the upper connection
	Size upper_connect_atom() const;

	// Branches / Non-polymer
	/// @brief  Return a list of indices of atoms at non-polymer connections.
	utility::vector1< uint > branch_connect_atoms() const;

	/// @brief  Return a list of names of atoms at non-polymer connections.
	utility::vector1< std::string > branch_connect_atom_names() const;

	/// @brief number of ResidueConnections, counting polymeric residue connections
	Size n_possible_residue_connections() const;

	Size
	n_polymeric_residue_connections() const {
		return n_polymeric_residue_connections_;
	}

	Size
	n_non_polymeric_residue_connections() const {
		return n_non_polymeric_residue_connections_;
	}

	/// @brief Get a ResidueConection.
	ResidueConnection const & residue_connection( Size const i ) const;

	/// @brief Does an atom form any inter-residue chemical bonds?
	bool
	atom_forms_residue_connection( Size const atomid ) const {
		return atom_2_residue_connection_map_.at( atomid ).size() != 0;
	}

	/// @brief How many inter-residue chemical bonds does a particular atom form?
	Size
	n_residue_connections_for_atom( Size const atomid ) const {
		return atom_2_residue_connection_map_[ atomid ].size();
	}

	/// @brief Convenience access function for the residue connection
	/// at a particular atom; requires that there is exactly one residue
	/// connection at this atom.
	Size
	residue_connection_id_for_atom( Size const atomid ) const {
		debug_assert( atom_2_residue_connection_map_[ atomid ].size() == 1 );
		return atom_2_residue_connection_map_[ atomid ][ 1 ];
	}

	//// @brief Accessor for the full complement of residue connections for a single atom.
	utility::vector1< Size > const &
	residue_connections_for_atom( Size const atomid ) const {
		return atom_2_residue_connection_map_[ atomid ];
	}

	bool
	residue_connection_is_polymeric( Size const resconn_id ) const {
		return ( resconn_id == lower_connect_id_ || resconn_id == upper_connect_id_ );
	}

	Size
	residue_connect_atom_index( Size const resconn_id ) const;

	/// @brief Does an atom with a given index have an icoor that depends, directly or indirectly, on the lower polymeric connection?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool atom_depends_on_lower_polymeric_connection( core::Size const atom_index ) const;

	/// @brief Does an atom with a given index have an icoor that depends, directly or indirectly, on the upper polymeric connection?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool atom_depends_on_upper_polymeric_connection( core::Size const atom_index ) const;

	/// @brief Does an atom with a given index have an icoor that depends, directly or indirectly, on the upper or lower polymeric connection?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool atom_depends_on_polymeric_connection( core::Size const atom_index ) const;

	/// @brief Does an atom with a given index have an icoor that depends, directly or indirectly, on a particular connection ID?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool atom_depends_on_connection( core::Size const atom_index, core::Size const connection_id ) const;

	//////////////////////////////////////////////////////////////////////
	///////////////// Actcoord                       /////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief require actcoord?
	bool
	requires_actcoord() const;

	/// @brief get indices for atoms used to define actcoord
	AtomIndices const &
	actcoord_atoms() const
	{
		return actcoord_atoms_indices_;

	}

	//////////////////////////////////////////////////////////////////////
	///////////////// properties /////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief is thiol?
	bool is_sidechain_thiol() const { return properties().has_property( SIDECHAIN_THIOL ); }

	/// @brief is disulfide?
	bool is_disulfide_bonded() const { return properties().has_property( DISULFIDE_BONDED ); }

	/// @brief is sidechain amine?
	bool is_sidechain_amine() const { return properties().has_property( SIDECHAIN_AMINE ); }

	/// @brief Is this an alpha-amino acid?
	bool is_alpha_aa() const { return properties().has_property( ALPHA_AA ); }

	/// @brief Is this a beta-amino acid?
	bool is_beta_aa() const { return properties().has_property( BETA_AA ); }

	/// @brief Is this a gamma-amino acid?
	bool is_gamma_aa() const { return properties().has_property( GAMMA_AA ); }

	/// @brief is this a water residue type?
	bool is_water() const { return properties().has_property( WATER ); }

	/// @brief Is this a residue type that can be virtualized or devirtualized by the packer?
	bool is_virtualizable_by_packer() const { return properties().has_property( VIRTUALIZABLE_BY_PACKER ); }

	/// @brief Is this an oligourea?
	bool is_oligourea() const { return properties().has_property( OLIGOUREA ); }

	/// @brief Is this an aramid?
	bool is_aramid() const { return properties().has_property( ARAMID ); }

	/// @brief Is this an ortho-aramid?
	bool is_ortho_aramid() const { return properties().has_property( ORTHO_ARAMID ); }

	/// @brief Is this a meta-aramid?
	bool is_meta_aramid() const { return properties().has_property( META_ARAMID ); }

	/// @brief Is this a para-aramid?
	bool is_para_aramid() const { return properties().has_property( PARA_ARAMID ); }

	/// @brief Is this an ortho-aramid?
	bool is_pre_methylene_ortho_aramid() const { return properties().has_property( PRE_METHYLENE_ORTHO_ARAMID ); }

	/// @brief Is this a meta-aramid?
	bool is_pre_methylene_meta_aramid() const { return properties().has_property( PRE_METHYLENE_META_ARAMID ); }

	/// @brief Is this a para-aramid?
	bool is_pre_methylene_para_aramid() const { return properties().has_property( PRE_METHYLENE_PARA_ARAMID ); }

	/// @brief Is this an ortho-aramid?
	bool is_post_methylene_ortho_aramid() const { return properties().has_property( POST_METHYLENE_ORTHO_ARAMID ); }

	/// @brief Is this a meta-aramid?
	bool is_post_methylene_meta_aramid() const { return properties().has_property( POST_METHYLENE_META_ARAMID ); }

	/// @brief Is this a para-aramid?
	bool is_post_methylene_para_aramid() const { return properties().has_property( POST_METHYLENE_PARA_ARAMID ); }

	/// @brief Is this an ortho-aramid?
	bool is_pre_methylene_post_methylene_ortho_aramid() const { return properties().has_property( PRE_METHYLENE_POST_METHYLENE_ORTHO_ARAMID ); }

	/// @brief Is this a meta-aramid?
	bool is_pre_methylene_post_methylene_meta_aramid() const { return properties().has_property( PRE_METHYLENE_POST_METHYLENE_META_ARAMID ); }

	/// @brief Is this a para-aramid?
	bool is_pre_methylene_post_methylene_para_aramid() const { return properties().has_property( PRE_METHYLENE_POST_METHYLENE_PARA_ARAMID ); }

	/// @brief Is this one of SRI's special heteropolymer building blocks?
	bool is_sri() const;

	/// @brief Is this a triazolemer?
	bool is_triazolemer() const;

	/// @brief Is this a D_AA, R_PEPTOID, or L_RNA?
	/// @details Convenience function to avoid quering is_d_aa(), is_r_peptoid(), and is_l_rna() repeatedly.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_mirrored_type() const { return ( is_d_aa() || is_r_peptoid() || is_l_rna() ); }

	/// @brief Is this residue N-methylated?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_n_methylated() const { return properties().has_property( N_METHYLATED ); }

	/// @brief is TNA?
	bool is_TNA() const { return properties().has_property( TNA ); }

	/// @brief is PNA?
	bool is_PNA() const { return properties().has_property( PNA ); }

	/// @brief is coarse?
	bool is_coarse() const;

	/// @brief is Nucleic Acid?
	bool is_NA() const { return is_DNA() || is_RNA() || is_TNA(); }

	/// @brief is purine?
	bool is_purine() const;

	/// @brief is pyrimidine?
	bool is_pyrimidine() const;

	/// @brief Is this a solvent molecule (SOLVENT property)?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_solvent() const;

	/// @brief Is this a canonical nucleic acid (CANONICAL_NUCLEIC property)?
	/// @details Only the standard nucliec acid types (dA, dC, dG, dT, A, C, G, U) are canonical.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_canonical_nucleic() const;

	/// @brief Is this a canonical amino acid (CANONICAL_AA property)?
	/// @details Only the standard amino acid types (ACDEFGHIKLMNPQRSTVWY) are canonical.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_canonical_aa() const;

	/// @brief Is this a canonical residue type (nucleic acid or amino acid)?
	/// @details Calls is_canonical_aa() and is_canonical_nucleic().
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_canonical() const;

	/// @brief is carbohydrate?
	bool is_carbohydrate() const;

	/// @brief is ligand?
	bool is_ligand() const;

	/// @brief is lipid?
	bool is_lipid() const;

	/// @brief Return true if this residue type is a metal ion, false otherwise.
	bool is_metal() const;

	/// @brief Return true if this residue type is a type that can bind to a metal ion (e.g. His, Asp, Cys, etc.),
	/// false otherwise.
	bool is_metalbinding() const;

	/// @brief is membrane?
	bool is_membrane() const;

	/// @brief is surface? (e.g. enamel)
	bool is_surface() const;

	/// @brief does this residue have sidechain orbitals?
	bool has_sc_orbitals() const;

	/// @brief is polar?
	bool is_polar() const;

	/// @brief is charged?
	bool is_charged() const;

	/// @brief is aromatic?
	bool is_aromatic() const;

	/// @brief is cyclic?
	bool is_cyclic() const;

	/// @brief is terminus?
	bool is_terminus() const { return is_upper_terminus() || is_lower_terminus() || properties().has_property( TERMINUS ); }

	/// @brief is lower terminus?
	bool is_lower_terminus() const;

	/// @brief is upper terminus?
	bool is_upper_terminus() const;

	/// @brief is a branch-point residue?
	bool is_branch_point() const;

	/// @brief is acetylated n terminus
	bool is_acetylated_nterminus() const;

	/// @brief is methylated c terminus
	bool is_methylated_cterminus() const;

	/// @brief  Check if residue is 'VIRTUAL_RESIDUE'
	///  This ONLY checks the VIRTUAL_RESIDUE PROPERTY!
	bool is_virtual_residue() const;

	/// @brief  Check if atom is an inverted virtual
	bool is_inverted_virtual_residue( ) const;

	/// @brief  Check if residue is 'VRT1'
	bool is_VRT1() const;

	/// @brief  Check if residue is a TP3 water.
	bool is_TP3() const;

	/// @brief is an adduct-modified residue?
	bool is_adduct() const;

	/// @brief Does this type have groups (not just single atoms) that are polymer-bond dependent?
	/// @details Always returns false for non-polymeric residue types.
	bool has_polymer_dependent_groups() const {
		return is_polymer() && has_polymer_dependent_groups_;
	}

	///////////////////////////////////////////////////////////////

	/// @brief Is this ResidueTypeBase a base type?
	bool is_base_type() const override;

	/// @brief Get a pointer to this ResidueTypeBase's base ResidueTypeBase.
	/// @details Returns the base_type_cop_ pointer if not null, self pointer if null.
	ResidueTypeCOP get_base_type_cop() const override;

	///////////////////////////////////////////////////////////////
	core::chemical::rna::RNA_Info const & RNA_info() const;

	/// @brief  Return the CarbohydrateInfo object containing sugar-specific properties for this residue.
	core::chemical::carbohydrates::CarbohydrateInfoCOP carbohydrate_info() const;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Dihedral Methods   /////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Return the indices for the set of atoms that define a particular
	/// intraresidue dihedral
	dihedral_atom_set const &
	dihedral( Size const dihe ) const
	{
		return dihedral_atom_sets_[ dihe ];
	}

	/// @brief Returns the list of all of the indices of all the intraresidue
	/// dihedrals a particular atom is involved in.
	/// Useful for calculating the derivatives for an atom.
	utility::vector1< Size > const &
	dihedrals_for_atom( Size atomno ) const
	{
		return dihedrals_for_atom_[ atomno ];
	}

	/// @brief Return the indices for the set of atoms that define a particular
	/// intraresidue angle
	bondangle_atom_set const &
	bondangle( Size const bondang ) const
	{
		return bondangle_atom_sets_[ bondang ];
	}

	/// @brief Returns the list of all of the indices of all the intraresidue
	/// bond angles a particular atom is involved in.
	/// Useful for calculating the derivatives for an atom.
	utility::vector1< Size > const &
	bondangles_for_atom( Size atomno ) const
	{
		return bondangles_for_atom_[ atomno ];
	}

	/// @brief get number of intraresidue bond angles
	Size
	num_bondangles() const
	{
		return bondangle_atom_sets_.size();
	}

	/// @brief print intraresidue bond angles to standard out
	void
	print_bondangles() const; // for debug

	/// @brief print chemical-bond path distances to standard out
	void
	print_pretty_path_distances() const; // for debug

	/// @brief Returns a list of those atoms within one bond of a residue connection.  For residue connection i,
	/// its position in this array is a list of pairs of atom-id's, the first of which is always the id
	/// for the atom forming residue connection i.
	utility::vector1< two_atom_set > const &
	atoms_within_one_bond_of_a_residue_connection( Size resconn ) const
	{
		return atoms_within_one_bond_of_a_residue_connection_[ resconn ];
	}

	/// @brief Returns a list of pairs for atom# atomid where
	/// first == the residue_connection id that lists atomid as being within one bond
	/// of a residue connection, and
	/// second == the index of the entry containing this atom in the
	/// atoms_within_one_bond_of_a_residue_connection_[ first ] array.
	/// Useful for calculating the derivatives for an atom.
	utility::vector1< std::pair< Size, Size > > const &
	within1bonds_sets_for_atom( Size atomid ) const
	{
		return within1bonds_sets_for_atom_[ atomid ];
	}

	/// @brief Returns the list of those atoms within two bonds of
	/// residue connection # resconn.  Each entry in this list is
	/// a triple of atom-id's, the first of which is always the id for
	/// the atom forming residue connection resconn.
	utility::vector1< three_atom_set > const &
	atoms_within_two_bonds_of_a_residue_connection( Size resconn ) const
	{
		return atoms_within_two_bonds_of_a_residue_connection_[ resconn ];
	}

	/// @brief Returns a list of pairs for atom # atomid where
	/// first == the residue_connection id that lists this atom as being within two bonds
	/// of a residue connection, and
	/// second == the index of the entry containing this atom in the
	/// atoms_within_two_bonds_of_a_residue_connection_[ first ] array.
	/// Useful for calculating the derivatives for an atom.
	utility::vector1< std::pair< Size, Size > >
	within2bonds_sets_for_atom( Size atomid ) const {
		return within2bonds_sets_for_atom_[ atomid ];
	}


	/////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	/// @brief Selects three atoms for orienting this residue type
	void
	select_orient_atoms( Size & center, Size & nbr1, Size & nbr2 ) const;

	/// @brief Pick atoms to use for orienting one Residue onto another, using standard logic.
	/// @details Standard logic applies to case in which (a) the residue has backbone atoms, and (b) the residue
	/// has sidechain atoms, and (c) the orient mode has not been set explicitly to force_nbr_atom_orient.  We loop through
	/// all backbone atoms and find the first atom that is bonded to a sidechain atom AND two other backbone atoms.  The
	/// first such atom becomes "center", and its two backbone neighbors become "nbr1" and "nbr2".
	/// @note If ignore_virtuals is true, none of the atoms involved can be virtuals.  If false, they can be.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void select_orient_atoms_standard_logic( Size & center, Size & nbr1, Size & nbr2, bool const ignore_virtuals ) const;

	/// @brief Selects three atoms for orienting this residue type
	/// @note Returns tuple of form [center, nbr1, nbr2]
	std::tuple<Size, Size, Size> select_orient_atoms() const;

	/// @brief  Generate string representation of ResidueType for debugging purposes.
	void show( std::ostream & output=std::cout, bool output_atomic_details=false ) const override;


	////////////////////////////////////////////////////////////////////////////
	// destruction observer methods
public:

	/// @brief attach RestypeDestructionEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( RestypeDestructionEvent const & )
	/// @param ptr **RAW** pointer to observer object
	/// @return Link that can be used to manage the connection.
	/// @remarks RestypeDestructionEvent observers will only be notified upon destruction of the ResidueType
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_destruction_obs( MemFn fn, Ptr ptr ) const {
#ifdef MULTI_THREADED
		std::lock_guard<std::mutex> lock( destruction_obs_mutex_ );
#endif
		return destruction_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach RestypeDestructionEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( RestypeDestructionEvent const & )
	/// @param ptr **RAW** pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	/// @remarks RestypeDestructionEvent observers will only be notified upon destruction of the ResidueType
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_destruction_obs( MemFn fn, Ptr ptr ) const {
#ifdef MULTI_THREADED
		std::lock_guard<std::mutex> lock( destruction_obs_mutex_ );
#endif
		return destruction_obs_hub_.disconnect( fn, ptr );
	}


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//// Private methods for MutableResidueType -> ResidueType conversion ////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
private:

	/// @brief Determine how to map the numeric indices of MutableResidueType to this ResidueType.
	/// ( As the MutableResidueType is const, we can assume a decent mapping here. )
	utility::vector1< core::Size >
	setup_atom_ordering( MutableResidueType const & source );

	/// @brief Iterate through the atoms, filling out the various atom info vectors
	void
	copy_atom_info( MutableResidueType const & source, utility::vector1< core::Size > const & old_to_new );

	/// @brief Fill in various other information which is being copied over from the MutableResidueType
	void
	copy_other_info( MutableResidueType const & source, utility::vector1< core::Size > const & old_to_new );

	/// @brief Precalculate other derived data from data already present
	void
	initialize_derived_data();

	/// @brief Do other pre-calculation checks, not for otherwise uninitialized variables,
	/// but for data consistency or otherwise "incorrectly" set values.
	void
	update_derived_data();

	/// @brief  If a ring has been added but no nu angles have been defined, automatically assign them.
	void auto_assign_nu_atoms();

	/// @brief If a polymer without explicitly set mainchain atoms,
	/// determine a list of main chain atoms by shortest path from LOWER to UPPER.
	/// @details Should only be called from copy_other_data() function.
	void define_mainchain_atoms();

	/// @brief Update the ring conformer sets derived information.
	void update_ring_conformer_sets();

	/// @brief compute the last controlling chi for an atom during the update_derived_data() stage.
	/// The last controlling chi is the furthest chi downstream of the mainchain which can change
	/// the location of an atom.
	/// @details Should only be called from initialize_derived_data() function.
	void
	update_last_controlling_chi();

	/// @brief Recursive subroutine invoked by update_last_controlling_chi().
	void
	note_chi_controls_atom( Size chi, Size atomno );

	/// @brief Determine which atoms are polymer bond-dependent.
	/// @details Should only be called from initialize_derived_data() function.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void update_polymer_dependent_groups();

	/// @brief Determine which atoms are nonpolymer bond-dependent.
	/// @details Should only be called from initialize_derived_data() function.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void update_nonpolymer_dependent_groups();

	/// @brief If there is an NCAARotamerLibrarySpecification, ensure that the rotamer backbone dependicies have been set.
	/// If they have not, set them to all mainchain torsions except omega (the final, inter-residue torsion).
	/// @details Should only be called from update_derived_data() function.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void update_ncaa_rotamer_library_specification_if_present();

	/// @brief Do updates which need a valid get_self_ptr()
	/// (And thus can't be called from the constructor.)
	void
	self_pointer_updates();

	/// @brief Final check of ResidueType data, called by ResidueType::make().
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void
	perform_checks() const;

private:

	//////////////////////////////////////////////////////////////////////
	// Numeric ResidueType properties.

	Size natoms_;
	Size nbonds_;

	/// @brief The number of heavy atoms
	Size nheavyatoms_;

	/// @brief The number of hbond_acceptors
	Size n_hbond_acceptors_;

	/// @brief number of hbond_donors
	Size n_hbond_donors_;

	/// @brief number of backbone heavy atoms
	Size n_backbone_heavyatoms_;

	/// @brief the index of first sidechain hydrogen atom
	Size first_sidechain_hydrogen_;

	//////////////////////////////////////////////////////////////////////
	// per-atom properties

	utility::vector1< std::string > atom_names_;
	std::map< std::string, core::Size > atom_name_to_index_;

	/// @brief The types for each atom.
	utility::vector1< AtomType > vec_atom_types_;
	utility::vector1< core::Size > vec_atom_type_indexes_;

	utility::vector1< ElementCOP > elements_;

	utility::vector1< core::Size > mm_atom_type_index_;
	utility::vector1< core::Size > gasteiger_atom_type_index_;

	utility::vector1< int > formal_charge_;
	utility::vector1< core::Real > partial_charge_;
	utility::vector1< core::Vector > ideal_xyz_;
	utility::vector1< AtomPropertiesOP > atomic_properties_;

	/// @brief Internal coordinates on how to build the given atom
	utility::vector1< AtomICoor > icoor_;

	/// @brief For each atom, list the other atom indices it's bonded to.
	utility::vector1< AtomIndices > bonded_neighbor_;

	/// @brief What is the bond type of the bond connecting the equivalent entry in bonded_neighbor_?
	utility::vector1<utility::vector1<BondName> > bonded_neighbor_type_;

	/// @brief What is the ringness type of the bond connecting the equivalent entry in bonded_neighbor_?
	utility::vector1<utility::vector1<BondRingness> > bonded_neighbor_ringness_;

	/// @brief indices of each heavyatom's first attached hydrogen
	utility::vector1< Size        > attached_H_begin_;

	/// @brief indices of each heavyatom's last attached hydrogen
	utility::vector1< Size        > attached_H_end_;

	/// @brief Data for the mm potentials.
	/// List all of the intra-residue dihedral angles and bond angles.
	/// vector of sets of atoms that make up dihedral angles in the residue
	/// Data for the mm potentials
	utility::vector1< dihedral_atom_set > dihedral_atom_sets_;

	/// @brief all intra-residue dihedral angles that each atom "participates" in
	utility::vector1< utility::vector1< Size > > dihedrals_for_atom_;

	/// @brief vector of sets of atoms that make up bond angles in the residue
	utility::vector1< bondangle_atom_set > bondangle_atom_sets_;

	/// @brief ???
	utility::vector1< utility::vector1< Size > > bondangles_for_atom_;

	/// @brief Data for controlling chi.
	/// Computed in update_last_controlling_chi() for each atom
	/// 0 means an atom whose location is not determined by any chi.
	utility::vector1< Size > last_controlling_chi_;

	/// @brief for chi i, the list of atoms last controlled by i
	/// E.g. chi2 on LEU list cd1, 1hd1, 1hd2, 1hd3, cd2, 2hd1, 2hd2, 2hd3, and hg1
	utility::vector1< AtomIndices > atoms_last_controlled_by_chi_;

	//////////////////////////////////////////////////////////////////////////
	// vectors of indices

	/// @brief For each atom, which orbitals are bonded to it.
	utility::vector1< utility::vector1<Size> > bonded_orbitals_;

	/// @brief indices of atoms with orbitals
	AtomIndices atoms_with_orb_index_;

	/// @brief indices of haro hydrogens
	AtomIndices Haro_index_;

	/// @brief indices of hpolar hydrogens
	AtomIndices Hpol_index_;

	/// @brief indices of Hbond acceptor positions
	AtomIndices accpt_pos_;

	/// @brief indices of polar Hydrogens for Hbond donors
	AtomIndices Hpos_polar_;

	/// @brief indices of apolar hydrogens
	AtomIndices Hpos_apolar_;

	/// @brief indices of Hbond acceptor positions that are part of the sidechain
	/// must be a subset of the atoms listed in the accpt_pos_ array
	AtomIndices accpt_pos_sc_;

	/// @brief indices of polar Hydrogens for Hbond donors that are part of the sidechain
	/// must be a subset of the atoms listed in the Hpos_polar_ array
	AtomIndices Hpos_polar_sc_;

	/// @brief Indices of all backbone atoms, hydrogens and heavyatoms
	AtomIndices all_bb_atoms_;

	/// @brief Indices of all sidechain atoms, hydrogens and heavyatoms
	AtomIndices all_sc_atoms_;

	/// @brief Is this atom a heavyatom with polar hydrogens attached?
	utility::vector1< bool > heavyatom_has_polar_hydrogens_;

	///////////////////////////////////////////////////////////////////////////
	// vectors of vectors of indices

	///// @brief the four atoms to build each chi angle
	utility::vector1< utility::vector1< core::Size > > chi_atoms_;
	/// @brief Is the corresponding chi in chi_atoms_ a proton chi?
	utility::vector1< bool        > is_proton_chi_;
	/// @brief Indices of the chi_atoms_ vector for proton chis
	utility::vector1< Size        > proton_chis_;
	/// @brief A "map" of chi indices to proteon_chi indices
	utility::vector1< Size        > chi_2_proton_chi_;
	/// @brief For a proton chi, the primary samples to diversify the rotamer library with
	utility::vector1< utility::vector1< Real > > proton_chi_samples_;
	/// @brief For a proton chi, how to handle extra ex_ levels
	utility::vector1< utility::vector1< Real > > proton_chi_extra_samples_;

	/// @brief Additional non-Dunbrack rotamer bins
	///
	///    pair<Real,Real>  ==>  mean,sdev
	///    for each chi angle i and rotamer j: chi_rotamers_[i][j]
	///
	utility::vector1< utility::vector1< std::pair< Real, Real > > > chi_rotamers_;

	/// @brief AtomIndex of four atoms to build each nu angle
	utility::vector1< AtomIndices > nu_atoms_;

	/// @brief AtomIndex of all ring atoms, not counting virtual atoms
	utility::vector1< AtomIndices > ring_atoms_;  // indexed by ring number

	/// @brief The saturation type of each ring in this residue
	utility::vector1< core::chemical::rings::RingSaturationType > ring_saturation_types_;  // indexed by ring number

	/// @brief The sets of all possible ring conformers, one per ring
	utility::vector1< rings::RingConformerSetOP > ring_conformer_sets_;

	/// @brief   Lowest-energy ring conformer for each ring
	/// @details used for setting up the RingConformerSets
	utility::vector1< std::string > lowest_ring_conformer_;  // indexed by ring number

	/// @brief   Low-energy ring conformers for each ring
	/// @details used for setting up the RingConformerSets
	utility::vector1< utility::vector1< std::string > > low_ring_conformers_;  // indexed by ring number

	///////////////////////////////////////////////////////////////////////////

	/// @brief number of bonds separated between any pair of atoms in this residue
	utility::vector1< utility::vector1< int > > path_distance_;


	///////////////////////////////////////////////////////////////////////////
	// features

	/// @brief Atom at the nominal root of the ICOOR tree
	core::Size root_atom_;

	///// @brief atom used for calculating residue-level neighbors
	core::Size nbr_atom_;
	/// @brief radius cutoff to define neighbors
	/// @details Should be the maximum distance from the nbr_atom_
	/// to any heavy atom in any valid rotamer.
	Real nbr_radius_;

	///// The isotopically averaged mass of the residue
	Real mass_;

	/// @brief  Vector of inter-residue connections expected for this residuetype
	/// NOW includes the polymer connections, as well as disulf-type connections
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	/// @note  Much of the code assumes at most one residue connection per atom.
	///        The pseudobond code has been partially written to handle cases where multiple connections
	///        to a single atom can exist, but much of the standard residue connection code assumes a
	///        simple correspondence between atoms and residue connections.  That code will have to be
	///        updated to support single-atom "backbones."
	utility::vector1< ResidueConnection > residue_connections_;
	/// @brief Mapping of atom indicies to residue connections
	utility::vector1< utility::vector1< Size > > atom_2_residue_connection_map_;

	/// @brief For calculating inter-residue bond angle and bond torsion energies,
	/// it is useful to have a list of those atoms within one bond of a residue connection atom.
	/// For residue connection i, its position in this array is a list of pairs of atom-id's,
	/// the first of which is always the id for the atom forming residue connection i.
	utility::vector1< utility::vector1< two_atom_set > > atoms_within_one_bond_of_a_residue_connection_;

	/// @brief For atom i, its position in this vector is a list of pairs where
	/// first == the residue_connection id that lists this atom as being within one bond
	/// of a residue connection, and second == the index of the entry containing this atom in the
	/// atoms_within_one_bond_of_a_residue_connection_[ first ] array.
	utility::vector1< utility::vector1< std::pair< Size, Size > > > within1bonds_sets_for_atom_;

	/// @brief For calculating inter-residue bond torsion energies, it is useful to have a list of those
	/// atoms within two bonds of a residue connection atom.  For residue connection i,
	/// its position in this array is a list of triples of atom-id's, the first of which is always the id for
	/// the atom forming residue connection i.
	utility::vector1< utility::vector1< three_atom_set > > atoms_within_two_bonds_of_a_residue_connection_;

	/// @brief For atom i, its position in this vector is a list of pairs where
	/// first == the residue_connection id that lists this atom as being within two bonds
	/// of a residue connection, and second == the index of the entry containing this atom in the
	/// atoms_within_two_bonds_of_a_residue_connection_[ first ] array.
	utility::vector1< utility::vector1< std::pair< Size, Size > > > within2bonds_sets_for_atom_;


	/// @brief Polymer lower connections
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	Size lower_connect_id_; // which connection is the lower connection?

	/// @brief  Polymer upper connections
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	Size upper_connect_id_; // which connection is the upper connection?

	/// @brief Number of non-polymeric residue connections
	Size n_non_polymeric_residue_connections_;
	/// @brief Number of polymeric residue connections
	Size n_polymeric_residue_connections_;

	/// @brief Const-access owning pointer to the base residue type.
	/// @details NULL if this ResidueType is a base type.  If this ResidueType is cloned or copied from a base type,
	/// this will point to the base type.  If this ResidueType is cloned or copied from a non-base type, this will point
	/// to the non-base type's base type.  The test for whether something is a base type is whether this owning pointer is
	/// NULL.
	/// @note The base type must not hold an owning pointer to itself!
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	ResidueTypeCOP base_type_cop_;

	/// @brief A container for properties unique to RNA.
	core::chemical::rna::RNA_InfoOP rna_info_;

	/// @brief A container for residue properties unique to carbohydrates.
	core::chemical::carbohydrates::CarbohydrateInfoOP carbohydrate_info_;

	/// @brief Index version of atom_base_
	utility::vector1<Size> atom_base_indices_;
	/// @brief Index version of abase2_
	utility::vector1<Size> abase2_indices_;
	/// @brief Index version of mainchain_atoms_
	AtomIndices mainchain_atoms_indices_;
	/// @brief Index version of actcoord_atoms_
	AtomIndices actcoord_atoms_indices_;
	/// @brief Index version of cut_bond_neighbor_
	utility::vector1< AtomIndices > cut_bond_neighbor_indices_;
	/// @brief Index version of atom_shadowed_
	utility::vector1< Size > atom_shadowed_indices_;

	/// @brief Does this ResidueType have groups with icoors that depend, directly or indirectly, on polymeric bonds?
	bool has_polymer_dependent_groups_;

	/// @brief vector of length corresponding to number of atoms indicating whether each atom is dependent on the lower
	/// polymer bond.
	/// @details empty if this is a ligand instead of a polymer bond.
	utility::vector1 < bool > atom_depends_on_lower_polymeric_connection_;

	/// @brief vector of length corresponding to number of atoms indicating whether each atom is dependent on the upper
	/// polymer bond.
	/// @details empty if this is a ligand instead of a polymer bond.
	utility::vector1 < bool > atom_depends_on_upper_polymeric_connection_;

	/// @brief An outer vector indexed to connection ID, and an inner vector indexed to atom ids.  "True" means that the atom id for the
	/// inner vector depends on the connection ID for the outer.
	utility::vector1< utility::vector1< bool > > atom_depends_on_connection_;

	/// @brief Who needs to be told if this ResidueType is destroyed?
	mutable utility::signals::BufferedSignalHub < void, RestypeDestructionEvent > destruction_obs_hub_;

#ifdef MULTI_THREADED
	/// @brief Mutex to control access to the destruction_obs_hub_ in multithreaded environment
	mutable std::mutex destruction_obs_mutex_;
#endif

public:

#ifdef    SERIALIZATION
public:
	friend class cereal::access;
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};  // class ResidueType

// Insertion operator (overloaded so that ResidueType can be "printed" in PyRosetta).
std::ostream & operator<<(std::ostream & output, ResidueType const & object_to_output);


} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_ResidueType )

#include <core/chemical/ResidueType.srlz.hh>
SPECIAL_COP_SERIALIZATION_HANDLING( core::chemical::ResidueType, core::chemical::serialize_residue_type, core::chemical::deserialize_residue_type )
#endif // SERIALIZATION

#endif // INCLUDED_core_chemical_Residues_HH
