// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   Residue.hh
/// @brief  Method declarations and simple accessor definitions for the Residue class
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_Residue_hh
#define INCLUDED_core_conformation_Residue_hh


// Unit headers
#include <core/conformation/Residue.fwd.hh>

// Package headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Orbital.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueType.hh> // also defines AtomIndices == vector1< Size >
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>

// No fwd.hh
namespace core { namespace chemical { namespace rings { struct RingConformer; } } }

#include <core/chemical/carbohydrates/CarbohydrateInfo.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzMatrix.fwd.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/access.hpp>
#endif

// C++ headers
#include <map>
#include <limits>



namespace core {
namespace conformation {

/// @brief  Instance Residue class, used for placed residues and rotamers
/**
@details This class is designed to be lightweight.  It holds a const-reference ("rsd_type_")
to a ResidueType object for access to information common to all instances
of a single type, e.g., Alanine or Thymine.  Residue stores any data unique
to a placed residue or rotamer, currently:

- a vector1 of Atoms, which holds the positions (and also the atom-types for
fast access during scoring);

- the sequence position and chain, both integers

- the backbone, side-chain, and internal ring (if applicable) torsion angles (of course backbone torsions are
not unique to a rotamer, and the chi angles are derivable from the coordinates,
but storing them in the residue is convenient for scoring purposes).

- the coordinates of an interaction center or centroid, used e.g., in the
knowledge-based full-atom pair term ("actcoord_").  Maybe this will also
hold the centroid position for centroid-mode scoring??
**/
class Residue : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Residue >
{
public:
	typedef chemical::AtomType       AtomType;
	typedef chemical::ResidueType    ResidueType;
	typedef chemical::ResidueTypeCOP ResidueTypeCOP;
	typedef chemical::AtomIndices    AtomIndices;

public:

	/// @brief constructor with a ResidueTypeCOP.  The dummy arg is there to prevent
	/// implicit casts from a ResidueTypeCOP to a Residue.
	Residue( ResidueTypeCOP rsd_type_in, bool const dummy_arg );

	/// @brief constructor with a ResidueType const &.  The dummy arg is there to prevent
	/// implicit casts from a ResidueType to a Residue.
	Residue( ResidueType const & rsd_type_in, bool const dummy_arg );

	// this is for boost serialize
	//Residue( ResidueType const & rsd_type_in, bool const /*dummy_arg*/, bool const /*dummy_arg2*/ ) :
	// utility::pointer::ReferenceCount(), rsd_type_(rsd_type_in) {}

	/// @brief  Rotamer-style constructor; orients ideal coords onto backbone of current_rsd
	Residue(
		ResidueType const & rsd_type_in,
		Residue const & current_rsd,
		Conformation const & conformation,
		bool preserve_c_beta = false
	);

	/// @brief Copy constructor.
	///
	Residue( Residue const & src );

	/// @brief Copy constructor that preserves everything EXCEPT flips ResidueType to its mirror counterpart.
	/// @details The flip_chirality option governs whether this residue keeps the old type or gets the mirror-image type.  All other data
	/// (including xyz coordinates) are preserved, so this ONLY switches the type, not the geometry.
	/// @author Vikram K. Mulligan.
	Residue( Residue const & src, bool const flip_chirality );

	/// @brief Destructor.
	///
	~Residue() override;

private:

	/// @brief Function called by both copy constructors, to avoid code duplication.
	/// @details As private member variables are added, add them to this to copy them.
	/// @author Vikram K. Mulligan.
	void init_residue_from_other( Residue const &src );

public:

	/// @brief Copy this residue( allocate actual memory for it )
	///
	ResidueOP clone() const;

	/// @brief Copy this residue( allocate actual memory for it ), keeping everything the same EXCEPT the type.
	/// @details Switches the ResidueType to the mirror type (D->L or L->D).  Preserves it for achiral residues.
	/// @note This function is the best way to convert a D-residue to its L-counterpart, or an L-residue to its D-counterpart.
	/// It assumes that you've already mirrored all of the coordinates, and just allows you to generate a replacement residue
	/// of the mirror type that preserves all other Residue information (connections, seqpos, xyz coordinates of all atoms,
	/// variant types, etc.).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	ResidueOP clone_flipping_chirality() const;


	/// @brief  Generate string representation of Residue for debugging purposes.
	void show( std::ostream & output=std::cout, bool output_atomic_details=false ) const;

	/// self pointers
	inline ResidueCOP get_self_ptr() const { return shared_from_this(); }
	inline ResidueOP get_self_ptr() { return shared_from_this(); }
	inline ResidueCAP get_self_weak_ptr() const { return ResidueCAP( shared_from_this() ); }
	inline ResidueAP get_self_weak_ptr() { return ResidueAP( shared_from_this() ); }

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Atom Functions              ////////////////////////
	//////////////////Atom Functions             /////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Returns the AtomType of this residue's atom with index number  <atomno>
	///
	/// example(s):
	///     residue.atom_type(3)
	/// See also:
	///     Residue
	///     Residue.atom_index
	///     AtomType
	///     Pose
	AtomType const &
	atom_type( int const atomno ) const
	{
		return rsd_type_.atom_type( atomno );
	}

	/// @brief Returns the AtomTypeSet of this residue
	///
	/// example(s):
	///     residue.atom_type_set()
	/// See also:
	///     Residue
	///     Residue.atom_type_index
	///     AtomType
	///     Pose
	chemical::AtomTypeSet const &
	atom_type_set() const
	{
		return rsd_type_.atom_type_set();
	}

	/// @brief Returns the atom_type_index of this residue's atom with index number  <atomno>
	/// atom_type_index is used to query this atom's AtomType from an AtomTypeSet,
	/// example: AtomTypeSet[atom_type_index] = AtomType
	///
	/// example(s):
	///     residue.atom_type_index(3)
	/// See also:
	///     Residue
	///     Residue.atom_index
	///     AtomType
	///     Pose
	Size
	atom_type_index( Size const atomno ) const;

	/// @brief Returns the atom charge of this residue's atom with index number  <atomno>
	///
	/// example(s):
	///     residue.atomic_charge(3)
	/// See also:
	///     Residue
	///     Residue.atom_index
	///     Pose
	Real
	atomic_charge( int const atomno ) const;

	/// @brief  Check if atom is virtual.
	bool
	is_virtual( Size const & atomno ) const;

	/// @brief  Check if atom is repulsive.
	bool
	is_repulsive( Size const & atomno ) const;

	/// @brief  Check if residue is virtual.
	bool
	is_virtual_residue() const
	{
		return rsd_type_.is_virtual_residue();
	}

	/// @brief Returns the index number of the  <atm>  in this residue
	/// example: residue.atom_index("CA") returns 2 for a normal amino acid
	///
	/// example(s):
	///     residue.atom_index("CA")
	/// See also:
	///     Residue
	///     AtomType
	///     Pose
	Size
	atom_index( std::string const & atm ) const
	{
		return rsd_type_.atom_index( atm );
	}

	/// @brief Returns the number of atoms in this residue
	///
	/// example(s):
	///     residue.natoms()
	/// See also:
	///     Residue
	///     Pose
	Size
	natoms() const
	{
		return rsd_type_.natoms();
	}

	/// @brief number of hbond_donors
	Size
	n_hbond_acceptors() const
	{
		return rsd_type_.n_hbond_acceptors();
	}

	/// @brief number of hbond_donors
	Size
	n_hbond_donors() const
	{
		return rsd_type_.n_hbond_donors();
	}

	/// @brief Returns the number of heavyatoms in this residue
	///
	/// example(s):
	///     residue.nheavyatoms()
	/// See also:
	///     Residue
	///     Pose

	Size
	nheavyatoms() const
	{
		return rsd_type_.nheavyatoms();
	}


	// NOTE: AtomIndices == vector1< Size >

	/// @brief Returns the AtomIndices of this residue's polar hydrogens
	/// @note: AtomIndices == vector1< Size >
	///
	/// example(s):
	///     residue.Hpos_polar()
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Residue.Hpol_index
	///     Residue.Hpos_apolar
	///     Pose
	AtomIndices const &
	Hpos_polar() const
	{
		return rsd_type_.Hpos_polar();
	}


	/// @brief Returns the AtomIndices of this residue's backbone atoms
	/// @note: heavyatoms and hydrogens, AtomIndices == vector1< Size >
	///
	/// example(s):
	///     residue.all_bb_atoms()
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Pose
	AtomIndices const &
	all_bb_atoms() const {
		return rsd_type_.all_bb_atoms();
	}

	/// @brief Returns the AtomIndices of this residue's aromatic hydrogens
	/// @note: AtomIndices == vector1< Size >
	///
	/// example(s):
	///     residue.Haro_index()
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Residue.Hpol_index
	///     Pose
	AtomIndices const &
	Haro_index() const{
		return rsd_type_.Haro_index();
	}

	/// @brief Returns the AtomIndices of this residue's polar hydrogens
	/// @note: AtomIndices == vector1< Size >
	///
	/// example(s):
	///     residue.Hpol_index()
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Residue.Hpol_index
	///     Residue.Hpos_apolar
	///     Residue.Hpos_polar
	///     Pose
	AtomIndices const &
	Hpol_index() const{
		return rsd_type_.Hpol_index();
	}

	/// @brief Returns the AtomIndices of this residue's apolar hydrogens
	/// @note: AtomIndices == vector1< Size >
	///
	/// example(s):
	///     residue.Hpos_apolar()
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Residue.Hpol_index()
	///     Residue.Hpos_polar
	///     Pose
	AtomIndices const &
	Hpos_apolar() const
	{
		return rsd_type_.Hpos_apolar();
	}

	/// @brief Returns the AtomIndices of this residue's polar sidechain hydrogens
	/// @note: AtomIndices == vector1< Size >
	///
	/// example(s):
	///     residue.Hpos_polar_sc()
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Residue.Hpol_index
	///     Residue.Hpos_polar
	///     Pose
	AtomIndices const &
	Hpos_polar_sc() const
	{
		return rsd_type_.Hpos_polar_sc();
	}


	/// @brief Returns the AtomIndices of this residue's h-bond acceptor atoms
	/// @note: AtomIndices == vector1< Size >
	///
	/// example(s):
	///     residue.accpt_pos()
	/// See also:
	///     Residue
	///     Residue.accpt_pos_sc
	///     Residue.atoms
	///     Pose
	AtomIndices const &
	accpt_pos() const
	{
		return rsd_type_.accpt_pos();
	}

	/// @brief Returns the AtomIndices of this residue's sidechain h-bond acceptor atoms
	/// @note: AtomIndices == vector1< Size >
	///
	/// example(s):
	///     residue.accpt_pos_sc()
	/// See also:
	///     Residue
	///     Residue.accpt_pos
	///     Residue.atoms
	///     Pose
	AtomIndices const &
	accpt_pos_sc() const
	{
		return rsd_type_.accpt_pos_sc();
	}

	/// @brief Gets the AtomIndices of the atoms in this residue that can bind to metals
	/// @note: AtomIndices == vector1< Size >
	/// @author: Vikram K. Mulligan (vmullig@uw.edu)
	void
	get_metal_binding_atoms( AtomIndices &metal_binding_indices ) const
	{
		rsd_type_.get_metal_binding_atoms(metal_binding_indices);
		return;
	}

	/// @brief Returns the number of virtual atoms in this residue.
	/// @details This calls the function with the same name in ResidueType, which counts virts on the fly (memory-efficient, performance-poor).
	/// This being the case, don't call this function repeatedly!  Call it once, and store the return value!
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Size
	n_virtual_atoms() const
	{
		return rsd_type_.n_virtual_atoms();
	}

	/// @brief Is a particular atom a heavy atom with chemically bound polar hydrogens? (i.e. a donor heavy atom)
	bool
	heavyatom_has_polar_hydrogens( Size ind ) const {
		return rsd_type_.heavyatom_has_polar_hydrogens( ind );
	}

	/// @brief Is a particular atom a heavy atom acceptor?
	bool
	heavyatom_is_an_acceptor( Size ind ) const {
		return rsd_type_.heavyatom_is_an_acceptor( ind );
	}

	/// @brief Is a particular atom a polar hydrogen?
	bool
	atom_is_polar_hydrogen( Size ind ) const {
		return rsd_type_.atom_is_polar_hydrogen( ind );
	}


	/// @brief Returns this residue's Atoms (const), a vector1 of Atom objects
	///
	/// example(s):
	///     residue.atoms()
	/// See also:
	///     Residue
	///     Pose
	Atoms const &
	atoms() const
	{
		return atoms_;
	}

	/// @brief Returns this residue's Atoms (non-const), a vector1 of Atom objects
	///
	/// example(s):
	///     residue.atoms()
	/// See also:
	///     Residue
	///     Pose
	Atoms &
	atoms()
	{
		return atoms_;
	}

	/// @brief begin interator, to iterate over atoms
	Atoms::iterator atom_begin() { return atoms_.begin(); }
	/// @brief end interator, to iterate over atoms
	Atoms::iterator atom_end  () { return atoms_.end  (); }

	Atoms::const_iterator atom_begin() const { return atoms_.begin(); }
	Atoms::const_iterator atom_end  () const { return atoms_.end  (); }

	/// @brief should be safe, given the atom ordering rules?
	Atoms::const_iterator sidechainAtoms_begin() const
	{
		return atoms_.begin() + first_sidechain_atom() - 1;
	}
	Atoms::const_iterator heavyAtoms_end() const
	{
		return atoms_.begin() + nheavyatoms();
	}

	/// @brief Returns this residue's Atom with index number  <atm_index>  (const)
	/// @note: Atom object is xyz and atom_type
	///
	/// example(s):
	///     residue.atom(3)
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Pose
	Atom const &
	atom( Size const atm_index ) const
	{
		return atoms_[ atm_index ];
	}

	/// @brief Returns this residue's Atom with index number  <atm_index>  (non-const)
	/// @note: Atom object is xyz and atom_type
	///
	/// example(s):
	///     residue.atom(3)
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Pose
	Atom &
	atom( Size const atm_index )
	{
		return atoms_[ atm_index ];
	}

	/// @brief Returns this residue's Atom with name  <atm_name>  (const)
	/// @note: Atom object is xyz and atom_type, slower but safer than hard-coding an integer index in code where you
	/// need a specific atom
	///
	/// example(s):
	///     residue.atom(3)
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Pose
	Atom const &
	atom( std::string const & atm_name ) const
	{
		return atoms_[ atom_index( atm_name ) ];
	}

	/// @brief Returns this residue's Atom with name  <atm_name>  (non-const)
	/// @note: Atom object is xyz and atom_type, slower but safer than hard-coding an integer index in code where you
	/// need a specific atom
	///
	/// example(s):
	///     residue.atom(3)
	/// See also:
	///     Residue
	///     Residue.atoms
	///     Pose
	Atom &
	atom( std::string const & atm_name )
	{
		return atoms_[ atom_index( atm_name ) ];
	}

	/// @brief Returns the position of this residue's atom with index number  <atm_index>
	/// @note: position is a Vector
	///
	/// example(s):
	///     residue.xyz(3)
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.set_xyz
	///     Pose
	Vector const &
	xyz( Size const atm_index ) const;

	/// @brief Returns the position of this residue's atom with name  <atm_name>
	/// @note: position is a Vector
	///
	/// example(s):
	///     residue.xyz("CA")
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.set_xyz
	///     Pose
	Vector const &
	xyz( std::string const & atm_name ) const;


	/// @brief Sets the position of this residue's atom with index number  <atm_index>
	///
	/// example(s):
	///
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.xyz
	///     Pose
	void
	set_xyz( core::Size const atm_index, Vector const & xyz_in );

	/// @brief Sets the position of this residue's atom with name  <atm_name>
	///
	/// example(s):
	///
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.xyz
	///     Pose
	void
	set_xyz( std::string const & atm_name, Vector const & xyz_in );

	/// @brief Returns the index number of the last backbone heavyatom
	/// @note  The heavyatoms come first in atom ordering,
	/// first backbone then sidechain, hydrogens follow the order
	/// of their attached heavyatom.
	///
	/// example(s):
	///     residue.last_backbone_atom()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.first_sidechain_atom
	///     Pose
	Size
	last_backbone_atom() const
	{
		return rsd_type_.last_backbone_atom();
	}

	/// @brief Returns the index number of the first sidechain heavyatom
	///
	/// example(s):
	///     residue.first_sidechain_atom()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.last_backbone_atom
	///     Pose
	Size
	first_sidechain_atom() const
	{
		return rsd_type_.first_sidechain_atom();
	}

	/// @brief Returns the index number of the first sidechain hydrogen
	///
	/// example(s):
	///     residue.first_sidechain_hydrogen()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.first_sidechain_atom
	///     Pose
	Size
	first_sidechain_hydrogen() const
	{
		return rsd_type_.first_sidechain_hydrogen();
	}


	/// @brief Returns the index number of the first hydrogen attached to the atom
	/// with index number  <atom>
	///
	/// example(s):
	///     residue.attached_H_begin()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.attached_H_end
	///     Pose
	Size
	attached_H_begin( int const atom ) const
	{
		return rsd_type_.attached_H_begin( atom );
	}

	/// @brief Returns the index number of the last hydrogen attached to the atom
	/// with index number  <atom>
	///
	/// example(s):
	///     residue.attached_H_end()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.attached_H_begin
	///     Pose
	Size
	attached_H_end( int const atom ) const
	{
		return rsd_type_.attached_H_end( atom );
	}

	/// @brief Returns the AtomIndices of the first hydrogen attached to each heavyatom
	///
	/// example(s):
	///     residue.attached_H_begin()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.attached_H_end
	///     Residue.nheavyatoms
	///     Pose
	AtomIndices const &
	attached_H_begin() const
	{
		return rsd_type_.attached_H_begin();
	}

	/// @brief Returns the AtomIndices of the last hydrogen attached to each heavyatom
	///
	/// example(s):
	///     residue.attached_H_end()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.attached_H_begin
	///     Residue.nheavyatoms
	///     Pose
	AtomIndices const &
	attached_H_end() const
	{
		return rsd_type_.attached_H_end();
	}

	/// @brief Returns true if residue has a valid lower connection.
	///
	/// example(s):
	///     residue.has_lower_connect()
	bool
	has_lower_connect() const
	{
		return rsd_type_.is_polymer() && (rsd_type_.lower_connect_id() != 0);
	}

	/// @brief Returns the index number of this residue's atom which connects to
	/// the residue before it in sequence
	/// @note: polymers only, example: for an amino acid, residue.lower_connect_atom() = atom_index("N")
	///
	/// example(s):
	///     residue.lower_connect_atom()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.upper_connect_atom
	///     Pose
	Size
	lower_connect_atom() const
	{
		return rsd_type_.lower_connect_atom();
	}

	/// @brief Returns true if residue has a valid upper connection.
	///
	/// example(s):
	///     residue.has_upper_connect()
	bool
	has_upper_connect() const
	{
		return rsd_type_.is_polymer() && (rsd_type_.upper_connect_id() != 0);
	}

	/// @brief Returns the index number of this residue's atom which connects to
	/// the residue after it in sequence
	/// @note: polymers only, example: for an amino acid, residue.upper_connect_atom() = atom_index("C")
	///
	/// example(s):
	///     residue.upper_connect_atom()
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.lower_connect_atom
	///     Pose
	Size
	upper_connect_atom() const
	{
		return rsd_type_.upper_connect_atom();
	}

	/// @brief Returns the index number of this residue's atom connected to the  <other>  Residue
	///
	/// example(s):
	///
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.lower_connect_atom
	///     Residue.upper_connect_atom
	///     Pose
	Size
	connect_atom( Residue const & other ) const;

	/// @brief  Get a list of heavy atoms connected to a given atom.
	utility::vector1< uint > get_adjacent_heavy_atoms( uint const atom_index ) const;

	/// @brief  Scan through the list of atoms connected to a given atom and return the 1st heavy atom found.
	uint first_adjacent_heavy_atom( uint const atom_index ) const;


	/// @brief Returns the shortest path distance from  <atom>  to any other atom in this residue
	///
	/// example(s):
	///
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Pose
	utility::vector1< int > const &
	path_distance( int atom ) const
	{
		return rsd_type_.path_distance( atom );
	}

	/// @brief Returns the shortest path distance for any atom pair in this residue
	/// example: path_distances()[atom1][atom2]
	///
	/// example(s):
	///
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.path_distance
	///     Pose
	utility::vector1< utility::vector1< int > > const &
	path_distances() const
	{
		return rsd_type_.path_distances();
	}
	/// @brief Returns the number of bonds separating atom  <at1>  from  <at2>
	///
	/// example(s):
	///
	/// See also:
	///     Residue
	///     Residue.atom
	///     Residue.atoms
	///     Residue.path_distance
	///     Pose
	int
	path_distance( int at1, int at2 ) const
	{
		return rsd_type_.path_distance( at1, at2 );
	}

	/// @brief Returns true if this residue's atom with index number  <atomno>  is a backbone atom
	///
	/// example(s):
	///     residue.atom_is_backbone(3)
	/// See also:
	///     Residue
	///     Residue.all_bb_atoms
	///     Residue.atom
	///     Residue.atoms
	///     Pose
	bool
	atom_is_backbone( int const atomno ) const
	{
		return rsd_type_.atom_is_backbone( atomno );
	}

	/// @brief Returns true if this residue's atom with index number  <atomno>  is a hydrogen
	///
	/// example(s):
	///     residue.atom_is_backbone(3)
	/// See also:
	///     Residue
	///     Residue.all_bb_atoms
	///     Residue.atom
	///     Residue.atoms
	///     Pose
	bool
	atom_is_hydrogen( Size const atomno ) const
	{
		return rsd_type_.atom_is_hydrogen( atomno );
	}


	/// @brief Returns the atom index of the  <atomno>  atom's base atom
	Size
	atom_base( int const atomno ) const
	{
		return rsd_type_.atom_base( atomno );
	}

	/// @brief Returns the atom index of the  <atomno>  atom's second base atom
	/// note: abase2 is this atom's first bonded neighbor other than
	/// this atom's base atom (unless it has only one neighbor)
	Size
	abase2( int const atomno ) const
	{
		return rsd_type_.abase2( atomno );
	}

	/// @brief Returns the AtomIndices for all bonded neighbor atoms of  <atm>
	AtomIndices const &
	bonded_neighbor( int const atm ) const
	{
		return rsd_type_.bonded_neighbor( atm );
	}


	/// @brief Returns the AtomIndices of this residue's mainchain atoms
	AtomIndices const &
	mainchain_atoms() const
	{
		return rsd_type_.mainchain_atoms();
	}

	/// @brief Returns the atom index of the residue's ith mainchain atom
	Size
	mainchain_atom( Size const i ) const
	{
		return rsd_type_.mainchain_atom( i );
	}

	/// @brief Returns the number of the residue's mainchain atoms
	Size
	n_mainchain_atoms() const
	{
		return rsd_type_.mainchain_atoms().size();
	}


	/// @brief Returns the AtomIndices of atoms that will be used to define this residue's actcoord.
	AtomIndices const &
	actcoord_atoms() const
	{
		return rsd_type_.actcoord_atoms();
	}


	/// @brief Return coordinates for building an atom from ideal internal coordinates,
	/// used for building missing atoms
	Vector
	build_atom_ideal(
		int const atomno,
		Conformation const & conformation // necessary for context, eg the C of the preceding residue for HN
	) const
	{
		return icoor( atomno ).build( *this, conformation );
	}

	/// @brief Returns the index number of this residue's atom used as a center for neighbor definition
	/// example: C-beta atom for some amino acids
	Size
	nbr_atom() const
	{
		return rsd_type_.nbr_atom();
	}

	/// @brief Returns the distance cutoff value used as a radius for neighbor definition
	Real
	nbr_radius() const
	{
		return rsd_type_.nbr_radius();
	}


	Vector const &
	nbr_atom_xyz() const
	{
		return atoms_[ rsd_type_.nbr_atom() ].xyz();
	}


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	/////////////////         Orbital Functions     //////////////////////
	////////////////          Orbital Functions     //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief
	AtomIndices const &
	atoms_with_orb_index() const
	{
		return rsd_type_.atoms_with_orb_index();
	}

	//Vector const
	//orbital_xyz(Size const orbital_index) const;


	Vector
	build_orbital_xyz( Size const orbital_index ) const
	{
		/*  core::chemical::orbitals::ICoorOrbitalData orb_icoor(rsd_type_.orbital_icoor_data(orbital_index));
		Vector stub1_xyz(this->atom(orb_icoor.stub1()).xyz());
		Vector stub2_xyz(this->atom(orb_icoor.stub2()).xyz());
		Vector stub3_xyz(this->atom(orb_icoor.stub3()).xyz());
		Vector orbital_vector(orb_icoor.build(stub1_xyz, stub2_xyz, stub3_xyz));

		//return orbitals_[ orbital_index ].xyz();
		return orbital_vector;*/

		core::chemical::orbitals::ICoorOrbitalData orb_icoor(rsd_type_.new_orbital_icoor_data(orbital_index));
		Vector stub1_xyz(this->atom(orb_icoor.get_stub1()).xyz());
		Vector stub2_xyz(this->atom(orb_icoor.get_stub2()).xyz());
		Vector stub3_xyz(this->atom(orb_icoor.get_stub3()).xyz());

		Vector orbital_vector(orb_icoor.build(stub1_xyz, stub2_xyz, stub3_xyz));
		return orbital_vector;
	}


	Vector const &
	orbital_xyz( Size const orbital_index ) const
	{
		return orbitals_[orbital_index].xyz();
	}

	void
	set_orbital_xyz( core::Size const orbital_index, Vector const & xyz_in )
	{
		orbitals_[ orbital_index ].xyz( xyz_in );
		orbitals_[orbital_index].type(rsd_type_.orbital(orbital_index).orbital_type_index() );
	}


	/// @brief Returns the number of orbitals in this residue
	Size
	n_orbitals() const
	{
		return rsd_type_.n_orbitals();
	}


	utility::vector1<core::Size> const &
	bonded_orbitals(int const atm) const{
		return rsd_type_.bonded_orbitals(atm);
	}


	std::string const &
	orbital_name(int const orbital_index) const{
		return rsd_type_.orbital(orbital_index).name();
	}

	chemical::orbitals::OrbitalType const &
	orbital_type(int const orbital_index) const{
		return rsd_type_.orbital_type(orbital_index);
	}

	Size
	orbital_type_index( Size const orbital_index ) const
	{
		return orbitals_[ orbital_index ].type();
	}

	void
	update_orbital_coords() {
		for (
				auto
				atoms_with_orb_index = rsd_type_.atoms_with_orb_index().begin(),
				atoms_with_orb_index_end = rsd_type_.atoms_with_orb_index().end();
				atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index
				) {
			utility::vector1<core::Size> orbital_indices(rsd_type_.bonded_orbitals(*atoms_with_orb_index));
			for (
					utility::vector1< core::Size >::const_iterator
					orbital_index = orbital_indices.begin(),
					orbital_index_end = orbital_indices.end();
					orbital_index != orbital_index_end; ++orbital_index
					) {
				Vector orb_xyz(this->build_orbital_xyz(*orbital_index));
				this->set_orbital_xyz(*orbital_index, orb_xyz );
			}
		}
	}


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	/////////////////         Residue Functions     //////////////////////
	////////////////          Residue Functions     //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////


	/// @brief Returns a ResidueOP for creating a rotamer of this residue
	/// Temporary hack until Residue hierarchy is worked out
	ResidueOP
	create_rotamer() const
	{
		return clone();
	}

	/// @brief Returns a ResidueOP for creating a copy of residue, same as clone()
	/// Temporary hack until Residue hierarchy is worked out
	ResidueOP
	create_residue() const
	{
		return clone();
	}

	/// @brief Returns this residue's ResidueType
	///
	/// example(s):
	///     residue.type()
	/// See also:
	///     Residue
	///     Residue.atom_type
	ResidueType const &
	type() const
	{
		return rsd_type_;
	}

	/// @brief Returns this residue's ResidueTypeSet
	chemical::ResidueTypeSetCOP
	residue_type_set() const
	{
		return rsd_type_.residue_type_set();
	}


	/// @brief Returns this residue's upper_connection
	/// a ResidueConnection has internal coords info
	/// on how to build the atom in the next residue which
	/// connects to this residue
	chemical::ResidueConnection const &
	upper_connect() const
	{
		return rsd_type_.upper_connect();
	}

	/// @brief Returns this residue's lower_connection
	/// a ResidueConnection has internal coords info
	/// on how to build the atom in the previous residue which
	/// connects to this residue
	chemical::ResidueConnection const &
	lower_connect() const
	{
		return rsd_type_.lower_connect();
	}

	/// Returns true if ???
	bool connections_match( Residue const & other ) const;

	/// @brief Returns the number of ResidueConnections on this residue
	/// including polymeric residue connections.
	///  This is the total number of possible connections from the ResidueType
	Size
	n_possible_residue_connections() const
	{
		return rsd_type_.n_possible_residue_connections();
	}

	/// @brief Returns the number of ResidueConnections on this residue
	/// including polymeric residue connections.
	///  This is the total number of actual connections to other residues.
	///  The index here does not nessessarily match with the connection index as n_possible_residue_connections does!!!!
	Size
	n_current_residue_connections() const;

	/// @brief Returns the number of polymeric ResidueConnections on this residue
	Size
	n_polymeric_residue_connections() const {
		return rsd_type_.n_polymeric_residue_connections();
	}

	/// @brief Returns the number of non-polymeric ResidueConnections on this residue
	Size
	n_non_polymeric_residue_connections() const {
		return rsd_type_.n_non_polymeric_residue_connections();
	}


	/// @brief Returns this residue's ResidueConnection
	/// a ResidueConnection has internal coords info
	/// on how to build the atom in a different residue which
	/// connects to this residue
	chemical::ResidueConnection const &
	residue_connection( int const resconn_index ) const
	{
		return rsd_type_.residue_connection( resconn_index );
	}

	Size
	residue_connect_atom_index( Size const resconn_id ) const {
		return rsd_type_.residue_connect_atom_index( resconn_id );
	}


	chemical::ResConnID
	connect_map( Size resconn_index ) const {
		return connect_map_[ resconn_index ];
	}

	/// @brief Returns the size (number of elements) in the connect_map_ object.
	///
	core::Size
	connect_map_size() const { return connect_map_.size(); }

	void
	clear_residue_connections();

	/// @brief Search through the other residue for connections to this residue, and
	/// ensure that this residue's connect_map is up to date with that residue's
	/// connection indices (and residue number).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void
	update_connections_to_other_residue( Residue const &other_rsd);

	void
	copy_residue_connections_from( Residue const & src );

	bool
	has_incomplete_connection() const;

	/// @brief Returns true is  <atomno>  has complete connectivity?
	bool
	has_incomplete_connection(
		core::Size const atomno
	) const;

	bool
	connection_incomplete( Size resconnid ) const;

	/// @brief Get the residue connection.
	///  How is this different than residue_connection?
	///
	chemical::ResConnID
	actual_residue_connection( Size resconnid ) const {
		return connect_map_[ resconnid ];
	}

	/// @brief Returns the residue number of a residue connected to this residue
	/// at this residue's connection resconn_index.
	/// @details  For example, in a normally-connected pose made of alpha-amino
	/// acids, calling residue_connection_partner(1) on residue 6 should return
	/// 5, since residue 6 is connected to residue 5 at the first connection of
	/// residue 6.
	/// Exactly the same as residue_connection_partner
	Size
	connected_residue_at_resconn( Size const resconn_index ) const {
		return connect_map_[ resconn_index ].resid();
	}

	/// @brief Returns the residue number of a residue connected to this residue
	/// at this residue's connection resconn_index.
	/// @details  For example, in a normally-connected pose made of alpha-amino
	/// acids, calling residue_connection_partner(1) on residue 6 should return
	/// 5, since residue 6 is connected to residue 5 at the first connection of
	/// residue 6.
	Size
	residue_connection_partner( Size const resconn_index ) const
	{
		return connect_map_[ resconn_index ].resid();
	}

	/// attempt to take residue connection info from src_rsd
	void
	copy_residue_connections( Residue const & src_rsd );


	/// @brief Returns the connection id on the OTHER residue connected to this residue at this residue's
	/// connection resconn_index.
	/// @details For example, if this were a normally-connected alpha amino acid, residue_connection_conn_id(1)
	/// would return 2, since the first connection in this residue is attached to the second connection in the
	/// previous residue.
	Size
	residue_connection_conn_id( Size const resconn_index ) const
	{
		return connect_map_[ resconn_index ].connid();
	}


	/// @brief set a connection to this residue by adding its partner's residue number
	void
	residue_connection_partner(
		Size const resconn_index, // ie, our connid
		Size const otherres,
		Size const other_connid
	);

	/// @brief Distance between a potential residue connection match and the position of the expected atom
	Distance
	connection_distance(
		conformation::Conformation const & conf,
		Size const resconn_index,
		Vector const & matchpoint
	) const;

	/// @brief Am I bonded to other?
	/// Meaningful for arbitrary topologies (e.g. circular peptides, disulfides)
	bool
	is_bonded( Residue const & other ) const;

	/// @brief Do I have any pseudobonds to other?
	bool
	is_pseudo_bonded( Residue const & other ) const
	{
		return is_pseudo_bonded( other.seqpos() );
	}

	/// @brief Am I bonded to other?
	/// Looks at all residue connections as opposed to doing arithmetic
	bool
	is_bonded( Size const other_index ) const;

	/// @brief Do I have any pseudobonds to other?
	bool
	is_pseudo_bonded( Size const other_index ) const
	{
		return pseudobonds_.find( other_index ) != pseudobonds_.end();
	}

	/// @brief Am I polymer bonded to other?
	bool
	is_polymer_bonded( Residue const & other ) const;

	/// @brief Am I polymer-bonded to other? checks lower and upper connections
	bool
	is_polymer_bonded( Size const other_index ) const;


	/// @brief Returns the vector1 of resconn ids that connect this residue to other
	utility::vector1< Size > const &
	connections_to_residue( Residue const & other ) const
	{
		return connections_to_residue( other.seqpos() );
	}

	/// @brief Returns the vector1 of resconn ids that connect this residue to other
	utility::vector1< Size > const &
	connections_to_residue( Size const other_resid ) const
	{
		debug_assert( connections_to_residues_.find( other_resid ) != connections_to_residues_.end() );
		return connections_to_residues_.find( other_resid )->second;
	}

	PseudoBondCollectionCOP
	get_pseudobonds_to_residue( Size resid ) const;

	std::map< Size, PseudoBondCollectionCOP > const &
	pseudobonds() const {
		return pseudobonds_;
	}


	void
	set_pseudobonds_to_residue( Size resid, PseudoBondCollectionCOP pbs );


	/// @brief Returns the chi rotamers available for this residue's chi angle  <chino>
	utility::vector1< std::pair< Real, Real > > const &
	chi_rotamers( Size const chino ) const
	{
		return rsd_type_.chi_rotamers( chino );
	}


	/// @brief atom indices for bonded neighbors to which atom-tree connections are disallowed.
	AtomIndices const &
	cut_bond_neighbor( int const atm ) const
	{
		return rsd_type_.cut_bond_neighbor( atm );
	}


	/// @brief Returns the number of atoms bonded to  <atomno>  in all residues?
	core::Size
	n_bonded_neighbor_all_res(
		core::Size const atomno,
		bool virt = false
	) const;


	/// @brief Convenience synonym for bonded_neighbor
	AtomIndices const &
	nbrs( int const atm ) const
	{
		return rsd_type_.bonded_neighbor( atm );
	}


	/// @brief Returns the mainchain torsion angles of this residue (const)
	///
	/// example(s):
	///     residue.mainchain_torsions()
	/// See also:
	///     Residue
	///     Pose
	///     Pose.omega
	///     Pose.phi
	///     Pose.psi
	utility::vector1< Real > const &
	mainchain_torsions() const
	{
		return mainchain_torsions_;
	}

	/// @brief Returns the mainchain torsion angles of this residue (non-const)
	utility::vector1< Real > &
	mainchain_torsions()
	{
		return mainchain_torsions_;
	}

	/// @brief Sets the mainchain torsion angles of this residue to  <torsions>
	///
	/// example(s):
	///     residue.mainchain_torsions()
	/// See also:
	///     Residue
	///     Pose
	///     Pose.set_omega
	///     Pose.set_phi
	///     Pose.set_psi
	void
	mainchain_torsions( utility::vector1< Real > const & torsions ) {
		mainchain_torsions_ = torsions;
	}

	/// @brief Returns a specific mainchain torsion angle for this residue
	/// example: mainchain_torsion(2) will be the psi angle for an amino acid
	///
	/// example(s):
	///     residue.mainchain_torsion(2)
	/// See also:
	///     Residue
	///     Pose
	///     Pose.omega
	///     Pose.phi
	///     Pose.psi
	Real
	mainchain_torsion( Size const torsion ) const
	{
		return mainchain_torsions_[ torsion ];
	}


	/// @brief Returns the number of chi angles this residue has
	///
	/// example(s):
	///     residue.nchi()
	/// See also:
	///     Residue
	///     Pose
	///     Pose.chi
	///     Pose.set_chi
	Size
	nchi() const
	{
		return rsd_type_.nchi();
	}


	/// @brief Returns the chi torsion angles of this residue (const)
	///
	/// example(s):
	///     residue.chi()
	/// See also:
	///     Residue
	///     Residue.nchi
	///     Pose
	///     Pose.chi
	utility::vector1< Real > const &
	chi() const
	{
		return chi_;
	}

	/// @brief Returns the chi torsion angles of this residue (non-const)
	utility::vector1< Real > &
	chi()
	{
		return chi_;
	}

	/// @brief Sets the chi torsion angles of this residue
	///
	/// CAUTION: This function does not cause updating to any internal coordinate data.
	/// See Residue::set_chi() and Residue::set_all_chi() functions for
	/// versions which handle coordinate updates.
	///
	/// example(s):
	///
	/// See also:
	///     Residue
	///     Pose
	///     Pose.set_chi
	void
	chi( utility::vector1< Real > const & chis ) {
		chi_ = chis;
	}

	/// @brief get a specific chi torsion angle
	///
	/// example(s):
	///     residue.chi(1)
	/// See also:
	///     Residue
	///     Pose
	///     Pose.chi
	Real
	chi( Size const chino ) const
	{
		return chi_[ chino ];
	}

	/// @brief Returns the AtomIndices of each set of four atoms defining each chi angle
	utility::vector1< AtomIndices > const &
	chi_atoms() const
	{
		return rsd_type_.chi_atoms();
	}

	/// @brief Returns the AtomIndices of the four atoms defining this residue's  <chino>  chi angle
	AtomIndices const &
	chi_atoms( int const chino ) const
	{
		return rsd_type_.chi_atoms( chino );
	}


	/// @brief   Return the number of nu (internal ring) angles this residue has.
	/// @details Example: residue.n_nus()\n
	/// See also:\n
	///     Residue\n
	///     Residue.nu\n
	///     Residue.nus\n
	///     Pose.set_ring_conformation\n
	///     Pose\n
	///     Pose.nu
	Size
	n_nus() const
	{
		return rsd_type_.n_nus();
	}


	/// @brief   Return the nu (internal ring) torsion angles of this residue.
	/// @details Example: residue.nus()\n
	/// See also:\n
	///     Residue\n
	///     Residue.nu\n
	///     Residue.n_nus\n
	///     Pose.set_ring_conformation\n
	///     Pose\n
	///     Pose.nu
	utility::vector1<core::Angle> const &
	nus() const
	{
		return nus_;
	}

	/// @brief   Return the nu (internal ring) torsion angles of this residue. (non-const)
	utility::vector1<core::Angle> &
	nus()
	{
		return nus_;
	}

	/// @brief   Get a specific nu (internal ring) torsion angle by index.
	/// @details Example: residue.nu(1)\n
	/// See also:\n
	///     Residue\n
	///     Residue.nus\n
	///     Residue.n_nus\n
	///     Pose.set_ring_conformation\n
	///     Pose\n
	///     Pose.nu
	/// @note    It is intentional that there is no set_nu() function; nu angles should only be modified together.
	/// Use Pose.set_ring_conformation() instead.
	core::Angle
	nu(core::uint const index) const
	{
		return nus_[index];
	}

	/// @brief Return the AtomIndices of each set of four atoms defining each nu angle.
	utility::vector1<AtomIndices> const
	nu_atoms() const
	{
		return rsd_type_.nu_atoms();
	}

	/// @brief Return the AtomIndices of the four atoms defining the specified nu angle.
	AtomIndices const
	nu_atoms(core::uint const index) const
	{
		return rsd_type_.nu_atoms(index);
	}


	/// @brief  Return the current RingConformer of this residue's nth ring.
	chemical::rings::RingConformer const & ring_conformer(
		core::uint const ring_num,
		core::Real limit=90.0) const;


	/// @brief Returns the sequence position of this residue
	Size
	seqpos() const
	{
		return seqpos_;
	}

	/// @brief Sets this residue's sequence position to  <setting>
	void
	seqpos( Size const setting )
	{
		seqpos_ = setting;
	}


	/// @brief Returns the sequence separation distance between this residue and  <other>
	/// @note: magnitude of distance only
	Size
	polymeric_sequence_distance( Residue const & other ) const
	{
		if ( ! is_polymer() ||  ! other.is_polymer() ) return Size( -1 );
		return ( chain_ == other.chain() ?
			( seqpos_ <= other.seqpos_ ? other.seqpos_ - seqpos_ : seqpos_ - other.seqpos_ ) : Size( -1 ) );
	}

	/// @brief Returns the sequence separation distance between this residue and  <other>
	/// positive if the other residue is downstream in sequence
	int
	polymeric_oriented_sequence_distance( Residue const & other ) const
	{
		if ( ! is_polymer() || ! other.is_polymer() ) return std::numeric_limits<int>::max();
		return ( chain_ == other.chain() ? other.seqpos_ - seqpos_ : std::numeric_limits<int>::max() );
	}


	/// @brief Returns this residue's chain id
	core::Size
	chain() const
	{
		return chain_;
	}

	/// @brief Sets this residue's chain id
	void
	chain( int const setting )
	{
		chain_ = setting;
	}


	/// @brief does this residue require an actcoord?
	bool
	requires_actcoord() const
	{
		return rsd_type_.requires_actcoord();
	}

	/// @brief Updates actcoord for this residue
	void
	update_actcoord();

	/// @brief Returns the coordinates used for pairE calculations (amino acids only)
	Vector const &
	actcoord() const
	{
		return actcoord_;
	}

	/// @brief Returns the coordinates used for pairE calculations (amino acids only)
	Vector &
	actcoord()
	{
		return actcoord_;
	}

	/// @brief Updates the sequence numbers for this residue and the numbers
	/// stored about its non-polymer connections
	/// called by our owning conformation when the sequence numbers are remapped
	void
	update_sequence_numbering( utility::vector1< Size > const & old2new );

	/// @brief Calculate the set of nu dihedrals from the xyz coordinates and store them.
	void update_nus();


	/////////////
	// properties

	/// @brief Returns true if this residue is a polymer
	bool
	is_polymer() const
	{
		return rsd_type_.is_polymer();
	}

	/// @brief Returns true if this residue is an amino acid
	bool
	is_protein() const
	{
		return rsd_type_.is_protein();
	}

	/// @brief Returns if this residue is a peptoid
	bool
	is_peptoid() const
	{
		return rsd_type_.is_peptoid();
	}

	/// @brief Is this one of SRI's special heteropolymer building blocks?
	bool is_sri() const
	{
		return rsd_type_.is_sri();
	}

	/// @brief Returns true if and only if this residue is a triazolemer
	///
	bool is_triazolemer() const
	{
		return rsd_type_.is_triazolemer();
	}

	/// @brief Returns true if this residue is a DNA residue
	bool
	is_DNA() const
	{
		return rsd_type_.is_DNA();
	}

	/// @brief Returns true if this residue is a RNA residue
	bool
	is_RNA() const
	{
		return rsd_type_.is_RNA();
	}

	/// @brief Returns true if this residue is a nucleic acid
	bool
	is_NA() const
	{
		return rsd_type_.is_NA();
	}

	/// @brief Returns true if this residue is a purine
	bool
	is_purine() const
	{
		return rsd_type_.is_purine();
	}

	/// @brief Returns true if this residue is a pyrimidine
	bool
	is_pyrimidine() const
	{
		return rsd_type_.is_pyrimidine();
	}

	/// @brief Returns true if this residue is a carbohydrate
	bool
	is_carbohydrate() const
	{
		return rsd_type_.is_carbohydrate();
	}

	/// @brief Returns true if this residue is a ligand
	bool
	is_ligand() const
	{
		return rsd_type_.is_ligand();
	}

	/// @brief Returns true if this residue is a metal ion, false otherwise.  The METAL property is specified in the params file under PROPERTIES.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool
	is_metal() const
	{
		return rsd_type_.is_metal();
	}

	/// @brief Returns true if this residue is a type capable of binding to a metal ion (e.g. His, Cys, etc.), false otherwise.  The METALBINDING property is specified in the params file under PROPERTIES.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool
	is_metalbinding() const
	{
		return rsd_type_.is_metalbinding();
	}

	/// @brief Returns true if this residue is a surface residue
	bool
	is_surface() const
	{
		return rsd_type_.is_surface();
	}


	/// @brief Returns true if the residue has side chain orbitals
	bool
	has_sc_orbitals() const
	{
		return rsd_type_.has_sc_orbitals();
	}


	/// @brief Returns true if the residue is polar
	bool
	is_polar() const
	{
		return rsd_type_.is_polar();
	}

	/// @briefReturns true if the residue is apolar
	/// @note: apolar is classified as NOT polar, aromatic, or charged
	bool
	is_apolar() const
	{
		if ( rsd_type_.is_polar() || rsd_type_.is_aromatic() || rsd_type_.is_charged() ) {
			return false;
		} else {
			return true;
		}
	}


	/// @brief Returns true if the residue is charged
	bool
	is_charged() const
	{
		return rsd_type_.is_charged();
	}

	/// @brief Returns true if the residue is aromatic
	bool
	is_aromatic() const
	{
		return rsd_type_.is_aromatic();
	}

	/// @brief residue is coarse (used for RNA right now)
	bool
	is_coarse() const
	{
		return rsd_type_.is_coarse();
	}

	/// @brief Returns true if the residue has a terminus variant
	bool
	is_terminus() const
	{
		return rsd_type_.is_terminus();
	}

	/// @brief Return true if the residue has an upper terminus variant
	bool
	is_upper_terminus() const
	{
		return rsd_type_.is_upper_terminus();
	}

	/// @brief Returns true if the residue has a lower terminus variant
	bool
	is_lower_terminus() const
	{
		return rsd_type_.is_lower_terminus();
	}

	/// @brief Return true if the residue is a branch point variant.
	bool
	is_branch_point() const
	{
		return rsd_type_.is_branch_point();
	}

	/// @brief Return true if the residue is a branch lower terminus variant.
	bool
	is_branch_lower_terminus() const
	{
		return rsd_type_.is_branch_lower_terminus();
	}

	/// @brief Returns true if the chi angles of another residue all fall within 5 deg
	bool
	is_similar_rotamer( Residue const & other ) const;

	/// @brief Returns true if the aa residue types are the same
	bool
	is_similar_aa( Residue const & other ) const
	{

		if ( rsd_type_.aa() != other.aa() ) {
			return false;
		} else {
			return true;
		}
	}

	/// @brief Return true if and only if this residue's type has shadow atoms.
	///
	inline bool has_shadow_atoms() const { return rsd_type_.has_shadow_atoms(); }

	/// @brief Return true if the residue has <property>.
	bool
	has_property( std::string const & property ) const
	{
		return rsd_type_.has_property( property );
	}

	/// @brief  Generic variant access
	bool
	has_variant_type( chemical::VariantType const variant_type ) const
	{
		return rsd_type_.has_variant_type( variant_type );
	}


	/////////////////////////////
	/// @brief Returns the name of this residue's atom with index number  <atm>
	std::string const &
	atom_name( int const atm ) const
	{
		return rsd_type_.atom_name( atm );
	}

	/// @brief Returns the mm_atom_name of this residue's atom with index number  <atom>
	std::string const &
	mm_atom_name(int const atom) const
	{
		return rsd_type_.atom(atom).mm_name();
	}

	/// @brief Returns this residue's ResidueType name
	/// @note: for proteins, this will be the amino acid type and variant type
	std::string const &
	name() const
	{
		return rsd_type_.name();
	}

	/// @brief Returns this residue's 3-letter representation
	/// @note: for proteins, this will be the 3-letter amino acid code
	std::string const &
	name3() const
	{
		return rsd_type_.name3();
	}

	/// @brief Returns this residue's 1-letter representation
	/// @note: for proteins, this will be the 1-letter amino acid code
	char
	name1() const
	{
		return rsd_type_.name1();
	}

	/// @brief Returns this residue's AA type, if any
	/// Used for knowledge-based scores, dunbrack, etc. could be "aa_unk"
	/// AA is enumeration
	chemical::AA const &
	aa() const
	{
		return rsd_type_.aa();
	}

	/// @brief Returns this residue's backbone_aa type, if any.
	/// @details  This is used for noncanonical alpha-amino acids that are templated on canonicals.
	/// For example, 4,5-dihydroxyisoleucine uses the ramachandran and p_aa_pp scoring of isoleucine.
	chemical::AA const &
	backbone_aa() const
	{
		return rsd_type_.backbone_aa();
	}

	/// @brief Returns the internal coordinates of this residue's atom with index number  <atm>
	chemical::AtomICoor const &
	icoor( int const atm ) const
	{
		return rsd_type_.icoor( atm );
	}

	/// @brief Does this atom depend on the LOWER_CONNECT?
	/// @details If recursive is set to true (the default), the dependencies are traced back and the function returns
	/// true if ANY atom that this atom depends on is dependent on LOWER_CONNECT.  If false, only the immediate stub1, stub2, and
	/// stub3 atoms are checked.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool atom_depends_on_lower( core::Size const atom_index, bool const recursive=true ) const;

	/// @brief Does this atom depend on the UPPER_CONNECT?
	/// @details If recursive is set to true (the default), the dependencies are traced back and the function returns
	/// true if ANY atom that this atom depends on is dependent on UPPER_CONNECT.  If false, only the immediate stub1, stub2, and
	/// stub3 atoms are checked.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool atom_depends_on_upper( core::Size const atom_index, bool const recursive=true ) const;

	/// @brief  bondlength analog to set_chi
	///    like set_chi, assumes changes propagate to atomtree
	///    keyed off of chi#, so we only allow distances corresponding to chi angles to refine
	///    distance corresponds to the distance between atoms 3 and 4 defining the chi
	///    chino==0 ==> CA-CB distance, which allows us to refine ALA CB position for example
	/// @author fpd
	void
	set_d( int const chino, Real const setting );

	/// @brief  bondangle analog to set_chi
	///    same idea as set_d
	/// @author fpd
	void
	set_theta( int const chino, Real const setting );

	/// @brief  set_theta for nus
	void
	set_tau( Size const nuno, Real const setting );

	void
	set_all_nu( utility::vector1< Real > const & nus, utility::vector1< Real > const & taus );

	void
	set_all_ring_nu( Size start, Size end, utility::vector1< Real > const & nus, utility::vector1< Real > const & taus );

	/// @brief Sets this residue's chi angle  <chino>  to  <setting>
	/// assuming that changes propagate according to the atom_base tree
	void
	set_chi( int const chino, Real const setting );

	/// @brief Sets all of this residue's chi angles using the set_chi function
	/// (wrapper function)
	void
	set_all_chi( utility::vector1< Real > const & chis );


	// Note: It is intentional that there is no set_nu() function; nu angles should only be modified together.
	// Do not add one.  Use Pose::set_ring_conformation() instead. ~Labonte


	/// @brief Returns true if this residue has an atom named  <atm>
	bool
	has( std::string const & atm ) const
	{
		return rsd_type_.has( atm );
	}

	/// @brief Builds coordinates for atoms missing from this residue
	/// assuming ideal internal coordinates
	void
	fill_missing_atoms(
		utility::vector1< bool > missing,
		Conformation const & conformation
	);

	/// @brief Selects three atoms for orienting this residue
	void
	select_orient_atoms(
		Size & center,
		Size & nbr1,
		Size & nbr2
	) const;

	/// @brief Orient our coords onto those of  <src>, using the atoms from select_orient_atoms
	void
	orient_onto_residue( Residue const & src );

	/// @brief Orient our coords onto those of  <src>, using the three atom pairs specified in the input
	/// @param atom_pairs
	//  Atom pairs used for alignment of the form:
	//  { src_center : center, src_nbr1 : nbr1, src_nbr2 : nbr1 }
	void
	orient_onto_residue(
		Residue const & src,
		utility::vector1< std::pair< std::string, std::string > > const & atom_pairs
	);

	/// @brief Orient our coords onto those of src, uses hard coded atom names (yuck) and will work for peptoid on to peptoid/peptide
	void
	orient_onto_residue_peptoid( Residue const & src, Conformation const & conformation );

	/// @brief Place this rotamer at the sequence position occupied by  <src>
	/// by reorienting the ideal side chain coordinates to match
	void
	place(
		Residue const & src,
		Conformation const & conformation,
		bool preserve_c_beta = false
	);

	/// @brief Applies a transform of the form Rx + v, where R is a rotation
	/// matrix, V is a vector, and x is the original position in xyz space
	void
	apply_transform_Rx_plus_v(
		numeric::xyzMatrix< Real > R,
		Vector v
	);

	/// @brief Return the RNA_residueType object. This is RNA specific.
	core::chemical::rna::RNA_Info const &
	RNA_type() const{
		return rsd_type_.RNA_type();
	}

	/// @brief  Return the CarbohydrateInfo object containing sugar-specific properties for this residue.
	core::chemical::carbohydrates::CarbohydrateInfoCOP carbohydrate_info() const;


	/// @brief BasicDataCache indexed by enum in residue_datacache.hh
	basic::datacache::BasicDataCacheCOP
	data_ptr() const;

	/// @brief BasicDataCache indexed by enum in residue_datacache.hh
	basic::datacache::BasicDataCacheOP
	nonconst_data_ptr();

	/// @brief Is this residue mirrored relative to its coordinates in the ResidueType?
	/// @details Only used for achiral residues like glycine in the context of mirror symmetry and other mirror-image structures.
	inline bool mirrored_relative_to_type() const { return mirrored_relative_to_type_; }

	/// @brief Set whether this residue is mirrored relative to its coordinates in the ResidueType.
	/// @details Only used for achiral residues like glycine in the context of mirror symmetry and other mirror-image structures.
	inline void set_mirrored_relative_to_type( bool const setting ) { mirrored_relative_to_type_ = setting; }

#ifdef SERIALIZATION
	/// @brief Serialize this object
	template < class Archive >
	void
	save( Archive & arc ) const;

	/// @brief Deserialize this object using a non-default-constructor
	/// since the Residue class requires its source ResidueType to be
	/// present before the constructor can be called (since a Residue
	/// holds a ResidueType reference).
	template < class Archive >
	static
	void
	load_and_construct(
		Archive & arc,
		cereal::construct< Residue > & construct
	);
#endif

	/////////////////////////////////////////////////////////////////////////////
	// private methods
	/////////////////////////////////////////////////////////////////////////////

private:

	/// @brief  Updates connections_to_residues_ using connect_map_
	void
	update_connections_to_residues();

	/// @brief apply transform of rotation R and translation V for all atoms downstream
	/// @note note this is not for general atom tree folding. only used in set_chi in which
	/// changes for a chi angle is fast propagated within one residue and not to invoke
	/// folding the whole atom tree.
	void
	apply_transform_downstream(
		int const atomno,
		numeric::xyzMatrix< Real > const & R,
		Vector const & v
	);

	void
	determine_nonstandard_polymer_status();

	/// @brief Assignment operator does not work for class Residue.
	/// This function is intentionally unimplemented and private.
	Residue const &
	operator = ( Residue const & rhs );

	/// @brief Orient coords onto those of <src>, using the specified atoms
	void orient_onto_residue(
		Residue const & src,
		Size center, Size nbr1, Size nbr2,
		Size src_center, Size src_nbr1, Size src_nbr2
	);

	/// @brief Determine the set of orbitals that are needed and compute their
	/// coordinates; invoked from within the constructors.
	void
	assign_orbitals();

	/// @brief Does this atom depend on the LOWER_CONNECT?  Checks recursively.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool atom_depends_on_lower( core::Size const atom_index, utility::vector1< core::Size > visited_atoms ) const;

	/// @brief Does this atom depend on the UPPER_CONNECT?  Checks recursively.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool atom_depends_on_upper( core::Size const atom_index, utility::vector1< core::Size > visited_atoms ) const;

	/// @brief Is a value in a list?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool in_list( core::Size const val, utility::vector1< core::Size > const &vect ) const;


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	/// @brief Pointer to the ResidueType that this Residue derives its chemical
	/// identity from.  The ResidueType is pointed to through both a pointer and
	/// a reference; the pointer for the sake of serialization, and the reference
	/// for access speed.
	ResidueTypeCOP rsd_type_ptr_;

	/// @brief our Residue type
	//-- should be CAP, perhaps?
	ResidueType const & rsd_type_;

	/// @brief our conformation atoms (not kinematic atom pointers) with xyz positions and atom type
	Atoms atoms_;

	utility::vector1< orbitals::OrbitalXYZCoords > orbitals_;


	/// @brief the sequence position
	Size seqpos_;

	/// @brief Is this residue mirrored relative to its coordinates in the ResidueType?
	/// @details Only used for achiral residues like glycine in the context of mirror symmetry and other mirror-image structures.
	/// False by default.
	bool mirrored_relative_to_type_;

	/// @brief the chain id number, starting from 1
	core::Size chain_;

	/// @brief our chi (side chain) torsion angles
	utility::vector1< Real > chi_;

	/// @brief our nu (internal ring) torsion angles
	utility::vector1<Angle> nus_;

	/// @brief our (possibly empty) backbone torsion angles
	utility::vector1< Real > mainchain_torsions_;

	/// @brief the action coordinate, an interaction centroid for knowledge-based terms like fa-pair
	/// in fact, only for fa-pair
	Vector actcoord_;

	basic::datacache::BasicDataCacheOP data_cache_;

	/////////////////////////////////
	/// Inter-residue connection data

	/// @brief true if is_polymer() and either upper_connect or lower_connect (if they exist)
	/// do not connect to seqpos()+1 or seqpos()-1
	bool nonstandard_polymer_;

	/// @brief map between connection ids on this residue and the
	/// connection points on other residues to which its bonded
	utility::vector1< chemical::ResConnID > connect_map_;

	/// @brief lists for each connected residue of the connection points on this residue
	/// that connect the pair.
	std::map< Size, utility::vector1< Size > > connections_to_residues_;

	/// @brief other residues within 4 bonds (connected through PseudoBonds)
	/// may include this residue (intra-residue pseudo-bonds)
	std::map< Size, PseudoBondCollectionCOP > pseudobonds_;

};

std::ostream & operator << ( std::ostream & os, Residue const & res );

} // conformation
} // core

#endif
