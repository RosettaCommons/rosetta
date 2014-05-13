// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ResidueType.hh
/// @brief Method declarations and simple accessors/getters for ResidueType
/// @author
/// Phil Bradley
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to lipids, carbohydrates, and other non-AAs)


#ifndef INCLUDED_core_chemical_ResidueType_hh
#define INCLUDED_core_chemical_ResidueType_hh

// Unit headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>
// Package headers
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.fwd.hh>
#ifdef WIN32
#include <core/chemical/Orbital.hh>
#include <core/chemical/ResidueConnection.hh>
#else
#include <core/chemical/Orbital.fwd.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#endif
#include <core/chemical/sdf/MolData.hh>
#include <core/chemical/rna/RNA_ResidueType.fwd.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/RingConformerSet.fwd.hh>

#include <core/types.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <map>
#include <functional>
#include <string>
#include <core/chemical/VariantType.fwd.hh>
#include <utility/vector1.hh>

// External headers
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>


namespace core {
namespace chemical {

typedef utility::keys::Key2Tuple< Size, Size > two_atom_set;
typedef utility::keys::Key3Tuple< Size, Size, Size > three_atom_set;
typedef utility::keys::Key3Tuple< Size, Size, Size > bondangle_atom_set;
typedef utility::keys::Key4Tuple< Size, Size, Size, Size > dihedral_atom_set;
typedef std::map< VD, std::string > ElementMap;

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
/// The data structure of a ResidueType is based on a boost::graph implementation.  Vertex descriptors (VD, yeah, I
/// know, the name is kind of bad) are the atoms, while the edge descriptors (ED, yet another bad name) are the bonds.
/// Initially, when a ResidueType is constructed, the following primary data are set:
///
///	atom_base_;
///	chi_atoms_;
///	nu_atoms_;
///	mainchain_atoms_;
///	nbr_atom_;
///	actcoord_atoms_;
///	cut_bond_neighbor_;
///	atom_shadowed_;
///
/// When this data is set, it is set based on vertex descriptors.  Because vertex descriptors never change, like atom
/// indices, there is no need to reorder this primary data; however, because Rosetta relies heavily on atom indices
/// to access data, atom indices for the above data have to be generated.  To do this, when finalized is called, a
/// function specifically designed to generate the atom indices for the primary data is called: generate_atom_indices.
/// This function iterates over the vertex descriptors assigned in the primary data and creates private data based on
/// atom indices.  For example, atom_base_ has atom_base_indices_.  When the function atom_base(atomno) is called, the
/// private datum atom_base_indices_ is called.  This allows for the external interface of ResidueType to be accessed by
/// atom indices while the internal functions in ResidueType work off of vertex descriptors.  This also removes the need
/// to have the former old2new reordering scheme.
///
///
/// Atoms: Setting of atoms includes indexing the atoms into vectors, saving their names into vectors/maps, saving the
/// associated mm_atom_type into a vector, saving bond connections into vectors, etc, etc.  On any given residue, the
/// heavy atoms are put into the vector first, (their indices are first,) and hydrogens are put in last.
///
/// Properties: Properties of a residue include things like DNA, PROTEIN, CHARGED, etc.  These properties indicate the
/// type of residue it is and what properties are associated with the residue.  They are set when read in.
/// Several lines of code must be modified to get them to work, all found in ResidueType.cc.
///
/// Orbitals: Orbitals are indexed separately from atoms.  They function much the same way as atoms, except for some
/// key differences.  To find atoms bonded to orbitals, you must provide the atom index, not the orbital index.  (I
/// haven't figured out how to get the reverse to work because of the separate indices.)  Orbital xyz coordinates are
/// not updated when atom coordinates are.  This is to keep speed consistent with just having atoms.  To output the
/// orbitals, use the flag -output_orbitals.
class ResidueType : public utility::pointer::ReferenceCount {

public:
	/// @brief destructor
	virtual ~ResidueType();

	/// @brief constructor
	/// @details We use the AtomTypeSet object to assign atom_types to atoms inside add_atom,
	/// and to identify (polar) hydrogens, acceptors, etc.
	ResidueType(
			AtomTypeSetCAP atom_types,
			ElementSetCAP element_types,
			MMAtomTypeSetCAP mm_atom_types,
			orbitals::OrbitalTypeSetCAP orbital_types//,
	);

	ResidueType(ResidueType const & residue_type);

	/// @brief make a copy
	ResidueTypeOP clone() const;


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Atom Functions              ////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief access by reference the atomset for which this residue is constructed
	AtomTypeSet const &
	atom_type_set() const
	{
		return *atom_types_;
	}

	/// @brief access by reference the atomset for which this residue is constructed
	ElementSet const &
	element_set() const
	{
		return *elements_;
	}

	/// @brief access by const pointer the atomset for which this residue is constructed
	AtomTypeSetCAP
	atom_type_set_ptr() const
	{
		return atom_types_;
	}

	//private: // For refactoring
	Atom & atom(Size const atom_index);
	//public:
	Atom const & atom(Size const atom_index) const;
	Atom & atom(VD const vd);
	Atom const & atom(VD const vd) const;
	Atom & atom(std::string const & atom_name);
	Atom const & atom(std::string const & atom_name) const;

	Orbital const & orbital(Size const orbital_index) const;
	Orbital const & orbital(std::string const & orbital_name) const;

	Bond & bond(ED const ed);
	Bond const & bond(ED const ed) const;

	/// @brief Get the chemical atom_type for this atom by it index number in this residue
	AtomType const &
	atom_type( Size const atomno ) const;

	/// @brief Get the chemical atom_type for this atom by it index number in this residue
	AtomType const &
	atom_type( VD const vd) const;

	/// @brief Get the chemical atom_type index number for this atom by its index number in this residue
	//	int
	//	atom_type_index( Size const atomno ) const;

	/// @brief number of atoms
	Size
	natoms() const
	{
		return graph_.num_vertices();
	}

	/// @brief number of heavy atoms
	Size
	nheavyatoms() const
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
	nbonds() const
	{
		return nbonds_;
	}

	/// @brief number of bonds for given atom
	Size
	nbonds( Size atom ) const
	{
		return boost::out_degree( ordered_atoms_[ atom ] , graph_);
	}

	/// @brief number of bonds for given atom
	Size
	nbonds( VD atom ) const
	{
		return boost::out_degree( atom , graph_ );
	}

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

	///@brief Counts the number of virtual atoms and returns the count.
	///@details The virtual count is not stored in the resiude type.  This count is performed on the fly, and
	///can hurt performance if reapeatedly carried out.  Not intended for use in large loops -- instead, call
	///once and store the value.
	///@author Vikram K. Mulligan (vmullig@uw.edu)
	Size
	n_virtual_atoms () const;

	///@brief indicates how many proton bonded neighbors an atom has
	Size
	number_bonded_hydrogens( Size const atomno ) const
	{
		if( attached_H_end_[ atomno ] == 0 ) return 0;
		else return attached_H_end_[ atomno ] - attached_H_begin_[ atomno ] + 1;
	}

	///@brief indicates how many heavyatom bonded neighbors an atom has
	Size
	number_bonded_heavyatoms( Size const atomno ) const;

	AtomIndices const &
	bonded_neighbor( Size const atomno ) const;

	AdjacentIterPair
	bonded_neighbor_iterators( VD const & atom ) const;


	utility::vector1<BondName> const & bonded_neighbor_types(Size const atomno) const;

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

	/// @brief indices of the atoms which are used to define a given chi angle (chino)
	AtomIndices const &
	chi_atoms( Size const chino ) const
	{
		assert(chino <= chi_atoms_indices_.size());
		return chi_atoms_indices_[ chino ];

	}

	/// @brief indices of the atoms which are used to define all the chi angles
	utility::vector1< AtomIndices > const &
	chi_atoms() const
	{
		return chi_atoms_indices_;
	}

	/// @brief Return indices of the atoms used to define a given nu (internal ring) angle.
	AtomIndices const &
	nu_atoms(core::uint const nu_index) const
	{
		assert(nu_index <= nu_atoms_indices_.size());
		return nu_atoms_indices_[nu_index];
	}

	/// @brief Return list of indices of the atoms used to define all the nu (internal ring) angles.
	utility::vector1< AtomIndices > const &
	nu_atoms() const
	{
		return nu_atoms_indices_;
	}

	/// @brief Gets indices of all atoms that can form bonds to metals
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void
	get_metal_binding_atoms( AtomIndices &metal_binding_indices ) const {
		metal_binding_indices.clear();

		for(core::Size i=1, imax=metal_binding_atoms_.size(); i<=imax; ++i) {
			if( has(metal_binding_atoms_[i]) ) {
				metal_binding_indices.push_back( atom_index( metal_binding_atoms_[i] )  );
			}
		}

		return;
	}

	///@brief Indices of all backbone atoms, hydrogens and heavyatoms
	AtomIndices const &
	all_bb_atoms() const {
		return all_bb_atoms_;
	}

	///@brief Indices of all sidechain atoms, hydrogens and heavyatoms
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

	bool
	heavyatom_has_polar_hydrogens( Size atomno ) const {
		assert( finalized_ );
		return graph_[ordered_atoms_[atomno]].heavyatom_has_polar_hydrogens();
	}

	bool
	heavyatom_is_an_acceptor( Size atomno ) const {
		assert( finalized_ );
		return graph_[ordered_atoms_[atomno]].is_acceptor();
	}

	bool
	atom_is_polar_hydrogen( Size atomno ) const {
		assert( finalized_ );
		return graph_[ordered_atoms_[atomno]].is_polar_hydrogen();
	}


	/// @brief indices of all mainchain atoms
	AtomIndices const &
	mainchain_atoms() const
	{
		return mainchain_atoms_indices_;
	}

	/// @brief index of mainchain atom
	Size
	mainchain_atom( Size const atm ) const
	{
		return mainchain_atoms_indices_[atm];
	}

	/// @brief set indices of all mainchain atoms
	void set_mainchain_atoms( AtomIndices const & mainchain );

	/// @brief is this atom present in this residue?
	bool
	has( std::string const & atom_name ) const
	{
		return atom_name_to_vd_.find(atom_name) != atom_name_to_vd_.end();
	}

	/// @brief get index of an atom's base atom
	Size
	atom_base( Size const atomno ) const;

	/// @brief get index of an atom's second base atom
	Size
	abase2( Size const atomno ) const;

	/// @brief get atom name by index
	std::string const &
	atom_name( Size const index ) const;


	/// @brief get atom index by name
	Size
	atom_index( std::string const & name ) const;

#ifdef WIN32
	// Fixes incorrect cast on WIN32 where atom_index("string") actually calls atom_index( VD )
	inline
	Size
	atom_index( const char *name ) const { return atom_index( std::string(name) ); }
#endif

	Size
	atom_index( VD const & vd) const;

	/// @brief get the vertex descriptor from the name of the atom.
	VD
	vd_from_name( std::string const & name) const;

	/// @brief Get the vertex descriptor from the atom index.
	VD
	vd_from_index(Size const & atomno) const;

	/// @brief Constant access to the underlying graph.
	ResidueGraph const &
	graph() const { return graph_; }

	void
	show_all_atom_names( std::ostream & out ) const;


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
		assert( finalized_ );
		return first_sidechain_hydrogen_;
	}

	VIterPair atom_iterators() const {
		return boost::vertices(graph_);
	}

	EIterPair bond_iterators() const {
		return boost::edges(graph_);
	}

	OutEdgeIterPair bond_iterators( VD const & atom ) const {
		return boost::out_edges(atom, graph_);
	}

	/// @brief is a backbone atom (heavy or hydrogen)?
	bool
	atom_is_backbone( Size const atomno ) const
	{
		assert( finalized_ );
		assert( atomno <= natoms() );
		return ( ( atomno <= n_backbone_heavyatoms_ ) ||
				( atomno > nheavyatoms_ && atomno < first_sidechain_hydrogen_ ) );
	}

	/// @brief quick lookup: is the atom with the given index a hydrogen or not?
	/// Atoms are sorted so that heavy atoms come first and hydrogen atoms come last.
	bool
	atom_is_hydrogen( Size const atomno ) const
	{
		assert( finalized_ );
		assert( atomno <= natoms() );
		return atomno > nheavyatoms_;
	}

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

	/// @brief get indices for atoms used to define actcoord
	AtomIndices const &
	actcoord_atoms() const
	{
		return actcoord_atoms_indices_;

	}

	VD const & get_vertex(Size atomno) const
	{
		return ordered_atoms_[atomno];
	}

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// MMAtom Functions              //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Get the MM atom_type for this atom by its index number in this residue
	MMAtomType const &
	mm_atom_type( Size const atomno ) const;


	/// @brief Get the MM atom_type index number for this atom by its index number in this residue

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Gasteiger Atom Type Functions              //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Get the MM atom_type for this atom by its index number in this residue
	gasteiger::GasteigerAtomTypeDataCOP
	gasteiger_atom_type( Size const atomno ) const;


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////          Orbital Functions     //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	chemical::orbitals::OrbitalType const &
	orbital_type(int const orbital_index) const;

	//	core::Size
	//	orbital_type_index( Size const orb_index ) const;


	/// @brief number of orbitals
	Size
	n_orbitals() const
	{
		return n_orbitals_;
	}

	/// @brief indices of the orbitals bonded to an atom
	utility::vector1<core::Size> const &
	bonded_orbitals(Size const atomno)const
	{
		return graph_[ ordered_atoms_[atomno] ].bonded_orbitals();
	}


	/// @brief is this orbital present in this residue?
	bool
	has_orbital( std::string const & orbital_name ) const
	{
		return ( orbitals_index_.find( orbital_name ) != orbitals_index_.end() );
	}

	//@brief indices of atoms with orbitals
	AtomIndices const &
	atoms_with_orb_index() const
	{
		return atoms_with_orb_index_;
	}


	/// @brief get orbital index by name
	core::Size
	orbital_index( std::string const & name ) const;


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////  Ring Conformer Set Functions  //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief    Return a pointer to the object containing the set of ring
	/// conformers possible for this cyclic residue.
	core::chemical::RingConformerSetCOP ring_conformer_set() const;


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////          Residue Functions     //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	ResidueTypeSet const &
	residue_type_set() const;

	///@brief set the residue type set of origin.
	void
	residue_type_set( ResidueTypeSetCAP set_in );

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


	// Connections
	// Lower
	ResidueConnection const & lower_connect() const;

	Size
	lower_connect_id() const
	{
		return lower_connect_id_;
	}

	/// @brief index number of the atom which connects to the lower connection
	Size
	lower_connect_atom() const;

	/// @brief set the atom which connects to the lower connection
	void
	set_lower_connect_atom( std::string const & atm_name );

	// Upper
	ResidueConnection const & upper_connect() const;

	Size
	upper_connect_id() const
	{
		return upper_connect_id_;
	}

	/// @brief index number of the atom which connects to the upper connection
	Size
	upper_connect_atom() const;

	/// @brief set the atom which connects to the upper connection
	void
	set_upper_connect_atom( std::string const & atm_name );


	/// @brief number of ResidueConnections, counting polymeric residue connections
	Size
	n_residue_connections() const;

	Size
	n_polymeric_residue_connections() const {
		return n_polymeric_residue_connections_;
	}

	Size
	n_non_polymeric_residue_connections() const {
		return n_non_polymeric_residue_connections_;
	}

	/// @brief get a ResidueConection
	ResidueConnection const & residue_connection( Size const i ) const;

	ResidueConnection & residue_connection( Size const i );


	/// @brief Does an atom form any inter-residue chemical bonds?
	bool
	atom_forms_residue_connection( Size const atomid ) {
		return atom_2_residue_connection_map_[ atomid ].size() != 0;
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
		assert( atom_2_residue_connection_map_[ atomid ].size() == 1 );
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


	/// @brief require actcoord?
	inline
	bool
	requires_actcoord() const { return is_protein_ && ( is_polar_ || is_aromatic_ ) && n_actcoord_atoms_ != 0; }

	/// @brief update actcoord
	void
	update_actcoord( conformation::Residue & rot ) const;


	//////////////////////////////////////////////////////////////////////
	/////////////////////////atoms////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief add an atom into this residue
	void
	add_atom(
			std::string const & atom_name,
			std::string const & atom_type_name,
			std::string const & mm_atom_type_name,
			Real const charge /*
		std::string const & count_pair_lower_name,
		std::string const & count_pair_upper_name,
		std::string const & count_pair_special*/
	);

	// Undefined, commenting out to fix PyRosetta build  void add_atom(Atom const & atom);

	/// @brief flag an atom for deletion by adding its index to the delete_atom_ list
	void
	delete_atom( std::string const & name );

	void
	delete_atom( Size const index );

	/// @brief set atom type
	void
	set_atom_type(
			std::string const & atom_name,
			std::string const & atom_type_name
	);


	/// @brief set mm atom type
	void
	set_mm_atom_type(
			std::string const & atom_name,
			std::string const & mm_atom_type_name
	);

	/// @brief Manually set the gasteiger typeset - will use the default set otherwise
	void set_gasteiger_typeset( gasteiger::GasteigerAtomTypeSetCOP gasteiger_atom_types );

	/// @brief set gasteiger atom type
	void
	set_gasteiger_atom_type(
		std::string const & atom_name,
		std::string const & gasteiger_atom_type_name
	);

	/// @brief Add an atom to the list of atoms that can potentially form a bond to a metal ion.
	/// Note that the atom must exist in the residue type (the function checks for this at runtime).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void
	add_metalbinding_atom ( std::string const atom_name );

	/// @brief add a bond between atom1 and atom2, if bond type is not specified, default to SingleBond
	void
	add_bond(
			std::string const & atom_name1,
			std::string const & atom_name2
	);

	// @brief add a bond between atom1 and atom2, specifying a bond type (SingleBond, DoubleBond, TripleBond, AromaticBond)
	void add_bond(std::string const & atom_name1, std::string const & atom_name2, BondName bondLabel);

	/// @brief add a bond between atom1 and atom2, if bond type is not specified, default to a SingleBond
	void
	add_cut_bond(
			std::string const & atom_name1,
			std::string const & atom_name2
	);


	// nbr_atom and nbr_radius are used for rsd-rsd neighbor calculation

	/// @brief set nbr_atom used to define residue-level neighbors
	void
	nbr_atom( std::string const & atom_name )
	{
		nbr_atom_ = ordered_atoms_[atom_index( atom_name )];
		nbr_atom_indices_ = atom_index( atom_name );
	}

	/// @brief get nbr_atom used to define residue-level neighbors
	Size
	nbr_atom() const
	{
		return nbr_atom_indices_;
	}

	/// @brief set nbr_radius_ used to define residue-level neighbors
	void
	nbr_radius( Real const radius )
	{
		nbr_radius_ = radius;
	}

	/// @brief get nbr_radius_ used to define residue-level neighbors
	Real
	nbr_radius() const
	{
		return nbr_radius_;
	}

	/// @brief get the molecular weight of this residue
	core::Real const &
	mass() const //mass
	{
		return mass_;
	}

	/// @brief sets atom_base[ atom1 ] = atom2
	void
	set_atom_base(
			std::string const & atom_name1,
			std::string const & atom_name2
	);


	/// @brief set an atom as backbone heavy atom
	/**
		 backbone stuff is a little tricky if we want to allow newly added atoms,
		 eg in patching, to be backbone atoms. We move any exsiting backbone heavy
		 atoms back into force_bb_ list and add the new one. Afterwards, the new
		 backbone heavy atom list will be generated in finalize() using info from
		 force_bb_.
	 **/
	void
	set_backbone_heavyatom( std::string const & name );

	/// @brief Dump out atomnames and icoor values
	void
	debug_dump_icoor();

	/// @brief AtomICoord of an atom
	AtomICoor const &
	icoor( Size const atm ) const;

	/// @brief set AtomICoor for an atom
	void
	set_icoor(
			std::string const & atm,
			Real const phi,
			Real const theta,
			Real const d,
			std::string const & stub_atom1,
			std::string const & stub_atom2,
			std::string const & stub_atom3,
			bool const update_xyz = false
	);

	/// @brief set AtomICoor for an atom
	void
	set_icoor(
			Size const & index,
			std::string const & atm,
			Real const phi,
			Real const theta,
			Real const d,
			std::string const & stub_atom1,
			std::string const & stub_atom2,
			std::string const & stub_atom3,
			bool const update_xyz = false
	);

	void assign_neighbor_atom();

	void assign_internal_coordinates();

	/// @brief calculate AtomICoor for an atom and set it
	void calculate_icoor(
			std::string const & child,
			std::string const & stub_atom1,
			std::string const & stub_atom2,
			std::string const & stub_atom3
	);

	void
	set_ideal_xyz(
			std::string const & atm,
			Vector const & xyz_in
	);

	void
	set_ideal_xyz(
			Size index,
			Vector const & xyz_in
	);

	void
	set_shadowing_atom(
			std::string const & atom,
			std::string const & atom_being_shadowed
	);

	/// @brief Reassign Rosetta atom types based on the current heuristics.
	/// emap is a map of VD->element strings. If an atom is not present in the element map,
	/// attempt to get the element string from the current type (it's an error if it doesn't have one.)
	/// If preserve is true, only retype those atoms which have an atom_type_index of zero.
	void
	retype_atoms(ElementMap const & emap, bool preserve=false);

	//////////////////////////////////////////////////////////////////////
	/////////////////////////orbitals/////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	///@brief add an orbital onto a residue based upon atom
	void
	add_orbital(
			std::string & orbital_name,
			std::string & orbital_type_name
	);

	///@brief add an orbital bond between an atom and an orbital.
	///@note NOTE!!!!! This is indexed based upon atoms, not orbitals. That means that in your params file
	/// you must have the atom as the first and orbital as the second.
	void
	add_orbital_bond(
			std::string const & atom_name1,
			std::string const & orbital_name
	);

	orbitals::ICoorOrbitalData const &
	orbital_icoor_data(Size const orbital_index) const;

	orbitals::ICoorOrbitalData const &
	new_orbital_icoor_data(Size const orbital_index) const;

	///@brief set OrbitalICoor for an orbital
	void
	set_orbital_icoor_id(
			std::string const & orbital,
			Real const phi,
			Real const theta,
			Real const d,
			std::string const & stub_atom1,
			std::string const & stub_atom2,
			std::string const & stub_atom3
	);

	//////////////////////////////////////////////////////////////////////
	/////////////////////////GRAPHS/////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	const HeavyAtomGraph
	heavy_atoms();

	const AcceptorAtomGraph
	acceptor_atoms();

	const HeavyAtomWithPolarHydrogensGraph
	heavy_atom_with_polar_hydrogens();

	const HeavyAtomWithHydrogensGraph
	heavy_atom_with_hydrogens();

	const HydrogenAtomGraph
	hydrogens();

	const PolarHydrogenGraph
	polar_hydrogens();

	const APolarHydrogenGraph
	apolar_hydrogens();

	const AromaticAtomGraph
	aromatic_atoms();


	//////////////////////////////////////////////////////////////////////
	/////////////////////////residues/////////////////////////////////////
	//////////////////////////////////////////////////////////////////////


	/// @brief Add a chi (side-chain) angle defined by four atoms.
	void
	add_chi(
			Size const chino,
			std::string const & atom_name1,
			std::string const & atom_name2,
			std::string const & atom_name3,
			std::string const & atom_name4
	);

	/// @brief Add a chi (side-chain) angle defined by four atoms to the end of the list of chis.
	void add_chi(
			std::string const & atom_name1,
			std::string const & atom_name2,
			std::string const & atom_name3,
			std::string const & atom_name4);

	/// @brief Add a nu (internal cyclic) angle defined by four atoms.
	void add_nu(
			core::uint const nu_index,
			std::string const & atom_name1,
			std::string const & atom_name2,
			std::string const & atom_name3,
			std::string const & atom_name4);

	/// @brief redefine a chi angle based on four atoms
	//    Added by Andy M. Chen in June 2009
	//    This is needed for certain PTM's
	void
	redefine_chi(
			Size const chino,
			std::string const & atom_name1,
			std::string const & atom_name2,
			std::string const & atom_name3,
			std::string const & atom_name4
	);

	void
	set_proton_chi(
			Size chino,
			utility::vector1< Real > dihedral_samples,
			utility::vector1< Real > extra_samples
	);

	/// @brief Add a rotamer bin for a given chi.
	void
	add_chi_rotamer(
			Size const chino,
			Real const mean,
			Real const sdev
	);

	/// @brief Adds a chi rotamer bin to the highest-indexed chi in the list of chis for this ResidueType.
	void add_chi_rotamer_to_last_chi(core::Angle const mean, core::Angle const sdev);


	/// @brief recalculate derived data, potentially reordering atom-indices
	void
	finalize();

	/// @brief an assertion function to ensure an ResidueType has been finalized
	inline
	void
	require_final() const
	{
		if ( !finalized_ ) {
			utility_exit_with_message( "need to finalize the residue first" );
		}
	}


	/// @brief  add a ResidueConnection
	/// Doesn't set the ideal geometry -- maybe it should?
	Size
	add_residue_connection( std::string const & atom_name );


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// properties /////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief acess specified properties
	utility::vector1< std::string > const &
	properties() const;

	/// @brief add properties
	void
	add_property( std::string const & property );

	///@brief add a numeric property
	void
	add_numeric_property(std::string const & tag,core::Real value);

	///@brief add a string property
	void
	add_string_property(std::string const & tag, std::string value);

	/// @brief delete properties
	//    Added by Andy M. Chen in June 2009
	//    This is needed for deleting properties, which occurs in certain PTM's
	void
	delete_property( std::string const & property );


	/// @brief add an atom to the list for calculating actcoord center
	void add_actcoord_atom( std::string const & atom );


	/// @brief is polymer?
	bool is_polymer() const { return is_polymer_; }

	/// @brief is protein?
	bool is_protein() const { return is_protein_; }

	/// @brief is this an alpha amino acid?
	bool is_alpha_aa() const { return is_alpha_aa_; }

	/// @brief is this a beta amino acid?
	bool is_beta_aa() const { return is_beta_aa_; }

	/// @brief is this a d-amino acid?
	bool is_d_aa() const { return is_d_aa_; }

  /// @brief is this an l-amino acid?
	bool is_l_aa() const { return is_l_aa_; }

	/// @brief is DNA?
	bool is_DNA() const{ return is_DNA_; }

	/// @brief is RNA?
	bool is_RNA() const{ return is_RNA_; }

	/// @brief is coarse?
	bool is_coarse() const{ return is_coarse_; }

	/// @brief is Nucleic Acid?
	bool is_NA() const{ return is_NA_; }

	/// @brief is carbohydrate?
	bool is_carbohydrate() const { return is_carbohydrate_; }

	/// @brief is lipid?
	bool is_lipid() const { return is_lipid_; }

	bool is_ligand() const { return is_ligand_; }

	/// @brief Returns true if this residue type is a metal ion, false otherwise.  The METAL property is specified in the params file under PROPERTIES.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool is_metal() const { return is_metal_; }

	/// @brief Returns true if this residue type is a type that can bind to a metal ion (e.g. His, Asp, Cys, etc.), false otherwise.  The METALBINDING property is specified in the params file under PROPERTIES.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool is_metalbinding() const { return is_metalbinding_; }

	/// @brief is surface? (e.g. enamel)
	bool is_surface() const { return is_surface_; }

	///@brief does this residue have sidechain orbitals?
	bool has_sc_orbitals() const { return has_sc_orbitals_; }

	/// @brief is polar?
	bool is_polar() const { return is_polar_; }

	/// @brief is charged?
	bool is_charged() const { return is_charged_; }

	/// @brief is aromatic?
	bool is_aromatic() const { return is_aromatic_; }

	/// @brief is cyclic?
	bool is_cyclic() const { return is_cyclic_; }

	/// @brief is terminus?
	bool is_terminus() const { return is_terminus_; }

	/// @brief is lower terminus?
	bool is_lower_terminus() const { return is_lower_terminus_; }

	/// @brief is upper terminus?
	bool is_upper_terminus() const { return is_upper_terminus_; }

	/// @brief is lower terminus of a branch?
	bool is_branch_lower_terminus() const { return is_branch_lower_terminus_; }

	/// @brief is acetylated n terminus
	bool is_acetylated_nterminus() const { return is_acetylated_nterminus_; }

	/// @brief is methylated c terminus
	bool is_methylated_cterminus() const { return is_methylated_cterminus_; }

	/// @brief  Check if atom is virtual.
	bool is_virtual( Size const & atomno ) const;

	/// @brief  Check if residue is 'VIRTUAL_RESIDUE'
	bool is_virtual_residue() const { return is_virtual_residue_; }

	/// @brief is an adduct-modified residue?
	bool is_adduct() const { return is_adduct_; }

	void set_adduct_flag( bool adduct_in ) { is_adduct_ = adduct_in; }

	// this probably isnt all that slow unless you have lots of properties, if its ever a problem
	// there is no reason not make properties_ an STL set, since this is basically what sets exist for
	/// @brief  Generic property access -- SLOW!!!!!
	bool
	has_property( std::string const & property ) const
	{
		return ( std::find( properties_.begin(), properties_.end(), property ) != properties_.end() );
	}

	///@brief get a numeric property
	core::Real
	get_numeric_property(std::string const & tag) const
	{
		std::map<std::string, core::Real>::const_iterator property_it(numeric_properties_.find(tag));
		if(property_it == numeric_properties_.end())
		{
			throw utility::excn::EXCN_KeyError(tag + " does not exist in ResidueType with name " + name3_);
			return 0.0; //keep compilers happy
		}

		return property_it->second;
	}

	///@brief get a numeric property
	std::string
	get_string_property(std::string const & tag) const
	{
		std::map<std::string, std::string>::const_iterator property_it(string_properties_.find(tag));
		if(property_it == string_properties_.end())
		{
			throw utility::excn::EXCN_KeyError(tag + " does not exist in ResidueType with name " + name3_);
			return "";
		}
		return property_it->second;
	}

	/// @brief  Generic variant access -- SLOW!!!!!
	bool
	has_variant_type( VariantType const & variant_type ) const
	{
		return ( std::find( variant_types_.begin(), variant_types_.end(), variant_type ) != variant_types_.end() );
	}

	/// @brief get all the variant types for this ResidueType
	utility::vector1< VariantType > const &
	variant_types() const
	{
		return variant_types_;
	}

	/// @brief  Does this residue have exactly the same set of properties as residue other?
	/// phil -- this code does not look correct to me
	/// should probably be other.has_variant_type not other.has_property
	///
	bool
	variants_match( ResidueType const & other ) const;

	/// @brief similar to variants_match(), but allows different
	/// adduct-modified states
	bool
	nonadduct_variants_match( ResidueType const & other ) const;


	/// @brief add one more variant type to this ResidueType
	void
	add_variant_type( VariantType const & variant_type )
	{
		if ( !has_variant_type( variant_type ) ) {
			variant_types_.push_back( variant_type );
		}
	}

	/// @brief set our aa-type (could be "UNK")
	void
	aa( std::string const & type )
	{
		aa_ = aa_from_name( type );
	}
	/// @brief AA to use for rotamer library
	void
	rotamer_aa( std::string const & type )
	{
		rotamer_aa_ = aa_from_name( type );
	}

	/// @brief AA to use for backbone scoring
	void
	backbone_aa( std::string const & type )
	{
		backbone_aa_ = aa_from_name( type );
	}


	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Names    //////////////////////////////////////////
	///////////////// Names     /////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////


	/// @brief get our (unique) residue name
	std::string const &
	name() const
	{
		return name_;
	}

	/// @brief set our (unique) residue name
	void
	name( std::string const & name_in )
	{
		name_ = name_in;
	}


	/// @brief get our 3letter code.  This is set in the
	/// ResidueType .params file through the IO_STRING
	/// tag along with the name1 string
	std::string const &
	name3() const
	{
		return name3_;
	}

	/// @brief set our 3letter code
	void
	name3( std::string const & name_in )
	{
		name3_ = name_in;
	}


	/// @brief get our 1letter code.  This is set in the
	/// ResidueType .params file through the IO_STRING
	/// tag along with the name3 string.
	char
	name1() const
	{
		return name1_;
	}


	/// @brief set our 1letter code
	void
	name1( char const code )
	{
		name1_ = code;
	}

	/// @brief get our interchangeability-group id.  Used to
	/// determine if two residue types are equivalent, except
	/// for their variant status.  E.g. ResidueTypes ALA and
	/// ALA_Nterm would be part of the same interchangeability
	/// group.  This has a degree of subjectivity; are TYR and
	/// pTYR in the same interchangeability group?  Probably not.
	/// This data can be specified in the ResidueTypes .params
	/// file with the INTERCHANGEABILITY_GROUP tag.
	std::string
	interchangeability_group() const {
		return interchangeability_group_;
	}

	/// @brief set our interchangeability-group id
	void
	interchangeability_group( std::string setting ) {
		interchangeability_group_ = setting;
	}

	///@brief set the MolData object
	void set_mol_data(sdf::MolData const & mol_data)
	{
		mol_data_ = mol_data;
	}

	sdf::MolData get_mol_data() const
	{
		return mol_data_;
	}


	/// @brief our traditional residue type, if any
	///
	/// @details Used for knowledge-based scores, dunbrack, etc.
	/// could be "aa_unk".
	///
	/// AA is an enum.  There are values for the 20 standard amino acids, the
	/// 19 canonical D-amino acids, common beta-amino acids and nucleic acids,
	/// and aa_unk as a general catch-all.
	AA const &
	aa() const
	{
		return aa_;
	}

	AA const &
	rotamer_aa() const
	{
		if(rotamer_aa_==aa_unk) return aa_;
		return rotamer_aa_;
	}

	/// @brief Returns the amino acid type to be used for backbone scoring (rama and p_aa_pp).
	AA const &
	backbone_aa() const
	{
		if(backbone_aa_==aa_unk) return aa_;
		return backbone_aa_;
	}


	void set_RotamerLibraryName( std::string const & filename );

	/// @brief A residue parameter file can refer to a set of "pdb rotamers" that can be
	/// superimposed onto a starting position for use in the packer.  These rotamers
	/// are loaded into the pack::dunbrack::RotamerLibrary at the time of their first use.
	std::string get_RotamerLibraryName() const;

	////////////////////////////////////////////////////////////////////////////
	/// dihedral methods
public:

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

	/// @brief Return the number of intraresidue dihedrals.  This covers all pairs
	/// of atoms that are separated by four bonds, and all pairs of intervening
	/// atoms.
	Size
	ndihe() const
	{
		return ndihe_;
	}

	///
	void
	print_dihedrals() const; // for debugging

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

	/// @brief Return the index of the atom that the "atom_shadowing"
	/// atom is shadowing; returns zero if the "atom_shadowing" atom is
	/// not shadowing anyone.
	Size
	atom_being_shadowed( Size atom_shadowing ) const {
		return atom_shadowed_indices_[ atom_shadowing ];
	}

	/// @brief print intraresidue bond angles to standard out
	void
	print_bondangles() const; // for debug

	/// @brief print chemical-bond path distances to standard out
	void
	print_pretty_path_distances() const; // for debug

	/// @brief Sets the path for the NCAA rotlib for the ResidueType
	void
	set_ncaa_rotlib_path( std::string const & path )
	{
		ncaa_rotlib_path_ = path;
	}

	/// @brief Returns the path to the NCAA rotlib for the residue type
	std::string const &
	get_ncaa_rotlib_path() const
	{
		return ncaa_rotlib_path_;
	}

	/// @brief Sets whether we are using a NCAA rotlib for the residue type
	void
	set_use_ncaa_rotlib( bool flag )
	{
		use_ncaa_rotlib_ = flag;
	}

	/// @brief Returns whether we are using a NCAA rotlib for the residue type
	bool
	get_use_ncaa_rotlib() const
	{
		return use_ncaa_rotlib_;
	}

	/// @brief Sets the number of rotatable bonds described by the NCAA rotlib (not nesesarily equal to nchi)
	void
	set_ncaa_rotlib_n_rotameric_bins( Size n_rots )
	{
		ncaa_rotlib_n_rots_ = n_rots;
	}

	/// @brief Returns the number of rotatable bonds described by the NCAA rotlib  (not nesesarily equal to nchi)
	Size
	set_ncaa_rotlib_n_rotameric_bins() const
	{
		return ncaa_rotlib_n_rots_;
	}

	/// @brief Sets the number of rotamers for each rotatable bond described by the NCAA rotlib
	void
	set_ncaa_rotlib_n_bin_per_rot( utility::vector1<Size> n_bins_per_rot );

	/// @brief Returns the number of rotamers for each rotatable bond described by the NCAA rotlib for a single bond
	Size
	get_ncaa_rotlib_n_bin_per_rot( Size n_rot )
	{
		return ncaa_rotlib_n_bins_per_rot_[ n_rot ];
	}

	/// @brief Returns the number of rotamers for each rotatable bond described by the NCAA rotlib for all bonds
	utility::vector1<Size> const &
	get_ncaa_rotlib_n_bin_per_rot() const
	{
		return ncaa_rotlib_n_bins_per_rot_;
	}

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

	///////////////////////////////////////////////////////////////
	core::chemical::rna::RNA_ResidueType const &	RNA_type() const;

	/// @brief  Return the CarbohydrateInfo object containing sugar-specific properties for this residue.
	core::chemical::carbohydrates::CarbohydrateInfoCOP carbohydrate_info() const;

	/// @brief Set force_nbr_atom_orient_, used to control orient atoms selected by select_orient_atoms
	void
	force_nbr_atom_orient( bool force_orient )
	{
		force_nbr_atom_orient_ = force_orient;
	}

	/// @brief Return force_nbr_atom_orient_, used to control orient atoms selected by select_orient_atoms
	bool force_nbr_atom_orient() const
	{
		return force_nbr_atom_orient_;
	}

	/// @brief Selects three atoms for orienting this residue type
	void
	select_orient_atoms(
			Size & center,
			Size & nbr1,
			Size & nbr2
	) const;

	/// @brief A graph-based function to determine the size of the smallest ring that involves a given atom.
	core::Size
	smallest_ring_size( VD const & atom, core::Size const & max_size = 999999) const;

	std::list<utility::vector1<ED> > const & rings(){
		return rings_and_their_edges_;
	}

	/// @brief  Generate string representation of ResidueType for debugging purposes.
	void show( std::ostream & output=std::cout ) const;


public:
	////////////////////////////////////////////////////////////////////////////
	// adduct methods
	/// @brief get the adducts defined for this residue
	utility::vector1< Adduct > const &
	defined_adducts() const
	{
		return defined_adducts_;
	}

	void
	add_adduct( Adduct & adduct_in )
	{
		defined_adducts_.push_back( adduct_in );
	}

	void
	report_adducts();

	////////////////////////////////////////////////////////////////////////////
	// private methods
private:

	/// set up atom ordering map old2new, called by finalize()
	void
	setup_atom_ordering();

	/// GRAPH FUNCTION to provide backward compatibility ////////
	void order_atoms();

	/// reorder primary data in ResidueType given the old2new map, called by finalize()
	void
	generate_atom_indices();

	/// update derived data in ResidueType, called by finalize()
	void
	update_derived_data();

	/// @brief Final check of ResidueType data, called by finalize().
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void
	perform_checks();

	void
	update_residue_connection_mapping();

	/// @brief compute the last controlling chi for an atom during the update_derived_data() stage.
	/// The last controlling chi is the furthest chi downstream of the mainchain which can change
	/// the location of an atom.
	void
	update_last_controlling_chi();

	/// @brief Recursive subroutine invoked by update_last_controlling_chi().
	void
	note_chi_controls_atom( Size chi, Size atomno );

	///@ recursive function to assign internal coordinates
	void assign_internal_coordinates(std::string const & current_atom);



	/////////////////////////////////////////////////////////////////////////////
	// data, see WARNING below
	/////////////////////////////////////////////////////////////////////////////

private:

	// AtomTypeSet Object
	/*
			used to define the set of allowed atomtypes for this residue and
			their properties
	*/
	AtomTypeSetCAP atom_types_;
	ElementSetCAP elements_;

	// MMAtomTypeSet
	MMAtomTypeSetCAP mm_atom_types_;

	// GasteigerAtomTypeSet
	gasteiger::GasteigerAtomTypeSetCOP gasteiger_atom_types_;

	// Orbital types
	orbitals::OrbitalTypeSetCAP orbital_types_;

	RingConformerSetOP conformer_set_;  // set of all possible ring conformers

	ResidueTypeSetCAP residue_type_set_;

	// Graph structures for residuetype
	ResidueGraph graph_; // Stores Atoms and Bonds as Nodes and Edges. First as duplicate material, then on its own.
	std::map<VD, Size> vd_to_index_; //vertex descriptor to an atom index
	std::map<VD, VD> atom_base_; //the atom base
	std::map<VD, VD> abase2_; //the base of the atom base


	// vector of atoms:
	/*
		note not pointers but Atom objects
		currently each Atom holds coords, atom_type, count_pair array index

		Atom order rules:
			(1) heavyatoms before hydrogens
			(2) backbone heavyatoms before sidechain heavyatoms
			(3) hydrogens are grouped by the heavyatom they are attached to
					 and come in the order of those heavyatoms
			(4) as a consequence of (2)+(3) --> backbone hydrogens come before
					 sidechain hydrogens
			(5) atom order in the residue file is preserved subject to rules 1-4
					see finalize() for the logic to determine the atom order

		WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
		WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
		WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

		If you add new properties that are associated with atoms, You need
		to make sure that you iterate over the VDs! Then in the generate
		indices function, associate your VDs to atom ordering.

		There is no more old2new ordering!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
		WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
		WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

	*/

	utility::vector1< Orbital > orbitals_;

	//////////////////////////////////////////////////////////////////
	// ints -- see the WARNING above if these are atom indices
	// number of heavy atoms
	Size nheavyatoms_;

	// number of hbond_acceptors
	Size n_hbond_acceptors_;

	// number of hbond_donors
	Size n_hbond_donors_;

	// number of orbitals
	Size n_orbitals_;

	// number of backbone heavy atoms
	Size n_backbone_heavyatoms_;

	// the index of first sidechain hydrogen atom
	Size first_sidechain_hydrogen_;

	// number of dihedral angle atom sets
	Size ndihe_;

	// number of bonds
	Size nbonds_;


	//////////////////////////////////////////////////////////////////////
	// per-atom properties
	// indices of the atoms psuedo bonded atoms. Used in orbital code
	utility::vector1< AtomIndices > bonded_neighbor_;
	utility::vector1<utility::vector1<BondName> > bonded_neighbor_type_;
	std::map<VD, utility::vector1<VD> > cut_bond_neighbor_;

	// indices of each heavyatom's first attached hydrogen
	utility::vector1< Size        > attached_H_begin_;

	// indices of each heavyatom's last attached hydrogen
	utility::vector1< Size        > attached_H_end_;
	utility::vector1<VD> parents_;
	std::map< VD, AtomICoor > icoor_;

	// Data for the mm potentials.  List all of the intra-residue dihedral angles and bond angles.
	// vector of sets of atoms that make up dihedral angles in the residue
	utility::vector1< dihedral_atom_set > dihedral_atom_sets_;

	// all intra-residue dihedral angles that each atom "participates" in.
	utility::vector1< utility::vector1< Size > > dihedrals_for_atom_;

	// vector of sets of atoms that make up bond angles in the residue
	utility::vector1< bondangle_atom_set > bondangle_atom_sets_;
	utility::vector1< utility::vector1< Size > > bondangles_for_atom_;

	// Data to describe virtual atoms that should shadow other atoms for the sake
	// of keeping intra-residue cycles closed when working with an atom tree, e.g.
	// NV shadows N on proline.  For each atom, the following vector lists the index
	// of the atom it is shadowing.
	std::map<VD, VD> atom_shadowed_;

	// Data for controlling chi.  Dependent data, computed in update_last_controlling_chi()
	// for each atom, the last controlling chi angle for that atom.  a chi of 0 represents
	// an atom whose location is not determined by any chi.
	utility::vector1< Size > last_controlling_chi_;

	// for chi i, the list of atoms last controlled by i.  E.g. chi2 on LEU
	// list cd1, 1hd1, 1hd2, 1hd3, cd2, 2hd1, 2hd2, 2hd3, and hg1
	utility::vector1< AtomIndices > atoms_last_controlled_by_chi_;

	//////////////////////////////////////////////////////////////////////////
	// vectors of indices
	// indices of Hbond acceptor positions

	//indices of atoms with orbitals
	AtomIndices atoms_with_orb_index_;

	//indices of haro hydrogens
	AtomIndices Haro_index_;

	//indices of hpolar hydrogens
	AtomIndices Hpol_index_;
	AtomIndices accpt_pos_;

	// indices of polar Hydrogens for Hbond donors
	AtomIndices Hpos_polar_;

	// indices of apolar hydrogens
	AtomIndices Hpos_apolar_;

	// indices of Hbond acceptor positions that are part of the sidechain
	// must be a subset of the atoms listed in the accpt_pos_ array
	AtomIndices accpt_pos_sc_;

	// indices of polar Hydrogens for Hbond donors that are part of the sidechain
	// must be a subset of the atoms listed in the Hpos_polar_ array
	AtomIndices Hpos_polar_sc_;

	// Indices of all backbone atoms, hydrogens and heavyatoms
	AtomIndices all_bb_atoms_;

	// Indices of all sidechain atoms, hydrogens and heavyatoms
	AtomIndices all_sc_atoms_;

	/// @brief Names of all of the atoms that are able to make a bond to a metal, for metal-binding residue types
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1 < std::string > metal_binding_atoms_;


	/// @brief indices of all mainchain atoms
	/// @details mainchain_atoms are those atoms on a path from polymer lower_connect
	/// to upper_connect. For protein, this will be N, CA and C.
	utility::vector1<VD> mainchain_atoms_;

	/// @brief indices of action coordinate centers
	/// @details the geometric center of the atoms listed defined the residue's "action coordinate"
	utility::vector1< VD > actcoord_atoms_;


	///////////////////////////////////////////////////////////////////////////
	// vectors of vectors of indices
	// indices of four atoms to build each chi angle

	utility::vector1<utility::vector1<VD> > chi_atoms_;
	utility::vector1< bool        > is_proton_chi_;
	utility::vector1< Size        > proton_chis_;
	utility::vector1< Size        > chi_2_proton_chi_;
	utility::vector1< utility::vector1< Real > > proton_chi_samples_;
	utility::vector1< utility::vector1< Real > > proton_chi_extra_samples_;

	// indices of four atoms to build each nu angle
	utility::vector1<utility::vector1<VD> > nu_atoms_;

	// number of bonds separated between any pair of atoms in this residue
	utility::vector1< utility::vector1< int > > path_distance_;

	// atom index lookup by atom name string
	std::map< std::string, VD > atom_name_to_vd_; //atom_graph_index_;

	// Legacy/backward compatibility device holds an ordered list of nodes' indices
	VDs ordered_atoms_; // Position in the vector represents Atom in "ordered arrangement"

	// atom index lookup by atom name string
	std::map< std::string, int > orbitals_index_;

	// Additional non-Dunbrack rotamer bins
	/*
			pair<Real,Real>  ==>  mean,sdev
			for each chi angle i and rotamer j: chi_rotamers_[i][j]
	 */
	utility::vector1< utility::vector1< std::pair< Real, Real > > > chi_rotamers_;

	std::string rotamer_library_name_;

	// NCAA rotlib stuff some of this is hardcoded elsewhere for the CAAs
	// whether or not we should use the NCAA rotlib if it exists
	bool use_ncaa_rotlib_;

	// path to the NCAA rotlib
	std::string ncaa_rotlib_path_;

	// the number of non-hydrogen chi angles in the NCAA rotlib
	Size ncaa_rotlib_n_rots_;

	// the number of rotamer bins for each chi angle in the NCAA rotlib
	utility::vector1< Size > ncaa_rotlib_n_bins_per_rot_;

	/////////////////////////////////////
	// properties -- some of these may be deducible from AA?
	//
	// if you add new things -- do they need to be initialized in the c-tor?
	// gcc debug does not seem to initialize bools to false!

	// residue properties as defined in the residue param files
	utility::vector1< std::string > properties_;
	bool is_polymer_;
	bool is_protein_;
	bool is_alpha_aa_;
	bool is_beta_aa_;
	bool is_l_aa_;
	bool is_d_aa_;
	bool is_charged_;
	bool is_polar_;
	bool has_sc_orbitals_;
	bool is_aromatic_;
	bool is_cyclic_;
	bool is_DNA_;
	bool is_RNA_;
	bool is_NA_;
	bool is_carbohydrate_;
	bool is_lipid_;
	bool is_ligand_;
	bool is_metal_; //Is this residue type a metal ion?
	bool is_metalbinding_; //Is this residue type a type capable of binding to a metal ion?
	bool is_surface_;
	bool is_terminus_; // last or first residue in a chain; set to TRUE during terminus patching
	bool is_lower_terminus_; // first residue in a chain; set to TRUE during terminus patching
	bool is_upper_terminus_; // last residue in a chain; set to TRUE during terminus patching
	bool is_branch_lower_terminus_;
	bool is_phosphonate_; // amino phosphonic acid instead of amino carboxylic acid
	bool is_phosphonate_upper_;
	bool is_acetylated_nterminus_;
	bool is_methylated_cterminus_;
	bool is_coarse_; //currently for coarse_RNA only
	bool is_adduct_;
	bool is_virtual_residue_;
	// etc., etc.

	// here we store the patch operations/variant types that describe this residue
	utility::vector1< VariantType > variant_types_;

	// Here we store arbitrary numeric properties with string names
	std::map<std::string,core::Real> numeric_properties_;

	// Here we store arbitrary string properties with string names
	std::map<std::string,std::string> string_properties_;

	//////////////////////////////////////////////////
	// features

	// standard rosetta aa-type for knowledge-based potentials, may be aa_unk
	// aa_ = THIS residue's aa-type; rotamer_aa_ = the aa-type on which rotamers will be based; backbone_aa_ = the aa-type on which the backbone scoring (rama, p_aa_pp) will be based.
	AA aa_, rotamer_aa_, backbone_aa_;

	// unique residue type id
	std::string name_;

	// pdb-file id, need not be unique
	std::string name3_;

	// one-letter code, also not necessarily unique
	char name1_;

	// interchangeability group lets a ResidueType claim to be functionally
	// interchangeable with any other ResidueType in the same group.  This
	// is used by the packer to decide which ResidueType from a desired group
	// has the right set of variants to be placed at a particular position.
	// E.g. if the interchangeability group is "ALA" and the packer is building
	// rotamers for residue 1, (the N-terminal residue) then, the packer will
	// select the "ALA:NTermProteinFull" ResidueType and build rotamers for it.
	std::string interchangeability_group_;

	// for rsd-rsd neighbor calculations
	// atom used for calculating residue-level neighbors
	VD nbr_atom_;

	// radius cutoff to define neighors
	Real nbr_radius_;

	// Controls which atoms are selected by "select_orient_atoms",
	// used to overlay residues during packing.
	bool force_nbr_atom_orient_;

	//Real mass_;
	Real mass_;

	// number of actcoord atoms
	/*
			the geometric center of the atoms listed defined the residue's "action coordinate"
	*/
	Size n_actcoord_atoms_;

	// the unprocessed metadata
	sdf::MolData mol_data_;

	/// @brief  Vector of inter-residue connections expected for this residuetype
	/// NOW includes the polymer connections, as well as disulf-type connections
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	/// @note  Much of the code assumes at most one residue connection per atom.
	///        The pseudobond code has been partially written to handle cases where multiple connections
	///        to a single atom can exist, but much of the standard residue connection code assumes a
	///        simple correspondence between atoms and residue connections.  That code will have to be
	///        updated to support single-atom "backbones."
	utility::vector1< ResidueConnection > residue_connections_;
	utility::vector1< utility::vector1< Size > > atom_2_residue_connection_map_;

	/// @brief For calculating inter-residue bond angle and bond torsion energies, it is useful to have a list of those
	/// atoms within one bond of a residue connection atom.  For residue connection i,
	/// its position in this array is a list of pairs of atom-id's, the first of which is always the id
	/// for the atom forming residue connection i.
	utility::vector1< utility::vector1< two_atom_set > > atoms_within_one_bond_of_a_residue_connection_;

	/// @brief For atom i, its position in this vector is a list of pairs where
	/// first == the residue_connection id that lists this atom as being within one bond
	/// of a residue connection, and
	/// second == the index of the entry containing this atom in the
	/// atoms_within_one_bond_of_a_residue_connection_[ first ] array.
	utility::vector1< utility::vector1< std::pair< Size, Size > > > within1bonds_sets_for_atom_;

	/// @brief For calculating inter-residue bond torsion energies, it is useful to have a list of those
	/// atoms within two bonds of a residue connection atom.  For residue connection i,
	/// its position in this array is a list of triples of atom-id's, the first of which is always the id for
	/// the atom forming residue connection i.
	utility::vector1< utility::vector1< three_atom_set > > atoms_within_two_bonds_of_a_residue_connection_;

	/// @brief For atom i, its position in this vector is a list of pairs where
	/// first == the residue_connection id that lists this atom as being within two bonds
	/// of a residue connection, and
	/// second == the index of the entry containing this atom in the
	/// atoms_within_two_bonds_of_a_residue_connection_[ first ] array.
	utility::vector1< utility::vector1< std::pair< Size, Size > > > within2bonds_sets_for_atom_;


	/// @brief  Polymer lower connections
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	// ResidueConnection lower_connect_; // deprecated
	Size lower_connect_id_; // which connection is the lower connection?

	/// @brief  Polymer upper connections
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	// ResidueConnection upper_connect_; // deprecated
	Size upper_connect_id_; // which connection is the upper connection?

	Size n_non_polymeric_residue_connections_;
	Size n_polymeric_residue_connections_;

	///////////////////////////////////////////////////////////////////////
	// These arrays are temporary, will be cleared in finalize():
	/// a list of atom indices to be deleted
	/** in the next call to finalize(), used in delete_atom which is called during patching**/
	/// atom indices forced to be considered backbone
	utility::vector1<VD> force_bb_;

	////////////////
	core::chemical::rna::RNA_ResidueTypeOP rna_residue_type_;

	// A container for residue properties unique to carbohydrates.
	core::chemical::carbohydrates::CarbohydrateInfoOP carbohydrate_info_;


	//All the rings and the edges. Defaulted to null until defined
	std::list<utility::vector1< ED > > rings_and_their_edges_;

	//ALL THE Indexed data
	utility::vector1<Size> atom_base_indices_;
	utility::vector1<Size> abase2_indices_;
	utility::vector1< AtomIndices > chi_atoms_indices_;
	utility::vector1<AtomIndices> nu_atoms_indices_;
	AtomIndices mainchain_atoms_indices_;
	Size nbr_atom_indices_;
	AtomIndices actcoord_atoms_indices_;
	utility::vector1< AtomIndices > cut_bond_neighbor_indices_;
	utility::vector1< Size > atom_shadowed_indices_;

	////////////////
	/// status
	bool finalized_;

	// Adducts defined for this residue
	utility::vector1< Adduct > defined_adducts_;


	// boost serialize stuff
#ifdef USEBOOSTSERIALIZE
	template<class Archive> friend void core::conformation::save_construct_data( Archive & ar, const core::conformation::Residue * t, const unsigned int file_version);
	template<class Archive> friend void core::conformation::load_construct_data( Archive & ar, core::conformation::Residue * t, const unsigned int file_version);
#endif

	bool nondefault_;
	std::string base_restype_name_;
public:
	void nondefault(bool in) { nondefault_ = in;}
	void base_restype_name(std::string const & in) { base_restype_name_ = in;}
	std::string base_restype_name() const {return base_restype_name_;}

	// this is a total hack, I'm tired
	mutable bool serialized_;
	// end hack?
};  // class ResidueType

// Insertion operator (overloaded so that ResidueType can be "printed" in PyRosetta).
std::ostream & operator<<(std::ostream & output, ResidueType const & object_to_output);

} // chemical
} // core

#endif // INCLUDED_core_chemical_Residues_HH
