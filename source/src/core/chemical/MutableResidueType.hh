// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MutableResidueType.hh
/// @brief Method declarations and simple accessors/getters for MutableResidueType
/// @author
/// Phil Bradley
/// Rocco Moretti (rmorettiase@gmail.com)
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to properties, rings, lipids, carbohydrates, and other non-AAs)


#ifndef INCLUDED_core_chemical_MutableResidueType_hh
#define INCLUDED_core_chemical_MutableResidueType_hh

// Unit headers
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/ResidueTypeBase.hh>
#include <core/chemical/ResidueType.fwd.hh> // For conversion constructor only

#include <core/chemical/ResidueGraphTypes.hh>

// Package headers
#include <core/chemical/MutableChiRecord.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.fwd.hh>
#include <core/chemical/rings/RingConformerSet.fwd.hh>
#include <core/chemical/rings/RingSaturationType.hh>

#ifdef WIN32
#include <core/chemical/MutableResidueConnection.hh>
#else
#include <core/chemical/MutableResidueConnection.fwd.hh>
#endif

// Project headers
#include <core/types.hh>

// Utility headers

// C++ headers

// External headers

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

/// @brief A class for defining a type of residue, modifiable version
/// @details This class contains the "chemical" information for residues.
/// A MutableResidueType in Rosetta can be a ligand, DNA, amino acid, or basically anything.
/// MutableResidueTypes are normally generated through .params files, which are read from the database chemical/residue_types.
/// However, there are several other ways of generating this class, and the MutableResidueTypes can even be modified during the run (hence their name).
///
/// The MutableResidueType differs from a plain ResidueType in that it's constructed to be modifiable.
/// It's also not intented to be used itself (aside from ResidueType modification).
/// Instead, the typical usage is to convert a MutableResidueType into a plain ResidueType.
///
/// Another MutableResidueType/ResidueType distinction is that in a MutableResidueType
/// the atom information is encoded in a molecular graph which should make it easier to add/remove/modify atoms.
/// These atoms are referred to primarily by their "vertex desciptor" in the graph, which should be invariant to insertion/deletion/etc.
/// While MutableResidueType does has a minimal sense of atom indexing, this is not the primary way to identify atoms,
/// and will not be robust to addition/deletion of atoms.
///
/// To actually use a MutableResidueType for a simulation, you need to create a plain ResidueType from it.
/// At this point all of the derived data will be calculated and converted to a more efficient struct-of-arrays format.

class MutableResidueType : public ResidueTypeBase
{
public:

private:
	MutableResidueType(); // private, not deleted, as serialization needs it.

public:

	/// @brief constructor
	/// @details We use the AtomTypeSet object to assign atom_types to atoms inside add_atom,
	/// and to identify (polar) hydrogens, acceptors, etc.
	MutableResidueType(
		AtomTypeSetCOP atom_types,
		ElementSetCOP element_types,
		MMAtomTypeSetCOP mm_atom_types,
		orbitals::OrbitalTypeSetCOP orbital_types//,
	);

	/// @brief Make a copy of a MutableResidueType.
	MutableResidueType( MutableResidueType const & residue_type );

	/// @brief Make a modifiable residue type from a standard one.
	explicit MutableResidueType( ResidueType const & residue_type );

	/// @brief make a copy
	MutableResidueTypeOP clone() const;

	/// @brief make a copy
	MutableResidueTypeOP placeholder_clone() const;

	/// @brief Move constructor, for rule of 5
	MutableResidueType( MutableResidueType && residue_type ) = default;

	/// @brief destructor, for rule of 5
	~MutableResidueType() override;

	/// @brief Copies  <src>  into the MutableResidueType
	MutableResidueType &
	operator=( MutableResidueType const & src );

	/// @brief Move assignment, for rule of 5
	MutableResidueType &
	operator=( MutableResidueType && ) = default;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Summary Functions               ////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief number of atoms
	/// Note: Don't use this to iterate over atoms -- use the vector returned by all_atoms() instead.
	Size
	natoms() const override
	{
		return graph_.num_vertices();
	}

	/// @brief number of heavy atoms
	Size
	nheavyatoms() const override;

	/// @brief number of bonds
	Size
	nbonds() const override
	{
		return graph_.num_edges();
	}

	/// @brief number of bonds for given atom
	Size
	nbonds( VD atom ) const
	{
		return boost::out_degree( atom , graph_ );
	}

	/// @brief Counts the number of virtual atoms and returns the count.
	/// @details The virtual count is not stored in the resiude type.  This count is performed on the fly, and
	///can hurt performance if reapeatedly carried out.  Not intended for use in large loops -- instead, call
	///once and store the value.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	Size
	n_virtual_atoms() const override;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Graph Functions             ////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Constant access to the underlying graph.
	ResidueGraph const &
	graph() const { return graph_; }

	VIterPair atom_iterators() const {
		return boost::vertices(graph_);
	}

	EIterPair bond_iterators() const {
		return boost::edges(graph_);
	}

	OutEdgeIterPair bond_iterators( VD const & atom ) const {
		return boost::out_edges(atom, graph_);
	}

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Atom Functions              ////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Convenience function for iterating over all the atoms
	/// e.g. `for ( VA atm: restype->all_atoms() ) {`
	utility::vector1< VD > const &
	all_atoms() const {
		return ordered_atoms_;
	}

	//Atom & atom(Size const atom_index);
	//Atom const & atom(Size const atom_index) const;
	Atom & atom(VD const atom_vd);
	Atom const & atom(VD const atom_vd) const;
	Atom & atom(std::string const & atom_name);
	Atom const & atom(std::string const & atom_name) const;

	using ResidueTypeBase::has;

	/// @brief is this atom present in this residue?
	bool
	has( std::string const & atom_name ) const override
	{
		return atom_name_to_vd_.find(atom_name) != atom_name_to_vd_.end();
	}

	/// @brief is this vertex descriptor present in this residue?
	bool
	has( VD const vd ) const
	{
		return core::chemical::has(graph_, vd);
	}

	/// @brief get atom name by vertex descriptor
	/// IMPORTANT NOTE: Atoms don't necessarily have names!
	std::string const &
	atom_name( VD const vd ) const;

	/// @brief get the vertex descriptor from the name of the atom.
	VD
	atom_vertex( std::string const & name ) const;

	/// @brief Get an integer value for the atom (the index for all_atoms())
	/// @details This is provided purely for a convenience in matching up external vector entries
	/// Use all_atoms() directly for normal iteration over atoms
	/// Will return 0 if the vd isn't in the MutableResidueType.
	Size
	atom_index( VD const & vd ) const;

	/// @brief Get the chemical atom_type for this atom
	AtomType const &
	atom_type( VD const vd ) const;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Bond Functions              ////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	Bond & bond(ED const ed);
	Bond const & bond(ED const ed) const;
	Bond & bond( VD vd1, VD vd2);
	Bond const & bond( VD vd1, VD vd2) const;
	Bond & bond(std::string const & atom1, std::string const & atom2);
	Bond const & bond(std::string const & atom1, std::string const & atom2) const;

	bool
	atoms_are_bonded( std::string const & atom1, std::string const & atom2 ) const;

	AdjacentIterPair
	bonded_neighbor_iterators( VD const & atom ) const;

	/// @brief Get the atoms bonded to the specified atom.
	/// @details This is recalculated each time the function is called.
	utility::vector1< VD >
	bonded_neighbors( VD const & atom ) const;

	///@brief indicates how many proton bonded neighbors an atom has
	Size
	number_bonded_hydrogens( VD atomvd ) const;

	/// @brief Get the atoms bonded to the specified atom, if they're annotated as hydrogens.
	/// @details This is recalculated each time the function is called.
	utility::vector1< VD >
	bonded_hydrogens( VD const & atom ) const;

	///@brief indicates how many heavyatom bonded neighbors an atom has, graph version
	Size
	number_bonded_heavyatoms( VD atomvd ) const;

	/// @brief Get the atoms bonded to the specified atom, if they're annotated as heavyatoms.
	/// @details This is recalculated each time the function is called.
	utility::vector1< VD >
	bonded_heavyatoms( VD const & atom ) const;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Chi Functions               ////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief number of chi angles
	Size
	nchi() const
	{
		return chis_.size();
	}

	/// @brief number of proton chis
	Size
	n_proton_chi() const;

	/// @brief number of proton chis
	bool
	is_proton_chi( Size const chino ) const
	{
		debug_assert( chis_[chino] != nullptr );
		return chis_[chino]->is_proton_chi();
	}

	/// @brief Will return if the chi is currently valid.
	bool
	chi_valid( Size const chino ) const;

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

	/// @brief VDs of the atoms which are used to define a given chi angle (chino)
	VDs const &
	chi_atom_vds( Size const chino ) const
	{
		debug_assert( chis_[ chino ] != nullptr );
		return chis_[ chino ]->chi_atoms();
	}

	/// @brief The proton chi samples
	/// Important - this takes the full chi numbering (not proton chi numbering)
	utility::vector1< Real > const &
	proton_chi_samples_for_chi( Size chino ) const
	{
		debug_assert( chis_[ chino ] != nullptr );
		return chis_[chino]->proton_chi_samples();
	}

	/// @brief The proton chi extra samples
	/// Important - this takes the full chi numbering (not proton chi numbering)
	utility::vector1< Real > const &
	proton_chi_extra_samples_for_chi( Size chino ) const
	{
		debug_assert( chis_[ chino ] != nullptr );
		return chis_[chino]->proton_chi_extra_samples();
	}

	/// @brief all rotamers bins (mean, std) for a given chi angle
	utility::vector1< std::pair< Real, Real > > const &
	chi_rotamers( Size const chino ) const
	{
		debug_assert( chis_[ chino ] != nullptr );
		return chis_[ chino ]->chi_rotamers();
	}

	utility::vector1< VD > const &
	nu_atoms( core::uint const nu_index ) const {
		return nu_atoms_[ nu_index ];
	}

	utility::vector1< VD > const &
	ring_atoms( core::uint const ring_num ) const {
		return ring_atoms_[ ring_num ];
	}

	/// @brief  Return the saturation level of this residue's nth cycle.
	core::chemical::rings::RingSaturationType
	ring_saturation_type( uint const ring_num ) const
	{
		return ring_saturation_types_[ ring_num ];
	}

	/// @brief  Low-energy ring conformers for the given ring
	utility::vector1< std::string > const &
	low_ring_conformers( core::uint const ring_num ) const {
		return low_ring_conformers_[ ring_num ];
	}

	/// @brief   Lowest-energy ring conformer for the given ring
	std::string const &
	lowest_ring_conformer( core::uint const ring_num ) const {
		return lowest_ring_conformer_[ ring_num ];
	}

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Connection Functions        ////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief number of ResidueConnections, counting polymeric residue connections
	Size n_possible_residue_connections() const;

	/// @brief Get a ResidueConection.
	MutableResidueConnection const & residue_connection( Size const i ) const;

	/// @brief Get a ResidueConection.
	MutableResidueConnection & residue_connection( Size const i );

	VD
	residue_connect_atom( Size const resconn_id ) const;

	///// @brief How many inter-residue chemical bonds does a particular atom form?
	Size
	n_residue_connections_for_atom( VD const atomid ) const;

	/// @brief  add a *non-polymeric* ResidueConnection
	/// @details For polymeric connections, see set_lower_connect() and set_upper_connect()
	/// Doesn't set the ideal geometry -- maybe it should?
	Size
	add_residue_connection( std::string const & atom_name );

	/// @brief Remove all connections (both polymeric and otherwise) for an atom
	void
	delete_residue_connection( VD atm );

	//////////////
	// Polymeric

	/// @brief The number of polymeric residue connections.
	Size
	n_polymeric_residue_connections() const {
		// Take advantage of `true = 1` to sum the number of non-zero connections
		return ( lower_connect_id_ != 0 ) + ( upper_connect_id_ != 0 );
	}

	bool
	residue_connection_is_polymeric( Size const resconn_id ) const {
		return ( resconn_id == lower_connect_id_ || resconn_id == upper_connect_id_ );
	}

	// Lower
	MutableResidueConnection const & lower_connect() const;

	Size
	lower_connect_id() const
	{
		return lower_connect_id_;
	}

	/// @brief index number of the atom which connects to the lower connection
	VD lower_connect_atom() const;

	/// @brief set the atom which connects to the lower connection
	void set_lower_connect_atom( std::string const & atm_name );

	// Upper
	MutableResidueConnection const & upper_connect() const;

	Size
	upper_connect_id() const
	{
		return upper_connect_id_;
	}

	/// @brief index number of the atom which connects to the upper connection
	VD upper_connect_atom() const;

	/// @brief set the atom which connects to the upper connection
	void set_upper_connect_atom( std::string const & atm_name );

	//////////////
	// other

	void
	add_metapatch_connect( std::string const & atom );

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////    Other Residue Functions     //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	//// nbr_atom and nbr_radius are used for rsd-rsd neighbor calculation

	/// @brief set nbr_atom used to define residue-level neighbors
	void
	nbr_atom( std::string const & atom_name ) {
		nbr_atom_ = atom_vertex( atom_name );
	}

	/// @brief set nbr_atom used to define residue-level neighbors
	void
	nbr_atom( VD vertex )
	{
		nbr_atom_ = vertex;
	}

	/// @brief get VD  used to define residue-level neighbors
	VD
	nbr_vertex() const
	{
		return nbr_atom_;
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

	/// @brief set an atom as backbone heavy atom
	/// @details When converted to a plain ResidueType, backbone atoms will be placed first.
	/// This is differenct from the mainchain atoms below, as the mainchain is a subset of the backbone
	void
	set_backbone_heavyatom( std::string const & name );

	bool
	is_backbone_heavyatom( VD atom ) const;

	/// @brief indices of all mainchain atoms
	utility::vector1< VD > const &
	mainchain_atoms() const { return mainchain_atoms_; }

	/// @brief set indices of all mainchain atoms
	/// TODO Should we make sure that mainchain atoms are also backbones?
	void set_mainchain_atoms( utility::vector1< VD > const & mainchain ) {
		mainchain_atoms_ = mainchain;
	}

	/// @brief get descriptors for atoms used to define actcoord
	utility::vector1< VD >
	actcoord_atoms() const;

	/// @brief add an atom to the list for calculating actcoord center
	void add_actcoord_atom( std::string const & atom );

	/// @brief Remove an atom from the list of act coord atoms
	/// (used in patching when it kills the valence that is thus used)
	/// @author Andrew Watkins (amw579@nyu.edu)
	void
	delete_actcoord_atom( std::string const & atom_name );

	/// @brief Add an alias name for an atom.
	void
	add_atom_alias( std::string const & rosetta_atom, std::string const & alias );

	/// @brief store canonical to alias mapping
	void
	add_canonical_atom_alias( std::string const & rosetta_atom, std::string const & alias );

	/// @brief Remove a given alias name for an atom.
	/// @details If error is true, raise error if the alias can't be found
	void
	delete_atom_alias( std::string const & alias, bool error=true );

	void
	set_shadowing_atom(
		std::string const & atom,
		std::string const & atom_being_shadowed
	);

	std::map< VD, VD > const &
	shadow_atoms() const {
		return atom_shadowed_;
	}

	/// @brief Add an atom to the list of atoms that can potentially form a bond to a metal ion.
	/// Note that the atom must exist in the residue type (the function checks for this at runtime).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void
	add_metalbinding_atom( std::string const & atom_name );

	/// @brief Remove an atom from the list of atoms that can potentially form a bond to a metal ion
	/// (used in patching when it kills the valence that is thus used)
	/// @author Andrew Watkins (amw579@nyu.edu)
	void
	delete_metalbinding_atom( std::string const & atom_name );

	/// @brief Is this ResidueType a base type?
	bool is_base_type() const override;

	/// @brief Get a pointer to this ResidueType's base ResidueType.
	/// @details Returns NULL if base_type_cop_ is NULL. (Difference from ResidueType behavior!)
	ResidueTypeCOP get_base_type_cop() const override;

	/// @brief Reset the base type COP to be null. This implies that this ResidueTypeBase is a base type.
	void reset_base_type_cop();

	/// @brief Set the base type COP.  This implies that this ResidueTypeBase is NOT a base type.
	void set_base_type_cop( ResidueTypeCOP new_base_type );

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////  Modification Functions        //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////
	/////////////////////////atoms////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief add an atom into this residue
	/// Will return the vertex descriptor of the added atom.
	VD
	add_atom(
		std::string const & atom_name,
		std::string const & atom_type_name,
		std::string const & mm_atom_type_name,
		Real const charge
	);

	/// @brief add an atom into this residue, with just the name.
	/// Will return the vertex descriptor of the added atom.
	VD
	add_atom( std::string const & atom_name = "" );

	VD
	add_atom(Atom const & atom, MutableICoorRecord const & icoor);

	/// @brief Remove the atom from the ResidueType
	void
	delete_atom( std::string const & name );

	/// @brief Remove the atom from the ResidueType
	void
	delete_atom( VD atom );

	/// @brief Rename the atom, updating the ResidueType-internal data mapping
	void
	rename_atom( VD atom, std::string const & name );

	/// @brief set atom type
	void
	set_atom_type(
		std::string const & name,
		std::string const & atom_type_name
	) { set_atom_type( atom_vertex( name ), atom_type_name); }

	/// @details Set atom type, correctly updating internal state.
	/// @note    This method also sets/updates the atom's element and the residue's mass.
	void
	set_atom_type( VD atom, std::string const & atom_type_name );

	/// @brief set gasteiger atom type
	void
	set_gasteiger_atom_type(
		std::string const & atom_name,
		std::string const & gasteiger_atom_type_name
	);

	/// @brief set gasteiger atom type
	void
	set_gasteiger_atom_type(
		VD atom,
		std::string const & gasteiger_atom_type_name
	);

	void
	delete_child_proton( std::string const & atom );

	//////////////////////////////////////////////////////////////////////
	/////////////////////////bonds////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief add a bond between atom1 and atom2, specifying a bond type (SingleBond, DoubleBond, TripleBond, AromaticBond)
	void add_bond(std::string const & atom_name1, std::string const & atom_name2, BondName bondLabel = SingleBond);

	/// @brief add a bond between atom1 and atom2, specifying a bond type (SingleBond, DoubleBond, TripleBond, AromaticBond)
	void add_bond(VD atom1, VD atom2, BondName bondLabel = SingleBond);

	/// @brief Delete a bond between the two atoms.
	/// @details Note that this might leave dangling atoms.
	void
	delete_bond(VD atom1, VD atom2);

	/// @brief  Change the bond type of the given bond from one type to another.
	/// @details Acts like a bond deletion + bond add, rather than a change.
	void change_bond_type(
		std::string const & atom_name1,
		std::string const & atom_name2,
		BondName const new_bond_label );

	/// @brief add a bond between atom1 and atom2, if bond type is not specified, default to a SingleBond
	void
	add_cut_bond(
		std::string const & atom_name1,
		std::string const & atom_name2
	);

	/// @brief Reset the bond distance to an atom whose internal coordinates have already been set.
	void reset_bond_distance_to_atom( std::string const & atm, core::Distance const d );

	//////////////////////////////////////////////////////////////////////
	/////////////////////////geometry/////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief get root_atom used as the base of the icoor tree.
	VD
	root_atom() const { return root_atom_; }

	/// @brief MutableICoorRecord of an atom
	MutableICoorRecordCOP
	icoor( VD const atm ) const;

	/// @brief set MutableICoorRecord for an atom
	/// @details phi and theta are in radians
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

	/// @brief set MutableICoorRecord for an atom, vertex descriptor version
	/// @details phi and theta are in radians
	void
	set_icoor(
		VD const & atm,
		Real const phi,
		Real const theta,
		Real const d,
		VD const & stub_atom1,
		VD const & stub_atom2,
		VD const & stub_atom3,
		bool const update_xyz = false
	);

	/// @brief set MutableICoorRecord for an atom, vertex descriptor version
	/// @details phi and theta are in radians
	void
	set_icoor(
		VD const & atm,
		Real const phi,
		Real const theta,
		Real const d,
		std::string const & stub_atom1,
		std::string const & stub_atom2,
		std::string const & stub_atom3,
		bool const update_xyz = false
	);

	// @brief Reset all the icoord records
	void
	clear_icoor();

	/// @brief The atom base is the distance atom the atom is bonded to.
	/// If the specified atom doesn't have valid ICOOR record yet, this function will return
	/// a null vertex
	VD
	atom_base( VD const atm ) const;

	void
	set_ideal_xyz(
		std::string const & atm,
		Vector const & xyz_in
	);

	void
	set_ideal_xyz(
		VD atm,
		Vector const & xyz_in
	);

	//////////////////////////////////////////////////////////////////////
	///////////////////////// chis   /////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Add a chi (side-chain) angle defined by four atoms.
	void
	add_chi(Size const chino,
		VD atom1,
		VD atom2,
		VD atom3,
		VD atom4);

	/// @brief Add a chi (side-chain) angle defined by four atoms.
	void
	add_chi(VD atom1,
		VD atom2,
		VD atom3,
		VD atom4);


	/// @brief Add a chi (side-chain) angle defined by four atoms.
	void
	add_chi(Size const chino,
		std::string const & atom_name1,
		std::string const & atom_name2,
		std::string const & atom_name3,
		std::string const & atom_name4);

	/// @brief Add a chi (side-chain) angle defined by four atoms to the end of the list of chis.
	void
	add_chi(std::string const & atom_name1,
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

	/// @brief Annotate a given chi as a proton chi, and set the sampling behavior
	/// If the chi is already listed as a proton chi, change the sampling behavior
	void
	set_proton_chi(
		Size chino,
		utility::vector1< Real > const & dihedral_samples,
		utility::vector1< Real > const & extra_samples
	);

	/// @brief Add a rotamer bin for a given chi.
	void
	add_chi_rotamer(
		Size const chino,
		Real const mean,
		Real const sdev
	);

	/// @brief Adds a chi rotamer bin to the highest-indexed chi in the list of chis for this MutableResidueType.
	void add_chi_rotamer_to_last_chi(core::Angle const mean, core::Angle const sdev);

	/// @brief Delete all of the chi rotamer bins from the specified chi for this MutableResidueType.
	void clear_chi_rotamers( core::uint const chi_no );

	/// @brief delete terminal chi angle
	void
	delete_terminal_chi();

	/// @brief Deletes a given chi, potentially renumbering the remaining chis
	void delete_chi( Size const chino );

	/// @brief Add a nu (internal cyclic) angle defined by four atoms.
	void add_nu(
		core::uint const nu_index,
		std::string const & atom_name1,
		std::string const & atom_name2,
		std::string const & atom_name3,
		std::string const & atom_name4);

	/// @brief Deletes a given nu, potentially renumbering the remaining nus
	void delete_nu( core::uint const nu_index );

	/// @brief  Add a ring definition.
	void add_ring( core::uint const ring_num,
		utility::vector1< std::string > const & ring_atoms,
		core::chemical::rings::RingSaturationType const saturation_type=core::chemical::rings::ALIPHATIC );

	/// @brief  Set this cyclic residue's lowest-energy ring conformer for the nth ring by IUPAC name.
	void set_lowest_energy_ring_conformer( core::uint const ring_num, std::string const & conformer );

	/// @brief  Set this cyclic residue's low-energy ring conformers for the nth ring by IUPAC name.
	void set_low_energy_ring_conformers( core::uint const ring_num, utility::vector1< std::string > const & conformers );

	/// @brief Deletes a given ring, potentially renumbering the remaining ring
	void delete_ring( core::uint const ring_num );

	//////////////////////////////////////////////////////////////////////
	/////////////////////////orbitals/////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief clear orbitals
	void
	clear_orbitals();

	/// @brief add an orbital onto a residue based upon atom
	void
	add_orbital(
		std::string & orbital_name,
		std::string & orbital_type_name
	);

	/// @brief add an orbital bond between an atom and an orbital.
	/// @note NOTE!!!!! This is indexed based upon atoms, not orbitals. That means that in your params file
	/// you must have the atom as the first and orbital as the second.
	void
	add_orbital_bond(
		std::string const & atom_name1,
		std::string const & orbital_name
	);

	/// @brief set OrbitalICoor for an orbital
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
	//////////////////////////////////////////////////////////////////////
	////////////////  Show Functions                //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief  Generate string representation of MutableResidueType for debugging purposes.
	void show( std::ostream & output=std::cout, bool output_atomic_details=false ) const override;

	void
	dump_vd_info() const;

	void
	show_all_atom_names( std::ostream & out ) const override;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////// Large-scale recalculation functions    //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	// RM: These should probably be moved out of the MutableResidueType itself,
	// into utility functions.

	/// @brief Figure out what the "center" atom of the residue is from current coordinates
	/// @details Assumes that all the ideal_xyz values have been set
	void assign_neighbor_atom();

	/// @brief Assign internal coordinates from the set ideal xyz coordinates.
	/// Note that it currently does not obey mainchain designations or cut bonds.
	void assign_internal_coordinates();

	/// @brief Function to assign internal coordinates from bonding patterns
	/// Note that it currently does not work well with polymers.
	void assign_internal_coordinates(core::chemical::VD new_root );

	/// @brief Assign ideal_xyz coordinates from icoor records.
	void
	fill_ideal_xyz_from_icoor();

	/// @brief Regenerate the rotatable chi bonds from the internal graph structure.
	/// If the number of proton chi samples would exceed max_proton_chi_samples, don't add extra sampling to proton chis.
	/// As a special case, if this is zero don't add any proton chi sampling at all.
	///
	/// Requires that Icoor and atom base records are up-to-date, and that ring bonds have been annotated.
	void
	autodetermine_chi_bonds( core::Size max_proton_chi_samples = 500 );

public:

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////  Special Functions             //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Run checks on this MutableResidueType, checking to make sure it's okay.
	/// Return true on success and false on failure
	/// (Intended usage is within something like a debug_assert())
	bool
	validate_residue_type() const;

	/// @brief Get access to the atom_name-to-vd map.
	/// Provided only for MutableResidueType->ResidueType conversion
	std::map< std::string, VD > const &
	get_atom_name_to_vd_map() const { return atom_name_to_vd_; }

	/// @brief Change which atom type set this MutableResidueType points to
	/// WARNING - While this tries to switch over the atom types, it will null out
	/// any that don't have exact name correspondences.
	/// You NEED to go through and manually reset the types on all atoms.
	/// Exposed for black-magic use only.
	void
	update_atom_type_set( AtomTypeSetCOP setting );

private:

	/// @brief After deleting a connection, update the icoors of any atoms dependent on HIGHER-numbered
	/// connections.
	/// @details Will not attempt to adjust any ICOOR referencing a delete UPPER/LOWER
	void update_icoors_after_connection_deletion( core::Size const conn_id_deleted );

	/// @brief Copy over the atom information from the ResidueType.
	/// Used for the conversion constructor
	/// @details The ordered_atoms_ vector matches the ResidueType order, so we don't need an old-to-new.
	void
	copy_atom_info( ResidueType const & source );

	/// @brief Fill in various other information which is being copied over from the ResidueType
	/// Used for the conversion constructor
	/// @details The ordered_atoms_ vector matches the ResidueType order, so we don't need an old-to-new.
	void
	copy_other_info( ResidueType const & source );

private:

	/// @brief The Atoms and Bonds of the MutableResidueType, stored as Nodes and Edges.
	ResidueGraph graph_;

	/// @brief atom index lookup by atom name string
	/// @details Note that not all atoms will necessarily have names or unique names.
	/// This map will merely map those which do,
	std::map< std::string, VD > atom_name_to_vd_;

	/// @brief The atoms in order.
	/// This is primarily to keep track of the order which that atoms have been added, for consistency.
	/// Unlike the main ResidueType class, there is *no* guarantees regarding property ordering of atoms.
	VDs ordered_atoms_;

	//////////////////////////////////////////////////////////////////////
	// General properties

	/// @brief Const-access owning pointer to the base residue type.
	/// @details NULL if this MutableResidueType doesn't have a base type.
	ResidueTypeCOP base_type_cop_;

	//////////////////////////////////////////////////////////////////////
	// Atom annotations

	VD root_atom_ = MutableResidueType::null_vertex;

	/// @brief atom used for calculating residue-level neighbors
	VD nbr_atom_ = MutableResidueType::null_vertex;
	/// @brief radius cutoff to define neighbors
	/// @details Should be the maximum distance from the nbr_atom_
	/// to any heavy atom in any valid rotamer.
	Real nbr_radius_ = 0;

	//////////////////////////////////////////////////////////////////////
	// per-atom properties

	/// @brief Data to describe virtual atoms that should shadow other atoms for the sake
	/// of keeping intraresidue cycles closed when working with an atom tree, e.g.
	/// NV shadows N on proline. For each atom, the following vector lists the index
	/// of the atom it is shadowing.
	std::map<VD, VD> atom_shadowed_;

	/// @brief Verticies of all mainchain atoms
	/// @details mainchain_atoms are those atoms on a path from polymer lower_connect
	/// to upper_connect. For protein, this will be N, CA and C.
	utility::vector1<VD> mainchain_atoms_;

	///////////////////////////////////////////////////////////////////////////
	// Chi information

	// The bundled information about the chis on this ResidueType
	utility::vector1< MutableChiRecordOP > chis_;

	/// @brief VDs of four atoms to build each nu angle
	utility::vector1< utility::vector1< VD > > nu_atoms_;

	/// @brief VDs of all ring atoms, not counting virtual atoms
	utility::vector1< utility::vector1< VD > > ring_atoms_;  // indexed by ring number

	/// @brief The saturation type of each ring in this residue
	utility::vector1< core::chemical::rings::RingSaturationType > ring_saturation_types_;  // indexed by ring number

	/// @brief   Low-energy ring conformers for each ring
	/// @details used for setting up the RingConformerSets
	utility::vector1< utility::vector1< std::string > > low_ring_conformers_;  // indexed by ring number

	/// @brief   Lowest-energy ring conformer for each ring
	/// @details used for setting up the RingConformerSets
	utility::vector1< std::string > lowest_ring_conformer_;  // indexed by ring number

	/// @brief  Vector of inter-residue connections expected for this residuetype
	/// NOW includes the polymer connections, as well as disulf-type connections
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	/// @note  Much of the code assumes at most one residue connection per atom.
	///        The pseudobond code has been partially written to handle cases where multiple connections
	///        to a single atom can exist, but much of the standard residue connection code assumes a
	///        simple correspondence between atoms and residue connections.  That code will have to be
	///        updated to support single-atom "backbones."
	utility::vector1< MutableResidueConnection > residue_connections_;

	/// @brief Polymer lower connections
	/// Updated in set_lower_connect_atom()
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	Size lower_connect_id_ = 0; // which connection is the lower connection?

	/// @brief  Polymer upper connections
	/// Updated in set_upper_connect_atom()
	/// @note  ResidueConnection objects store the ideal internal coordinates for the connected atom
	Size upper_connect_id_ = 0; // which connection is the upper connection?

public:

	static VD const null_vertex;

#ifdef    SERIALIZATION
public:
	friend class cereal::access;
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};  // class MutableResidueType

// Insertion operator (overloaded so that MutableResidueType can be "printed" in PyRosetta).
std::ostream & operator<<(std::ostream & output, MutableResidueType const & object_to_output);


} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_MutableResidueType )
#endif // SERIALIZATION

#endif // INCLUDED_core_chemical_Residues_HH
