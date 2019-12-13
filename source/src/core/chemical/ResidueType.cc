// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResidueType.cc
/// @brief Method definitions for ResidueType
/// @author
/// Phil Bradley
/// Rocco Moretti (rmorettiase@gmail.com)
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to rings, properties, rings, lipids, carbohydrates, and other non-AAs)

// Unit headers
#include <core/chemical/ResidueType.hh>

// Package Headers

// Project Headers
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomProperties.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/Orbital.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/RestypeDestructionEvent.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>

//// For conversion
#include <core/chemical/MutableResidueType.hh> // For conversion functions
#include <core/chemical/MutableResidueConnection.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// External headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzVector.io.hh>

// C++ headers
#include <algorithm>

#ifdef    SERIALIZATION
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Utility serialization headers
#include <numeric/xyz.serialization.hh>
#include <utility/vector1.srlz.hh>
#include <utility/keys/Key2Tuple.srlz.hh>
#include <utility/keys/Key3Tuple.srlz.hh>
#include <utility/keys/Key4Tuple.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.ResidueType" );

ResidueType::ResidueType() = default; // private, not deleted because of serialization

/// @details The construction of a ResidueType from a MutableResidueType happens in several stages
/// 1. The canonical order of atoms (backbone heavy/sidechain heavy/hydrogens) is established.
/// 2. Atom and bond information is converted from the graph-based to the list-base representation
/// 3. The other information in the MutableResidueType is converted.
/// 4. Derived data can be pre-computed from the initial ResidueType state.
/// 5. Subobjects which need back-OPs to the ResidueType can be updated (done by ResidueType::make(), not the constructor).
ResidueType::ResidueType( MutableResidueType const & mrt ):
	ResidueTypeBase( mrt ),
	base_type_cop_( mrt.get_base_type_cop() )
{
	if ( ! mrt.validate_residue_type() ) {
		utility_exit_with_message( "Residue type " + mrt.name() + " is not appropriately configured." );
	}

	utility::vector1< core::Size > old_to_new = setup_atom_ordering( mrt );
	copy_atom_info( mrt, old_to_new );
	copy_other_info( mrt, old_to_new );

	initialize_derived_data();
	update_derived_data();
}

ResidueTypeCOP
ResidueType::make( MutableResidueType const & mrt ) {
	ResidueTypeOP new_restype( new ResidueType( mrt ) );
	new_restype->self_pointer_updates();

	// The ResidueType should be finished now.
	// Just double check if we have noticable issues.
	new_restype->perform_checks();

	return new_restype;
}

ResidueType::~ResidueType()
{
	destruction_obs_hub_( RestypeDestructionEvent( this ) ); // Notify the destruction observers that this residue type is being destroyed.
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////////////// Summary Functions               ////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/// @brief Counts the number of virtual atoms and returns the count.
/// @details The virtual count is not stored in the residue type.  This count is performed on the fly, and
///can hurt performance if repeatedly carried out.  Not intended for use in large loops -- instead, call
///once and store the value.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
Size
ResidueType::n_virtual_atoms() const
{
	core::Size virtcount = 0;
	// RM: Do we want to only look at the atom type, or do we need to also look at atom properties?
	for ( auto const & type: vec_atom_types_ ) {
		if ( type.is_virtual() ) {
			++virtcount;
		}
	}
	return virtcount;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////////// Atom Functions              ////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

core::chemical::element::Elements
ResidueType::element(core::Size atomno) const {
	ElementCOP element( element_type(atomno) );
	if ( ! element ) {
		return core::chemical::element::UnknownElement;
	} else {
		return element->element();
	}
}

bool
ResidueType::heavyatom_has_polar_hydrogens( Size atomno ) const {
	PyAssert((atomno > 0) && (atomno <= heavyatom_has_polar_hydrogens_.size()), "ResidueType::heavyatom_has_polar_hydrogens(): atomno is not in this ResidueType!");
	return heavyatom_has_polar_hydrogens_[ atomno ];
}

/// @brief  Check if atom is virtual.
bool
ResidueType::is_virtual( Size const atomno ) const
{
	return ( atom_type( atomno ).is_virtual() || atom_has_property( atomno, VIRTUAL_ATOM ) );
}

/// @brief  Check if atom is repulsive.
bool
ResidueType::is_repulsive( Size const atomno ) const
{
	return ( atom_type( atomno ).is_repulsive() );
}


bool
ResidueType::atom_has_property( Size const atomno, AtomProperty const property ) const {
	PyAssert((atomno > 0) && (atomno <= atomic_properties_.size()), "ResidueType::atom_has_property(): atomno is not in this ResidueType!");
	if ( atomic_properties_[ atomno ] == nullptr ) {
		return false; // Can't have the property if you don't have properties
	} else {
		return atomic_properties_[ atomno ]->has_property( property );
	}
}

void
ResidueType::show_all_atom_names( std::ostream & out ) const {
	for ( core::Size ii(1); ii <= atom_names_.size(); ++ii ) {
		out << ii << "  '" << atom_names_[ii] << "'" << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////////// Bond Functions              ////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

Size
ResidueType::number_bonded_heavyatoms( Size const atomno ) const
{
	return bonded_neighbor(atomno).size() - number_bonded_hydrogens( atomno );
}

/// @brief indices of the bonded neighbors for an atom
AtomIndices const &
ResidueType::bonded_neighbor( Size const atomno ) const
{
	return bonded_neighbor_[atomno];
}

utility::vector1< std::pair< core::Size, core::Size > >
ResidueType::bonds() const {
	utility::vector1< std::pair< core::Size, core::Size > > bondlist;
	for ( core::Size ii(1); ii <= bonded_neighbor_.size(); ++ii ) {
		for ( core::Size bn: bonded_neighbor_[ii] ) {
			if ( ii <= bn ) {
				bondlist.push_back( std::make_pair( ii, bn ) );
			}
		}
	}
	debug_assert( bondlist.size() == nbonds_ );
	return bondlist;
}

/// @brief Indicates whether or not two atom indices have a chemical bond linking them.
/// @details Note that this assumes that the Rosetta machinery is set up so that if
/// atom 1 is bonded to atom 2, atom 2 is bonded to atom 1.  This function breaks if
/// that assumption breaks.
/// @author Vikram K. Mulligan
bool
ResidueType::atoms_are_bonded(
	core::Size const atom_index1,
	core::Size const atom_index2
) const {
	AtomIndices const atom1partners = nbrs(atom_index1);
	for ( core::Size i=1, imax=atom1partners.size(); i<=imax; ++i ) {
		if ( atom1partners[i] == atom_index2 ) return true;
	}
	return false;
}

utility::vector1<BondName> const &
ResidueType::bonded_neighbor_types(Size const atomno) const
{
	return bonded_neighbor_type_[atomno];
}

utility::vector1<BondRingness> const &
ResidueType::bonded_neighbor_ringnesses(Size const atomno) const
{
	return bonded_neighbor_ringness_[atomno];
}

BondName
ResidueType::bond_type(core::Size const atom_index1, core::Size const atom_index2 ) const {
	AtomIndices const & atom1partners = bonded_neighbor_[atom_index1];
	for ( core::Size ii=1, iimax=atom1partners.size(); ii<=iimax; ++ii ) {
		if ( atom1partners[ii] == atom_index2 ) {
			return bonded_neighbor_type_[atom_index1][ii];
		}
	}
	return UnknownBond;
}

BondRingness
ResidueType::bond_ringness(core::Size const atom_index1, core::Size const atom_index2 ) const {
	AtomIndices const & atom1partners = bonded_neighbor_[atom_index1];
	for ( core::Size ii=1, iimax=atom1partners.size(); ii<=iimax; ++ii ) {
		if ( atom1partners[ii] == atom_index2 ) {
			return bonded_neighbor_ringness_[atom_index1][ii];
		}
	}
	return UnknownRingness;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////////////// Chi Functions               ////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/// @brief Return whether this atom is in a particular ring
bool
ResidueType::is_ring_atom( uint const ring_num, uint const atom_id ) const{
	debug_assert( ring_num <= ring_atoms_.size() );
	return std::find(ring_atoms_[ ring_num ].begin(), ring_atoms_[ ring_num ].end(), atom_id) != ring_atoms_[ ring_num ].end();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////////////// Icoord Functions               //////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/// @brief AtomICoord of an atom
AtomICoor const &
ResidueType::icoor( Size const atm ) const
{
	debug_assert( 1 <= atm && atm <= icoor_.size() );
	return icoor_[ atm ];
}

/// @brief get index of an atom's base atom
Size
ResidueType::atom_base( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <= atom_base_indices_.size()), "ResidueType::atom_base( Size const atomno ): atomno is not in this ResidueType!");
	return atom_base_indices_[atomno];
}

/// @brief get index of an atom's second base atom
Size
ResidueType::abase2( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <= atom_base_indices_.size()), "ResidueType::abase2( Size const atomno ): atomno is not in this ResidueType!");
	return abase2_indices_[atomno];
}

void
ResidueType::debug_dump_icoor() const
{

	TR.Debug << "ICoor for " << name3() << std::endl;
	for ( Size ii = 1 ; ii <= natoms() ; ++ii ) {
		TR.Debug << " Atom " << ii << " atom name: " << atom_name( ii ) << " ideal xyz " << ideal_xyz(ii) << std::endl;
	}
	pretty_print_atomicoor(TR.Debug, *this);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////////////// MMAtom Functions              //////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/// @brief Get the MM atom_type for this atom by its index number in this residue
MMAtomType const &
ResidueType::mm_atom_type( Size const atomno ) const
{
	MMAtomTypeSetCOP mm_atom_types( mm_atom_types_ptr() );
	debug_assert( mm_atom_types != nullptr );
	return ( *mm_atom_types )[ mm_atom_type_index( atomno ) ];
}

std::string const &
ResidueType::mm_name( Size const atomno ) const
{
	MMAtomTypeSetCOP mm_atom_types( mm_atom_types_ptr() );
	debug_assert( mm_atom_types != nullptr );
	return ( *mm_atom_types )[ mm_atom_type_index( atomno ) ].name();
}

gasteiger::GasteigerAtomTypeDataCOP
ResidueType::gasteiger_atom_type(core::Size atomno) const
{
	if ( gasteiger_atom_type_index_[atomno] == 0 ) { return nullptr; }
	gasteiger::GasteigerAtomTypeSetCOP gast_set( gasteiger_atom_typeset() );
	if ( gast_set == nullptr ) { return nullptr; }
	return (*gast_set)[ gasteiger_atom_type_index_[atomno] ];
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////////////////          Orbital Functions     //////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

orbitals::ICoorOrbitalData const &
ResidueType::orbital_icoor_data(Size const orbital_index) const{
	return orbital( orbital_index ).icoor();
}


orbitals::ICoorOrbitalData const &
ResidueType::new_orbital_icoor_data(Size const orbital_index) const{
	return orbital( orbital_index ).new_icoor();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////////////////  Ring Conformer Set Functions  //////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

// Return a pointer to the object containing the set of ring conformers possible for this residue's nth cycle.
core::chemical::rings::RingConformerSetCOP
ResidueType::ring_conformer_set( core::uint ring_num ) const {
	if ( ring_num <= ring_conformer_sets_.size() ) {
		return ring_conformer_sets_[ ring_num ];
	} else {
		return nullptr;
	}
}

core::Size
ResidueType::n_ring_conformer_sets() const
{
	return ring_conformer_sets_.size();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////////////////      Connection Functions      //////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

// Connections ////////////////////////////////////////////////////////////////
// Lower
ResidueConnection const &
ResidueType::lower_connect() const
{
	debug_assert( properties().has_property( POLYMER ) );
	debug_assert( lower_connect_id_ != 0 );
	return residue_connections_[ lower_connect_id_ ];
}

Size
ResidueType::lower_connect_atom() const {
	debug_assert( properties().has_property( POLYMER ) );
	debug_assert( lower_connect_id_ != 0 );
	return residue_connections_[ lower_connect_id_ ].atomno();
}

// Upper
ResidueConnection const &
ResidueType::upper_connect() const
{
	debug_assert( properties().has_property( POLYMER ) );
	debug_assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ];
}

Size
ResidueType::upper_connect_atom() const
{
	debug_assert( properties().has_property( POLYMER ) );
	debug_assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ].atomno();
}

// Branches / Non-polymer
// Return a list of indices of atoms at non-polymer connections.
utility::vector1< uint >
ResidueType::branch_connect_atoms() const
{
	utility::vector1< uint > atoms;
	Size const n_connections( n_possible_residue_connections() );
	for ( uint i( 1 ); i <= n_connections; ++i ) {
		if ( i == lower_connect_id_ || i == upper_connect_id_ ) { continue; }
		atoms.push_back( residue_connect_atom_index( i ) );
	}
	debug_assert( atoms.size() == n_non_polymeric_residue_connections_ );

	return atoms;
}

// Return a list of names of atoms at non-polymer connections.
utility::vector1< std::string >
ResidueType::branch_connect_atom_names() const
{
	utility::vector1< uint > const atoms( branch_connect_atoms() );
	Size const n_atoms( atoms.size() );
	utility::vector1< std::string > names( n_atoms );
	for ( uint i( 1 ); i <= n_atoms; ++i ) {
		names[ i ] = atom_name( atoms[ i ] );
	}
	return names;
}

Size
ResidueType::n_possible_residue_connections() const {
	return residue_connections_.size();
}

ResidueConnection const &
ResidueType::residue_connection( Size const i ) const
{
	return residue_connections_[ i ];
}

Size
ResidueType::residue_connect_atom_index( Size const resconn_id ) const {
	return residue_connections_[ resconn_id ].atomno();
}

/// @brief Does an atom with a given index have an icoor that depends, directly or indirectly, on the lower polymeric connection?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::atom_depends_on_lower_polymeric_connection( core::Size const atom_index ) const {
	debug_assert( atom_index > 0 && atom_index <= natoms() );
	if ( !is_polymer() ) return false;
	return atom_depends_on_lower_polymeric_connection_[atom_index];
}

/// @brief Does an atom with a given index have an icoor that depends, directly or indirectly, on the upper polymeric connection?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::atom_depends_on_upper_polymeric_connection( core::Size const atom_index ) const {
	debug_assert( atom_index > 0 && atom_index <= natoms() );
	if ( !is_polymer() ) return false;
	return atom_depends_on_upper_polymeric_connection_[atom_index];
}

/// @brief Does an atom with a given index have an icoor that depends, directly or indirectly, on the upper or lower polymeric connection?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::atom_depends_on_polymeric_connection( core::Size const atom_index ) const {
	debug_assert( atom_index > 0 && atom_index <= natoms() );
	if ( !is_polymer() ) return false;
	return atom_depends_on_lower_polymeric_connection_[atom_index] || atom_depends_on_upper_polymeric_connection_[atom_index];
}

/// @brief Does an atom with a given index have an icoor that depends, directly or indirectly, on a particular connection ID?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::atom_depends_on_connection( core::Size const atom_index, core::Size const connection_id ) const {
	debug_assert( atom_index > 0 && atom_index <= natoms() );
	debug_assert( connection_id > 0 && connection_id <= n_possible_residue_connections() );
	return atom_depends_on_connection_[connection_id][atom_index];
}

//////////////////////////////////////////////////////////////////////
///////////////// Actcoord                       /////////////////////
//////////////////////////////////////////////////////////////////////

bool
ResidueType::requires_actcoord() const
{
	return properties().has_property( PROTEIN ) &&
		( properties().has_property( POLAR ) || properties().has_property( AROMATIC ) ) &&
		actcoord_atoms_indices_.size() != 0;
}

//////////////////////////////////////////////////////////////////////
///////////////// properties /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/// @brief Is this one of SRI's special heteropolymer building blocks?
///
bool
ResidueType::is_sri() const {
	return properties().has_property( SRI );
}

/// @brief Is this a triazolemer?
///
bool
ResidueType::is_triazolemer() const
{
	return properties().has_property( TRIAZOLE_LINKER );
}

bool
ResidueType::is_coarse() const
{
	return properties().has_property( COARSE );
}

bool
ResidueType::is_purine() const
{
	return properties().has_property( PURINE );
}

bool
ResidueType::is_pyrimidine() const
{
	return properties().has_property( PYRIMIDINE );
}

/// @brief Is this a solvent molecule (SOLVENT property)?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::is_solvent() const
{
	return properties().has_property( SOLVENT );
}

/// @brief Is this a canonical nucleic acid (CANONICAL_NUCLEIC property)?
/// @details Only the standard nucliec acid types (dA, dC, dG, dT, A, C, G, U) are canonical.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::is_canonical_nucleic() const
{
	return properties().has_property( CANONICAL_NUCLEIC );
}

/// @brief Is this a canonical amino acid (CANONICAL_AA property)?
/// @details Only the standard amino acid types (ACDEFGHIKLMNPQRSTVWY) are canonical.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::is_canonical_aa() const
{
	return properties().has_property( CANONICAL_AA );
}

/// @brief Is this a canonical residue type (nucleic acid or amino acid)?
/// @details Calls is_canonical_aa() and is_canonical_nucleic().
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::is_canonical() const
{
	return is_canonical_aa() || is_canonical_nucleic();
}

bool
ResidueType::is_carbohydrate() const
{
	return properties().has_property( CARBOHYDRATE );
}

bool
ResidueType::is_ligand() const
{
	return properties().has_property( LIGAND );
}

bool
ResidueType::is_lipid() const
{
	return properties().has_property( LIPID );
}

/// @details The METAL property is specified in the params file under PROPERTIES.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
ResidueType::is_metal() const
{
	return properties().has_property( METAL );
}

/// @details The METALBINDING property is specified in the params file under PROPERTIES.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
ResidueType::is_metalbinding() const
{
	return properties().has_property( METALBINDING );
}

bool
ResidueType::is_membrane() const
{
	return properties().has_property( MEMBRANE );
}

bool
ResidueType::is_surface() const
{
	return properties().has_property( SURFACE );
}

bool
ResidueType::has_sc_orbitals() const
{
	return properties().has_property( SC_ORBITALS );
}

bool
ResidueType::is_polar() const
{
	return properties().has_property( POLAR );
}

bool
ResidueType::is_charged() const
{
	return properties().has_property( CHARGED );
}

bool
ResidueType::is_aromatic() const
{
	return properties().has_property( AROMATIC );
}

bool
ResidueType::is_cyclic() const
{
	return properties().has_property( CYCLIC );
}

bool
ResidueType::is_lower_terminus() const
{
	return properties().has_property( LOWER_TERMINUS );
}

bool
ResidueType::is_upper_terminus() const
{
	return properties().has_property( UPPER_TERMINUS );
}

bool
ResidueType::is_branch_point() const
{
	return properties().has_property( BRANCH_POINT );
}

bool
ResidueType::is_acetylated_nterminus() const
{
	return properties().has_property( ACETYLATED_NTERMINUS );
}

bool
ResidueType::is_methylated_cterminus() const
{
	return properties().has_property( METHYLATED_CTERMINUS );
}

bool
ResidueType::is_virtual_residue() const
{
	return properties().has_property( VIRTUAL_RESIDUE );
}

/// @brief  Check if residue is 'INVERTING_VIRTUAL_RESIDUE'
/// @details Used by the symmetry machinery for mirror symmetry operations.
bool
ResidueType::is_inverted_virtual_residue() const
{
	return properties().has_property( INVERTED_VIRTUAL_RESIDUE );
}

/// @brief  Check if residue is 'VRT1'
bool
ResidueType::is_VRT1() const
{
	return properties().has_property( VRT1 );
}

/// @brief  Check if residue is a TP3 water.
bool
ResidueType::is_TP3() const
{
	return properties().has_property( TP3 );
}


bool
ResidueType::is_adduct() const
{
	return properties().has_property( ADDUCT );
}

/// @brief Is this ResidueTypeBase a base type?
/// @details Checks the base_type_cop_ pointer.  If it's null, this is assumed to be a base type.
bool
ResidueType::is_base_type() const {
	return base_type_cop_ == nullptr;
}

/// @brief Get a pointer to this ResidueTypeBase's base ResidueType.
/// @details Returns the base_type_cop_ pointer if not null, self pointer if null.
ResidueTypeCOP
ResidueType::get_base_type_cop() const {
	if ( base_type_cop_ ) {
		return base_type_cop_;
	}
	return get_self_ptr();
}

core::chemical::rna::RNA_Info const &
ResidueType::RNA_info() const{
	return ( *rna_info_ );
}

// Return the CarbohydrateInfo object containing sugar-specific properties for this residue.
core::chemical::carbohydrates::CarbohydrateInfoCOP
ResidueType::carbohydrate_info() const
{
	return carbohydrate_info_;
}

//////////////////////////////////////////////////////////////////////
///////////////// Dihedral Methods   /////////////////////////////////
//////////////////////////////////////////////////////////////////////

void
ResidueType::print_bondangles() const
{
	TR.Debug << "START BOND ANGLES ATOM NAMES" << std::endl;
	TR.Debug << "Number of bond angles: " << bondangle_atom_sets_.size() << std::endl;
	for ( Size i = 1; i <= bondangle_atom_sets_.size(); ++i ) {
		AtomType at1 = atom_type( bondangle_atom_sets_[ i ].key1() );
		AtomType at2 = atom_type( bondangle_atom_sets_[ i ].key2() );
		AtomType at3 = atom_type( bondangle_atom_sets_[ i ].key3() );
		MMAtomType at5 = mm_atom_type( bondangle_atom_sets_[ i ].key1() );
		MMAtomType at6 = mm_atom_type( bondangle_atom_sets_[ i ].key2() );
		MMAtomType at7 = mm_atom_type( bondangle_atom_sets_[ i ].key3() );

		TR.Debug << "PDB:" << "\t"
			<< atom_name( bondangle_atom_sets_[ i ].key1() ) << "\t"
			<< atom_name( bondangle_atom_sets_[ i ].key2() ) << "\t"
			<< atom_name( bondangle_atom_sets_[ i ].key3() ) << "\t"
			<< "MM:" << "\t"
			<< mm_name( bondangle_atom_sets_[ i ].key1() ) << "\t"
			<< mm_name( bondangle_atom_sets_[ i ].key2() ) << "\t"
			<< mm_name( bondangle_atom_sets_[ i ].key3() ) << "\t"
			<< "MM2:" << "\t"
			<< at5.name() << "\t"
			<< at6.name() << "\t"
			<< at7.name() << "\t"
			<< "ROS:" << "\t"
			<< at1.name() << "\t"
			<< at2.name() << "\t"
			<< at3.name() << "\t"
			<< std::endl;
	}
	TR.Debug << "END BOND ANGLES ATOM NAMES" << std::endl;
}

void
ResidueType::print_pretty_path_distances() const
{
	TR.Debug << "START PATH DISTANCES" << std::endl;
	// print header line
	for ( Size i = 1; i <= natoms(); ++i ) {
		TR.Debug << "\t" << atom_name( i );
	}
	TR.Debug << std::endl;

	for ( Size j = 1; j <= natoms(); ++j ) {
		TR.Debug << atom_name( j ) << "\t";
		for ( Size k = 1; k <= natoms(); ++k ) {
			TR.Debug << path_distance_[j][k] << "\t";
		}
		TR.Debug << std::endl;
	}
	TR.Debug << "END PATH DISTANCES" << std::endl;
}

/////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void
ResidueType::select_orient_atoms(
	Size & center,
	Size & nbr1,
	Size & nbr2
) const
{
	center = 0;
	nbr1 = 0;
	nbr2 = 0;

	// No backbone atoms, all backbone atoms, or orient mode explicitly set to nbr_atom
	if ( first_sidechain_atom() == 1 || first_sidechain_atom() > natoms() || force_nbr_atom_orient() ) {
		// If no backbone atoms (or all bb atoms), assume nbr_atom will be close to center-of-mass.
		center = nbr_atom();

		// If is hydrogen or too few neighbors, try trekking up the atom tree
		while ( center > nheavyatoms() || bonded_neighbor(center).size() < 2 ) {
			center = atom_base(center);
		}
		AtomIndices const & nbrs( bonded_neighbor(center) );
		// First try to find two neighbors that are heavyatoms
		for ( Size j=1; j<= nbrs.size(); ++j ) {
			Size const nbr( nbrs[j] );
			if ( nbr <= nheavyatoms() ) {
				if ( nbr1 ) nbr2 = nbr;
				else nbr1 = nbr;
			}
		}
		// Failing that, just try for two neighbors!
		if ( !( center && nbr1 && nbr2 ) ) {
			for ( Size j=1; j<= nbrs.size(); ++j ) {
				Size const nbr( nbrs[j] );
				if ( nbr1 ) nbr2 = nbr;
				else nbr1 = nbr;
			}
		}
		if ( !( center && nbr1 && nbr2 ) ) {
			//debug_assert() isn't enough for these cases b/c they're typically ligands
			// and thus depend on user input -- need to be caught even in release mode.
			utility_exit_with_message( "Cannot superimpose residues of type " + name() );
		}
	} else {
		// look for a backbone atom, one of whose neighbors is a sidechain atom
		// center will be this atom
		// nbr1 and nbr2 will be the backbone heavyatom nbrs of this atom
		// eg center = CA, nbr1 = N. nbr2 = C in the protein case
		select_orient_atoms_standard_logic( center, nbr1, nbr2, true ); //First try ignoring virtuals.
		if ( !( center && nbr1 && nbr2 ) ) {
			center = 0; nbr1 = 0; nbr2 = 0;
			select_orient_atoms_standard_logic( center, nbr1, nbr2, false ); //If we fail, try again allowing virtuals.
		}
	}
	//TR.Debug << "Superimposing on " << atom_name(center) << " " << atom_name(nbr1) << " " << atom_name(nbr2) << "\n";
	runtime_assert_string_msg( center && nbr1 && nbr2, "Error in core::chemical::ResidueType::select_orient_atoms(): Could not find three atoms to use for orienting residue type " + name() + "." );
}

/// @brief Pick atoms to use for orienting one Residue onto another, using standard logic.
/// @details Standard logic applies to case in which (a) the residue has backbone atoms, and (b) the residue
/// has sidechain atoms, and (c) the orient mode has not been set explicitly to force_nbr_atom_orient.  We loop through
/// all backbone atoms and find the first atom that is bonded to a sidechain atom AND two other backbone atoms.  The
/// first such atom becomes "center", and its two backbone neighbors become "nbr1" and "nbr2".
/// @note If ignore_virtuals is true, none of the atoms involved can be virtuals.  If false, they can be.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
ResidueType::select_orient_atoms_standard_logic(
	Size & center,
	Size & nbr1,
	Size & nbr2,
	bool const ignore_virtuals
) const {
	// look for a backbone atom, one of whose neighbors is a sidechain atom
	// center will be this atom
	// nbr1 and nbr2 will be the backbone heavyatom nbrs of this atom
	// eg center = CA, nbr1 = N. nbr2 = C in the protein case
	for ( Size atom_index( 1 ); atom_index <= natoms(); ++atom_index ) {
		if ( atom_is_backbone( atom_index ) ) {
			AtomIndices const & nbrs( bonded_neighbor( atom_index ) );
			center = 0; nbr1 = 0; nbr2 = 0;
			for ( Size nbr_index( 1 ); nbr_index <= nbrs.size(); ++nbr_index ) {
				Size const nbr( nbrs[ nbr_index ] );
				if ( ignore_virtuals && is_virtual(nbr) ) continue;
				if ( !atom_is_backbone( nbr ) && atom_base( nbr ) == atom_index ) {
					// nbr is a sidechain atom that branches from the atom at atom_index.
					// It is also not a virtual atom, (which is important for avoiding vector normalization errors.)
					center = atom_index;
				} else if ( atom_is_backbone( nbr ) && nbr <= nheavyatoms() ) {
					// nbr is a (real) backbone heavy atom neighbor of the atom at atom_index
					if ( nbr1 ) { nbr2 = nbr; }
					else { nbr1 = nbr; }
				}
			}
		} // atom_index is backbone
		if ( center && nbr1 && nbr2 ) break;
	} // atom_index
}

std::tuple<Size, Size, Size>
ResidueType::select_orient_atoms() const
{
	Size center, nbr1, nbr2;
	select_orient_atoms(center, nbr1, nbr2);

	return std::tuple<Size, Size, Size>(center, nbr1, nbr2);
}

///////////////////////////////////////////////////////////////
/// @author Labonte <JWLabonte@jhu.edu>
void
ResidueType::show( std::ostream & output, bool output_atomic_details ) const
{
	using namespace std;
	using namespace utility;

	output << name() << " (" << name3() << ", " << name1() << "):" << endl;

	output << "Base: " << base_name() << std::endl;

	properties().show( output );

	output << " Main-chain atoms:";
	Size const n_mainchain_atoms(mainchain_atoms_indices_.size() );
	for ( uint i = 1; i <= n_mainchain_atoms; ++i ) {
		output << ' ' << atom_name( mainchain_atoms_indices_[ i ] );
	}
	output << endl;

	output << " Backbone atoms:  ";
	Size const n_bb_atoms( all_bb_atoms_.size() );
	for ( uint i = 1; i <= n_bb_atoms; ++i ) {
		output << ' ' << atom_name( all_bb_atoms_[ i ] );
	}
	output << endl;

	if ( is_cyclic() ) {
		for ( uint i( 1 ); i <= n_rings(); ++i ) {
			output << " Ring atoms:  ";
			Size const n_ring_atoms( ring_atoms_[ i ].size() );
			for ( uint j = 1; j <= n_ring_atoms; ++j ) {
				output << ' ' << atom_name( ring_atoms_[ i ][ j ] );
			}
			output << endl;
		}
	}

	output << " Side-chain atoms:";
	Size const n_sc_atoms( all_sc_atoms_.size() );
	for ( uint i = 1; i <= n_sc_atoms; ++i ) {
		output << ' ' << atom_name( all_sc_atoms_[ i ] );
	}
	output << endl;

	if ( is_branch_point() ) {
		output << " Branch-point atoms:";
		vector1< string > const atom_names( branch_connect_atom_names() );
		Size const n_atoms( atom_names.size() );
		for ( uint i( 1 ); i <= n_atoms; ++i ) {
			output << ' ' << atom_names[ i ];
		}
		output << endl;
	}

	if ( properties().has_property( CARBOHYDRATE ) ) {
		carbohydrate_info_->show( output );
	}

	if ( output_atomic_details ) {
		output << " Atomic Details:" << endl;
		Size const n_atoms( natoms() );
		for ( uint i( 1 ); i <= n_atoms; ++i ) {
			output << "  Atom " << i << ": " << atom_name( i );
			if ( ! is_virtual(i) ) {
				output << " (" << element_type(i)->get_chemical_symbol() << ')';
			} else {
				output << " (virtual)";
			}
			output << endl;

			output << "   Types (type set indices): ";
			output << "Rosetta: " /*<< type_name_*/ << " (" << atom_type_index(i) << "); ";
			output << "CHARMm: " << mm_name(i) << " (" << mm_atom_type_index(i) << "); ";
			output << "Gasteiger: ";
			if ( gasteiger_atom_typeset() ) {
				output << (*gasteiger_atom_typeset())[ gasteiger_atom_type_index_[i] ]->get_name();
			} else {
				output << "None";
			}
			output << endl;

			output << "   Charge: " << "partial: " << atom_charge(i) << "; formal: ";
			if ( formal_charge(i) > 0 ) {
				output << '+';
			}
			output << formal_charge(i) << endl;

			output << "   Properties: ";
			if ( heavyatom_has_polar_hydrogens(i) ) {
				output << "H-bond donor, ";
			}
			if ( heavyatom_is_an_acceptor(i) ) {
				output << "H-bond acceptor, ";
			}
			if ( atom_is_hydrogen(i) ) {
				if ( atom_is_polar_hydrogen(i) ) {
					output << "polar hydrogen, ";
				} else if ( atom_is_aro_hydrogen(i) ) {
					output << "aromatic hydrogen, ";
				} else {
					output << "non-polar hydrogen, ";
				}
			}
			output << std::endl;
		}
	}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//// Private methods for MutableResidueType -> ResidueType conversion ////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


///@details Because there isn't necessarily any consistent order in the MutableResidueType
/// (particularly because patching often adds and removes atoms)
///we need to reorder the atoms. The atoms are reordered with the backbone
///atoms coming first, sidechain atoms second, and finally hydrogens last.
///
/// We attempt to keep the atom ordering more-or-less the same as the implicit ordering in the MutableResidueType.
//
// This will also set the following data on the ResidueType
// * natoms_
// * nheavyatoms_
// * n_backbone_heavyatoms_
// * first_sidechain_hydrogen_
// * attached_H_begin_
// * attached_H_end_
//
// The return type is the MutableResidueType indices in the new order
// (This is by index rather than VD to keep ResidueGraphTypes from
// being needed in the ResidueType.hh header.)
utility::vector1< core::Size >
ResidueType::setup_atom_ordering( MutableResidueType const & mrt )
{

	utility::vector1<VD> const & all_atoms( mrt.all_atoms() );

	natoms_ = all_atoms.size();

	// Now separate out the atoms by their type.
	utility::vector1< core::Size > bb_atoms, sidechain_atoms;
	bb_atoms.reserve( all_atoms.size() );
	sidechain_atoms.reserve( all_atoms.size() );
	for ( core::Size ii(1); ii <= all_atoms.size(); ++ii ) {
		Atom const & atom( mrt.atom( all_atoms[ii] ) );
		if ( atom.is_hydrogen() ) {
			continue; // Deal with this below
		} else if ( mrt.is_backbone_heavyatom( all_atoms[ii] ) ) {
			bb_atoms.push_back( ii );
		} else {
			sidechain_atoms.push_back( ii );
		}
	}
	n_backbone_heavyatoms_ = bb_atoms.size();
	nheavyatoms_ = bb_atoms.size() + sidechain_atoms.size();

	utility::vector1< core::Size > atom_order( bb_atoms );
	atom_order.reserve( all_atoms.size() );
	atom_order.append( sidechain_atoms );

	attached_H_begin_.reserve( nheavyatoms_ ); // Only for the heavy atoms
	attached_H_end_.reserve( nheavyatoms_ );
	for ( Size ii(1); ii <= nheavyatoms_; ++ii ) {
		utility::vector1< VD > bonded_h( mrt.bonded_hydrogens( all_atoms[ atom_order[ii] ] ) );
		attached_H_begin_.push_back( atom_order.size()+1 ); // One past the current (even if it puts it off the end
		if ( bonded_h.empty() ) { // Catch no hydrogens case
			attached_H_end_.push_back( 0 );
		} else {
			for ( VD bh: bonded_h ) {
				atom_order.push_back( mrt.atom_index(bh) );
			}
			attached_H_end_.push_back( atom_order.size() );
		}
	}
	if ( sidechain_atoms.empty() ) {
		first_sidechain_hydrogen_ = natoms_ + 1; // One off the end
	} else {
		first_sidechain_hydrogen_ = attached_H_begin_[n_backbone_heavyatoms_+1];
	}

	return atom_order;
}

/// @details This function updates most of the atomistic data in the ResidueType, by iterating over the atoms in the proper order
void
ResidueType::copy_atom_info(
	MutableResidueType const & mrt,
	utility::vector1< core::Size > const & old_to_new
) {

	debug_assert( atom_type_set_ptr() );
	debug_assert( mm_atom_types_ptr() );

	///////////////////////////////////////////////
	// Initialize data members for in-place updates
	nbonds_ = 0;
	n_hbond_acceptors_ = 0;
	n_hbond_donors_ = 0;
	mass_ = 0.0;
	atom_names_.reserve( natoms_ );
	vec_atom_types_.reserve( natoms_ );
	vec_atom_type_indexes_.reserve( natoms_ );
	elements_.reserve( natoms_ );
	mm_atom_type_index_.reserve( natoms_ );
	gasteiger_atom_type_index_.reserve( natoms_ );
	formal_charge_.reserve( natoms_ );
	partial_charge_.reserve( natoms_ );
	ideal_xyz_.reserve( natoms_ );
	atomic_properties_.reserve( natoms_ );
	bonded_orbitals_.reserve( natoms_ );
	bonded_neighbor_.resize( natoms_ ); // resize as we modify in-place
	bonded_neighbor_type_.resize( natoms_ ); // resize as we modify in-place
	bonded_neighbor_ringness_.resize( natoms_ ); // resize as we modify in-place
	cut_bond_neighbor_indices_.resize( natoms_ ); // Defaults empty
	heavyatom_has_polar_hydrogens_.reserve( natoms_ );
	// Others aren't reserved here, as they don't necessarily have the full length.

	//////////////////////////////////////////////////////////////////////////////////
	// Iterate over all the atoms in the MRT's graph, copying over the old information.

	for ( core::Size new_index(1); new_index <= old_to_new.size(); ++new_index ) {
		Size old_index = old_to_new[ new_index ];
		VD old_vd = mrt.all_atoms()[old_index];
		Atom const & atom( mrt.atom(old_vd ) );

		///////////////////////////////////////////////
		// Copy over the information in the Atom object
		atom_names_.push_back( atom.name() );
		atom_name_to_index_[ atom.name() ] = new_index;
		atom_name_to_index_[ utility::stripped_whitespace( atom.name() ) ] = new_index;

		vec_atom_type_indexes_.push_back( atom.atom_type_index() );
		AtomType const & atype( atom_type_set()[ atom.atom_type_index() ] );
		vec_atom_types_.push_back( atype );
		elements_.push_back( atom.element_type() );
		if ( atom.element_type() != nullptr ) {
			mass_ += atom.element_type()->weight();
		}
		if ( atom.mm_name().empty() ) {
			TR << "Atom " << atom.name() << " on residue " << ResidueTypeBase::name() << " is mising the MM type. Typing as VIRT." << std::endl;
			mm_atom_type_index_.push_back( mm_atom_types_ptr()->atom_type_index( "VIRT" ) ); // Type as VIRT if missing??
		} else {
			mm_atom_type_index_.push_back( mm_atom_types_ptr()->atom_type_index( atom.mm_name() ) );
		}
		if ( gasteiger_atom_typeset() == nullptr || atom.gasteiger_atom_type() == nullptr ) {
			gasteiger_atom_type_index_.push_back( 0 );
		} else {
			gasteiger_atom_type_index_.push_back( gasteiger_atom_typeset()->atom_type_index( atom.gasteiger_atom_type()->get_name() ) );
		}

		formal_charge_.push_back( atom.formal_charge() );
		partial_charge_.push_back( atom.charge() );
		ideal_xyz_.push_back( atom.ideal_xyz() );
		// Clone because we will modify it later (and don't want any additional modifications to propagate).
		atomic_properties_.push_back( utility::pointer::make_shared< AtomProperties >( atom.properties() ) );
		bonded_orbitals_.push_back( atom.bonded_orbitals() );

		// We defer on icoor, as we need the name->index mapping

		//////////////////////////////////////////////////////////
		// Update atom property annotations based on current atom
		if ( atype.is_acceptor() ) { ++n_hbond_acceptors_; }
		if ( atype.is_donor() ) { ++n_hbond_donors_; }
		if ( atype.atom_has_orbital() ) { atoms_with_orb_index_.push_back(new_index); }
		if ( atype.is_haro() ) { Haro_index_.push_back(new_index); }
		if ( atype.is_polar_hydrogen() ) { Hpol_index_.push_back( new_index ); }
		if ( ! atom.is_virtual() && ! atype.is_virtual() ) {
			if ( atype.is_acceptor() ) {
				accpt_pos_.push_back( new_index );
				if ( new_index > n_backbone_heavyatoms_ ) {
					accpt_pos_sc_.push_back( new_index );
				}
			}
			if ( atype.is_polar_hydrogen() ) {
				Hpos_polar_.push_back( new_index );
				if ( new_index > first_sidechain_hydrogen_ ) {
					Hpos_polar_sc_.push_back( new_index );
				}
			}
		} // endif is not virtual
		if ( atype.is_hydrogen() && ! atype.is_polar_hydrogen() ) {
			Hpos_apolar_.push_back( new_index );
		}

		if ( atype.is_hydrogen() ) {
			if ( new_index < first_sidechain_hydrogen_ ) {
				all_bb_atoms_.push_back( new_index );
			} else {
				all_sc_atoms_.push_back( new_index );
			}
		} else {
			if ( new_index <= n_backbone_heavyatoms_ ) {
				all_bb_atoms_.push_back( new_index );
			} else {
				all_sc_atoms_.push_back( new_index );
			}
		}

		////////////////////////////////////////////////////////////////////////////
		// Copy over the information for Bonds that are attached to the current Atom

		utility::vector1< VD > bond_nbr_vd( mrt.bonded_neighbors( old_vd ) );
		utility::vector1< Size > & bond_nbr_index( bonded_neighbor_[ new_index ] );
		utility::vector1<BondName> & bond_nbr_type( bonded_neighbor_type_[ new_index ] );
		utility::vector1<BondRingness> & bond_nbr_ringness( bonded_neighbor_ringness_[ new_index ] );
		nbonds_ += bond_nbr_vd.size(); // This double-counts (once for each atom) -- we'll correct that later

		bool has_polar_H = false;
		// We want to sort the bonded neghbor entries in increasing order of (new) atom index.
		utility::vector1< std::pair< core::Size, VD > > pair_vd;
		for ( VD bnd_vd: bond_nbr_vd ) {
			pair_vd.push_back( std::make_pair( old_to_new.index( mrt.atom_index( bnd_vd ) ), bnd_vd ) );
		}
		std::sort( pair_vd.begin(), pair_vd.end() ); // Sorts lexicographically

		for ( auto const & bnd_entry: pair_vd ) {
			core::Size nbr_index = bnd_entry.first;
			VD bnd_vd = bnd_entry.second;

			bond_nbr_index.push_back( nbr_index );
			Bond const & bond( mrt.bond(old_vd, bnd_vd ) );
			bond_nbr_type.push_back( bond.bond_name() );
			bond_nbr_ringness.push_back( bond.ringness() );

			if ( bond.cut_bond() ) {
				cut_bond_neighbor_indices_[ new_index ].push_back( nbr_index );
			}

			AtomType const & bnd_atype( atom_type_set()[ mrt.atom(bnd_vd).atom_type_index() ] );
			if ( bnd_atype.is_polar_hydrogen() ) {
				has_polar_H = true;
			}
		}
		heavyatom_has_polar_hydrogens_.push_back( has_polar_H );

	} // Iteration over atoms

	nbonds_ /= 2; // We double counted this (one for each direction)
}

void
ResidueType::copy_other_info(
	MutableResidueType const & mrt,
	utility::vector1< core::Size > const & old_to_new
) {

	atom_2_residue_connection_map_.resize( natoms_ ); // Empty by default
	icoor_.reserve( natoms_ );
	atom_base_indices_.reserve( natoms_ );
	abase2_indices_.resize( natoms_, 0 ); // Defaults to zero
	chi_atoms_.resize( mrt.nchi() ); // Will be altered in-place
	is_proton_chi_.reserve( mrt.nchi() );
	proton_chis_.reserve( mrt.n_proton_chi() );
	chi_2_proton_chi_.reserve( mrt.nchi() );
	proton_chi_samples_.reserve( mrt.n_proton_chi() );
	proton_chi_extra_samples_.reserve( mrt.n_proton_chi() );
	chi_rotamers_.reserve( mrt.nchi() );
	nu_atoms_.resize( mrt.n_nus() ); // Will be altered in-place
	ring_atoms_.resize( mrt.n_rings() ); // Will be altered in-place
	ring_saturation_types_.reserve( mrt.n_rings() );
	low_ring_conformers_.reserve( mrt.n_rings() );
	lowest_ring_conformer_.reserve( mrt.n_rings() );
	mainchain_atoms_indices_.reserve( mrt.mainchain_atoms().size() );
	actcoord_atoms_indices_.reserve( mrt.actcoord_atoms().size() );
	atom_shadowed_indices_.resize( natoms_, 0 ); // Defaults to zero

	//////////////////////
	// Atom name mapping
	/////////////////////

	// Pull in the atom_aliases
	for ( auto const & pair: mrt.atom_aliases() ) {
		if ( atom_name_to_index_.count( pair.second ) == 1 ) { // If we have the target
			if ( atom_name_to_index_.count( pair.first ) == 0 ) { // And we don't already have the source
				atom_name_to_index_[ pair.first ] = atom_name_to_index_[ pair.second ];
			}
			// Also add the stripped version.
			std::string const & stripped = utility::stripped_whitespace( pair.first  );
			if ( atom_name_to_index_.count( stripped ) == 0 ) {
				atom_name_to_index_[ stripped ] = atom_name_to_index_[ pair.second ];
			}
		}
	}
	// TODO: This means that when we bounce back to the MutableResidueType we're "polluting" the atom_name_to_index_.

	// TODO: Should we also be treating the canonical atom aliases?
	// Or generally, what should the difference in handling between atom aliases and canonical ones be?

	// We've already built a partial list of name->index maps, but we need to add additional mappings
	// (these are residual ones from renames, etc.)
	for ( auto const & pair: mrt.get_atom_name_to_vd_map() ) {
		if ( atom_name_to_index_.count( pair.first ) == 0 ) { // Only take things which aren't in the map
			atom_name_to_index_[ pair.first ] = old_to_new.index( mrt.atom_index( pair.second ) );
		}
	}

	//////////////////////
	// Residue Connections
	//////////////////////

	// First set up residue connections (preserve order)
	lower_connect_id_ = mrt.lower_connect_id();
	upper_connect_id_ = mrt.upper_connect_id();
	n_polymeric_residue_connections_ = mrt.n_polymeric_residue_connections();
	n_non_polymeric_residue_connections_ = mrt.n_possible_residue_connections() - mrt.n_polymeric_residue_connections();

	for ( core::Size ii(1); ii <= mrt.n_possible_residue_connections(); ++ii ) {
		MutableResidueConnection const & mrc( mrt.residue_connection(ii) );
		Size atomno = old_to_new.index( mrt.atom_index( mrc.vertex() ) );
		MutableICoorRecord const & ic( mrc.icoor() );
		std::string conn_name;
		if ( ii == lower_connect_id_ ) {
			conn_name = "LOWER";
		} else if ( ii == upper_connect_id_ ) {
			conn_name = "UPPER";
		} else {
			conn_name = "CONN" + std::to_string( ii );
		}
		AtomICoor aicoor( conn_name, ic.phi(), ic.theta(), ic.d(), ic.stub_atom1(), ic.stub_atom2(), ic.stub_atom3(), *this );
		ResidueConnection rc( atomno, aicoor, ii );
		residue_connections_.push_back( rc );

		atom_2_residue_connection_map_[ atomno ].push_back( ii );
	}

	//////////////////////
	// ICoords
	//////////////////////

	// (As I understand it, for small vectors, linear search is faster than setting up a std::map)
	root_atom_ = old_to_new.index( mrt.atom_index( mrt.root_atom() ) );
	nbr_atom_ = old_to_new.index( mrt.atom_index( mrt.nbr_vertex() ) );
	nbr_radius_ = mrt.nbr_radius();

	for ( core::Size new_index(1); new_index <= old_to_new.size(); ++new_index ) {
		Size old_index = old_to_new[ new_index ];
		VD old_vd = mrt.all_atoms()[old_index];
		Atom const & atom( mrt.atom(old_vd ) );
		MutableICoorRecordCOP ic( atom.icoor() );
		if ( ic == nullptr ) {
			if ( mrt.natoms() == 1 ) {
				// Special case one-atom residue types (don't need fancy icoor for these.).
				ic = utility::pointer::make_shared< MutableICoorRecord >(0,0,0, atom_names_[ new_index ], atom_names_[ new_index ], atom_names_[ new_index ] );
			} else {
				utility_exit_with_message("When convering ResidueType '" + name() + "' the internal coordinates for atom '" + atom_names_[ new_index ] + "' aren't set.");
			}
		}
		// In order for this to work, this->atom_index(std::string) and this->n_possible_residue_connections() must be functioning
		AtomICoor aicoor( atom_names_[ new_index ], ic->phi(), ic->theta(), ic->d(), ic->stub_atom1(), ic->stub_atom2(), ic->stub_atom3(), *this );
		icoor_.push_back( aicoor );

		core::Size base( atom_index( ic->stub_atom1() ) );
		debug_assert( base != 0 );
		if ( base != new_index ) {
			atom_base_indices_.push_back( base );
		} else if ( natoms_ == 1 ) {
			// Single atom residue type
			atom_base_indices_.push_back( base );
		} else {
			// Root atom
			base = atom_index( ic->stub_atom2() );
			debug_assert( base != 0 );
			atom_base_indices_.push_back( base );
		}
	}

	// abase2, just for acceptor positions (others get 0)
	// needs base and bonded_neighbors_ set up
	for ( core::Size atm: accpt_pos_ ) {
		core::Size base = atom_base_indices_[ atm ];
		debug_assert( base != 0 );
		AtomIndices const & nbrs = bonded_neighbor_[atm];
		if ( nbrs.empty() ) {
			utility_exit_with_message( "failed to set abase2 for acceptor atom, it has no nbrs!" );
		} else if ( nbrs.size() == 1 ) {
			abase2_indices_[ atm ] = atom_base_indices_[ base ];
			//iwd  The first child of the root is root's atom_base.
			//iwd  But if that child has no children, it ends up as its own abase2.
			//iwd  So instead we use the second child of the parent,
			//iwd  which must exist if there are 3+ atoms in this tree.
			if ( abase2_indices_[ atm ] == atm ) {
				AtomIndices const & base_nbrs( bonded_neighbor_[base] );
				for ( core::Size bn: base_nbrs ) {
					if ( bn != atm ) {
						abase2_indices_[ atm ] = bn;
						break;
					}
				}
			}
		} else if ( nbrs[1] == base ) {
			abase2_indices_[ atm ] = nbrs[2];
		} else {
			abase2_indices_[ atm ] = nbrs[1];
		}
	}

	//////////////////////
	// Chis
	//////////////////////

	// Convert chi info
	for ( core::Size ii(1); ii <= mrt.nchi(); ++ii ) {
		utility::vector1< core::Size > & chi( chi_atoms_[ii] );
		debug_assert( mrt.chi_atom_vds(ii).size() == 4 );
		for ( VD atm: mrt.chi_atom_vds(ii) ) {
			Size atomno = old_to_new.index( mrt.atom_index( atm ) );
			chi.push_back( atomno );
		}
		is_proton_chi_.push_back( mrt.is_proton_chi(ii) );
		if ( mrt.is_proton_chi(ii) ) {
			proton_chis_.push_back( ii );
			proton_chi_samples_.push_back( mrt.proton_chi_samples_for_chi( ii ) );
			proton_chi_extra_samples_.push_back( mrt.proton_chi_extra_samples_for_chi( ii ) );
			debug_assert( proton_chi_samples_.size() == proton_chi_extra_samples_.size() );
			chi_2_proton_chi_.push_back( proton_chi_samples_.size() );
		} else {
			chi_2_proton_chi_.push_back( 0 );
		}
		chi_rotamers_.push_back( mrt.chi_rotamers(ii) );
	}

	//////////////////////
	// Nu & Ring info
	//////////////////////

	for ( core::Size ii(1); ii <= mrt.n_rings(); ++ii ) {
		utility::vector1< core::Size > & ring( ring_atoms_[ii] );
		for ( VD atm: mrt.ring_atoms(ii) ) {
			Size atomno = old_to_new.index( mrt.atom_index( atm ) );
			ring.push_back( atomno );
		}
		ring_saturation_types_.push_back( mrt.ring_saturation_type(ii) );
		low_ring_conformers_.push_back( mrt.low_ring_conformers(ii) );
		lowest_ring_conformer_.push_back( mrt.lowest_ring_conformer(ii) );
	}
	// ring_conformer_sets_ is set later on

	for ( core::Size ii(1); ii <= mrt.n_nus(); ++ii ) {
		utility::vector1< core::Size > & nu( nu_atoms_[ii] );
		for ( VD atm: mrt.nu_atoms(ii) ) {
			Size atomno = old_to_new.index( mrt.atom_index( atm ) );
			nu.push_back( atomno );
		}
	}
	auto_assign_nu_atoms(); // Autogenerate nus, if not already specified.

	//////////////////////
	// Mainchain
	//////////////////////

	if ( mrt.mainchain_atoms().empty() ) {
		define_mainchain_atoms();
	} else {
		for ( VD atm: mrt.mainchain_atoms() ) {
			Size atomno = old_to_new.index( mrt.atom_index( atm ) );
			mainchain_atoms_indices_.push_back( atomno );
		}
	}

	//////////////////////
	// Actcoord
	//////////////////////

	for ( VD atm: mrt.actcoord_atoms() ) {
		Size atomno = old_to_new.index( mrt.atom_index( atm ) );
		actcoord_atoms_indices_.push_back( atomno );
	}

	//////////////////////
	// Atom Shadow
	//////////////////////

	for ( auto const & pair: mrt.shadow_atoms() ) {
		Size atomno1 = old_to_new.index( mrt.atom_index( pair.first ) );
		Size atomno2 = old_to_new.index( mrt.atom_index( pair.second ) );
		atom_shadowed_indices_[ atomno1 ] = atomno2;
	}

}

void
ResidueType::initialize_derived_data() {

	core::Size natoms = ResidueType::natoms(); // Can't do virtual dispatch from constructor

	//////////////////////
	// Atom Path distances
	//////////////////////

	//Requires that nbrs()/bonded_neighbor_ be correct
	ObjexxFCL::FArray2D_int path_distances( get_residue_path_distances( *this ));
	path_distance_.resize( natoms );
	for ( Size ii = 1; ii <= natoms; ++ii ) {
		path_distance_[ ii ].resize( natoms );
		for ( Size jj = 1; jj <= natoms; ++jj ) {
			path_distance_[ ii ][ jj ] = path_distances( ii, jj );
		}
	}

	//////////////////////
	// Dihedral atom sets
	//////////////////////

	dihedral_atom_sets_.clear();
	dihedrals_for_atom_.resize( natoms );
	for ( Size ii = 1; ii <= natoms; ++ii ) dihedrals_for_atom_[ ii ].clear();

	// get for all pairs of atoms separated by 1 bond
	for ( Size central_atom1 = 1; central_atom1 < natoms; ++central_atom1 ) {
		for ( Size central_atom2 = central_atom1+1; central_atom2 <= natoms; ++central_atom2 ) {
			if ( path_distance_[ central_atom1 ][ central_atom2 ] == 1 ) {

				// get all atoms separated from central_atom1/2 by one bond that are not central_atom2/1
				utility::vector1< Size > ca1d1;
				utility::vector1< Size > ca2d1;

				// ca1
				for ( Size i = 1; i <= natoms; ++i ) {
					if ( ( path_distance_[ central_atom1 ][ i ] == 1 ) && ( i != central_atom2 ) ) {
						ca1d1.push_back( i );
					}
				}
				// ca2
				for ( Size i = 1; i <= natoms; ++i ) {
					if ( ( path_distance_[ central_atom2 ][ i ] == 1 ) && ( i != central_atom1 ) ) {
						ca2d1.push_back( i );
					}
				}

				// for each pair of dihedral angle start or end atoms create a dihedral angle using central atom
				for ( core::Size & terminal_atom1 : ca1d1 ) {
					for ( core::Size & terminal_atom2 : ca2d1 ) {
						dihedral_atom_set temp( terminal_atom1, central_atom1, central_atom2, terminal_atom2 );
						dihedral_atom_sets_.push_back( temp );
						Size const which_dihedral = dihedral_atom_sets_.size();
						dihedrals_for_atom_[ terminal_atom1 ].push_back( which_dihedral );
						dihedrals_for_atom_[   central_atom1 ].push_back( which_dihedral );
						dihedrals_for_atom_[   central_atom2 ].push_back( which_dihedral );
						dihedrals_for_atom_[ terminal_atom2 ].push_back( which_dihedral );
					}
				}

			}
		}
	}

	///////////////////////////////
	// Bond Angle atom sets
	///////////////////////////////

	bondangle_atom_sets_.clear();
	bondangles_for_atom_.resize( natoms );
	for ( Size ii = 1; ii <= natoms; ++ii ) bondangles_for_atom_[ ii ].clear();

	// iterate over all atoms that could be a central atom
	for ( Size central_atom = 1; central_atom <= natoms; ++central_atom ) {

		AtomIndices const & bonded_neighbors(bonded_neighbor(central_atom) );
		Size const num_bonded_neighbors( bonded_neighbors.size() );

		// create all possible combinations of branching atoms
		for ( Size i = 1; i < num_bonded_neighbors; ++i ) {
			for ( Size j = i+1; j <= num_bonded_neighbors; ++j ) {
				bondangle_atom_set temp( bonded_neighbors[i], central_atom, bonded_neighbors[j] );
				bondangle_atom_sets_.push_back( temp );
				Size const which_angle = bondangle_atom_sets_.size();
				bondangles_for_atom_[ bonded_neighbors[i] ].push_back( which_angle );
				bondangles_for_atom_[        central_atom ].push_back( which_angle );
				bondangles_for_atom_[ bonded_neighbors[j] ].push_back( which_angle );
			}
		}
	}

	///////////////////////////////
	// Atoms Close to Connections
	///////////////////////////////

	// Now for inter-residue connections, find the sets of atoms that are within one and within two bonds
	// of a residue connection point.  From these sets, all inter-residue bond angle and bond torsions may
	// be enumerated when evaluating residue pair energies.  Also compute the backwards mapping: a list for
	// each atom of the within-1-bond and within-2-bond sets that the atom is listed as being part of. These
	// lists are needed when evaluating atom derivatives wrt the bond dihedrals and angles.
	atoms_within_one_bond_of_a_residue_connection_.resize( residue_connections_.size() );
	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) atoms_within_one_bond_of_a_residue_connection_[ ii ].clear();

	within1bonds_sets_for_atom_.resize( natoms );
	for ( Size ii = 1; ii <= natoms; ++ii ) within1bonds_sets_for_atom_[ ii ].clear();

	atoms_within_two_bonds_of_a_residue_connection_.resize( residue_connections_.size() );
	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) atoms_within_two_bonds_of_a_residue_connection_[ ii ].clear();

	within2bonds_sets_for_atom_.resize( natoms );
	for ( Size ii = 1; ii <= natoms; ++ii ) within2bonds_sets_for_atom_[ ii ].clear();

	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) {
		Size const ii_resconn_atom = residue_connections_[ ii ].atomno();

		AtomIndices const & ii_bonded_neighbors(bonded_neighbor(ii_resconn_atom)  );
		Size const ii_num_bonded_neighbors( ii_bonded_neighbors.size() );

		for ( Size jj = 1; jj <= ii_num_bonded_neighbors; ++jj ) {
			Size const jj_atom = ii_bonded_neighbors[ jj ];

			// Record that ii_resconn_atom and jj_atom are within a single bond of residue connection ii.
			two_atom_set wi1( ii_resconn_atom, jj_atom );
			atoms_within_one_bond_of_a_residue_connection_[ ii ].push_back( wi1 );

			// For atoms ii_resconn_atom and jj_atom, mark residue connection ii as a
			// connection point the are within one bond of.
			Size const which_wi1 = atoms_within_one_bond_of_a_residue_connection_[ ii ].size();
			within1bonds_sets_for_atom_[ ii_resconn_atom ].push_back( std::make_pair( ii, which_wi1 ) );
			within1bonds_sets_for_atom_[ jj_atom ].push_back( std::make_pair( ii, which_wi1 ));

			AtomIndices const & jj_bonded_neighbors(bonded_neighbor(jj_atom)   );
			Size const jj_num_bonded_neighbors( jj_bonded_neighbors.size() );

			for ( Size kk = 1; kk <= jj_num_bonded_neighbors; ++kk ) {
				Size const kk_atom = jj_bonded_neighbors[ kk ];
				if ( kk_atom == ii_resconn_atom ) continue; // skip iiat->jjat->iiat

				three_atom_set wi2( ii_resconn_atom, jj_atom, kk_atom );
				atoms_within_two_bonds_of_a_residue_connection_[ ii ].push_back( wi2 );

				Size const which_wi2 = atoms_within_two_bonds_of_a_residue_connection_[ ii ].size();
				within2bonds_sets_for_atom_[ ii_resconn_atom ].push_back( std::make_pair( ii, which_wi2 ) );
				within2bonds_sets_for_atom_[ jj_atom ].push_back( std::make_pair( ii, which_wi2 ));
				within2bonds_sets_for_atom_[ kk_atom ].push_back( std::make_pair( ii, which_wi2 ));
			}
		}
	}

	////////////////////////////////////
	// Properties & Property-based info
	////////////////////////////////////

	// Ring conformer sets
	update_ring_conformer_sets();

	// Last controlling chi
	// Also update RNA_Info object, as that special-cases last_controlling_chi
	if ( properties().has_property( RNA ) || properties().has_property( TNA ) ) {
		rna_info_ = utility::pointer::make_shared< core::chemical::rna::RNA_Info >();
		//update_last_controlling_chi is treated separately for RNA case. Parin Sripakdeevong, June 26, 2011
		if ( nchi() >= 4 || properties().has_property( TNA ) ) {
			// safety against hypothetical RNA RTs without 4 chi AND vs.
			// premature finalize() calls.
			rna_info_->rna_update_last_controlling_chi( *this, last_controlling_chi_, atoms_last_controlled_by_chi_);
			rna_info_->update_derived_rna_data( *this );
		} else {
			TR.Debug << "ResidueType " << name() << " is an RNA with less than three chis - using fall-back behavior " << std::endl;
			update_last_controlling_chi(); // Fallback.
		}
	} else {
		update_last_controlling_chi();
	}

	////////////////////////////////////
	// Connection-dependent atoms
	////////////////////////////////////

	update_polymer_dependent_groups();
	update_nonpolymer_dependent_groups();
}

void
ResidueType::update_derived_data() {

	////////////////////////////////////
	// Properties Property-based info
	////////////////////////////////////

	// RM: Is this really needed, and if so, can/should it be generalized?
	// Set up some atom properties
	for ( Size ii = 1; ii <= ResidueType::natoms(); ++ii ) {
		debug_assert( atomic_properties_[ ii ] != nullptr );
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		atomic_properties_[ ii ]->set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( mm_name( ii ) == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( atom_is_aro_hydrogen( ii_bonded_neighbors[ jj ] ) ) {
					atomic_properties_[ ii ]->set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	////////////////////////////////////
	// Rama info
	////////////////////////////////////

	// Set the RamaPrePro potential name to be this residue's name.
	if ( ResidueType::is_base_type() ) {
		if ( !get_rama_prepro_map_file_name(false).empty() && get_rama_prepro_mainchain_torsion_potential_name(false).empty() ) {
			set_rama_prepro_mainchain_torsion_potential_name( name(), false );
		}
		if ( !get_rama_prepro_map_file_name(true).empty() && get_rama_prepro_mainchain_torsion_potential_name(true).empty() ) {
			set_rama_prepro_mainchain_torsion_potential_name( name(), true );
		}
	}

	////////////////////////////////////
	// Rotamer Libraries
	////////////////////////////////////

	update_ncaa_rotamer_library_specification_if_present();
}

/// If a ring has been added but no nu angles have been defined, automatically assign them.
/// @details  Don't even bother trying if there are multiple rings.
/// The creator of the .params file needs to add the definitions manually.
/// @author   Labonte <JWLabonte@jhu.edu>
void
ResidueType::auto_assign_nu_atoms()
{
	if ( n_rings() == 0 ) { return; }  // Nothing to do here.
	if ( n_nus() != 0 ) { return; }  // Nus are already set, hopefully correctly.
	if ( n_rings() > 1 ) {
		utility_exit_with_message( "ResidueType::auto_assign_nu_atoms: "
			"ResidueTypes with more than one ring must explicitly define their nu angles." );
	}

	AtomIndices const & ring_atoms( ring_atoms_[ 1 ] );
	Size const n_ring_atoms( ring_atoms.size() );
	nu_atoms_.resize( n_ring_atoms );

	AtomIndices nu_atoms( 4 );  //  exactly 4 atoms per torsion definition
	for ( core::uint i( 1 ); i <= n_ring_atoms; ++i ) {
		// In the conditional operations below, the value, if zero, becomes n_ring_atoms.
		// (Bummer, using GNU ?: conditional expressions throws a warning.  Now the math has to be done twice; lame.)
		nu_atoms[ 1 ] = ring_atoms[ ( i - 1 ) ? ( i - 1 ) : n_ring_atoms ];
		nu_atoms[ 2 ] = ring_atoms[ i ];
		nu_atoms[ 3 ] = ring_atoms[ ( ( i + 1 ) % n_ring_atoms ) ? ( ( i + 1 ) % n_ring_atoms ) : n_ring_atoms ];
		nu_atoms[ 4 ] = ring_atoms[ ( ( i + 2 ) % n_ring_atoms ) ? ( ( i + 2 ) % n_ring_atoms ) : n_ring_atoms ];
		nu_atoms_[ i ] = nu_atoms;
	}
}


/// @details Needs connections, bonded neighbors set
void
ResidueType::define_mainchain_atoms()
{
	if ( !is_polymer() ) {
		return; // Don't do anything if we're not a polymer
	}

	// Test that we're actually a decent polymer
	if ( !upper_connect_id() || !upper_connect_atom() ||
			!lower_connect_id() || !lower_connect_atom() ) {
		TR.Warning << "Residue " << name() << " claims it's a polymer, " <<
			"but it doesn't have the appropriate UPPER and LOWER connection points specified.  " <<
			"Set MAINCHAIN_ATOMS in the topology file to remove this warning." << std::endl;
		return;
	}

	Size upper_connect( upper_connect_atom() ), lower_connect( lower_connect_atom() );
	// Default main chain: defined by shortest path from LOWER to UPPER.
	// IMO, everyone should really explicitly define the main chain from the topology file.  ~Labonte
	ObjexxFCL::FArray2D_int path_dist_matrix( get_residue_path_distances( *this ) );
	uint atom( lower_connect );
	while ( atom != upper_connect ) {
		mainchain_atoms_indices_.push_back( atom );
		AtomIndices const & nbrs( bonded_neighbor( atom ) );
		int min_d( path_dist_matrix( atom, upper_connect ) );
		uint next_atom( atom );

		for ( uint i=1; i<= nbrs.size(); ++i ) {
			uint const nbr( nbrs[i] );
			if ( path_dist_matrix( nbr, upper_connect ) < min_d ) {
				min_d = path_dist_matrix( nbr, upper_connect );
				next_atom = nbr;
			}
		}
		debug_assert( next_atom != atom );
		atom = next_atom;
	}
	mainchain_atoms_indices_.push_back( upper_connect );

}

void
ResidueType::update_ring_conformer_sets() {
	// Assign (a) set(s) of possible ring conformations.
	if ( properties().has_property( CYCLIC ) ) {
		ring_conformer_sets_.resize( n_rings() );
		for ( uint i( 1 ); i <= n_rings(); ++i ) {
			ring_conformer_sets_[ i ] = utility::pointer::make_shared< rings::RingConformerSet >(
				ring_atoms_[ i ].size(), ring_saturation_types_[ i ],
				lowest_ring_conformer_[ i ], low_ring_conformers_[ i ] );
		}
	}
}

void
ResidueType::update_last_controlling_chi() {
	core::Size natoms = ResidueType::natoms(); // This can be called from the constructor, where we can't do virtual dispatch

	last_controlling_chi_.resize( natoms );
	std::fill( last_controlling_chi_.begin(), last_controlling_chi_.end(), 0 );

	/// 1. First we have to mark all the atoms who are direct descendants of the 3rd
	/// atom in each chi; this prevents the note_chi_controls_atom recursion from overwriting
	/// the last-controlling chi for atoms descending from a particular chi.
	for ( Size ii = 1; ii <= nchi(); ++ii ) {
		AtomIndices atoms(chi_atoms(ii));
		Size const iiat3 = atoms[3];//chi_atoms_[ ii ][ 3 ];
		// This may be unnecessary; I believe two atoms pair as each other's bases only at the mainchain.
		Size const iiat3base = atom_base(iiat3);
		AtomIndices const & ii_nbrs(bonded_neighbor(iiat3));
		for ( Size jj = 1; jj <= ii_nbrs.size(); ++jj ) {
			Size const jj_atom = ii_nbrs[ jj ];
			if ( atom_base(jj_atom) == iiat3 && iiat3base != jj_atom ) {
				last_controlling_chi_[ jj_atom ] = ii;
			}
		}
	}

	/// 2. Now, lets recurse through all the atoms that are not direct descendants
	/// of the 3rd atom in a chi.  E.g. chi2 in PHE controls several more atoms than
	/// just CD1 and CD2.
	for ( Size ii = nchi(); ii >= 1; --ii ) {
		/// Note children of atom 3 of chi_ii as being controlled by chi ii.
		AtomIndices atoms(chi_atoms(ii));
		Size const iiat3 = atoms[3];//chi_atoms_[ ii ][ 3 ];
		// This may be unnecessary; I believe two atoms pair as each other's bases only at the mainchain.
		Size const iiat3base = atom_base(iiat3);
		AtomIndices const & ii_nbrs(bonded_neighbor(iiat3)  );
		for ( Size jj = 1; jj <= ii_nbrs.size(); ++jj ) {
			Size const jj_atom = ii_nbrs[ jj ];
			if ( atom_base(jj_atom) == iiat3 && iiat3base != jj_atom ) {
				note_chi_controls_atom( ii, jj_atom );
			}
		}
	}

	/// Now compute the atoms_last_controlled_by_chi_ arrays.

	/// get ready to allocate space in the atoms_last_controlled_by_chi_ arrays
	utility::vector1< Size > natoms_for_chi( nchi(), 0 );
	for ( Size ii = 1; ii <= natoms; ++ii ) {
		if ( last_controlling_chi_[ ii ] != 0 ) {
			++natoms_for_chi[ last_controlling_chi_[ ii ] ];
		}
	}

	/// allocate space
	atoms_last_controlled_by_chi_.resize( nchi() );
	for ( Size ii = 1; ii <= nchi(); ++ii ) {
		atoms_last_controlled_by_chi_[ ii ].clear();
		atoms_last_controlled_by_chi_[ ii ].reserve( natoms_for_chi[ ii ] );
	}

	/// fill the arrays
	for ( Size ii = 1; ii <= natoms; ++ii ) {
		if ( last_controlling_chi_[ ii ] != 0 ) {
			atoms_last_controlled_by_chi_[ last_controlling_chi_[ ii ]].push_back( ii );
		}
	}

}

/// @details O(N) recursive algorithm for determining the last chi for each atom.
/// Each atom is visited at most twice.
void
ResidueType::note_chi_controls_atom( Size chi, Size atomno )
{
	/// This should never be called on the "root" atom or it will enter an infinite loop
	debug_assert(  atom_base(atomno) != atomno );

	/// End the recursion: this atom already has had it's last chi identified, and it's not
	/// the chi we're currently labeling atoms with.
	if ( last_controlling_chi_[ atomno ] != 0 && last_controlling_chi_[ atomno ] != chi ) return;

	last_controlling_chi_[ atomno ] = chi;

	AtomIndices const & nbrs(bonded_neighbor(atomno) );
	for ( Size ii = 1; ii <= nbrs.size(); ++ii ) {
		/// descend into atoms who list atomno as their parent;
		/// atom_base_ defines a tree except at the root, where
		/// atom_base_[ atom_base_[ ii ]] == ii
		if ( atom_base(nbrs[ii]) == atomno ) {
			note_chi_controls_atom( chi, nbrs[ ii ] );
		}
	}
}

/// @brief Determine which atoms are polymer bond-dependent.
/// @details Should only be called from update_derived_data() function.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ResidueType::update_polymer_dependent_groups() {
	core::Size const n_atom( ResidueType::natoms() ); // This can be called from the constructor, where we can't do virtual dispatch

	has_polymer_dependent_groups_ = false;
	atom_depends_on_lower_polymeric_connection_.resize( n_atom, false );
	atom_depends_on_upper_polymeric_connection_.resize( n_atom, false );

	if ( !is_polymer() ) return;

	// This logic scales with natoms^3, which isn't ideal, but it's only called once on ResidueType initialization.  Still, it prevents
	// the necessity of implementing a recursive function, which is what I'm trying to avoid:
	for ( core::Size attempts(1); attempts<=n_atom; ++attempts ) { //The maximum number of iterations needed should be n_atom, though it will be much lower in most cases.
		bool nothing_changed = true;
		for ( core::Size i(1); i<=n_atom; ++i ) {
			if ( icoor(i).depends_on_polymer_lower() || icoor(i).depends_on_a_true_index( atom_depends_on_lower_polymeric_connection_ ) ) {
				nothing_changed = false;
				has_polymer_dependent_groups_ = true;
				atom_depends_on_lower_polymeric_connection_[i] = true;
			}
			if ( icoor(i).depends_on_polymer_upper() || icoor(i).depends_on_a_true_index( atom_depends_on_upper_polymeric_connection_ ) ) {
				nothing_changed = false;
				has_polymer_dependent_groups_ = true;
				atom_depends_on_upper_polymeric_connection_[i] = true;
			}
		}
		if ( nothing_changed ) break;
	}
}

/// @brief Determine which atoms are nonpolymer bond-dependent.
/// @details Should only be called from update_derived_data() function.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ResidueType::update_nonpolymer_dependent_groups() {
	core::Size const n_atom( ResidueType::natoms() ); // This can be called from the constructor, where we can't do virtual dispatch
	atom_depends_on_connection_.resize( n_possible_residue_connections(), utility::vector1< bool >( n_atom, false ) );

	// This logic scales with natoms^3, which isn't ideal, but it's only called once on ResidueType initialization.  Still, it prevents
	// the necessity of implementing a recursive function, which is what I'm trying to avoid:
	for ( core::Size j(1), jmax(n_possible_residue_connections()); j<=jmax; ++j ) {
		for ( core::Size attempts(1); attempts<=n_atom; ++attempts ) { //The maximum number of iterations needed should be n_atom, though it will be much lower in most cases.
			bool nothing_changed = true;
			for ( core::Size i(1); i<=n_atom; ++i ) {
				if ( icoor(i).depends_on_residue_connection(j) || icoor(i).depends_on_a_true_index( atom_depends_on_connection_[j] ) ) {
					nothing_changed = false;
					atom_depends_on_connection_[j][i] = true;
				}
			}
			if ( nothing_changed ) break;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @brief If there is an NCAARotamerLibrarySpecification, ensure that the
/// rotamer backbone dependicies have been set. If they have not, set them to
/// all mainchain torsions except omega (the final, inter-residue torsion).
/// @author Vikram K. Mulligan (vmullig@uw.edu).
///////////////////////////////////////////////////////////////////////////////
void
ResidueType::update_ncaa_rotamer_library_specification_if_present() {
	if ( rotamer_library_specification() == nullptr ) return;

	using namespace core::chemical::rotamers;
	NCAARotamerLibrarySpecificationOP libspec(
		utility::pointer::dynamic_pointer_cast< NCAARotamerLibrarySpecification >(
		rotamer_library_specification_nonconst()
		)
	);

	if ( libspec == nullptr ) return; //This isn't an NCAARotamerLibrarySpecification

	if ( libspec->rotamer_bb_torsion_indices().size() != 0 ) return; //Nothing to do -- indices have already been set.

	//Add all mainchain torsion indices EXCEPT omega:
	for ( core::Size i(1), imax(mainchain_atoms_indices_.size()); i<imax; ++i ) {
		libspec->add_rotamer_bb_torsion_index(i);
	}
}

void
ResidueType::self_pointer_updates() {
	if ( properties().has_property( CARBOHYDRATE ) ) {
		carbohydrate_info_ =
			carbohydrates::CarbohydrateInfoOP( new carbohydrates::CarbohydrateInfo( get_self_weak_ptr() ) );
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @brief Final check of ResidueType data, called by constuctor
/// @details These checks are meant to be quick and low-expense, and are only
/// called on ResidueType creation, so they shouldn't generally add much to Rosetta
/// processing time.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
///////////////////////////////////////////////////////////////////////////////
void
ResidueType::perform_checks() const
{
	bool checkspass = true;
	std::stringstream msg;
	msg << "One or more internal errors have occurred in residue type setup for " << name() << " (" << name3() << ", " << name1() << ")"<< std::endl;

	for ( core::Size chino(1); chino <= nchi(); ++chino ) {
		// These are check which are made in Residue::set_chi
		AtomIndices const & chi_atms( chi_atoms( chino ) );
		if ( atom_base( chi_atms[3] ) != chi_atms[2] ) {
			msg << "In chi #" << chino << ", the base of the third atom (" << atom_name(chi_atms[3]) <<") is " << atom_name(atom_base( chi_atms[3] ));
			msg << ", rather than the second atom of the chi (" << atom_name(chi_atms[2]) << ")" << std::endl;
			checkspass=false;
		}
		if ( atom_base( chi_atms[4] ) != chi_atms[3]  ) {
			msg << "In chi #" << chino << ", the base of the fourth atom (" << atom_name(chi_atms[4]) <<") is " << atom_name(atom_base( chi_atms[4] ));
			msg << ", rather than the third atom of the chi (" << atom_name(chi_atms[3]) << ")" << std::endl;
			checkspass=false;
		}
	}

	if ( is_metal() && (1 > nheavyatoms_ || is_virtual(1) ) ) {
		msg << "A metal residue type " << name() << " has a non-metal atom as atom 1." << std::endl;
		checkspass=false;
	}

	if ( is_metalbinding() && get_metal_binding_atoms().size()==0 ) {
		msg << "A metal-binding residue " << name() << " has no metal binding atoms listed in its params file (PROPERTIES METALBINDING "
			"without METAL_BINDING_ATOMS list)." << std::endl;
		checkspass=false;
	} else if ( !is_metalbinding() && get_metal_binding_atoms().size()>0 ) {
		msg << "A residue " << name() << " that has not been declared as a metal-binding residue has metal binding atoms listed in its "
			"params file (METAL_BINDING_ATOMS list without PROPERTIES METALBINDING)." << std::endl;
		checkspass=false;
	}

	if ( properties().has_property( ALPHA_AA ) && properties().has_property( BETA_AA ) ) {
		msg << "Error!  A residue type " << name() << " specifies that it is both an alpha and a beta amino acid in its params file." <<
			std::endl;
		checkspass=false;
	}
	if ( properties().has_property( L_AA ) && properties().has_property( D_AA ) ) {
		msg << "Error!  A residue type " << name() << " specifies that it is both an L-amino acid and a D-amino acid in its params "
			"file." << std::endl;
		checkspass=false;
	}
	if ( (backbone_aa_raw() != core::chemical::aa_unk) && ! properties().has_property( ALPHA_AA ) ) {
		msg << "Error!  A residue type " << name() << " specifies a standard alpha amino acid to use as a template for backbone scoring"
			" (rama and p_aa_pp scoring functions) without specifying that it is itself an alpha amino acid "
			"(PROPERTIES ALPHA_AA)." << std::endl;
		checkspass=false;
	}
	if ( (na_analogue_raw() != core::chemical::aa_unp) && ( ! properties().has_property( RNA )  || properties().has_property( CANONICAL_NUCLEIC ) ) ) {

		msg << "Error!  A residue type " << name() << " specifies a standard nucleic acid to use as a fragment analogue"
			" but it is not itself an RNA residue OR it is a canonical RNA residue" << std::endl;
		checkspass=false;
	}

	for ( Size n = 1; n <= Hpol_index_.size(); n++ ) {
		if ( !Hpos_polar_.has_value( Hpol_index_[n] ) ) {
			msg << "Hpos_polar " << atom_name( Hpol_index()[n] ) << " not in Hpol_index!?" << std::endl;
			checkspass = false;
		}
	}

	for ( Size n = 1; n <= Hpos_polar_.size(); n++ ) {
		if ( !Hpol_index_.has_value( Hpos_polar_[n] ) ) {
			msg << "Hpol_index " << atom_name( Hpol_index()[n] ) << " not in Hpos_polar!?" << std::endl;
			checkspass = false;
		}
	}

	if ( vec_atom_types_.empty() ) {
		msg << "Atom types vector is empty!" << std::endl;
	}

	if ( !checkspass ) {
		utility_exit_with_message(msg.str());
	}

	return;
}

////////////////////////////////////////////////////////////////////
// Helper methods //////////////////////////////////////////////////

// Insertion operator (overloaded so that ResidueType can be "printed" in PyRosetta).
std::ostream &
operator<<(std::ostream & output, ResidueType const & object_to_output)
{
	object_to_output.show(output);
	return output;
}

} // chemical
} // core


#ifdef    SERIALIZATION

template< class Archive >
void
core::chemical::ResidueType::save( Archive & arc ) const {
	using namespace core::chemical;

	arc( cereal::base_class< core::chemical::ResidueTypeBase >( this ) );

	arc( CEREAL_NVP( natoms_ ) ); // Size
	arc( CEREAL_NVP( nbonds_ ) ); // Size
	arc( CEREAL_NVP( nheavyatoms_ ) ); // Size
	arc( CEREAL_NVP( n_hbond_acceptors_ ) ); // Size
	arc( CEREAL_NVP( n_hbond_donors_ ) ); // Size
	arc( CEREAL_NVP( n_backbone_heavyatoms_ ) ); // Size
	arc( CEREAL_NVP( first_sidechain_hydrogen_ ) ); // Size

	arc( CEREAL_NVP( atom_names_ ) ); // utility::vector1< std::string >
	arc( CEREAL_NVP( atom_name_to_index_ ) ); // std::map< std::string, core::Size >
	arc( CEREAL_NVP( vec_atom_types_ ) ); // utility::vector1< AtomType >
	arc( CEREAL_NVP( vec_atom_type_indexes_ ) ); // utility::vector1< core::Size >
	utility::vector1< core::Size > element_indexes;
	for ( ElementCOP const & ele: elements_ ) {
		if ( ele == nullptr ) {
			element_indexes.push_back(0);
		} else {
			core::Size index( ele->get_index() );
			debug_assert( index != 0 );
			debug_assert( ele == element_set()[ index ] );
			element_indexes.push_back( index );
		}
	}
	arc( element_indexes ); // EXEMPT elements_
	arc( CEREAL_NVP( mm_atom_type_index_ ) ); // utility::vector1< core::Size >
	arc( CEREAL_NVP( gasteiger_atom_type_index_ ) ); // utility::vector1< core::Size >
	arc( CEREAL_NVP( formal_charge_ ) ); // utility::vector1< int >
	arc( CEREAL_NVP( partial_charge_ ) ); // utility::vector1< core::Real >
	arc( CEREAL_NVP( ideal_xyz_ ) ); // utility::vector1< core::Vector >
	arc( CEREAL_NVP( atomic_properties_ ) ); // utility::vector1< AtomPropertiesOP >
	arc( CEREAL_NVP( icoor_ ) ); // utility::vector1< AtomICoor >

	arc( CEREAL_NVP( bonded_neighbor_ ) ); // utility::vector1<AtomIndices>
	arc( CEREAL_NVP( bonded_neighbor_type_ ) ); // utility::vector1<utility::vector1<BondName> >
	arc( CEREAL_NVP( bonded_neighbor_ringness_ ) ); // utility::vector1<utility::vector1<BondRingness> >
	arc( CEREAL_NVP( attached_H_begin_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( attached_H_end_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( dihedral_atom_sets_ ) ); // utility::vector1<dihedral_atom_set>
	arc( CEREAL_NVP( dihedrals_for_atom_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( bondangle_atom_sets_ ) ); // utility::vector1<bondangle_atom_set>
	arc( CEREAL_NVP( bondangles_for_atom_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( last_controlling_chi_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( atoms_last_controlled_by_chi_ ) ); // utility::vector1<AtomIndices>

	arc( CEREAL_NVP( bonded_orbitals_ ) ); // utility::vector1< utility::vector1<Size> >
	arc( CEREAL_NVP( atoms_with_orb_index_ ) ); // AtomIndices
	arc( CEREAL_NVP( Haro_index_ ) ); // AtomIndices
	arc( CEREAL_NVP( Hpol_index_ ) ); // AtomIndices
	arc( CEREAL_NVP( accpt_pos_ ) ); // AtomIndices
	arc( CEREAL_NVP( Hpos_polar_ ) ); // AtomIndices
	arc( CEREAL_NVP( Hpos_apolar_ ) ); // AtomIndices
	arc( CEREAL_NVP( accpt_pos_sc_ ) ); // AtomIndices
	arc( CEREAL_NVP( Hpos_polar_sc_ ) ); // AtomIndices
	arc( CEREAL_NVP( all_bb_atoms_ ) ); // AtomIndices
	arc( CEREAL_NVP( all_sc_atoms_ ) ); // AtomIndices
	arc( CEREAL_NVP( heavyatom_has_polar_hydrogens_ ) ); // utility::vector1< bool >

	arc( CEREAL_NVP( chi_atoms_ ) ); // utility::vector1< utility::vector1< core::Size > >
	arc( CEREAL_NVP( is_proton_chi_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( proton_chis_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( chi_2_proton_chi_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( proton_chi_samples_ ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( proton_chi_extra_samples_ ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( chi_rotamers_ ) ); // utility::vector1<utility::vector1<std::pair<Real, Real> > >
	arc( CEREAL_NVP( nu_atoms_ ) ); // utility::vector1< AtomIndices >
	arc( CEREAL_NVP( ring_atoms_ ) ); // utility::vector1< AtomIndices >
	arc( CEREAL_NVP( ring_saturation_types_ ) ); // utility::vector1< core::chemical::rings::RingSaturationType >
	arc( CEREAL_NVP( lowest_ring_conformer_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( low_ring_conformers_ ) ); // utility::vector1<utility::vector1<std::string> >
	arc( CEREAL_NVP( path_distance_ ) ); // utility::vector1<utility::vector1<int> >

	arc( CEREAL_NVP( root_atom_ ) ); // core::Size
	arc( CEREAL_NVP( nbr_atom_ ) ); // core::Size
	arc( CEREAL_NVP( nbr_radius_ ) ); // Real
	arc( CEREAL_NVP( mass_ ) ); // Real

	arc( CEREAL_NVP( residue_connections_ ) ); // utility::vector1< ResidueConnection >
	arc( CEREAL_NVP( atom_2_residue_connection_map_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( atoms_within_one_bond_of_a_residue_connection_ ) ); // utility::vector1<utility::vector1<two_atom_set> >
	arc( CEREAL_NVP( within1bonds_sets_for_atom_ ) ); // utility::vector1<utility::vector1<std::pair<Size, Size> > >
	arc( CEREAL_NVP( atoms_within_two_bonds_of_a_residue_connection_ ) ); // utility::vector1<utility::vector1<three_atom_set> >
	arc( CEREAL_NVP( within2bonds_sets_for_atom_ ) ); // utility::vector1<utility::vector1<std::pair<Size, Size> > >
	arc( CEREAL_NVP( lower_connect_id_ ) ); // Size
	arc( CEREAL_NVP( upper_connect_id_ ) ); // Size
	arc( CEREAL_NVP( n_non_polymeric_residue_connections_ ) ); // Size
	arc( CEREAL_NVP( n_polymeric_residue_connections_ ) ); // Size

	arc( CEREAL_NVP( base_type_cop_ ) ); // ResidueTypeCOP

	arc( CEREAL_NVP( atom_base_indices_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( abase2_indices_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( mainchain_atoms_indices_ ) ); // AtomIndices
	arc( CEREAL_NVP( actcoord_atoms_indices_ ) ); // AtomIndices
	arc( CEREAL_NVP( cut_bond_neighbor_indices_ ) ); // utility::vector1<AtomIndices>
	arc( CEREAL_NVP( atom_shadowed_indices_ ) ); // utility::vector1<Size>

	arc( CEREAL_NVP( has_polymer_dependent_groups_ ) ); //bool
	arc( CEREAL_NVP( atom_depends_on_lower_polymeric_connection_ ) ); //utility::vector1<bool>
	arc( CEREAL_NVP( atom_depends_on_upper_polymeric_connection_ ) ); //utility::vector1<bool>
	arc( CEREAL_NVP( atom_depends_on_connection_ ) ); //utility::vector1< utility::vector1< bool > >

	// Observers aren't being serialized - any observer on the remote side will have to reattach itself.
	// EXEMPT destruction_obs_hub_
	// EXEMPT destruction_obs_mutex_;

	// These are derived information which will be re-initialized from scratch on the other side.
	// EXEMPT rna_info_
	// EXEMPT carbohydrate_info_
	// EXEMPT ring_conformer_sets_
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ResidueType::load( Archive & arc ) {
	using namespace core::chemical;

	arc( cereal::base_class< core::chemical::ResidueTypeBase >( this ) );

	arc( natoms_ ); // Size
	arc( nbonds_ ); // Size
	arc( nheavyatoms_ ); // Size
	arc( n_hbond_acceptors_ ); // Size
	arc( n_hbond_donors_ ); // Size
	arc( n_backbone_heavyatoms_ ); // Size
	arc( first_sidechain_hydrogen_ ); // Size

	arc( atom_names_ ); // utility::vector1< std::string >
	arc( atom_name_to_index_ ); // std::map< std::string, core::Size >
	arc( vec_atom_types_ ); // utility::vector1< AtomType >
	arc( vec_atom_type_indexes_ ); // utility::vector1< core::Size >
	utility::vector1< core::Size > element_indexes;
	arc( element_indexes ); // EXEMPT elements_
	for ( core::Size index: element_indexes ) {
		if ( index == 0 ) {
			elements_.push_back( nullptr );
		} else {
			elements_.push_back( element_set()[ index ] );
		}
	}
	arc( mm_atom_type_index_ ); // utility::vector1< core::Size >
	arc( gasteiger_atom_type_index_ ); // utility::vector1< core::Size >
	arc( formal_charge_ ); // utility::vector1< int >
	arc( partial_charge_ ); // utility::vector1< core::Real >
	arc( ideal_xyz_ ); // utility::vector1< core::Vector >
	arc( atomic_properties_ ); // utility::vector1< AtomPropertiesOP >
	arc( icoor_ ); // utility::vector1< AtomICoor >

	arc( bonded_neighbor_ ); // utility::vector1<AtomIndices>
	arc( bonded_neighbor_type_ ); // utility::vector1<utility::vector1<BondName> >
	arc( bonded_neighbor_ringness_ ); // utility::vector1<utility::vector1<BondRingness> >
	arc( attached_H_begin_ ); // utility::vector1<Size>
	arc( attached_H_end_ ); // utility::vector1<Size>
	arc( dihedral_atom_sets_ ); // utility::vector1<dihedral_atom_set>
	arc( dihedrals_for_atom_ ); // utility::vector1<utility::vector1<Size> >
	arc( bondangle_atom_sets_ ); // utility::vector1<bondangle_atom_set>
	arc( bondangles_for_atom_ ); // utility::vector1<utility::vector1<Size> >
	arc( last_controlling_chi_ ); // utility::vector1<Size>
	arc( atoms_last_controlled_by_chi_ ); // utility::vector1<AtomIndices>

	arc( bonded_orbitals_ ); // utility::vector1< utility::vector1<Size> >
	arc( atoms_with_orb_index_ ); // AtomIndices
	arc( Haro_index_ ); // AtomIndices
	arc( Hpol_index_ ); // AtomIndices
	arc( accpt_pos_ ); // AtomIndices
	arc( Hpos_polar_ ); // AtomIndices
	arc( Hpos_apolar_ ); // AtomIndices
	arc( accpt_pos_sc_ ); // AtomIndices
	arc( Hpos_polar_sc_ ); // AtomIndices
	arc( all_bb_atoms_ ); // AtomIndices
	arc( all_sc_atoms_ ); // AtomIndices
	arc( heavyatom_has_polar_hydrogens_ ); // utility::vector1< bool >

	arc( chi_atoms_ ); // utility::vector1< utility::vector1< core::Size > >
	arc( is_proton_chi_ ); // utility::vector1<_Bool>
	arc( proton_chis_ ); // utility::vector1<Size>
	arc( chi_2_proton_chi_ ); // utility::vector1<Size>
	arc( proton_chi_samples_ ); // utility::vector1<utility::vector1<Real> >
	arc( proton_chi_extra_samples_ ); // utility::vector1<utility::vector1<Real> >
	arc( chi_rotamers_ ); // utility::vector1<utility::vector1<std::pair<Real, Real> > >
	arc( nu_atoms_ ); // utility::vector1< AtomIndices >
	arc( ring_atoms_ ); // utility::vector1< AtomIndices >
	arc( ring_saturation_types_ ); // utility::vector1< core::chemical::rings::RingSaturationType >
	arc( lowest_ring_conformer_ ); // utility::vector1<std::string>
	arc( low_ring_conformers_ ); // utility::vector1<utility::vector1<std::string> >
	arc( path_distance_ ); // utility::vector1<utility::vector1<int> >

	arc( root_atom_ ); // core::Size
	arc( nbr_atom_ ); // core::Size
	arc( nbr_radius_ ); // Real
	arc( mass_ ); // Real

	arc( residue_connections_ ); // utility::vector1< ResidueConnection >
	arc( atom_2_residue_connection_map_ ); // utility::vector1<utility::vector1<Size> >
	arc( atoms_within_one_bond_of_a_residue_connection_ ); // utility::vector1<utility::vector1<two_atom_set> >
	arc( within1bonds_sets_for_atom_ ); // utility::vector1<utility::vector1<std::pair<Size, Size> > >
	arc( atoms_within_two_bonds_of_a_residue_connection_ ); // utility::vector1<utility::vector1<three_atom_set> >
	arc( within2bonds_sets_for_atom_ ); // utility::vector1<utility::vector1<std::pair<Size, Size> > >
	arc( lower_connect_id_ ); // Size
	arc( upper_connect_id_ ); // Size
	arc( n_non_polymeric_residue_connections_ ); // Size
	arc( n_polymeric_residue_connections_ ); // Size

	arc( base_type_cop_ ); // ResidueTypeCOP

	arc( atom_base_indices_ ); // utility::vector1<Size>
	arc( abase2_indices_ ); // utility::vector1<Size>
	arc( mainchain_atoms_indices_ ); // AtomIndices
	arc( actcoord_atoms_indices_ ); // AtomIndices
	arc( cut_bond_neighbor_indices_ ); // utility::vector1<AtomIndices>
	arc( atom_shadowed_indices_ ); // utility::vector1<Size>

	arc( has_polymer_dependent_groups_ ); //bool
	arc( atom_depends_on_lower_polymeric_connection_ ); //utility::vector1<bool>
	arc( atom_depends_on_upper_polymeric_connection_ ); //utility::vector1<bool>
	arc( atom_depends_on_connection_ ); //utility::vector1< utility::vector1< bool > >

	// Observers aren't being serialized - any observer on the remote side will have to reattach itself.
	// EXEMPT destruction_obs_hub_
	// EXEMPT destruction_obs_mutex_;

	// Regenerate values we don't bother to pass across
	if ( properties().has_property( RNA ) || properties().has_property( TNA ) ) { // Matches logic in ResidueType::initialize_derived_data()
		rna_info_ = utility::pointer::make_shared< core::chemical::rna::RNA_Info >(); // EXEMPT rna_info_
		if ( nchi() >= 4 || properties().has_property( TNA ) ) {
			rna_info_->update_derived_rna_data( *this );
		}
	}
	update_ring_conformer_sets(); // EXEMPT ring_conformer_sets_
	self_pointer_updates(); // EXEMPT carbohydrate_info_
}
SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ResidueType );
CEREAL_REGISTER_TYPE( core::chemical::ResidueType )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_ResidueType )
#endif // SERIALIZATION
