// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MutableResidueType.cc
/// @brief Method definitions for MutableResidueType
/// @author
/// Phil Bradley
/// Rocco Moretti (rmorettiase@gmail.com)
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to rings, properties, lipids, carbohydrates, and other non-AAs)

// Unit headers
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/MutableResidueConnection.hh>

// Package Headers
#include <core/chemical/ResidueProperties.hh>

// Project Headers
#include <core/chemical/residue_support.hh>
#include <core/chemical/icoor_support.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/bond_support.hh>
#include <core/chemical/RestypeDestructionEvent.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>
#include <core/chemical/Orbital.hh>

#include <core/chemical/ResidueType.hh> // For conversion functions
#include <core/chemical/ResidueConnection.hh> // For conversion functions

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>
#include <utility/graph/ring_detection.hh>
#include <utility/pointer/memory.hh>

// External headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/conversions.hh>

// C++ headers
#include <algorithm>

#ifdef    SERIALIZATION
#include <core/chemical/ResidueGraphTypes.srlz.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Utility serialization headers
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

using namespace ObjexxFCL;

static basic::Tracer TR( "core.chemical.MutableResidueType" );

VD const MutableResidueType::null_vertex = boost::graph_traits<ResidueGraph>::null_vertex();

MutableResidueType::MutableResidueType() = default; // private, not deleted because of serialization

MutableResidueType::MutableResidueType(
	AtomTypeSetCOP atom_types,
	ElementSetCOP elements,
	MMAtomTypeSetCOP mm_atom_types,
	orbitals::OrbitalTypeSetCOP orbital_types
) : ResidueTypeBase(atom_types, elements, mm_atom_types, orbital_types)
{}

MutableResidueType::MutableResidueType( MutableResidueType const & residue_type ):
	ResidueTypeBase( residue_type ) // Needed for compiler warning
{
	*this = residue_type; // Also copies base class info.
}

MutableResidueType::MutableResidueType( ResidueType const & residue_type ):
	ResidueTypeBase( residue_type )
{
	copy_atom_info( residue_type );

	copy_other_info( residue_type );

	// If we're coming from a ResidueType, we should be valid, but double check this.
	debug_assert( validate_residue_type() );
}

MutableResidueType::~MutableResidueType() = default;

// Utility function for operator=
utility::vector1<VD>
remap_vect( utility::vector1<VD> const & vect, std::map<VD, VD> const & old_to_new ) {
	utility::vector1<VD> new_vect;
	for ( VD vd: vect ) {
		debug_assert( old_to_new.count( vd ) );
		new_vect.push_back( old_to_new.at( vd ) );
	}
	return new_vect;
}

/// @details This needs to be explicit, as we're copying over VDs.
MutableResidueType &
MutableResidueType::operator=( MutableResidueType const & residue_type )
{
	// Short-circuit self copies (I don't know how robust the below is to such things).
	if ( this == & residue_type ) { return *this; } // Short-circuit self copies

	ResidueTypeBase::operator=( residue_type ); // Copy base class data

	nbr_radius_ = residue_type.nbr_radius_;
	ring_saturation_types_ = residue_type.ring_saturation_types_;
	lowest_ring_conformer_ = residue_type.lowest_ring_conformer_;
	low_ring_conformers_ = residue_type.low_ring_conformers_;
	lower_connect_id_ = residue_type.lower_connect_id_;
	upper_connect_id_ = residue_type.upper_connect_id_;


	// Copying the graph will change the VDs.
	// As such, any reference to the old VDs must be updated to the new VDs
	graph_ = residue_type.graph_;

	// Create a mapping from the old VDs to the new ones.
	std::map<VD, VD> old_to_new;
	old_to_new[ MutableResidueType::null_vertex ] = MutableResidueType::null_vertex; // Null verticies in one are null verticies in the other.
	for (
			VIterPair vp = boost::vertices(graph_), old_vp= boost::vertices(residue_type.graph_);
			vp.first != vp.second;
			++vp.first, ++old_vp.first
			) {
		auto v_iter= vp.first;
		auto old_v_iter= old_vp.first;
		VD vd = *v_iter;
		VD old_vd = *old_v_iter;
		// We assume that boost::graph copy preserves the iteration order of verticies
		debug_assert( graph_[ vd ] == residue_type.graph_[old_vd]);
		old_to_new[old_vd] = vd; //Assuming the boost::graph copy preserves ordering within the vertices list
	}

	// If ever we add VD/ED info to the Atoms/Bonds, we need to update them here.

	atom_name_to_vd_.clear();
	for ( auto const & pair: residue_type.atom_name_to_vd_ ) {
		debug_assert( old_to_new.count( pair.second ) );
		atom_name_to_vd_[ pair.first ] = old_to_new[ pair.second ];
	}

	ordered_atoms_ = remap_vect( residue_type.ordered_atoms_, old_to_new );

	debug_assert( old_to_new.count( residue_type.root_atom_ ) );
	root_atom_ = old_to_new[ residue_type.root_atom_ ];
	debug_assert( old_to_new.count( residue_type.nbr_atom_ ) );
	nbr_atom_ = old_to_new[ residue_type.nbr_atom_ ];

	atom_shadowed_.clear();
	for ( auto const & pair : residue_type.atom_shadowed_ ) {
		debug_assert( old_to_new.count( pair.first ) && old_to_new.count( pair.second ) );
		atom_shadowed_[ old_to_new[pair.first] ] = old_to_new[pair.second];
	}

	mainchain_atoms_ = remap_vect( residue_type.mainchain_atoms_, old_to_new );

	chis_.clear();
	for ( MutableChiRecordOP const & chi: residue_type.chis_ ) {
		if ( chi != nullptr ) {
			chis_.push_back( utility::pointer::make_shared< MutableChiRecord >( *chi ) );
			chis_.back()->remap_atom_vds( old_to_new );
		} else {
			chis_.push_back( nullptr );
		}
	}

	nu_atoms_.clear();
	for ( auto const & subvect: residue_type.nu_atoms_ ) {
		nu_atoms_.push_back( remap_vect( subvect, old_to_new ) );
	}

	ring_atoms_.clear();
	for ( auto const & subvect: residue_type.ring_atoms_ ) {
		ring_atoms_.push_back( remap_vect( subvect, old_to_new ) );
	}

	// Copy the connections, and update the VDs.
	residue_connections_ = residue_type.residue_connections_;
	for ( auto & conn: residue_connections_ ) {
		conn.remap_atom_vds( old_to_new );
	}

	return *this;
}

//////////////////////////////////////////////////////////////////////////////

/// @brief make a copy
MutableResidueTypeOP
MutableResidueType::clone() const
{
	MutableResidueTypeOP rsd_ptr( new MutableResidueType( *this ) );
	return rsd_ptr;
}


/// @brief make a copy, but only with all the stuff needed by patch selectors needs to be filled.
MutableResidueTypeOP
MutableResidueType::placeholder_clone() const
{
	MutableResidueTypeOP rsd( new MutableResidueType( atom_type_set_ptr(), element_set_ptr(), mm_atom_types_ptr(), orbital_types_ptr() ) );
	rsd->name ( name() );
	rsd->name1( name1() );
	rsd->name3( name3() );
	rsd->base_name( base_name() );
	rsd->aa( aa() );
	rsd->interchangeability_group( interchangeability_group() );
	rsd->set_properties( utility::pointer::make_shared<ResidueProperties>( properties() ) );
	return rsd;
}

//////////////////////////////////////////////////////////////////////////////

Size
MutableResidueType::nheavyatoms() const {
	core::Size nheavy(0);
	for ( VD atm: ordered_atoms_ ) {
		if ( ! atom(atm).is_hydrogen() ) {
			++nheavy;
		}
	}
	return nheavy;
}

/// @brief Counts the number of virtual atoms and returns the count.
/// @details The virtual count is not stored in the residue type.  This count is performed on the fly, and
///can hurt performance if repeatedly carried out.  Not intended for use in large loops -- instead, call
///once and store the value.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
Size
MutableResidueType::n_virtual_atoms() const
{
	core::Size virtcount = 0;
	for ( core::Size ia=1, iamax=natoms(); ia<=iamax; ++ia ) {
		// Must call this function instead of is_virtual because we use n_virtual_atoms() in patching.
		// If we start using this a lot in other contexts, write n_virtual_atoms_finalizedonly.
		if ( atom_type_set()[ graph_[ ordered_atoms_[ ia ] ].atom_type_index() ].is_virtual() ) ++virtcount;
	}
	return virtcount;
}

//////////////////////////////////////////////////////////////////////////////

Atom & MutableResidueType::atom(std::string const & atom_name) {
	return graph_[ atom_name_to_vd_[atom_name] ];
}
Atom const & MutableResidueType::atom(std::string const & atom_name) const {
	auto found = atom_name_to_vd_.find( atom_name );
	debug_assert( found != atom_name_to_vd_.end());
	return graph_[  found->second ];
}
Atom & MutableResidueType::atom(VD const atom_vd) {
	debug_assert( has(atom_vd) );
	return graph_[ atom_vd ];
}
Atom const & MutableResidueType::atom(VD const atom_vd) const {
	debug_assert( has(atom_vd) );
	return graph_[ atom_vd ];
}

/// @brief Get atom name by vertex descriptor
std::string const &
MutableResidueType::atom_name( VD vd ) const
{
	static std::string const unknown("<bad_vd>");
	if ( ! has(vd) ) { // Don't fail if the VD isn't present here.
		return unknown;
	}
	return graph_[ vd ].name();
}

VD
MutableResidueType::atom_vertex( std::string const & name ) const{
	auto atom_name_to_vd_iter( atom_name_to_vd_.find( name ) );
	if ( atom_name_to_vd_iter == atom_name_to_vd_.end() ) {
		TR.Error << "atom name : '" << name << "' not available in residue " << name3() << std::endl;
		show_all_atom_names( TR.Error );
		TR.Error << std::endl;
		utility_exit_with_message("Atom name not found in MutableResidueType.");
	}

	return atom_name_to_vd_iter->second;
}

Size
MutableResidueType::atom_index( VD const & vd ) const
{
	return ordered_atoms_.index( vd ); // Will return zero if not found
}


AtomType const &
MutableResidueType::atom_type( VD const vd ) const
{
	PyAssert( (vd != ResidueGraph::null_vertex() ) && has( vd ),
		"MutableResidueType::atom_type( VD const vd ): vd is not in this MutableResidueType!");
	debug_assert( (vd != ResidueGraph::null_vertex()) &&
		( std::find(ordered_atoms_.begin(), ordered_atoms_.end(), vd) != ordered_atoms_.end()) );
	return atom_type_set()[ graph_[ vd ].atom_type_index() ];
}

//////////////////////////////////////////////////////////////////////////////

Bond &
MutableResidueType::bond(ED const ed) {
	return graph_[ ed ];
}
Bond const &
MutableResidueType::bond(ED const ed) const {
	return graph_[ ed ];
}

Bond &
MutableResidueType::bond(VD vd1, VD vd2) {
	ED ed;
	bool found;
	boost::tie( ed, found ) = boost::edge( vd1, vd2 , graph_ );
	if ( ! found ) {
		utility_exit_with_message( "Cannot find bond between " + atom_name(vd1) + " and " + atom_name(vd2) + " in residue " + name() );
	}
	return graph_[ ed ];
}

Bond const &
MutableResidueType::bond(VD vd1, VD vd2) const {
	ED ed;
	bool found;
	boost::tie( ed, found ) = boost::edge( vd1, vd2 , graph_ );
	if ( ! found ) {
		utility_exit_with_message( "Cannot find bond between " + atom_name(vd1) + " and " + atom_name(vd2) + " in residue " + name() );
	}
	return graph_[ ed ];
}

Bond &
MutableResidueType::bond(std::string const & atom1, std::string const & atom2) {
	return bond( atom_vertex( atom1 ), atom_vertex( atom2 ) );
}

Bond const &
MutableResidueType::bond(std::string const & atom1, std::string const & atom2) const {
	return bond( atom_vertex( atom1 ), atom_vertex( atom2 ) );
}

bool
MutableResidueType::atoms_are_bonded( std::string const & atom1, std::string const & atom2 ) const {
	ED ed;
	bool found;
	boost::tie( ed, found ) = boost::edge( atom_vertex( atom1 ), atom_vertex( atom2 ), graph_ );
	return found;
}

AdjacentIterPair
MutableResidueType::bonded_neighbor_iterators( VD const & atom ) const {
	return boost::adjacent_vertices( atom, graph_ );
}

utility::vector1< VD >
MutableResidueType::bonded_neighbors( VD const & atom ) const {
	utility::vector1< VD > nbrs;
	AdjacentIter itr, itr_end;
	for ( boost::tie(itr, itr_end) = boost::adjacent_vertices( atom, graph_ ); itr != itr_end; ++itr ) {
		nbrs.push_back( *itr );
	}
	return nbrs;
}

Size
MutableResidueType::number_bonded_hydrogens( VD atomvd ) const {
	return bonded_hydrogens( atomvd ).size();
}

utility::vector1< VD >
MutableResidueType::bonded_hydrogens( VD const & atm ) const {
	utility::vector1< VD > nbrs;
	AdjacentIter itr, itr_end;
	for ( boost::tie(itr, itr_end) = boost::adjacent_vertices( atm, graph_ ); itr != itr_end; ++itr ) {
		if ( atom( *itr ).is_hydrogen() ) {
			nbrs.push_back( *itr );
		}
	}
	return nbrs;
}

Size
MutableResidueType::number_bonded_heavyatoms( VD atomvd ) const {
	return bonded_heavyatoms( atomvd ).size();
}

utility::vector1< VD >
MutableResidueType::bonded_heavyatoms( VD const & atm ) const {
	utility::vector1< VD > nbrs;
	AdjacentIter itr, itr_end;
	for ( boost::tie(itr, itr_end) = boost::adjacent_vertices( atm, graph_ ); itr != itr_end; ++itr ) {
		if ( ! atom( *itr ).is_hydrogen() ) {
			nbrs.push_back( *itr );
		}
	}
	return nbrs;
}

// Chis ////////////////////////////////////////////////////////////////

Size
MutableResidueType::n_proton_chi() const {
	Size npchi = 0;
	for ( MutableChiRecordOP chi: chis_ ) {
		if ( chi != nullptr && chi->is_proton_chi() ) {
			++npchi;
		}
	}
	return npchi;
}

/// @brief Will return if the chi is currently valid.
bool
MutableResidueType::chi_valid( Size const chino ) const {
	if ( chino == 0 || chino > chis_.size() ) { return false; }
	if ( chis_[ chino ] == nullptr ) { return false; }
	if ( chis_[ chino ]->chi_atoms().size() != 4 ) { return false; }
	if ( !has(chis_[ chino ]->chi_atoms()[1]) ||
			!has(chis_[ chino ]->chi_atoms()[2]) ||
			!has(chis_[ chino ]->chi_atoms()[3]) ||
			!has(chis_[ chino ]->chi_atoms()[4]) ) {
		return false;
	}
	return true;
}

// Connections ////////////////////////////////////////////////////////////////

Size
MutableResidueType::n_possible_residue_connections() const {
	return residue_connections_.size();
}

// Get a ResidueConection.
MutableResidueConnection const &
MutableResidueType::residue_connection( Size const i ) const
{
	return residue_connections_[ i ];
}

MutableResidueConnection &
MutableResidueType::residue_connection( Size const i )
{
	return residue_connections_[ i ];
}

VD
MutableResidueType::residue_connect_atom( Size const resconn_id ) const {
	return residue_connections_[ resconn_id ].vertex();
}

Size
MutableResidueType::n_residue_connections_for_atom( VD const atomid ) const {
	Size nconn = 0;
	for ( auto const & mrc: residue_connections_ ) {
		if ( mrc.vertex() == atomid ) {
			++nconn;
		}
	}
	return nconn;
}

Size
MutableResidueType::add_residue_connection( std::string const & atom_name )
{
	residue_connections_.push_back( MutableResidueConnection( atom_vertex( atom_name ) ) );
	return residue_connections_.size();
}

// Note: This needs to be robust to removing a VD that's not in the restype
void
MutableResidueType::delete_residue_connection( VD atm ) {
	// Want to delete connections in reverse order, as we're deleting while we're iterating
	for ( core::Size ii(residue_connections_.size()); ii >= 1; --ii ) {
		if ( residue_connections_[ii].vertex() == atm ) {
			if ( ii == lower_connect_id_ ) { lower_connect_id_ = 0; }
			if ( ii == upper_connect_id_ ) { upper_connect_id_ = 0; }
			if ( ii < lower_connect_id_ ) { --lower_connect_id_; }
			if ( ii < upper_connect_id_ ) { --upper_connect_id_; }
			// Also need to update any ICOOR records which reference connections greater than or equal to the deleted.
			update_icoors_after_connection_deletion( ii );
			residue_connections_.erase( residue_connections_.begin() + (ii-1) );
		}
	}
}

// Lower
MutableResidueConnection const &
MutableResidueType::lower_connect() const
{
	debug_assert( lower_connect_id_ != 0 );
	return residue_connections_[ lower_connect_id_ ];
}

VD
MutableResidueType::lower_connect_atom() const {
	debug_assert( lower_connect_id_ != 0 );
	return residue_connections_[ lower_connect_id_ ].vertex();
}

// Set the atom which connects to the lower connection.
void
MutableResidueType::set_lower_connect_atom( std::string const & atm_name )
{
	if ( atm_name == "NONE" ) {
		if ( lower_connect_id_ != 0 ) {
			TR.Debug << "ERASING LOWER_CONNECT: " << lower_connect_id_ << " lcid: " << upper_connect_id_ << std::endl;
			auto to_erase( residue_connections_.begin() );
			to_erase += lower_connect_id_ - 1;
			residue_connections_.erase( to_erase );
			update_icoors_after_connection_deletion( lower_connect_id_ );
			if ( lower_connect_id_ < upper_connect_id_ ) { --upper_connect_id_; }
			lower_connect_id_ = 0;

		}
	} else {
		if ( lower_connect_id_ == 0 ) {
			MutableResidueConnection rc( atom_vertex( atm_name ) );
			residue_connections_.push_back( rc );
			lower_connect_id_ = residue_connections_.size();
		} else {
			residue_connections_[ lower_connect_id_ ].vertex( atom_vertex( atm_name ) );
		}
	}
}

// Upper
MutableResidueConnection const &
MutableResidueType::upper_connect() const
{
	debug_assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ];
}

VD
MutableResidueType::upper_connect_atom() const
{
	debug_assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ].vertex();
}

// Set the atom which connects to the upper connection.
void
MutableResidueType::set_upper_connect_atom( std::string const & atm_name )
{
	if ( atm_name == "NONE" ) {
		if ( upper_connect_id_ != 0 ) {
			TR.Debug << "ERASING UPPER_CONNECT: " << upper_connect_id_ << " lcid: " << lower_connect_id_  << std::endl;
			auto to_erase( residue_connections_.begin() );
			to_erase += upper_connect_id_ - 1;
			residue_connections_.erase( to_erase );
			if ( upper_connect_id_ < lower_connect_id_ ) { --lower_connect_id_; }
			update_icoors_after_connection_deletion( upper_connect_id_ );
			upper_connect_id_ = 0;
		}
	} else {
		if ( upper_connect_id_ == 0 ) {
			MutableResidueConnection rc( atom_vertex( atm_name ) );
			residue_connections_.push_back( rc );
			upper_connect_id_ = residue_connections_.size();
		} else {
			residue_connections_[ upper_connect_id_ ].vertex( atom_vertex( atm_name ) );
		}
	}
}

void
MutableResidueType::add_metapatch_connect( std::string const & atom ) {
	// Provide unique variant name
	// We have to do this or connections get dropped--not all variants get put
	// back in. This is worse than you think--because they DON'T get dropped by
	// the metal!
	debug_assert( has(atom) );

	using namespace numeric::conversions;
	std::string varname_base( "MP-" + atom + "-CONNECT" );
	std::string res_varname( varname_base );

	// The below is needed just to ask has_variant_type on a custom variant.
	enable_custom_variant_types();

	Size count=0;
	while ( true ) {
		if ( count > 20 ) {
			utility_exit_with_message( "Could not find a new VariantType for MutableResidueType: " + name() );
		}
		++count;
		if ( count == 1 ) {
			if ( ! has_variant_type( res_varname ) ) break;
		} else {
			res_varname = varname_base + utility::to_string( count );
			if ( ! has_variant_type( res_varname ) ) break;
		}
	}
	add_variant_type( res_varname );

	VD atm = atom_vertex( atom );
	if ( bonded_hydrogens( atm ).size() == 0 ) {
		Size const connid( add_residue_connection( atom ) );
		MutableICoorRecordCOP aicoor = icoor( atm );
		debug_assert( aicoor != nullptr );

		// These coordinates are generic.
		set_icoor( "CONN"+ObjexxFCL::string_of( connid ),
			3.14159,
			70.600000*3.14159/180.000000,
			1.37,
			atom,
			aicoor->stub_atom1(),
			aicoor->stub_atom2() );
	} else {
		VD proton_index = bonded_hydrogens( atm )[1];
		MutableICoorRecordCOP aicoor = icoor( proton_index );
		debug_assert( aicoor != nullptr );

		Size const connid( add_residue_connection( atom ) );
		set_icoor( "CONN"+ObjexxFCL::string_of( connid ),
			aicoor->phi()+radians(180.0),
			aicoor->theta(),
			1.37,
			aicoor->stub_atom1(),
			aicoor->stub_atom2(),
			aicoor->stub_atom3() );
	}
}

// Other Residue Functions ////////////////////////////////////////////////////////////////

void
MutableResidueType::set_backbone_heavyatom( std::string const & name )
{
	if ( !has(name) ) {
		utility_exit_with_message("Trying to set bb atom that does not exist in residuetype");
	}
	graph_[ atom_name_to_vd_[name] ].is_backbone( true );
}

bool
MutableResidueType::is_backbone_heavyatom( VD atom ) const {
	return graph_[ atom ].is_backbone();
}

utility::vector1< VD >
MutableResidueType::actcoord_atoms() const {
	utility::vector1< VD > actcoords;
	for ( VD atm: all_atoms() ) {
		if ( graph_[ atm ].is_actcoord() ) {
			actcoords.push_back( atm );
		}
	}
	return actcoords;
}

/// @details add an atom to the list for calculating actcoord center
void
MutableResidueType::add_actcoord_atom( std::string const & atm )
{
	TR.Trace << "adding act coord atom: " << name() << ' ' << atm << std::endl;
	atom( atm ).is_actcoord( true );
}

/// @brief Remove an atom from the list of act coord atoms
/// (used in patching when it kills the valence that is thus used)
/// @author Andrew Watkins (amw579@nyu.edu)
void
MutableResidueType::delete_actcoord_atom(
	std::string const & atom_name
) {
	if ( !has( atom_name ) ) {
		std::string message = "Error in removing act coord atom from residue type " + name3() + ". Atom " + atom_name + " was not found.";
		utility_exit_with_message(message);
	}
	atom( atom_name ).is_actcoord( false );
}


/// @brief Add an alias name for an atom.
void
MutableResidueType::add_atom_alias( std::string const & rosetta_atom, std::string const & alias ) {
	if ( ! has( rosetta_atom ) ) {
		utility_exit_with_message( "Unable to add atom alias for non-existent atom " + rosetta_atom );
	}
	std::string const stripped_alias( utility::stripped_whitespace( alias ) );
	if ( stripped_alias.size() == 0 ) {
		utility_exit_with_message( "Cannot alias atom name to empty or all whitespace string." );
	}
	if ( atom_aliases().count( stripped_alias ) ) {
		utility_exit_with_message( "Cannot add atom alias; ResidueType " + name() +
			" already has an alias named " + stripped_alias + " mapped to " + atom_aliases()[ stripped_alias ] );
	}
	if ( has( stripped_alias ) ) {
		utility_exit_with_message( "Cannot add atom alias; ResidueType " + name() +
			" already has an atom named " + stripped_alias );
	}

	atom_aliases()[ alias ] = rosetta_atom;
	atom_aliases()[ stripped_alias ] = rosetta_atom;
}

/// @brief store map of canonical name to atom alias
void
MutableResidueType::add_canonical_atom_alias( std::string const & rosetta_atom, std::string const & alias ) {

	if ( canonical_atom_aliases().count( rosetta_atom ) != 0 ) {
		utility_exit_with_message( "Cannot add atom alias mapping, residue type already has an atom or alias named "+rosetta_atom );
	}
	canonical_atom_aliases()[ rosetta_atom ] = alias;
}

/// @brief Remove a given alias name for an atom.
void
MutableResidueType::delete_atom_alias( std::string const & alias, bool error ) {
	if ( atom_aliases().count(alias) ) {
		atom_aliases().erase( atom_aliases().find(alias) );
	} else if ( error ) {
		utility_exit_with_message( "Cannot remove atom alias " + alias + " as it does not exist as an alias." );
	}
	std::string const stripped_alias( utility::stripped_whitespace( alias ) );
	if ( atom_aliases().count(stripped_alias) ) { // Double check, as it might not be there, or alias==stripped_alias
		atom_aliases().erase( atom_aliases().find(stripped_alias) );
	}
}


void
MutableResidueType::set_shadowing_atom(
	std::string const & atom,
	std::string const & atom_being_shadowed
)
{
	atom_shadowed_[ atom_vertex( atom ) ] = atom_vertex( atom_being_shadowed );
}

/// @brief Add an atom to the list of atoms that can potentially form a bond to a metal ion.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// TODO Should this be an annotation on the Atom?
void
MutableResidueType::add_metalbinding_atom (
	std::string const & atom_name
) {
	if ( !has(atom_name) ) {
		std::string message = "Error in adding metal-binding atom to residue type " + name3() + ". Atom " + atom_name + " was not found.";
		utility_exit_with_message(message);
	}
	metal_binding_atoms().push_back( atom_name ); //Store names rather than indices, since indices might change.
	return;
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Remove an atom from the list of atoms that can potentially form a bond to a metal ion
/// (used in patching when it kills the valence that is thus used)
/// @author Andrew Watkins (amw579@nyu.edu)
/// TODO Should this be an annotation on the Atom?
void
MutableResidueType::delete_metalbinding_atom(
	std::string const & atom_name
) {
	if ( !has( atom_name ) ) {
		std::string message = "Error in removing metal-binding atom from residue type " + name3() + ". Atom " + atom_name + " was not found.";
		utility_exit_with_message(message);
	}
	metal_binding_atoms().erase( std::remove( metal_binding_atoms().begin(), metal_binding_atoms().end(), atom_name ),
		metal_binding_atoms().end() );

	return;
}

/// @brief Is this ResidueTypeBase a base type?
/// @details Checks the base_type_cop_ pointer.  If it's null, this is assumed to be a base type.
bool
MutableResidueType::is_base_type() const {
	return bool( !base_type_cop_ );
}

ResidueTypeCOP
MutableResidueType::get_base_type_cop() const {
	return base_type_cop_; // Direct access -- don't access self.
}

/// @brief Reset the base type COP to be null.  This implies that this ResidueTypeBase is a base type.
///
void
MutableResidueType::reset_base_type_cop() {
	base_type_cop_ = nullptr;
}

/// @brief Set the base type COP.  This implies that this ResidueTypeBase is NOT a base type.
///
void
MutableResidueType::set_base_type_cop(
	ResidueTypeCOP new_base_type
) {
	base_type_cop_ = new_base_type;
}

// Modification Functions ////////////////////////////////////////////////////////////////

// Atom Modification ////////////////////////////////////////////////////////////////

/// @note this does not set xyz coordinates for the added atom
VD
MutableResidueType::add_atom(
	std::string const & atom_name,
	std::string const & atom_type_name,
	std::string const & mm_atom_type_name,
	Real const charge
) {

	if ( atom_name_to_vd_.find(atom_name) != atom_name_to_vd_.end() ) {
		utility_exit_with_message("Can't add atom named `" + atom_name + "` to ResidueType " + name() + " as it already has one with that name.");
	}
	if ( atom_name_to_vd_.find( stripped_whitespace(atom_name)) != atom_name_to_vd_.end() ) {
		utility_exit_with_message("Can't add atom named `" + stripped_whitespace(atom_name) + "` to ResidueType " + name() + " as it already has one with that name.");
	}

	Atom atom(
		atom_name,
		mm_atom_type_name,
		nullptr, // Element information will be set by set_atom_type() below.
		charge,
		Vector(0.0)
	);

	VD v = graph_.add_vertex( atom );
	atom_name_to_vd_[ atom_name ] = v;
	atom_name_to_vd_[ stripped_whitespace( atom_name ) ] = v;
	ordered_atoms_.push_back(v);
	debug_assert( boost::num_vertices(graph_) == ordered_atoms_.size() );

	// Needed to reset the internal state.
	set_atom_type( v, atom_type_name );

	return v;
}

VD
MutableResidueType::add_atom(
	std::string const & atom_name /* = "" */
) {
	VD v = graph_.add_vertex( Atom() );

	if ( atom_name.size() ) {
		if ( atom_name_to_vd_.find(atom_name) != atom_name_to_vd_.end() ) {
			utility_exit_with_message("Can't add atom named `" + atom_name + "` to ResidueType " + name() + " as it already has one with that name.");
		}
		if ( atom_name_to_vd_.find( stripped_whitespace(atom_name)) != atom_name_to_vd_.end() ) {
			utility_exit_with_message("Can't add atom named `" + stripped_whitespace(atom_name) + "` to ResidueType " + name() + " as it already has one with that name.");
		}
		graph_[v].name(atom_name);
		atom_name_to_vd_[ atom_name ] = v;
		atom_name_to_vd_[ stripped_whitespace( atom_name ) ] = v;
	}

	ordered_atoms_.push_back(v);
	debug_assert( boost::num_vertices(graph_) == ordered_atoms_.size() );

	core::chemical::regenerate_graph_vertex_index(graph_);
	return v;
}

VD
MutableResidueType::add_atom(Atom const & atom, MutableICoorRecord const & icoor){

	VD v = graph_.add_vertex(atom);
	std::string const & atom_name( atom.name() );
	if ( atom_name.size() ) {
		if ( atom_name_to_vd_.find(atom_name) != atom_name_to_vd_.end() ) {
			utility_exit_with_message("Can't add atom named `" + atom_name + "` to ResidueType " + name() + " as it already has one with that name.");
		}
		if ( atom_name_to_vd_.find( stripped_whitespace(atom_name)) != atom_name_to_vd_.end() ) {
			utility_exit_with_message("Can't add atom named `" + stripped_whitespace(atom_name) + "` to ResidueType " + name() + " as it already has one with that name.");
		}
		atom_name_to_vd_[ atom_name ] = v;
		atom_name_to_vd_[ stripped_whitespace( atom_name ) ] = v;
	}

	ordered_atoms_.push_back(v);
	debug_assert( boost::num_vertices(graph_) <= ordered_atoms_.size() );

	graph_[v].icoor( utility::pointer::make_shared< MutableICoorRecord >( icoor ) ); // set icoor for this atom

	core::chemical::regenerate_graph_vertex_index(graph_);
	return v;
}

void
MutableResidueType::delete_atom( std::string const & name )
{
	debug_assert( has( name ) );
	delete_atom( atom_vertex(name) );
}

void
MutableResidueType::delete_atom( VD atm )
{
	runtime_assert( has( atm ) );
	// Delete any references to this atom
	// Do this first so we don't attempt to use any deleted atom information

	if ( root_atom_ == atm ) {
		root_atom_ = MutableResidueType::null_vertex;
	}
	if ( nbr_atom_ == atm ) {
		nbr_atom_ = MutableResidueType::null_vertex;
	}
	// Modifying map when iterating can be problematic.
	for ( auto it = atom_name_to_vd_.begin(); it != atom_name_to_vd_.end(); /*none*/ ) {
		if ( it->second == atm ) {
			it = atom_name_to_vd_.erase( it ); // erase gives the next valid iterator.
		} else {
			++it;
		}
	}
	atom_shadowed_.erase( atm ); // Robust to not being presenti
	for ( auto it = atom_shadowed_.begin(); it != atom_shadowed_.end(); /*none*/ ) {
		if ( it->second == atm ) {
			it = atom_shadowed_.erase( it ); // erase gives the next valid iterator.
		} else {
			++it;
		}
	}
	ordered_atoms_.pop( atm ); // Will renumber atoms, but noone calling delete_atom() should depend on ordering.
	for ( auto it = atom_name_to_vd_.begin(); it != atom_name_to_vd_.end(); /*none*/ ) {
		if ( it->second == atm ) {
			it = atom_name_to_vd_.erase( it );
		} else {
			++it;
		}
	}
	mainchain_atoms_.pop( atm );
	// Iterate backwards, so we delete in reverse order.
	for ( core::Size ii(chis_.size()); ii >= 1; --ii ) {
		// Patching unfortunately relies on absolute Chi numbers,
		// so we can't actually delete this Chi - we just empty it.
		// This can be funky if we never actually set it again.
		if ( chis_[ii] != nullptr && chis_[ii]->chi_atoms().contains( atm ) ) {
			// TODO: Fix patching such that this is an actual error.
			// (That is, make it such that patching doesn't rely on absolute chi numbers.)
			TR.Debug << "Deleting atom on " << name() << " invalidates chi " << ii << std::endl;
			chis_[ii] = nullptr;
		}
	}
	for ( core::Size ii(nu_atoms_.size()); ii >= 1; --ii ) {
		if ( nu_atoms_[ii].contains( atm ) ) {
			delete_nu( ii );
		}
	}
	for ( core::Size ii(ring_atoms_.size()); ii >= 1; --ii ) {
		if ( ring_atoms_[ii].contains( atm ) ) {
			delete_ring( ii );
		}
	}

	delete_residue_connection( atm ); // Should be robust to atm not being in a connection.

	// Note that we don't attempt to delete any name-based references.
	// (We keep those around in case they re-appear.)
	// Is this what we want?

	graph_.clear_vertex(atm);
	graph_.remove_vertex(atm);
}

/// @brief Rename the atom, updating the ResidueType-internal data mapping
void
MutableResidueType::rename_atom( VD atm, std::string const & atom_name ) {
	debug_assert( has( atm ) );
	debug_assert(atom_name.empty() || atom_name_to_vd_.count(atom_name) == 0);
	debug_assert(atom_name.empty() || atom_name_to_vd_.count(utility::stripped_whitespace(atom_name)) == 0);

	graph_[atm].name(atom_name);
	if ( ! atom_name.empty() ) { // If we're blanking the name, we don't need to update name mapping
		atom_name_to_vd_[ atom_name ] = atm;
		atom_name_to_vd_[ stripped_whitespace( atom_name ) ] = atm;
	}
	// Note that the old atom names are still valid.
}

void
MutableResidueType::set_atom_type( VD atom, std::string const & atom_type_name )
{
	debug_assert( has( atom ) );
	Atom & a = graph_[ atom ];

	// Get the new AtomType and its index.
	// (includes internal check for invalid type name)
	if ( atom_type_name.empty() ) {
		TR.Debug << "Setting atom type on atom " << atom_name( atom ) << " to null." << std::endl;
		a.atom_type_index( 0 ); // Null it out
		// We keep the element, as without it it's harder to subsequently retype the atom
		return; // Nothing more to do
	}

	core::uint const atom_type_index = atom_type_set().atom_type_index( atom_type_name );
	AtomType const & atype = atom_type_set()[atom_type_index];

	// Set/update AtomType index.
	a.atom_type_index( atom_type_index );

	// Set/update element and mass.
	if ( element_set_ptr() ) {  // Be robust if elements_ isn't defined.
		std::string const & element_name( atom_type_set()[ atom_type_index ].element() );
		core::uint const element_index = element_set().element_index( element_name );
		a.element_type( element_set()[ element_index ] );
	} else {
		TR.Warning << "Elements set undefined." << std::endl;
	}

	// TODO Kill this redundancy??
	a.is_virtual( atype.is_virtual() );
}

/// @brief set gasteiger atom type
void
MutableResidueType::set_gasteiger_atom_type(
	std::string const & atom_name,
	std::string const & gasteiger_atom_type_name
) {
	set_gasteiger_atom_type( atom_vertex( atom_name ) ,gasteiger_atom_type_name);
}

/// @brief set gasteiger atom type
void
MutableResidueType::set_gasteiger_atom_type(
	VD atom,
	std::string const & gasteiger_atom_type_name
) {
	gasteiger::GasteigerAtomTypeDataCOP gasteiger_type;
	if ( gasteiger_atom_type_name == "" ) {
		gasteiger_type = nullptr;
	} else {
		if ( ! gasteiger_atom_typeset() ) {
			set_gasteiger_atom_typeset( ChemicalManager::get_instance()->gasteiger_atom_type_set() );
		}
		gasteiger_type = gasteiger_atom_typeset()->atom_type( gasteiger_atom_type_name );
	}
	Atom & a = graph_[ atom ];
	a.gasteiger_atom_type( gasteiger_type );
}

// Bond Modification ////////////////////////////////////////////////////////////////

/// @details add a bond between atom1 and atom2 and add a BondType object referencing the bond using the specified bondName
void
MutableResidueType::add_bond(std::string const & atom_name1, std::string const & atom_name2, BondName bondLabel /*=SingleBond*/)
{
	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
		std::string message = "add_bond: atom " + atom_name1 + " and/or " + atom_name2 + " don't exist!";
		utility_exit_with_message( message );
	}

	/////// Standard Version /////////

	VD v1( atom_vertex( atom_name1 ) );
	VD v2( atom_vertex( atom_name2 ) );

	add_bond( v1, v2, bondLabel );
}

void
MutableResidueType::delete_child_proton( std::string const & atom ) {
	std::string res_varname( "MP-" + atom + "-PRUNEH" );
	Size count = 0;
	enable_custom_variant_types();
	while ( true ) {
		if ( count > 20 ) {
			utility_exit_with_message( "Could not find a new VariantType for MutableResidueType: " + name() );
		}
		++count;
		if ( count == 1 ) {
			if ( ! has_variant_type( res_varname ) ) break;
		} else {
			std::string const new_res_varname = res_varname + utility::to_string( count );
			if ( ! has_variant_type( new_res_varname ) ) {
				res_varname = new_res_varname;
				break;
			}
		}
	}
	add_variant_type( res_varname );

	// AMW: It seems like when we "delete" a proton, or fail to do so and virt
	// instead, it doesn't keep track of it...
	VD atm = atom_vertex( atom );
	utility::vector1< VD > bond_h = bonded_hydrogens( atm );
	core::Size nhydrogens = bond_h.size();
	if ( nhydrogens == 0 ) {
		TR.Trace << "No bonded hydrogens at " << atom << " in " << name() << std::endl;
	} else {
		// delete last proton
		VD proton_index = bond_h.back();

		// If this proton is part of any chis, we need to delete them.
		// (Do this before deleting, as atom deletion invalidates the chi.)
		// RM: Ideally, the atom deletion would take care of this, but that messes up patching ...
		for ( core::Size ii(chis_.size()); ii >= 1; --ii ) { // Iterate backwards, so we delete in reverse order.
			if ( chis_[ii] != nullptr && chis_[ii]->chi_atoms().contains( proton_index ) ) {
				delete_chi( ii ); // Actually erases, doesn't just void out.
			}
		}

		// Delete
		TR.Trace << "Removing " << atom_name( proton_index ) << std::endl;
		delete_atom( proton_index );

		// If there is more than one proton, allow the remain proton to occupy other positions
		if ( nhydrogens > 1 ) {
			VD alt_proton_index = bond_h.front();
			MutableICoorRecordCOP aicoor = icoor( alt_proton_index );
			debug_assert( aicoor != nullptr );

			for ( Size ii = 1; ii <= nchi(); ++ii ) {
				if ( chi_atom_vds( ii )[ 4 ] != proton_index )  continue;
				utility::vector1< Real > dihedral_samples;
				for ( Size jj = 0; jj<nhydrogens; ++jj ) {
					dihedral_samples.push_back( fmod( aicoor->phi() + jj*(360.0/nhydrogens), 360.0) );
				}
				set_proton_chi( ii, dihedral_samples, utility::vector1< Real >() );
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////

// @details add a bond between atom1 and atom2 and add a BondType object referencing the bond using the specified bondName
void
MutableResidueType::add_bond(VD atom1, VD atom2, BondName bondLabel /*=SingleBond*/)
{
	if ( !has( atom1 ) || !has( atom2 ) ) {
		utility_exit_with_message( "add_bond: atoms don't exist!" );
	}

	// check if bond already exists...
	if ( boost::edge(atom1, atom2, graph_).second ) {
		utility_exit_with_message( "Can't add bonds between " + atom_name(atom1) + " and " + atom_name(atom2) + " in residue " + name() + ": a bond already exists!." );
	}

	ResidueGraph::edge_descriptor e_added;
	bool added;
	boost::tie(e_added, added) = graph_.add_edge( atom1, atom2, Bond(-1, bondLabel)); /// -1 means Bond distance not set here. This will be fixed in the future
	debug_assert(added);
}

// Change the bond type of the given bond from one type to another.
/// @author   Labonte <JWLabonte@jhu.edu>
void
MutableResidueType::change_bond_type(
	std::string const & atom_name1,
	std::string const & atom_name2,
	BondName const new_bond_label )
{
	if ( ! ( has( atom_name1 ) || has( atom_name2 ) ) ) {
		utility_exit_with_message( "Change_bond_type: atoms " + atom_name1 + " and " + atom_name2 + " don't exist in residue " + name() );
	}

	VD const vd_source = atom_vertex( atom_name1 );
	VD const vd_target = atom_vertex( atom_name2 );
	graph_.remove_edge( vd_source, vd_target );
	graph_.add_edge( vd_source, vd_target, Bond( -1, new_bond_label ) ); /// -1 means Bond distance not set here.
}

/// @brief Delete a bond between the two atoms.
/// @details Note that this might leave dangling atoms.
void
MutableResidueType::delete_bond(VD atom1, VD atom2) {
	ED edge;
	bool found;
	boost::tie(edge, found) = boost::edge(atom1, atom2, graph_ );
	if ( !found ) {
		utility_exit_with_message("Cannot delete non-existent bond."); // TODO: add checking for atom presense.
	}
	boost::remove_edge( edge, graph_ );
}


/// @details add a cut_bond between atom1 and atom2, which disallows an atom-tree connection,
///            though the atoms are really bonded.
void
MutableResidueType::add_cut_bond(
	std::string const & atom_name1,
	std::string const & atom_name2
)
{
	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
		std::string message = "add_cut_bond: atoms " + atom_name1 + " and " + atom_name2 + " don't exist!";
		utility_exit_with_message( message  );
	}

	VD const vd_source = atom_vertex( atom_name1 );
	VD const vd_target = atom_vertex( atom_name2 );

	std::pair< ED, bool> ed_present_pair( boost::edge(vd_source,vd_target,graph_) );
	if ( ! ed_present_pair.second ) {
		utility_exit_with_message( "Cut bond set for bond which does not exist " + atom_name1 + " -- " + atom_name2 );
	}

	graph_[ ed_present_pair.first ].cut_bond( true );
}

// Reset the bond distance to an atom whose internal coordinates have already been set.
/// @details Looks up the internal coordinates to build the given atom and then resets the bond distance, updating
/// the xyz coordinates afterward.\n
/// One cannot currently reset the bond distance of residue connections using this method.
/// @param   <atm>: the string name of the atom
/// @param   <d>: the new bond distance in angstroms
/// @note    This method is primarily useful for patch operations, which may need to change the hybridization of an
/// atom and thus the bond length from the atom to which it is attached.
/// @author  Labonte <JWLabonte@jhu.edu>
void
MutableResidueType::reset_bond_distance_to_atom( std::string const & atm, core::Distance const d )
{
	MutableICoorRecordCOP const & atm_ic( icoor( atom_vertex( atm ) ) );
	// There used to be a "is_internal()" check here. That's not as easy to do in the new scheme,
	// and I'm not sure how necessary it is.
	set_icoor( atm, atm_ic->phi(), atm_ic->theta(), d,
		atm_ic->stub_atom1(), atm_ic->stub_atom2(), atm_ic->stub_atom3(), true );
}


// Geometry ////////////////////////////////////////////////////////////////

/// @brief AtomICoord of an atom
MutableICoorRecordCOP
MutableResidueType::icoor( VD const atm ) const
{
	debug_assert( has(atm) );
	return graph_[ atm ].icoor();
}

/// @details set AtomICoor for an atom
///
/// will update the xyz coords as well if desired, useful inside a patching operation where new
/// atoms are being added.
void
MutableResidueType::set_icoor(
	std::string const & atm,
	Real const phi,
	Real const theta,
	Real const d,
	std::string const & stub_atom1,
	std::string const & stub_atom2,
	std::string const & stub_atom3,
	bool const update_xyz // = false
)
{
	MutableICoorRecordOP ic( new MutableICoorRecord(phi, theta, d, stub_atom1, stub_atom2, stub_atom3 ) );

	switch ( string_to_icoord_type( atm ) ) {
	case ICoordAtomIDType::INTERNAL :
		debug_assert( has( atm ) );
		set_icoor( atom_vertex( atm ), phi, theta, d, stub_atom1, stub_atom2, stub_atom3, update_xyz );
		break;
	case ICoordAtomIDType::CONNECT :
		residue_connections_[ get_connection_number( atm ) ].icoor( *ic );
		break;
	case ICoordAtomIDType::POLYMER_LOWER :
		debug_assert( lower_connect_id_ != 0 );
		residue_connections_[ lower_connect_id_ ].icoor( *ic );
		break;
	case ICoordAtomIDType::POLYMER_UPPER :
		debug_assert( upper_connect_id_ != 0 );
		residue_connections_[ upper_connect_id_ ].icoor( *ic );
		break;
	default :
		utility_exit_with_message( "unrecognized stub atom id type!" );
		break; //to silence warning
	}
}

void
MutableResidueType::set_icoor(
	VD const & atm,
	Real const phi,
	Real const theta,
	Real const d,
	VD const & stub_atom1,
	VD const & stub_atom2,
	VD const & stub_atom3,
	bool const update_xyz /*= false*/
) {
	debug_assert( has(atm) && has(stub_atom1) && has(stub_atom2) && has(stub_atom3) );

	// Slighly icky assumption that we can get decent atom names with the stubs
	set_icoor( atm, phi, theta, d, atom_name(stub_atom1), atom_name(stub_atom2), atom_name(stub_atom3), update_xyz );
}

void
MutableResidueType::set_icoor(
	VD const & atm,
	Real const phi,
	Real const theta,
	Real const d,
	std::string const & stub_atom1,
	std::string const & stub_atom2,
	std::string const & stub_atom3,
	bool const update_xyz // = false
)
{
	debug_assert( has(atm) );

	MutableICoorRecordOP ic( new MutableICoorRecord(phi, theta, d, stub_atom1, stub_atom2, stub_atom3) );

	if ( has( stub_atom1 ) && atm == atom_vertex( stub_atom1 ) ) {
		//Root atom case
		if ( root_atom_ != MutableResidueType::null_vertex && root_atom_ != atm ) {
			TR.Error << "Can't reset root. Was " << atom_name( root_atom_ ) << " Resetting to " << atom_name( atm ) << std::endl;
			utility_exit_with_message( "Attempted to inappropriately reset ICOOR root atom." );
		}
		root_atom_ = atm;
		// Now we can continue as normal.
	}

	atom( atm ).icoor( ic );

	if ( update_xyz ) {
		set_ideal_xyz( atm, ic->build( *this ) );
	}
}

void
MutableResidueType::clear_icoor() {
	root_atom_ = MutableResidueType::null_vertex;
	for ( VD atm: ordered_atoms_ ) {
		atom(atm).icoor( nullptr );
	}
}

/// @brief get vd of an atom's base atom
VD
MutableResidueType::atom_base( VD const atomvd ) const
{
	PyAssert( has(atomvd), "MutableResidueType::atom_base( VD const atomvd ): atomvd is not in this MutableResidueType!");

	MutableICoorRecordCOP icoor = atom( atomvd ).icoor();
	if ( icoor == nullptr ) {
		return MutableResidueType::null_vertex;
	}
	std::string const & name = icoor->stub_atom1();
	if ( ! has( name ) ) {
		return MutableResidueType::null_vertex;
	}
	VD base = atom_vertex( name );
	if ( base == atomvd && natoms() >= 2 ) {
		// Root atom case. The base for the root atom is the angle atom.
		std::string const & name2 = icoor->stub_atom2();
		if ( ! has( name2 ) ) {
			return MutableResidueType::null_vertex;
		}
		return atom_vertex( name2 );
	} else {
		return base;
	}
}

void
MutableResidueType::set_ideal_xyz(
	std::string const & atm,
	Vector const & xyz_in
)
{
	set_ideal_xyz(atom_vertex(atm), xyz_in);
}

void
MutableResidueType::set_ideal_xyz(
	VD atm,
	Vector const & xyz_in
)
{
	if ( ! has(atm) ) {
		utility_exit_with_message("Cannot set ideal coordinates for non-existent atom.");
	}
	Atom & a = graph_[ atm ];
	a.ideal_xyz( xyz_in );
}

// Chi ////////////////////////////////////////////////////////////////

// Add a chi (side-chain) angle defined by four atoms.
void
MutableResidueType::add_chi( Size const chino,
	VD atom1,
	VD atom2,
	VD atom3,
	VD atom4)
{
	if ( !has( atom1 ) || !has( atom2 ) ||
			!has( atom3 ) || !has( atom4 ) ) {
		utility_exit_with_message("MutableResidueType::add_chi: atoms don't exist!" );
	}

	if ( chis_.size() < chino ) {
		if ( chino != chis_.size() +1 ) {
			// If you're adding chis out-of-order, something has probably gotten messed up.
			TR.Warning << "When working with " << name() << " (potentially when patching) attempted to add a chi (" << chino << ") which isn't the next chi in sequence. This will be an error if the intermediate chis aren't filled in." << std::endl;
		}
		chis_.resize( chino );
	}

	chis_[chino] = utility::pointer::make_shared< MutableChiRecord >( atom1, atom2, atom3, atom4 );
}  // add_chi

void
MutableResidueType::add_chi(VD atom1,
	VD atom2,
	VD atom3,
	VD atom4)
{
	add_chi(nchi() + 1, atom1, atom2, atom3, atom4);
}

// Add a chi (side-chain) angle defined by four atoms.
void
MutableResidueType::add_chi(
	Size const chino,
	std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4)
{
	add_chi( chino,
		atom_vertex(atom_name1),
		atom_vertex(atom_name2),
		atom_vertex(atom_name3),
		atom_vertex(atom_name4));
}

// Add a chi (side-chain) angle defined by four atoms to the end of the list of chis.
/// @details This method is needed for combinatorial patching of MutableResidueTypes for which the number of chis is variable.
/// Its primary purpose is to be used with add_chi_rotamer_to_last_chi() that adds rotamer bins to the last chi in the
/// list.  In this way, a new chi can be added by a patch file and its rotamer bins set without needing to designate a
/// chi index.
/// @note    See also add_chi_rotamer_to_last_chi().
/// @author  Labonte
void
MutableResidueType::add_chi(std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4)
{
	add_chi(nchi() + 1, atom_name1, atom_name2, atom_name3, atom_name4);
}

// redefine a chi angle based on four atoms
/// @details This function is almost an exact copy of the add_chi function except that vector resizing does NOT occur.
/// It is needed for certain PTMs that affect proton chis (e.g., phosphorylation and sulfation).
/// @author Andy M. Chen (June 2009)
void
MutableResidueType::redefine_chi(
	Size const chino,
	std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4
)
{
	if ( !has( atom_name1 ) || !has( atom_name2 ) ||
			!has( atom_name3 ) || !has( atom_name4 ) ) {
		utility_exit_with_message("MutableResidueType::redefine_chi: atoms dont exist!" );
	}

	// Right now add_chi() redefines the chi if it's in range.
	add_chi( chino, atom_vertex( atom_name1 ), atom_vertex( atom_name2 ), atom_vertex( atom_name3 ), atom_vertex( atom_name4 ) );
} // redefine_chi

/// @details Describe proton behavior for residue type; where should rotamer samples be considered,
/// and if expanded rotamers are desired, what deviations from the original rotamer samples
/// should be included.
/// E.g. dihedral_samples of 60, -60, and 180 could have an extra_sample of
/// 20 which would produce rotamers at 40 60 & 80, -40 -60 & -80, and -160, 180 & 160.
/// Extra_samples at 10 and 20 would produce 15 different rotamer samples.
void
MutableResidueType::set_proton_chi(
	Size chino,
	utility::vector1< Real > const & dihedral_samples,
	utility::vector1< Real > const & extra_samples
)
{
	if ( chino > chis_.size() ) {
		utility_exit_with_message("Error setting proton chi: Chi to set as proton chi does not exist.");
	}

	chis_[ chino ]->set_proton_chi( dihedral_samples, extra_samples );
}

///////////////////////////////////////////////////////////////////////////////

// Add a rotamer bin for a given chi.
/// @details A rotamer bin has the mean and standard deviation.
void
MutableResidueType::add_chi_rotamer(
	Size const chino,
	Real const mean,
	Real const sdev
)
{
	if ( chino > chis_.size() ) {
		utility_exit_with_message("Can't add chi rotamer to chi which does not exist.");
	}
	chis_[chino]->add_chi_rotamer(mean, sdev);
}

// Adds a chi rotamer bin to the highest-indexed chi in the list of chis for this MutableResidueType.
/// @details This method is needed for combinatorial patching of MutableResidueTypes for which the number of chis is variable.
/// Its primary purpose is to be used with the overloaded version of add_chi() that adds a new chi to the end of the
/// list.  In this way, a new chi can be added by a patch file and its rotamer bins set without needing to designate a
/// chi index.
/// @note    See also add_chi().
/// @author  Labonte
void
MutableResidueType::add_chi_rotamer_to_last_chi(core::Angle const mean, core::Angle const sdev)
{
	add_chi_rotamer(nchi(), mean, sdev);
}

// Delete all of the chi rotamer bins from the specified chi for this MutableResidueType.
/// @details This method is useful if one has redefined a chi within a patch file such that the old rotamer bins need
/// to be regenerated.
/// @author  Labonte <JWLabonte@jhu.edu>
void
MutableResidueType::clear_chi_rotamers( core::uint const chi_no )
{
	chis_[ chi_no ]->clear_chi_rotamers();
}

// delete the terminal chi angle
/// @author Andrew M. Watkins (April 2015)
void
MutableResidueType::delete_terminal_chi(
)
{
	debug_assert( chis_.size() >= 1 );
	chis_.resize( chis_.size() -1 );
} // delete_terminal_chi

void
MutableResidueType::delete_chi( core::Size const chino ) {
	debug_assert( chino != 0 && chino <= chis_.size() );

	TR.Debug << "Deleting chi " << chino << " on ResidueType " << name() << std::endl;

	chis_.erase_index( chino );
}


// Add a nu (internal cyclic) angle defined by four atoms.
void
MutableResidueType::add_nu( core::uint const nu_index,
	std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4 )
{
	if ( ! has( atom_name1 ) || ! has( atom_name2 ) || ! has( atom_name3 ) || ! has( atom_name4 ) ) {
		utility_exit_with_message( "MutableResidueType::add_nu: Requested atoms don't exist in this MutableResidueType!" );
	}

	utility::vector1< VD > atoms;
	atoms.push_back( atom_vertex( atom_name1 ) );
	atoms.push_back( atom_vertex( atom_name2 ) );
	atoms.push_back( atom_vertex( atom_name3 ) );
	atoms.push_back( atom_vertex( atom_name4 ) );

	if ( nu_atoms_.size() < nu_index ) {
		nu_atoms_.resize( nu_index );
	}
	nu_atoms_[ nu_index ] = atoms;
}

void
MutableResidueType::delete_nu( core::uint const nu_index ) {
	debug_assert( nu_index != 0 && nu_index <= nu_atoms_.size() );
	TR.Debug << "Deleting nu " << nu_index << " on ResidueType " << name() << std::endl;
	nu_atoms_.erase_index( nu_index );
}

// Add a ring definition.
void
MutableResidueType::add_ring(
	core::uint const ring_num,
	utility::vector1< std::string > const & ring_atoms,
	core::chemical::rings::RingSaturationType const saturation_type
) {
	Size const ring_size( ring_atoms.size() );
	utility::vector1< VD > atoms( ring_size );
	for ( uint i( 1 ); i <= ring_size; ++i ) {
		if ( ! has( ring_atoms[ i ] ) ) {
			utility_exit_with_message( "MutableResidueType::add_ring: Requested atoms don't exist in this MutableResidueType!" );
		}
		atoms[ i ] = atom_vertex( ring_atoms[ i ] );
	}

	if ( ring_atoms_.size() < ring_num ) {
		ring_atoms_.resize( ring_num );
	}
	if ( ring_saturation_types_.size() < ring_num ) {
		ring_saturation_types_.resize( ring_num );
	}
	if ( low_ring_conformers_.size() < ring_num ) {
		low_ring_conformers_.resize( ring_num );
	}
	if ( lowest_ring_conformer_.size() < ring_num ) {
		lowest_ring_conformer_.resize( ring_num );
	}
	ring_atoms_[ ring_num ] = atoms;
	ring_saturation_types_[ ring_num ] = saturation_type;
}

// Set this cyclic residue's lowest-energy ring conformer for the nth ring by IUPAC name.
void
MutableResidueType::set_lowest_energy_ring_conformer( core::uint const ring_num, std::string const & conformer )
{
	lowest_ring_conformer_[ ring_num ] = conformer;
}

// Set this cyclic residue's low-energy ring conformers for the nth ring by IUPAC name.
void
MutableResidueType::set_low_energy_ring_conformers( core::uint const ring_num, utility::vector1< std::string > const & conformers )
{
	low_ring_conformers_[ ring_num ] = conformers;
}

void
MutableResidueType::delete_ring( core::uint const ring_index ) {
	debug_assert( ring_index != 0 && ring_index <= ring_atoms_.size() );
	debug_assert( ring_atoms_.size() == low_ring_conformers_.size() );
	debug_assert( ring_atoms_.size() == lowest_ring_conformer_.size() );
	debug_assert( ring_atoms_.size() == ring_saturation_types_.size() );
	TR.Debug << "Deleting ring " << ring_index << " on ResidueType " << name() << std::endl;
	ring_atoms_.erase_index( ring_index );
	low_ring_conformers_.erase_index( ring_index );
	lowest_ring_conformer_.erase_index( ring_index );
	ring_saturation_types_.erase_index( ring_index );
}

// Orbitals ////////////////////////////////////////////////////////////////

void
MutableResidueType::clear_orbitals()
{
	get_orbitals().clear();
	get_orbitals_index().clear();
	for ( VD atm: ordered_atoms_ ) {
		graph_[atm].set_bonded_orbitals( utility::vector1<core::Size>{} ); // Reset all the associated atom-orbital correspondence
	}
}

/// @note this does not set xyz coordinates for the added orbital but sets the index of the orbital and maps
/// it to the type of orbital.
void
MutableResidueType::add_orbital(
	std::string & orbital_name,
	std::string & orbital_type_name
) {

	// store the atom type
	// the next call will fail if the orbital type name is unrecognized
	orbitals::OrbitalTypeSetCOP orbital_types( orbital_types_ptr() );
	debug_assert( orbital_types != nullptr );
	Size type( orbital_types->orbital_type_index( orbital_type_name ) );

	// store the name
	get_orbitals().push_back(Orbital(orbital_name, type, Vector(0.0)));

	get_orbitals_index()[ orbital_name ] = n_orbitals();
	get_orbitals_index()[ utility::stripped_whitespace( orbital_name ) ] = n_orbitals();
}

///////////////////////////////////////////////////////////////////////////////
/// @brief add an orbital bond between an atom and an orbital.
/// @note NOTE!!!!! This is indexed based upon atoms, not orbitals. That means that in your params file
/// you must have the atom as the first and orbital as the second.
void
MutableResidueType::add_orbital_bond(
	std::string const & atom_name1,
	std::string const & orbital_name
)
{
	if ( !has( atom_name1 ) || !has_orbital( orbital_name ) ) {
		std::string message = "add_bond: atoms " + atom_name1 + " and " + orbital_name + " dont exist!";
		utility_exit_with_message( message  );
	}

	if ( atom_name_to_vd_.find(atom_name1) == atom_name_to_vd_.end() ) {
		utility_exit_with_message("atom_name: " + atom_name1 +" not found. Improper params file!");
	}

	graph_[ atom_vertex( atom_name1 ) ].add_bonded_orbital( orbital_index( orbital_name ) );
}

//set the orbital icoor data.
void
MutableResidueType::set_orbital_icoor_id(
	std::string const & orbital,
	Real const phi,
	Real const theta,
	Real const d,
	std::string const & s1,
	std::string const & s2,
	std::string const & s3
)
{
	Size orb_indx(orbital_index(orbital));
	orbitals::ICoorOrbitalData new_icoor( phi, theta, d, s1, s2, s3 );

	get_orbitals()[ orb_indx ].new_icoor( new_icoor );
}


// Show Functions /////////////////////////////////////////////////////////////////////

/// @author Labonte <JWLabonte@jhu.edu>
void
MutableResidueType::show( std::ostream & output, bool output_atomic_details ) const
{
	using namespace std;
	using namespace utility;

	output << name() << " (" << name3() << ", " << name1() << "):" << endl;

	output << "Base: " << base_name() << std::endl;

	properties().show( output );

	// For display purposes, we assume that all the atoms have a name

	output << " Main-chain atoms:";
	for ( VD atm: mainchain_atoms_ ) {
		output << ' ' << atom_name( atm );
	}
	output << endl;

	output << " Backbone atoms:  ";
	for ( VD atm: all_atoms() ) {
		if ( graph_[atm].is_backbone() ) {
			output << ' ' << atom_name( atm );
		}
	}
	output << endl;

	output << " Side-chain atoms:  ";
	for ( VD atm: all_atoms() ) {
		if ( ! graph_[atm].is_backbone() ) {
			output << ' ' << atom_name( atm );
		}
	}
	output << endl;

	if ( ! ring_atoms_.empty() ) {
		for ( auto const & ring: ring_atoms_ ) {
			output << " Ring atoms:  ";
			for ( VD atm: ring ) {
				output << ' ' << atom_name( atm );
			}
			output << endl;
		}
	}

	if ( output_atomic_details ) {
		output << " Atomic Details:" << endl;
		for ( VD atm: all_atoms() ) {
			output << "  Atom: ";
			atom( atm ).show( output );
		}
	}
}

void
MutableResidueType::dump_vd_info() const {
	TR << "Residue " << name() << std::endl;

	for ( core::Size ii(1); ii <= ordered_atoms_.size(); ++ii ) {
		TR << " atom index " << ii << ": " << ordered_atoms_[ii] << " In graph? " << has(ordered_atoms_[ii]) << std::endl;
	}
	TR << "-------------" << std::endl;
	VIterPair allverts( boost::vertices( graph_ ) );
	for ( auto iter(allverts.first); iter != allverts.second; ++iter ) {
		TR << " atom " << atom_name(*iter) << " vd: " << *iter << std::endl;
	}
	TR << "-------------" << std::endl;
}


void
MutableResidueType::show_all_atom_names( std::ostream & out ) const {

	out << "atom_name atom_index atom_vertex" << std::endl;
	for ( VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first ) {
		auto v_iter= vp.first;
		VD vd = *v_iter;
		Atom a = graph_[vd];
		// Can't use atom_index(), as if that fails it would call show_all_atom_names, leading to recursion.
		core::Size index(0);
		for ( core::Size ii(1); ii <= ordered_atoms_.size(); ++ii ) {
			if ( ordered_atoms_[ii] == vd ) {
				index = ii;
				break;
			}
		}
		out << "'" << a.name() << "' " << index << " " << vd << std::endl;
	}

}

// Large-scale recalculation functions /////////////////////////////////////////////////

void MutableResidueType::assign_neighbor_atom()
{
	if ( ordered_atoms_.empty() ) {
		utility_exit_with_message("Can't assign neighbor atom for " + name() + " -- has no atoms!" );
	}

	//calculate the geometric center of all atoms in the residue
	Vector total(0.0,0.0,0.0);
	for ( VD atm: ordered_atoms_ ) {
		total += graph_[atm].ideal_xyz();
	}

	Vector center = total/ordered_atoms_.size();

	//locate the atom which is closest to the center
	VD min_atm = MutableResidueType::null_vertex;
	core::Real min_distance = 50000.0;

	for ( VD atm: ordered_atoms_ ) {
		core::Real distance = center.distance(graph_[atm].ideal_xyz());
		if ( (distance < min_distance) && (!graph_[atm].is_hydrogen()) ) {
			min_distance = distance;
			min_atm = atm;
		}
	}
	debug_assert(min_atm != MutableResidueType::null_vertex);
	//set neighbor atom
	nbr_atom( min_atm );
}

void MutableResidueType::assign_internal_coordinates()
{
	// Reuse the existing root, or failing that, the neighbor atom
	// As a last resort, just use atom #1
	VD new_root( root_atom_ );
	if ( new_root == MutableResidueType::null_vertex || ! has(new_root) ) {
		new_root = nbr_atom_;
	}
	if ( new_root == MutableResidueType::null_vertex || ! has(new_root) ) {
		// Make sure we're not assigning a null vertex to the root.
		for ( VD oa: ordered_atoms_ ) {
			if ( oa != MutableResidueType::null_vertex && has(oa) ) {
				new_root = oa;
				break;
			}
		}
	}
	debug_assert( new_root != MutableResidueType::null_vertex );
	debug_assert( has( new_root ) );
	assign_internal_coordinates( new_root );
}

void MutableResidueType::assign_internal_coordinates(core::chemical::VD new_root)
{
	//%TODO: right now we're ignoring M FRAG lines and M SPLT lines in molfiles
	if ( n_possible_residue_connections() != 0 ) {
		TR.Error << "Residue " << name() << " has connections - can't assign internal coordinates.";
		utility_exit_with_message("Cannot currently assign internal coordinates for polymeric residue.");
	}
	debug_assert( new_root != MutableResidueType::null_vertex );
	// Reset the root atom so we can re-root the tree
	root_atom_ = MutableResidueType::null_vertex;
	reroot_restype(*this, graph_, new_root);
}

void
MutableResidueType::fill_ideal_xyz_from_icoor() {
	core::chemical::fill_ideal_xyz_from_icoor(*this, graph_);
}

/// @brief Regenerate the rotatable chi bonds from the internal graph structure.
/// If the number of proton chi samples would exceed max_proton_chi_samples, don't add extra sampling to proton chis.
/// As a special case, if this is zero don't add any proton chi sampling at all.
///
/// Requires that Icoor and atom base records are up-to-date, and that ring bonds have been annotated.
void
MutableResidueType::autodetermine_chi_bonds( core::Size max_proton_chi_samples ) {
	// First, clear off the current internal information on rotatable chis;
	chis_.clear();
	// TODO: Notate if we have a "non-standard" setting for any of the proton chis

	// As far as I can tell, there isn't any specified ordering for chis.
	// The canonical residues go from root out, but for ligands
	// there doesn't look to be any ordering guarantee.
	utility::vector1<VDs> found_chis( core::chemical::find_chi_bonds( *this ) );
	utility::vector1< core::Size > proton_chis; // Not the member variable as set_proton_chi modifies that.

	if ( is_protein() ) {

		utility::vector1<VDs> true_chis; // filtered and ordered from found_chis.
		// Note that this algorithm to get down to the 'real' chis is pretty
		// gross, but when N is < 10 most reasonable big-Os are fine, right?

		// Step 1. Get chi1 (it's the one with N as first or fourth)
		// Other criterion -- atom 3 can't be C (it'll find N CA C O)
		for ( VDs const & chi : found_chis ) {
			TR.Trace << "looking at found chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
			if ( atom_name( chi[ 1 ] ) == "N" && atom_name( chi[ 3 ] ) != "C" ) {
				true_chis.push_back( chi );

				break;
			}
		}

		// Step 2. Get remainder of chis by asking each one to start with the
		// second atom of the prior chi[s]. Note that this will potentially
		// confuse branches, but branched sidechains with lots of chis are treated
		// poorly by essentially any chi system.
		std::string target_first_atom;
		if ( true_chis.size() > 0 )  target_first_atom = atom_name( true_chis[ 1 ][ 2 ] );
		while ( true ) {

			// this extra loop is to future-proof a bit against branching: multiple
			// chis per pass may start with the target_first_atom and therefore
			// we don't want to update it right away.

			std::string candidate_new_atom = target_first_atom;
			for ( VDs const & found_chi : found_chis ) {
				TR.Trace << "looking at found chi: " << atom_name( found_chi[1] ) << " " << atom_name( found_chi[2] ) << " " << atom_name( found_chi[3] ) << " " << atom_name( found_chi[4] ) << std::endl;
				if ( atom_name( found_chi[ 1 ] ) == target_first_atom ) {
					true_chis.push_back( found_chi );
					candidate_new_atom = atom_name( found_chi[ 2 ] );
				}

			}
			if ( candidate_new_atom == target_first_atom ) break;

			// This may have to become a vector -- where we accumulated many
			// candidate_new_atom -- later.
			target_first_atom = candidate_new_atom;
		}
		for ( VDs const & true_chi : true_chis ) {
			debug_assert( true_chi.size() == 4 );
			TR.Debug << "looking at true chi: " << atom_name( true_chi[1] ) << " " << atom_name( true_chi[2] ) << " " << atom_name( true_chi[3] ) << " " << atom_name( true_chi[4] ) << std::endl;
			add_chi( true_chi[1], true_chi[2], true_chi[3], true_chi[4] );
			if ( atom( true_chi[4] ).element_type()->element() == core::chemical::element::H ) {
				// proton chi
				proton_chis.push_back( nchi() );
			}
		} // for all found chis

		// TODO: This is probably out of place here - we should move it to a more general location
		core::chemical::annotate_backbone( *this );

	} else if ( is_RNA() ) {
		utility::vector1<VDs> true_chis; // filtered and ordered from found_chis.

		//CHI 1 C2' C1' N9  C4
		//CHI 2 C4' C3' C2' C1'
		//CHI 3 C3' C2' C1' N9
		//CHI 4 C3' C2' O2' HO2'
		// First base atom is either N1 or N9

		// Will never actually remain this value -- but we need
		// to know that we're not using uninitialized, and we
		// need a value we can always initialize to for RNA rsd.
		// (This is tough just because of thenature of VDs.)
		// So we use a runtime_assert after the loop.
		VD first_base_atom = atom_vertex( "P" );
		for ( VDs const & chi : found_chis ) {
			TR.Trace << "looking at found chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
			if ( atom_name( chi[ 1 ] ) == "C2'" && atom_name( chi[ 2 ] ) != "O2'" ) {
				true_chis.push_back( chi );
				first_base_atom = chi[3];
				// Third atom of chi1 is first base atom... but not first SC atom (formally)
				// setting it that way is HELL for hbond types
				break;
			}
		}
		runtime_assert( true_chis.size() == 0 || first_base_atom != atom_vertex( "P" ) );

		// Step 2. Hard-fix three chis: two rings, and proton chi for HO2'.
		VDs new_chi{atom_vertex("C4'"), atom_vertex("C3'"), atom_vertex("C2'"), atom_vertex("C1'")};
		true_chis.emplace_back( new_chi );
		// Skip this chi for N
		if ( first_base_atom != atom_vertex( "P" ) ) {
			new_chi = VDs{ atom_vertex("C3'"), atom_vertex("C2'"), atom_vertex("C1'"), first_base_atom };
		}
		//true_chis.emplace_back( { atom_vertex("C4'"), atom_vertex("C3'"), atom_vertex("C2'"), atom_vertex("C1'") } );
		true_chis.emplace_back( new_chi );
		// What to do absent HO2'?
		// answer: whatever else O2' is bonded to that's not C2'
		if ( has( "HO2'" ) ) {
			//chi = ;
			true_chis.emplace_back( VDs{ atom_vertex("C3'"), atom_vertex("C2'"), atom_vertex("O2'"), atom_vertex("HO2'") } );

		} else {
			// actually essential for pdb_T38 for example.
			// Oh: we also *must* figure out why atom ordering matters so much for the backbone
			// sidechain distinction (rather than reordering after a Correct grouping) and how to
			// do the latter. Right now I have to relabel...
			for ( VDs const & found_chi : found_chis ) {
				if ( atom_name( found_chi[ 3 ] ) == "O2'" && atom_name( found_chi[ 2 ] ) == "C2'" ) {
					true_chis.emplace_back( found_chi );
				}
			}
		}

		// Theoretical final step (AMW TODO): add all chis that are children of the base
		// or of O2', in a protein-y way, as chis 5+.
		std::string target_first_atom;
		if ( true_chis.size() > 0 ) target_first_atom = atom_name( true_chis[ 1 ][ 2 ] );

		while ( true ) {

			// this extra loop is to future-proof a bit against branching: multiple
			// chis per pass may start with the target_first_atom and therefore
			// we don't want to update it right away.

			std::string candidate_new_atom = target_first_atom;
			for ( VDs const & chi : found_chis ) {
				TR.Trace << "looking at found chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
				if ( atom_name( chi[ 1 ] ) == target_first_atom ) {
					true_chis.push_back( chi );
					candidate_new_atom = atom_name( chi[ 2 ] );
				}
			}
			if ( candidate_new_atom == target_first_atom ) break;

			// This may have to become a vector -- where we accumulated many
			// candidate_new_atom -- later.
			target_first_atom = candidate_new_atom;
		}

		for ( VDs const & chi : true_chis ) {
			TR.Debug << "looking at true chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
			debug_assert( chi.size() == 4 );
			add_chi( chi[1], chi[2], chi[3], chi[4] );
			if ( atom( chi[4] ).element_type()->element() == core::chemical::element::H ) {
				// proton chi
				proton_chis.push_back( nchi() );
			}
		} // for all found chis

		// TODO: This is probably out of place here - we should move it to a more general location
		core::chemical::annotate_backbone( *this );

	} else {
		// ligand logic: far simpler.

		for ( VDs const & chi : found_chis ) {
			debug_assert( chi.size() == 4 );
			add_chi( chi[1], chi[2], chi[3], chi[4] );
			if ( atom( chi[4] ).element_type()->element() == core::chemical::element::H ) {
				// proton chi
				proton_chis.push_back( nchi() );
			}
		} // for all found chis
	}

	if ( max_proton_chi_samples == 0 ) {
		return;
	}

	core::Size num_H_confs(1);  // TODO: Need to have introspection about base conformers?
	for ( core::Size pchi(1); pchi <= proton_chis.size(); ++pchi ) {
		if ( core::chemical::is_sp2_proton_chi( proton_chis[ pchi ], *this ) ) {
			num_H_confs *= 2; // The base amount of proton chi samples
		} else {
			num_H_confs *= 3;
		}
	}
	utility::vector1< Real > extra_samples;
	if ( num_H_confs > max_proton_chi_samples ) {
		TR.Warning << "Number of base proton chi samples (" << num_H_confs << ") for " << name() << " exceeds requested number of samples" << std::endl;
	} else if ( 3* num_H_confs > max_proton_chi_samples ) {
		TR << "Skipping extra samples for proton chis on " << name() << "; would give " << num_H_confs << " conformers." << std::endl;
	} else {
		extra_samples.push_back( 20 );
	}
	for ( core::Size pchi(1); pchi <= proton_chis.size(); ++pchi ) {
		if ( core::chemical::is_sp2_proton_chi( proton_chis[ pchi ], *this ) ) {
			utility::vector1< Real > sp2_sampling; //C++x11 static initializers would be nice here
			sp2_sampling.push_back( 0 );
			sp2_sampling.push_back( 180 );
			set_proton_chi( proton_chis[ pchi ], sp2_sampling, extra_samples );
		} else {
			utility::vector1< Real > sp3_sampling; //C++x11 static initializers would be nice here
			sp3_sampling.push_back( 60 );
			sp3_sampling.push_back( -60 );
			sp3_sampling.push_back( 180 );
			set_proton_chi( proton_chis[ pchi ], sp3_sampling, extra_samples );
		}
	}
}

// Special Functions ///////////////////////////////////////////////////////////

bool
MutableResidueType::validate_residue_type() const {
	// Unless we trigger a test, we've passed.
	bool success = true;

	// Check if the all_atoms() vector is messed up.
	if ( MutableResidueType::natoms() != ordered_atoms_.size() ) {
		success = false;
		TR.Error << "Residue " << name() << " has an inconsistent atom list: " << ordered_atoms_.size() << " vs. " << natoms() << std::endl;
	}

	// Pass through a set to check for duplicates
	std::set< VD > all_atom_set( ordered_atoms_.begin(), ordered_atoms_.end() );
	if ( all_atom_set.size() != ordered_atoms_.size() ) {
		success = false;
		TR.Error << "Residue " << name() << " has an atom list with " << ordered_atoms_.size() - all_atom_set.size() << " duplicates." << std::endl;
	}

	// Check that we don't have any unknown VDs in the atom list.
	for ( VD atm: all_atom_set ) {
		if ( ! has( atm ) ) {
			success = false;
			TR.Error << "Residue " << name() << " has an unknown vertex descriptor in its atom list." << std::endl;
		}
	}

	// Check that our chi_atoms_ has been correctly filled out.
	for ( core::Size ii(1); ii <= chis_.size(); ++ii ) {
		if ( chis_[ii] == nullptr ) {
			success = false;
			TR.Error << "Residue " << name() << " chi number " << ii << " is not properly set." << std::endl;
			continue; // No info to do the remaining checks
		}
		for ( VD vd: chis_[ii]->chi_atoms() ) {
			if ( all_atom_set.count( vd ) == 0 ) {
				success = false;
				TR.Error << "Residue " << name() << " has an atom in chi " << ii << " which isn't actually in the residue." << std::endl;
			}
		}
	}

	// Check that the icoor records have all been set.
	for ( VD atm: all_atom_set ) {
		if ( atom(atm).icoor() == nullptr ) {
			success = false;
			TR.Error << "Residue " << name() << " does not have an appropriate ICOOR record for atom '" << atom(atm).name() << "'" << std::endl;
		} else if ( ! atom(atm).icoor()->buildable(*this, /*verbose=*/true ) ) {
			success = false;
			TR.Error << "Residue " << name() << " does not have an buildable ICOOR record for atom '" << atom(atm).name() << "'" << std::endl;
		}
	}

	return success;
}


void
MutableResidueType::update_atom_type_set( AtomTypeSetCOP setting ) {
	debug_assert( setting );

	utility::vector1< std::string > old_types;
	if ( atom_type_set_ptr() ) {
		for ( VD atm: ordered_atoms_ ) {
			std::string const & old_name( atom_type_set()[ graph_[ atm ].atom_type_index() ].name() );
			old_types.push_back( old_name );
			set_atom_type( atm, "" ); // Null it out under the current atom type set.
		}
	}

	atom_type_set( setting );

	for ( core::Size ii(1); ii <= ordered_atoms_.size(); ++ii ) {
		if ( old_types.size() >= ii && setting->has_atom( old_types[ii] ) ) {
			set_atom_type( ordered_atoms_[ii], old_types[ii] );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

/// @brief A short utility function for update_icoors_after_connection_deletion()
std::string
correct_conn(std::string const & old, core::Size const conn_id_deleted ) {
	if ( string_to_icoord_type(old) != ICoordAtomIDType::CONNECT ) {
		return old;
	}
	core::Size conn_num = get_connection_number( old );
	if ( conn_num > conn_id_deleted ) {
		return "CONN" + std::to_string( conn_num - 1 );
	} else {
		return old;
	}
}


void
MutableResidueType::update_icoors_after_connection_deletion( core::Size const conn_id_deleted ) {
	for ( VD atm: ordered_atoms_ ) {
		MutableICoorRecordCOP const & ic = graph_[atm].icoor();
		if ( ic == nullptr ) { continue; }
		bool needs_update = false;
		// Do a two-pass version to prevent building new ICoor records if we don't need to.
		for ( core::Size ii(1); ii <= 3; ++ii ) {
			if ( ic->stub_type(ii) == ICoordAtomIDType::CONNECT ) {
				core::Size conn_num = get_connection_number( ic->stub_atom(ii) );
				if ( conn_num == conn_id_deleted ) {
					graph_[atm].icoor( nullptr ); // This icoor is now invalid.
					needs_update = false;
					break;
				} else if ( conn_num > conn_id_deleted ) {
					needs_update = true;
				}
			}
		}
		if ( needs_update ) {
			std::string stub1 = correct_conn( ic->stub_atom1(), conn_id_deleted );
			std::string stub2 = correct_conn( ic->stub_atom2(), conn_id_deleted );
			std::string stub3 = correct_conn( ic->stub_atom3(), conn_id_deleted );
			MutableICoorRecordOP new_icoor( new MutableICoorRecord(ic->phi(), ic->theta(), ic->d(), stub1, stub2, stub3 ) );
			graph_[atm].icoor( new_icoor );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void
MutableResidueType::copy_atom_info( ResidueType const & rt ) {

	// Iterate through all the atoms in the ResidueType, replicating them on our end.
	for ( core::Size ii(1); ii <= rt.natoms(); ++ii ) {
		VD new_vd = add_atom( rt.atom_name(ii), rt.atom_type(ii).name(), rt.mm_name(ii), rt.atom_charge(ii) );
		Atom & atom = graph_[ new_vd ];
		atom.gasteiger_atom_type( rt.gasteiger_atom_type(ii) );
		atom.formal_charge( rt.formal_charge(ii) );
		atom.ideal_xyz( rt.ideal_xyz(ii) );
		AtomICoor const & aic = rt.icoor(ii);
		MutableICoorRecordCOP ic = utility::pointer::make_shared< MutableICoorRecord const >( aic.phi(), aic.theta(), aic.d(),
			aic.stub_atom1().name( rt ),
			aic.stub_atom2().name( rt ),
			aic.stub_atom3().name( rt ) );
		atom.icoor( ic );
		atom.set_bonded_orbitals( rt.bonded_orbitals(ii) );
		atom.reset_all_properies( rt.atom_properties(ii) );
		// !!! Absolute Stereochemistry is not round-tripped
		// !!! Greek distance is not round-tripped

		if ( ii < rt.first_sidechain_atom() ) {
			atom.is_backbone( true );
		}
	}

	for ( auto const & pair: rt.bonds() ) {
		VD atm1 = ordered_atoms_[ pair.first ];
		VD atm2 = ordered_atoms_[ pair.second ];
		add_bond( atm1, atm2, rt.bond_type( pair.first, pair.second ) );
		bond(atm1, atm2).ringness( rt.bond_ringness( pair.first, pair.second ) );
		/// !!! The other attributes don't round-trip.
	}

	// Per atom info that needs to be deferred until all atoms are preset
	for ( core::Size ii(1); ii <= rt.natoms(); ++ii ) {
		//Cut bonds (defaults to false)
		for ( core::Size bnd: rt.cut_bond_neighbor(ii) ) {
			bond( ordered_atoms_[ii], ordered_atoms_[bnd] ).cut_bond(true);
		}
		//Shadowed atoms
		if ( rt.atom_being_shadowed(ii) != 0 ) {
			atom_shadowed_[ ordered_atoms_[ii] ] = ordered_atoms_[ rt.atom_being_shadowed(ii) ];
		}
	}

}

void
MutableResidueType::copy_other_info( ResidueType const & rt ) {

	if ( ! rt.is_base_type() ) { // Keep as nullptr if it is a base type
		base_type_cop_ = rt.get_base_type_cop();
	}

	///////////////////////
	// ICOOR
	//////////////////////

	root_atom_ = ordered_atoms_[ rt.root_atom() ];
	nbr_atom_ = ordered_atoms_[ rt.nbr_atom() ];
	nbr_radius_ = rt.nbr_radius();

	///////////////////////
	// CHI
	//////////////////////
	chis_.resize( rt.nchi() ); // Construct in-place

	for ( core::Size ii(1); ii <= rt.nchi(); ++ii ) {
		utility::vector1< VD > chi_vert;
		for ( core::Size atm: rt.chi_atoms(ii) ) {
			chi_vert.push_back( ordered_atoms_[atm] );
		}
		chis_[ii] = utility::pointer::make_shared< MutableChiRecord >( chi_vert );
		if ( rt.is_proton_chi(ii) ) {
			core::Size prot_chi = rt.chi_2_proton_chi(ii);
			chis_[ii]->set_proton_chi( rt.proton_chi_samples(prot_chi), rt.proton_chi_extra_samples(prot_chi) );
		}
		chis_[ii]->set_chi_rotamers( rt.chi_rotamers(ii) );
	}

	///////////////////////
	// NU
	//////////////////////
	nu_atoms_.resize( rt.n_nus() );

	for ( core::Size ii(1); ii <= rt.n_nus(); ++ii ) {
		for ( core::Size aa: rt.nu_atoms(ii) ) {
			nu_atoms_[ii].push_back( ordered_atoms_[aa] );
		}
	}

	///////////////////////
	// RINGS
	//////////////////////
	ring_atoms_.resize( rt.n_rings() );

	for ( core::Size ii(1); ii <= rt.n_rings(); ++ii ) {
		for ( core::Size aa: rt.ring_atoms(ii) ) {
			ring_atoms_[ii].push_back( ordered_atoms_[aa] );
		}
	}

	ring_saturation_types_ = rt.ring_saturation_types();
	lowest_ring_conformer_ = rt.lowest_ring_conformers();
	low_ring_conformers_ = rt.low_ring_conformers();

	///////////////////////
	// Mainchain atoms
	//////////////////////
	for ( core::Size aa: rt.mainchain_atoms() ) {
		mainchain_atoms_.push_back( ordered_atoms_[aa] );
	}

	///////////////////////
	// Actcoord atoms
	//////////////////////
	for ( core::Size aa: rt.actcoord_atoms() ) {
		graph_[ ordered_atoms_[aa] ].is_actcoord( true );
	}

	///////////////////////
	// Connections
	//////////////////////
	lower_connect_id_ = rt.lower_connect_id();
	upper_connect_id_ = rt.upper_connect_id();

	for ( core::Size ii(1); ii <= rt.n_possible_residue_connections(); ++ii ) {
		ResidueConnection const & old_conn = rt.residue_connection(ii);
		MutableResidueConnection new_conn( ordered_atoms_[ old_conn.atomno() ] );
		AtomICoor const & aic = old_conn.icoor();
		MutableICoorRecord ic( aic.phi(), aic.theta(), aic.d(),
			aic.stub_atom1().name( rt ),
			aic.stub_atom2().name( rt ),
			aic.stub_atom3().name( rt ) );
		new_conn.icoor( ic );
		residue_connections_.push_back( new_conn );
	}
}


///////////////////////////////////////////////////////////////////////////////
// Helper methods //////////////////////////////////////////////////////////////

// Insertion operator (overloaded so that MutableResidueType can be "printed" in PyRosetta).
std::ostream &
operator<<(std::ostream & output, MutableResidueType const & object_to_output)
{
	object_to_output.show(output);
	return output;
}

} // chemical
} // core


#ifdef    SERIALIZATION

template< class Archive >
void
core::chemical::MutableResidueType::save( Archive & arc ) const {
	using namespace core::chemical;

	arc( cereal::base_class< core::chemical::ResidueTypeBase >( this ) );

	// We can't serialize the Boost graph directly, and, besides, the VDs will change
	// when we reconstruct the graph. To circumvent this, we serialize atoms directly
	// in index order (as atom names might not be unique ... theoretically).
	arc( ::cereal::make_nvp( "natoms", ordered_atoms_.size() ) );
	for ( VD vd : ordered_atoms_ ) {
		SERIALIZE_VD( arc, vd );
		arc( graph_[ vd ] ); // EXEMPT graph_ ordered_atoms_
	}

	arc( ::cereal::make_nvp( "nbonds", boost::num_edges(graph_) ) );
	for ( EIterPair biter_pair( boost::edges(graph_) ); biter_pair.first != biter_pair.second; ++biter_pair.first ) {
		ED ed( *biter_pair.first );
		VD source( boost::source( ed, graph_ ) ), target( boost::target( ed, graph_ ) );
		SERIALIZE_VD( arc, source );
		SERIALIZE_VD( arc, target );
		arc( graph_[ ed ] );
	}

	SERIALIZE_T_VD_MAP( arc, atom_name_to_vd_ );

	arc( CEREAL_NVP( base_type_cop_ ) ); //ResidueTypeCOP

	SERIALIZE_VD( arc, root_atom_, "root_atom_" ); // EXEMPT root_atom_
	SERIALIZE_VD( arc, nbr_atom_, "nbr_atom_" ); // EXEMPT nbr_atom_
	arc( CEREAL_NVP( nbr_radius_ ) ); // Real

	SERIALIZE_VD_VD_MAP( arc, atom_shadowed_ ); // EXEMPT atom_shadowed_
	SERIALIZE_VD_VECTOR( arc, mainchain_atoms_ ); // EXEMPT mainchain_atoms_

	SERIALIZE_NESTED_VD_VECTOR( arc, nu_atoms_ ); // EXEMPT nu_atoms_
	SERIALIZE_NESTED_VD_VECTOR( arc, ring_atoms_ ); // EXEMPT ring_atoms_

	// Need to be VD-adjusted on the back end.
	arc( CEREAL_NVP( chis_ ) ); // utility::vector1< MutableChiRecordOP >
	arc( CEREAL_NVP( residue_connections_ ) ); // utility::vector1<MutableResidueConnection>

	arc( CEREAL_NVP( ring_saturation_types_ ) );  // utility::vector1< core::chemical::rings::RingSaturationType >
	arc( CEREAL_NVP( low_ring_conformers_ ) ); // utility::vector1<utility::vector1<std::string> >
	arc( CEREAL_NVP( lowest_ring_conformer_ ) ); // utility::vector1<std::string>

	arc( CEREAL_NVP( lower_connect_id_ ) ); // Size
	arc( CEREAL_NVP( upper_connect_id_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::MutableResidueType::load( Archive & arc ) {
	using namespace core::chemical;

	arc( cereal::base_class< core::chemical::ResidueTypeBase >( this ) );

	std::map< VD, VD > old_to_new;
	old_to_new[ MutableResidueType::null_vertex ] = MutableResidueType::null_vertex; // Null vertex self-converts.

	// EXEMPT graph_
	core::Size natoms; arc( natoms );
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		VD old_vd; DESERIALIZE_VD( arc, old_vd );
		Atom atom; arc( atom );
		atom.update_typesets( *this );
		VD new_vd = graph_.add_vertex(atom);
		old_to_new[ old_vd ] = new_vd;
		ordered_atoms_.push_back(new_vd); // EXEMPT ordered_atoms_
	}
	// Now we have a completely set old_to_nes

	core::Size nbonds; arc( nbonds );
	for ( core::Size ii(1); ii <= nbonds; ++ii ) {
		VD source, target;
		DESERIALIZE_VD( arc, source );
		source = old_to_new.at(source);
		DESERIALIZE_VD( arc, target );
		target = old_to_new.at(target);
		Bond bond; arc( bond );
		bool added; ED new_ed;
		boost::tie( new_ed, added ) = graph_.add_edge( source, target, bond );
		runtime_assert( added );
	}

	DESERIALIZE_T_VD_MAP( arc, atom_name_to_vd_, old_to_new ); // EXEMPT atom_name_to_vd_

	arc( base_type_cop_ ); // ResidueTypeCOP

	DESERIALIZE_VD( arc, root_atom_, old_to_new ); // EXEMPT root_atom_
	DESERIALIZE_VD( arc, nbr_atom_, old_to_new ); // EXEMPT nbr_atom_
	arc( nbr_radius_ ); // Real

	DESERIALIZE_VD_VD_MAP( arc, atom_shadowed_, old_to_new ); // EXEMPT atom_shadowed_
	DESERIALIZE_VD_VECTOR( arc, mainchain_atoms_, old_to_new ); // EXEMPT mainchain_atoms_

	DESERIALIZE_NESTED_VD_VECTOR( arc, nu_atoms_, old_to_new ); // EXEMPT nu_atoms_
	DESERIALIZE_NESTED_VD_VECTOR( arc, ring_atoms_, old_to_new ); // EXEMPT ring_atoms_

	arc( chis_ ); // utility::vector1< MutableChiRecordOP >
	for ( MutableChiRecordOP & chi: chis_ ) {
		if ( chi != nullptr ) {
			chi->remap_atom_vds( old_to_new );
		}
	}

	arc( residue_connections_ ); // utility::vector1<MutableResidueConnection>
	for ( MutableResidueConnection & rescon : residue_connections_ ) {
		rescon.remap_atom_vds( old_to_new );
	}

	arc( ring_saturation_types_ );  // utility::vector1< core::chemical::rings::RingSaturationType >
	arc( low_ring_conformers_ ); // utility::vector1<utility::vector1<std::string> >
	arc( lowest_ring_conformer_ ); // utility::vector1<std::string>

	arc( lower_connect_id_ ); // Size
	arc( upper_connect_id_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::MutableResidueType );
CEREAL_REGISTER_TYPE( core::chemical::MutableResidueType )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_MutableResidueType )
#endif // SERIALIZATION
