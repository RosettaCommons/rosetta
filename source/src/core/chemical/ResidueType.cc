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
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to rings, properties, lipids, carbohydrates, and other non-AAs)

// Unit headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>

// Package Headers
#include <core/chemical/ResidueProperties.hh>
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/chemical/residue_support.hh>
#include <core/chemical/icoor_support.hh>
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

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>
#include <utility/graph/ring_detection.hh>

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

static THREAD_LOCAL basic::Tracer tr( "core.chemical.ResidueType" );

// must be a better place for this, probably already exists!
inline
std::string
strip_whitespace( std::string const & name )
{
	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to dothis?
	return trimmed_name;
}

VD const ResidueType::null_vertex = boost::graph_traits<ResidueGraph>::null_vertex();

ResidueType::ResidueType() = default; // private, not deleted because of serialization

ResidueType::ResidueType(
	AtomTypeSetCOP atom_types,
	ElementSetCOP elements,
	MMAtomTypeSetCOP mm_atom_types,
	orbitals::OrbitalTypeSetCOP orbital_types
) : utility::pointer::ReferenceCount(),
	atom_types_( atom_types ),
	elements_( elements ),
	mm_atom_types_( mm_atom_types ),
	gasteiger_atom_types_(),
	orbital_types_( orbital_types ),
	mode_( INVALID_t ),
	graph_(),
	orbitals_(),
	nheavyatoms_(0),
	n_hbond_acceptors_(0),
	n_hbond_donors_(0),
	n_backbone_heavyatoms_(0),
	first_sidechain_hydrogen_( 0 ),
	disulfide_atom_name_( "NONE" ),
	rotamer_library_specification_( nullptr ),
	properties_( ResiduePropertiesOP( new ResidueProperties( this ) ) ),
	aa_( aa_unk ),
	rotamer_aa_( aa_unk ),
	backbone_aa_( aa_unk ),
	na_analogue_( aa_unp ),
	base_name_(""),
	base_type_cop_(), //Assumes that this is a base type by default.
	name_(""),
	name3_(""),
	name1_(' '),
	interchangeability_group_(),
	root_atom_( ResidueType::null_vertex ),
	nbr_atom_( ResidueType::null_vertex ),
	nbr_radius_( 0 ),
	force_nbr_atom_orient_(false),
	remap_pdb_atom_names_(false),
	mass_(0),
	lower_connect_id_( 0 ),
	upper_connect_id_( 0 ),
	n_non_polymeric_residue_connections_( 0 ),
	n_polymeric_residue_connections_( 0 ),
	carbohydrate_info_(/* NULL */),
	rings_and_their_edges_(),
	nbr_atom_indices_( 0 ),
	rama_prepro_mainchain_torsion_potential_name_(),
	rama_prepro_mainchain_torsion_potential_name_beforeproline_(),
	rama_prepro_map_file_name_(),
	rama_prepro_map_file_name_beforeproline_(),
	finalized_(false),
	nondefault_(false)
{
	if ( atom_types != nullptr ) {
		// For any actual ResidueType atom_types should be valid, but there's tricky RTS bootstraping logic
		// in Patch.cc that uses nullptr AtomTypeSets.
		mode_ = atom_types->mode();
	}
}

ResidueType::~ResidueType()
{
	destruction_obs_hub_( RestypeDestructionEvent( this ) ); // Notify the destruction observers that this residue type is being destroyed.
}

ResidueType::ResidueType( ResidueType const & residue_type ):
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< ResidueType >()
{
	*this = residue_type;
}

ResidueType &
ResidueType::operator=( ResidueType const & residue_type )
{
	// We're overwriting the contents of the ResidueType - it's effectively being destroyed.
	// Notify observers and then disconnect them.
	destruction_obs_hub_( RestypeDestructionEvent( this ) );
	destruction_obs_hub_.clear();

	atom_types_ = residue_type.atom_types_;
	elements_ = residue_type.elements_;
	mm_atom_types_ = residue_type.mm_atom_types_;
	gasteiger_atom_types_ = residue_type.gasteiger_atom_types_;
	orbital_types_ = residue_type.orbital_types_;
	conformer_sets_ = residue_type.conformer_sets_;
	mode_ = residue_type.mode_;
	graph_ = residue_type.graph_; // copying Will change the VDs
	vd_to_index_.clear(); // must be regenerated for new VDs
	atom_base_ = residue_type.atom_base_;
	abase2_ = residue_type.abase2_;
	ordered_atoms_.clear(); // This must be regenerated to hold the new vertex_descriptors
	orbitals_ = residue_type.orbitals_;
	nheavyatoms_ = residue_type.nheavyatoms_;
	n_hbond_acceptors_ = residue_type.n_hbond_acceptors_;
	n_hbond_donors_ = residue_type.n_hbond_donors_;
	n_backbone_heavyatoms_ = residue_type.n_backbone_heavyatoms_;
	first_sidechain_hydrogen_ = residue_type.first_sidechain_hydrogen_;
	bonded_neighbor_ = residue_type.bonded_neighbor_;
	bonded_neighbor_type_ = residue_type.bonded_neighbor_type_;
	cut_bond_neighbor_ = residue_type.cut_bond_neighbor_;
	attached_H_begin_ = residue_type.attached_H_begin_;
	attached_H_end_ = residue_type.attached_H_end_;
	icoor_ = residue_type.icoor_;
	dihedral_atom_sets_ = residue_type.dihedral_atom_sets_;
	dihedrals_for_atom_ = residue_type.dihedrals_for_atom_;
	improper_dihedral_atom_sets_ = residue_type.improper_dihedral_atom_sets_;
	improper_dihedrals_for_atom_ = residue_type.improper_dihedrals_for_atom_;
	bondangle_atom_sets_ = residue_type.bondangle_atom_sets_;
	bondangles_for_atom_ = residue_type.bondangles_for_atom_;
	atom_shadowed_ = residue_type.atom_shadowed_;
	last_controlling_chi_ = residue_type.last_controlling_chi_;
	atoms_last_controlled_by_chi_ = residue_type.atoms_last_controlled_by_chi_;
	atoms_with_orb_index_ = residue_type.atoms_with_orb_index_;
	Haro_index_ = residue_type.Haro_index_;
	Hpol_index_ = residue_type.Hpol_index_;
	accpt_pos_ = residue_type.accpt_pos_;
	Hpos_polar_ = residue_type.Hpos_polar_;
	Hpos_apolar_ = residue_type.Hpos_apolar_;
	accpt_pos_sc_ = residue_type.accpt_pos_sc_;
	Hpos_polar_sc_ = residue_type.Hpos_polar_sc_;
	all_bb_atoms_ = residue_type.all_bb_atoms_;
	all_sc_atoms_ = residue_type.all_sc_atoms_;
	metal_binding_atoms_ = residue_type.metal_binding_atoms_;
	disulfide_atom_name_ = residue_type.disulfide_atom_name_;
	mainchain_atoms_ = residue_type.mainchain_atoms_;
	actcoord_atoms_ = residue_type.actcoord_atoms_;
	chi_atoms_ = residue_type.chi_atoms_;
	is_proton_chi_ = residue_type.is_proton_chi_;
	proton_chis_ = residue_type.proton_chis_;
	chi_2_proton_chi_ = residue_type.chi_2_proton_chi_;
	proton_chi_samples_ = residue_type.proton_chi_samples_;
	proton_chi_extra_samples_ = residue_type.proton_chi_extra_samples_;
	nu_atoms_ = residue_type.nu_atoms_;
	ring_atoms_ = residue_type.ring_atoms_;
	path_distance_ = residue_type.path_distance_;
	atom_name_to_vd_.clear(); // This must be regenerated below to hold the new vertex_descriptors
	atom_aliases_ = residue_type.atom_aliases_;
	orbitals_index_ = residue_type.orbitals_index_;
	chi_rotamers_ = residue_type.chi_rotamers_;
	rotamer_library_specification_ = residue_type.rotamer_library_specification_;
	lowest_ring_conformer_ = residue_type.lowest_ring_conformer_;
	low_ring_conformers_ = residue_type.low_ring_conformers_;
	properties_ = ResiduePropertiesOP( new ResidueProperties( *residue_type.properties_, this ) );
	aa_ = residue_type.aa_;
	rotamer_aa_ = residue_type.rotamer_aa_;
	backbone_aa_ = residue_type.backbone_aa_;
	na_analogue_ = residue_type.na_analogue_;
	base_name_ = residue_type.base_name_;
	base_type_cop_ = residue_type.base_type_cop_; //If the residue type that we're copying has a base type, copy the base type pointer, too.
	runtime_assert( base_type_cop_.get() != this );
	name_ = residue_type.name_;
	name3_ = residue_type.name3_;
	name1_ = residue_type.name1_;
	interchangeability_group_ = residue_type.interchangeability_group_;
	root_atom_ = residue_type.root_atom_;
	nbr_atom_ = residue_type.nbr_atom_;
	nbr_radius_ = residue_type.nbr_radius_;
	force_nbr_atom_orient_ = residue_type.force_nbr_atom_orient_;
	remap_pdb_atom_names_ = residue_type.remap_pdb_atom_names_;
	mass_ = residue_type.mass_;
	residue_connections_ = residue_type.residue_connections_;
	atom_2_residue_connection_map_ = residue_type.atom_2_residue_connection_map_;
	atoms_within_one_bond_of_a_residue_connection_ = residue_type.atoms_within_one_bond_of_a_residue_connection_;
	within1bonds_sets_for_atom_ = residue_type.within1bonds_sets_for_atom_;
	atoms_within_two_bonds_of_a_residue_connection_ = residue_type.atoms_within_two_bonds_of_a_residue_connection_;
	within2bonds_sets_for_atom_ = residue_type.within2bonds_sets_for_atom_;
	lower_connect_id_ = residue_type.lower_connect_id_;
	upper_connect_id_ = residue_type.upper_connect_id_;
	n_non_polymeric_residue_connections_ = residue_type.n_non_polymeric_residue_connections_;
	n_polymeric_residue_connections_ = residue_type.n_polymeric_residue_connections_;
	force_bb_ = residue_type.force_bb_;
	rna_info_ = residue_type.rna_info_;
	// CarbohydrateInfo has a back-pointer to ResidueType and must be reset during finalize
	carbohydrate_info_ = nullptr; /* NULL */
	//rings_and_their_edges_=residue_type.rings_and_their_edges_; //Apparently updated below
	atom_base_indices_ = residue_type.atom_base_indices_;
	abase2_indices_ = residue_type.abase2_indices_;
	chi_atoms_indices_ = residue_type.chi_atoms_indices_;
	nu_atoms_indices_ = residue_type.nu_atoms_indices_;
	ring_atoms_indices_ = residue_type.ring_atoms_indices_;
	mainchain_atoms_indices_ = residue_type.mainchain_atoms_indices_;
	nbr_atom_indices_ = residue_type.nbr_atom_indices_;
	actcoord_atoms_indices_ = residue_type.actcoord_atoms_indices_;
	cut_bond_neighbor_indices_ = residue_type.cut_bond_neighbor_indices_;
	atom_shadowed_indices_ = residue_type.atom_shadowed_indices_;
	rama_prepro_mainchain_torsion_potential_name_ = residue_type.rama_prepro_mainchain_torsion_potential_name_;
	rama_prepro_mainchain_torsion_potential_name_beforeproline_ = residue_type.rama_prepro_mainchain_torsion_potential_name_beforeproline_;
	rama_prepro_map_file_name_ = residue_type.rama_prepro_map_file_name_;
	rama_prepro_map_file_name_beforeproline_ = residue_type.rama_prepro_map_file_name_beforeproline_;
	finalized_ = residue_type.finalized_;
	defined_adducts_ = residue_type.defined_adducts_;
	nondefault_ = residue_type.nondefault_;

	// When you copy vertex descriptors from cached data, the vertex descriptors are pointing to the old copied graph.
	// New vertices are assigned.  You have to map the old vertex to the new vertex.

	std::map<ED, ED> old_edges_to_new_edges;
	for (
			EIterPair ep = boost::edges(graph_), old_ep = boost::edges(residue_type.graph_);
			ep.first != ep.second;
			++ep.first, ++old_ep.first
			) {
		EIter e_iter = ep.first;
		EIter old_e_iter = old_ep.first;
		ED  ed = *e_iter;
		ED old_ed = *old_e_iter;
		old_edges_to_new_edges[old_ed] = ed;
	}

	std::list<utility::vector1<ED> >  old_rings(rings_and_their_edges_);
	rings_and_their_edges_.clear();
	for ( std::list<utility::vector1<ED> >::const_iterator it = old_rings.begin(); it != old_rings.end(); ++it ) {
		utility::vector1<ED> old_edges = *it;
		utility::vector1<ED> new_edges;
		for ( Size i = 1; i <= old_edges.size(); ++i ) {
			new_edges.push_back(old_edges_to_new_edges[old_edges[i]]);
		}
		rings_and_their_edges_.push_back(new_edges);
	}


	// Setup the atom_name_to_vd_ and map old to new...
	std::map<VD, VD> old_to_new;
	old_to_new[ ResidueType::null_vertex ] = ResidueType::null_vertex; // Null verticies in one are null verticies in the other.
	for (
			VIterPair vp = boost::vertices(graph_), old_vp= boost::vertices(residue_type.graph_);
			vp.first != vp.second;
			++vp.first, ++old_vp.first
			) {
		auto v_iter= vp.first;
		auto old_v_iter= old_vp.first;
		VD vd = *v_iter;
		VD old_vd = *old_v_iter;
		old_to_new[old_vd] = vd; //Assuming the boost::graph copy preserves ordering within the vertices list
		Atom & a = graph_[vd];
		debug_assert( a == residue_type.graph_[old_vd]);
#ifndef NDEBUG
		NameVDInserted const name_vd_inserted =
#endif
			atom_name_to_vd_.insert(NameVDPair(a.name(), vd));
		debug_assert(name_vd_inserted.second); // Don't add atoms with the same name
		atom_name_to_vd_.insert( NameVDPair( strip_whitespace( a.name() ), vd) );
		//debug_assert(strip_name_vd_inserted.second); // If this is 4 chars, than it will be the same as before.
	}
	// We also have to add in the aliased names
	for ( auto iter( residue_type.atom_aliases_.begin() ), iter_end( residue_type.atom_aliases_.end() );
			iter != iter_end; ++iter ) {
		debug_assert( !atom_name_to_vd_.count( iter->first ) );
		debug_assert( atom_name_to_vd_.count( iter->second ) );
		VD vd( atom_name_to_vd_[ iter->second ] );
		atom_name_to_vd_.insert( NameVDPair( iter->first, vd) );
	}
	// Setup the temporary ordered_atoms_ vector for refactor
	// VKM, 5 Sept 2015: I'm not sure the code below is working as desired.  Valgrind issues might be caused here.
	auto begin = residue_type.ordered_atoms_.begin();
	VDs::const_iterator const end = residue_type.ordered_atoms_.end();
	for ( ; begin != end; ++begin ) {
		VD old_vd = *begin;
		VD vd = old_to_new[old_vd];
		ordered_atoms_.push_back(vd);
		vd_to_index_[vd] = ordered_atoms_.size();
	}
	std::map<VD, VD> old_atom_base(atom_base_);
	atom_base_.clear();
	for ( auto & it : old_atom_base ) {
		VD old_key = it.first;
		VD old_value = it.second;
		VD new_key = old_to_new[old_key];
		VD new_value = old_to_new[old_value];
		atom_base_[new_key] = new_value;
	}

	std::map<VD, VD> old_atom_shadowed(atom_shadowed_);
	atom_shadowed_.clear();
	for ( auto & it : old_atom_shadowed ) {
		VD old_key = it.first;
		VD old_value = it.second;
		VD new_key = old_to_new[old_key];
		VD new_value = old_to_new[old_value];
		atom_shadowed_[new_key] = new_value;
	}

	utility::vector1<utility::vector1<VD> > old_chi_atoms(chi_atoms_);
	chi_atoms_.clear();
	//chi_atoms_.resize(old_chi_atoms.size());
	for ( auto const & old_vector : old_chi_atoms ) {
		debug_assert(old_vector.size() == 4);
		utility::vector1<VD> new_vector;
		for ( Size i= 1; i<= old_vector.size(); ++i ) {
			new_vector.push_back(old_to_new[old_vector[i]]);
		}
		chi_atoms_.push_back(new_vector);
	}

	utility::vector1< utility::vector1< VD > > old_nu_atoms( nu_atoms_ );
	nu_atoms_.clear();
	for ( auto const & old_vector : old_nu_atoms ) {
		debug_assert( old_vector.size() == 4 );
		utility::vector1< VD > new_vector;
		for ( Size i= 1; i<= old_vector.size(); ++i ) {
			new_vector.push_back( old_to_new[ old_vector[ i ] ] );
		}
		nu_atoms_.push_back( new_vector );
	}

	utility::vector1< utility::vector1< VD > > old_ring_atoms( ring_atoms_ );
	ring_atoms_.clear();
	for ( auto const & old_vector : old_ring_atoms ) {
		utility::vector1< VD > new_vector;
		for ( Size i( 1 ); i <= old_vector.size(); ++i ) {
			new_vector.push_back( old_to_new[ old_vector[ i ] ] );
		}
		ring_atoms_.push_back( new_vector );
	}

	utility::vector1<VD> old_mainchain( mainchain_atoms_ );
	mainchain_atoms_.clear();
	for ( Size i= 1; i <= old_mainchain.size(); ++i ) {
		mainchain_atoms_.push_back( old_to_new[ old_mainchain[i] ]);
	}

	// We map null_vertex to null_vertex
	VD old_root = root_atom_;
	root_atom_ = old_to_new[old_root];

	VD old_nbr = nbr_atom_;
	nbr_atom_ = old_to_new[old_nbr];

	utility::vector1<VD> old_act(actcoord_atoms_);
	actcoord_atoms_.clear();
	for ( Size i=1; i<= old_act.size(); ++i ) {
		actcoord_atoms_.push_back(old_to_new[old_act[i]]);
	}


	std::map<VD, utility::vector1<VD> > old_cut_bonds(cut_bond_neighbor_);
	cut_bond_neighbor_.clear();
	for ( auto const & elem : old_cut_bonds ) {
		VD old_key = elem.first;
		utility::vector1<VD> old_value = elem.second;
		VD new_key = old_to_new[old_key];
		utility::vector1<VD> new_value;
		for ( Size i=1; i<= old_value.size(); ++i ) {
			new_value.push_back(old_to_new[ old_value[i] ] );
		}
		cut_bond_neighbor_[ new_key ] = new_value;
	}

	std::map< VD, AtomICoor > old_icoor(icoor_);
	icoor_.clear();
	for ( auto const & elem : old_icoor ) {
		VD old_key = elem.first;
		VD new_key = old_to_new[old_key];
		AtomICoor old_icoor = elem.second; //now we have to change the vertex descriptors within icoor. They are pointing to an old vertex descriptor
		old_icoor.built_atom_vertex( old_to_new[ old_icoor.built_atom_vertex() ] );
		for ( Size i=1; i<= 3; ++i ) {
			ICoorAtomID & stub_atom( old_icoor.stub_atom(i) );
			if ( stub_atom.type() == ICoorAtomID::INTERNAL ) {
				VD old_vertex = stub_atom.vertex();
				VD new_vertex = old_to_new[old_vertex];
				stub_atom.vertex(new_vertex);
			}
		}
		icoor_[new_key] = old_icoor;
	}

	for ( auto & residue_connection : residue_connections_ ) {
		VD old_vertex = residue_connection.vertex();
		residue_connection.vertex(old_to_new[old_vertex]);
		residue_connection.atomno(vd_to_index_.find(residue_connection.vertex())->second);

		AtomICoor new_icoor = residue_connection.icoor();
		for ( Size j = 1; j <= 3; ++j ) {
			VD old_vd = new_icoor.stub_atom( j ).vertex();
			new_icoor.stub_atom( j ).vertex( old_to_new[old_vd] );
			new_icoor.stub_atom( j ).atomno( vd_to_index_.find(new_icoor.stub_atom(j).vertex())->second );
		}
		residue_connection.icoor( new_icoor );

	}

	for ( Size index=1; index<= natoms(); ++index ) {
		utility::vector1< core::Size > const orbs(graph_[ordered_atoms_[index]].bonded_orbitals());
		for ( core::Size orb : orbs ) {
			orbitals_[orb].new_icoor().vertex1( old_to_new[ orbitals_[orb].new_icoor().vertex1()  ] );
			orbitals_[orb].new_icoor().vertex2(old_to_new[ orbitals_[orb].new_icoor().vertex2()  ] );
			orbitals_[orb].new_icoor().vertex3( old_to_new[ orbitals_[orb].new_icoor().vertex3()  ] );
		}
	}

	utility::vector1<VD> old_bb(force_bb_);
	force_bb_.clear();
	for ( Size i=1; i<= old_bb.size(); ++i ) {
		force_bb_.push_back(old_to_new[ old_bb[i] ]);
	}

	return *this;
}

//////////////////////////////////////////////////////////////////////////////

/// @brief make a copy
ResidueTypeOP
ResidueType::clone() const
{
	ResidueTypeOP rsd_ptr( new ResidueType( *this ) );
	return rsd_ptr;
}


/// @brief make a copy, but only with all the stuff needed by patch selectors needs to be filled.
ResidueTypeOP
ResidueType::placeholder_clone() const
{
	ResidueTypeOP rsd( new ResidueType( atom_type_set_ptr(), element_set_ptr(), mm_atom_types_ptr(), orbital_types_ptr() ) );
	rsd->name ( name() );
	rsd->name1( name1() );
	rsd->name3( name3() );
	rsd->aa( aa() );
	rsd->interchangeability_group( interchangeability_group() );
	ResiduePropertiesOP properties( new ResidueProperties( *properties_, &( *rsd ) ) );
	rsd->set_properties( properties );
	// following would totally work, but needs atom_names -- might be a way to refactor?
	// if ( properties->has_property( CARBOHYDRATE ) ) {
	//  rsd->carbohydrate_info_ = carbohydrates::CarbohydrateInfoOP( new carbohydrates::CarbohydrateInfo( get_self_weak_ptr() ) );
	// }
	return rsd;
}

Atom & ResidueType::atom(Size const atom_index){
	return graph_[ ordered_atoms_[atom_index] ];
}
Atom const & ResidueType::atom(Size const atom_index) const{
	return graph_[ ordered_atoms_[atom_index] ];
}
Atom & ResidueType::atom(std::string const & atom_name){
	return graph_[ atom_name_to_vd_[atom_name] ];
}
Atom const & ResidueType::atom(std::string const & atom_name) const{
	auto found = atom_name_to_vd_.find( atom_name );
	debug_assert( found != atom_name_to_vd_.end());
	return graph_[  found->second ];
}
Atom & ResidueType::atom(VD const atom_vd){
	debug_assert( has(atom_vd) );
	return graph_[ atom_vd ];
}
Atom const & ResidueType::atom(VD const atom_vd) const{
	debug_assert( has(atom_vd) );
	return graph_[ atom_vd ];
}

Orbital const & ResidueType::orbital(Size const orbital_index) const{
	return orbitals_[orbital_index];
}
Orbital const & ResidueType::orbital(std::string const & orbital_name) const{
	return orbitals_[ orbital_index(orbital_name) ];
}

Bond &
ResidueType::bond(ED const ed){
	return graph_[ ed ];
}
Bond const &
ResidueType::bond(ED const ed) const{
	return graph_[ ed ];
}

Bond &
ResidueType::bond(VD vd1, VD vd2){
	ED ed;
	bool found;
	boost::tie( ed, found ) = boost::edge( vd1, vd2 , graph_ );
	if ( ! found ) {
		utility_exit_with_message( "Cannot find bond between " + atom_name(vd1) + " and " + atom_name(vd2) + " in residue " + name() );
	}
	return graph_[ ed ];
}

Bond const &
ResidueType::bond(VD vd1, VD vd2) const{
	ED ed;
	bool found;
	boost::tie( ed, found ) = boost::edge( vd1, vd2 , graph_ );
	if ( ! found ) {
		utility_exit_with_message( "Cannot find bond between " + atom_name(vd1) + " and " + atom_name(vd2) + " in residue " + name() );
	}
	return graph_[ ed ];
}

Bond &
ResidueType::bond(std::string const & atom1, std::string const & atom2) {
	return bond( atom_vertex( atom1 ), atom_vertex( atom2 ) );
}

Bond const &
ResidueType::bond(std::string const & atom1, std::string const & atom2) const {
	return bond( atom_vertex( atom1 ), atom_vertex( atom2 ) );
}

// Connections ////////////////////////////////////////////////////////////////
// Lower
ResidueConnection const &
ResidueType::lower_connect() const
{
	debug_assert( properties_->has_property( POLYMER ) );
	debug_assert( lower_connect_id_ != 0 );
	return residue_connections_[ lower_connect_id_ ];
}

Size
ResidueType::lower_connect_atom() const {
	debug_assert( properties_->has_property( POLYMER ) );
	debug_assert( lower_connect_id_ != 0 );
	return vd_to_index_.find(residue_connections_[ lower_connect_id_ ].vertex())->second;
}

// Set the atom which connects to the lower connection.
void
ResidueType::set_lower_connect_atom( std::string const & atm_name )
{
	finalized_ = false;
	if ( atm_name == "NONE" ) {
		if ( lower_connect_id_ != 0 ) {
			tr.Debug << "ERASING LOWER_CONNECT: " << lower_connect_id_ << " lcid: " << upper_connect_id_ << std::endl;
			auto to_erase( residue_connections_.begin() );
			to_erase += lower_connect_id_ - 1;
			residue_connections_.erase(  to_erase );
			update_residue_connection_mapping();
			debug_assert( n_polymeric_residue_connections_ != 0 );
			--n_polymeric_residue_connections_;
			if ( lower_connect_id_ < upper_connect_id_ ) { --upper_connect_id_; }
			lower_connect_id_ = 0;

		}
	} else {
		if ( lower_connect_id_ == 0 ) {
			ResidueConnection rc( atom_index( atm_name ), ordered_atoms_[atom_index( atm_name )] );
			residue_connections_.push_back( rc );
			lower_connect_id_ = residue_connections_.size();
			++n_polymeric_residue_connections_;
		} else {
			residue_connections_[ lower_connect_id_ ].atomno( atom_index( atm_name ) );
			residue_connections_[ lower_connect_id_ ].vertex( ordered_atoms_[atom_index( atm_name )] );
		}
	}
	update_residue_connection_mapping();
}


// Upper
ResidueConnection const &
ResidueType::upper_connect() const
{
	debug_assert( properties_->has_property( POLYMER ) );
	debug_assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ];
}

Size
ResidueType::upper_connect_atom() const
{
	debug_assert( properties_->has_property( POLYMER ) );
	debug_assert( upper_connect_id_ != 0 );
	return vd_to_index_.find(residue_connections_[ upper_connect_id_ ].vertex())->second;
}

// Set the atom which connects to the upper connection.
void
ResidueType::set_upper_connect_atom( std::string const & atm_name )
{
	finalized_ = false;
	if ( atm_name == "NONE" ) {
		if ( upper_connect_id_ != 0 ) {
			tr.Debug << "ERASING UPPER_CONNECT: " << upper_connect_id_ << " lcid: " << lower_connect_id_  << std::endl;
			auto to_erase( residue_connections_.begin() );
			to_erase += upper_connect_id_ - 1;
			residue_connections_.erase(  to_erase );
			debug_assert( n_polymeric_residue_connections_ != 0 );
			--n_polymeric_residue_connections_;
			if ( upper_connect_id_ < lower_connect_id_ ) { --lower_connect_id_; }
			upper_connect_id_ = 0;
		}
	} else {
		if ( upper_connect_id_ == 0 ) {
			ResidueConnection rc( atom_index( atm_name ), ordered_atoms_[atom_index( atm_name )] );
			residue_connections_.push_back( rc );
			upper_connect_id_ = residue_connections_.size();
			++n_polymeric_residue_connections_;
		} else {
			residue_connections_[ upper_connect_id_ ].atomno( atom_index( atm_name ) );
			residue_connections_[ upper_connect_id_ ].vertex( ordered_atoms_[atom_index( atm_name )] );
		}
	}
	update_residue_connection_mapping();
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

	// Branch lower connects should be treated like lower connects, so remove them from this list.
	if ( is_branch_lower_terminus() ) {
		// If this is a branch lower terminus, the branch lower connection SHOULD be at the first of the atoms found
		// above.  So return all but the first element of the vector.
		return utility::vector1< uint >( atoms.begin() + 1, atoms.end() );
	}

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


// General
// Number of ResidueConnections, counting polymeric residue connections
Size
ResidueType::n_possible_residue_connections() const
{
	return residue_connections_.size();
}

// Get a ResidueConection.
ResidueConnection const &
ResidueType::residue_connection( Size const i ) const
{
	return residue_connections_[ i ];
}

ResidueConnection &
ResidueType::residue_connection( Size const i )
{
	return residue_connections_[ i ];
}

Size
ResidueType::residue_connect_atom_index( Size const resconn_id ) const {
	return vd_to_index_.find( residue_connections_[ resconn_id ].vertex() )->second;
}


AtomType const &
ResidueType::atom_type( VD const vd ) const
{
	//TODO: Is there a better way of validating the VD instead of checking all the ordered atoms?
	PyAssert( (vd != ResidueGraph::null_vertex() ) &&
		( std::find(ordered_atoms_.begin(), ordered_atoms_.end(), vd) != ordered_atoms_.end()),
		"ResidueType::atom_type( VD const vd ): vd is not in this ResidueType!");
	debug_assert( (vd != ResidueGraph::null_vertex()) &&
		( std::find(ordered_atoms_.begin(), ordered_atoms_.end(), vd) != ordered_atoms_.end()) );
	return ( *atom_types_ )[ graph_[ vd ].atom_type_index() ];
}

/// @brief Get the atom name by index
std::string const &
ResidueType::atom_name( Size const index ) const
{
	PyAssert((index > 0) && (index <= ordered_atoms_.size()), "ResidueType::atom_name( Size const index ): index is not in this ResidueType!");
	debug_assert((index > 0) && (index <= ordered_atoms_.size()));
	return graph_[ ordered_atoms_[index] ].name();
}

/// @brief Get atom name by vertex descriptor
std::string const &
ResidueType::atom_name( VD vd ) const
{
	PyAssert( has(vd), "ResidueType::atom_name( VD vd ): vertex descriptor is not in this ResidueType!");
	debug_assert( has(vd) );
	return graph_[ vd ].name();
}

/// @brief get index of an atom's base atom
Size
ResidueType::atom_base( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <= ordered_atoms_.size()), "ResidueType::atom_base( Size const atomno ): atomno is not in this ResidueType!");
	return atom_base_indices_[atomno];
}

/// @brief get vd of an atom's base atom
VD
ResidueType::atom_base( VD const atomvd ) const
{
	PyAssert( has(atomvd), "ResidueType::atom_base( VD const atomvd ): atomvd is not in this ResidueType!");
	return atom_base_.find( atomvd )->second;
}

/// @brief get index of an atom's second base atom
Size
ResidueType::abase2( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <=  ordered_atoms_.size()), "ResidueType::abase2( Size const atomno ): atomno is not in this ResidueType!");
	return abase2_indices_[atomno];
}

/// @brief Counts the number of virtual atoms and returns the count.
/// @details The virtual count is not stored in the residue type.  This count is performed on the fly, and
///can hurt performance if repeatedly carried out.  Not intended for use in large loops -- instead, call
///once and store the value.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
Size
ResidueType::n_virtual_atoms () const
{
	core::Size virtcount = 0;
	for ( core::Size ia=1, iamax=natoms(); ia<=iamax; ++ia ) {
		if ( is_virtual(ia) ) ++virtcount;
	}
	return virtcount;
}

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

AdjacentIterPair
ResidueType::bonded_neighbor_iterators( VD const & atom ) const {
	return boost::adjacent_vertices( atom, graph_ );
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

/// @brief Return whether this atom is in a particular ring
bool
ResidueType::is_ring_atom( uint const ring_num, uint const atom_id ) const{
	debug_assert( ring_num <= ring_atoms_indices_.size() );
	return std::find(ring_atoms_indices_[ ring_num ].begin(), ring_atoms_indices_[ ring_num ].end(), atom_id) != ring_atoms_indices_[ ring_num ].end();
}

const HeavyAtomGraph
ResidueType::heavy_atoms(){
	HeavyAtomFilter filter(graph_, atom_types_);
	HeavyAtomGraph fg(graph_, boost::keep_all(), filter);
	return fg;
}

const AcceptorAtomGraph
ResidueType::acceptor_atoms(){
	AcceptorAtomFilter filter(graph_, atom_types_);
	AcceptorAtomGraph graph(graph_, boost::keep_all(), filter);
	return graph;
}

const HeavyAtomWithPolarHydrogensGraph
ResidueType::heavy_atom_with_polar_hydrogens(){
	HeavyAtomWithPolarHydrogensFilter filter(graph_, atom_types_);
	HeavyAtomWithPolarHydrogensGraph graph(graph_, boost::keep_all(), filter);
	return graph;
}

const HeavyAtomWithHydrogensGraph
ResidueType::heavy_atom_with_hydrogens(){
	HeavyAtomWithHydrogensFilter filter(graph_, atom_types_);
	HeavyAtomWithHydrogensGraph graph(graph_, boost::keep_all(), filter);
	return graph;
}


const HydrogenAtomGraph
ResidueType::hydrogens(){
	HydrogenAtomFilter filter(graph_, atom_types_);
	HydrogenAtomGraph fg(graph_, boost::keep_all(), filter);
	return fg;
}

const PolarHydrogenGraph
ResidueType::polar_hydrogens(){
	PolarHydrogenFilter filter(graph_, atom_types_);
	PolarHydrogenGraph fg(graph_, boost::keep_all(), filter);
	return fg;
}

const APolarHydrogenGraph
ResidueType::apolar_hydrogens(){
	APolarHydrogenFilter filter(graph_, atom_types_);
	APolarHydrogenGraph fg(graph_, boost::keep_all(), filter);
	return fg;
}

const AromaticAtomGraph
ResidueType::aromatic_atoms(){
	AromaticAtomFilter filter(graph_, atom_types_);
	AromaticAtomGraph fg(graph_, boost::keep_all(), filter);
	return fg;
}


/// @note this does not set xyz coordinates for the added atom
VD
ResidueType::add_atom(
	std::string const & atom_name,
	std::string const & atom_type_name,
	std::string const & mm_atom_type_name,
	Real const charge
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	debug_assert(atom_name_to_vd_.find(atom_name) == atom_name_to_vd_.end());
	debug_assert(atom_name_to_vd_.find( strip_whitespace(atom_name)) == atom_name_to_vd_.end());

	// the next calls will fail if the type name is unrecognized
	MMAtomTypeSetCOP mm_atom_types( mm_atom_types_ );
	Size const mm_type( mm_atom_types->atom_type_index( mm_atom_type_name ) );

	// Element information will be set by set_atom_type() below.
	ElementCOP element;

	Atom atom(
		atom_name,
		mm_atom_type_name,
		mm_type,
		element,
		charge,
		Vector(0.0)
	);

	VD v = graph_.add_vertex( atom );
	atom_name_to_vd_[ atom_name ] = v;
	atom_name_to_vd_[ strip_whitespace( atom_name ) ] = v;
	ordered_atoms_.push_back(v);

	// Less than here, as patching can remove atoms from the graph before it removes them from ordered_atoms_
	debug_assert( boost::num_vertices(graph_) <= ordered_atoms_.size() );

	// Atom const * graph_atom = &graph_[v];

	// Needed to reset the internal state.
	set_atom_type( v, atom_type_name );

	// allocate space for the new atom !!!!!!!!!!!!!!!!!!!!!!
	// eg, in the atom/resconn map
	debug_assert( atom_2_residue_connection_map_.size() == ordered_atoms_.size()-1 );

	atom_2_residue_connection_map_.resize( ordered_atoms_.size() );
	bonded_neighbor_.resize(ordered_atoms_.size());
	bonded_neighbor_type_.resize(ordered_atoms_.size());
	vd_to_index_[v] = ordered_atoms_.size();
	atom_base_[v] =v;  //graph_atom->atom_base(ordered_atoms_.size()); // base defaults to self
	icoor_[v] =  AtomICoor( atom_name, 0.0, 0.0, 0.0, atom_name, atom_name, atom_name, *this );

	return v;
}

VD
ResidueType::add_atom(
	std::string const & atom_name /* = "" */
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	VD v = graph_.add_vertex( Atom() );

	if ( atom_name.size() ) {
		debug_assert(atom_name_to_vd_.find(atom_name) == atom_name_to_vd_.end());
		debug_assert(atom_name_to_vd_.find( strip_whitespace(atom_name)) == atom_name_to_vd_.end());
		graph_[v].name( atom_name);
		atom_name_to_vd_[ atom_name ] = v;
		atom_name_to_vd_[ strip_whitespace( atom_name ) ] = v;
	}

	ordered_atoms_.push_back(v);

	// Less than here, as patching can remove atoms from the graph before it removes them from ordered_atoms_
	debug_assert( boost::num_vertices(graph_) <= ordered_atoms_.size() );

	// allocate space for the new atom !!!!!!!!!!!!!!!!!!!!!!
	// eg, in the atom/resconn map
	debug_assert( atom_2_residue_connection_map_.size() == ordered_atoms_.size()-1 );

	atom_2_residue_connection_map_.resize( ordered_atoms_.size() );
	bonded_neighbor_.resize(ordered_atoms_.size());
	bonded_neighbor_type_.resize(ordered_atoms_.size());
	vd_to_index_[v] = ordered_atoms_.size();
	atom_base_[v] =v;  // base defaults to self
	icoor_[v] =  AtomICoor( v, 0.0, 0.0, 0.0, v, v, v, *this );

	return v;
}

VD
ResidueType::add_atom(Atom const & atom, AtomICoor const & icoor){
	finalized_ = false;

	VD v = graph_.add_vertex(atom);
	std::string atom_name( atom.name());
	if ( atom_name.size() ) {
		debug_assert(atom_name_to_vd_.find(atom_name) == atom_name_to_vd_.end());
		debug_assert(atom_name_to_vd_.find( strip_whitespace(atom_name)) == atom_name_to_vd_.end());
		atom_name_to_vd_[ atom_name ] = v;
		atom_name_to_vd_[ strip_whitespace( atom_name ) ] = v;
	}

	ordered_atoms_.push_back(v);

	// allocate space for the new atom !!!!!!!!!!!!!!!!!!!!!!
	// eg, in the atom/resconn map
	debug_assert( atom_2_residue_connection_map_.size() == ordered_atoms_.size()-1 );

	atom_2_residue_connection_map_.resize( ordered_atoms_.size() );
	bonded_neighbor_.resize(ordered_atoms_.size());
	bonded_neighbor_type_.resize(ordered_atoms_.size());
	vd_to_index_[v] = ordered_atoms_.size();
	atom_base_[v] =v;  // base defaults to self
	icoor_[v] =  icoor;

	return v;

}


/// @brief flag an atom for deletion by adding its index to the delete_atom_ list
void
ResidueType::delete_atom( std::string const & name )
{
	debug_assert( has( name ) );
	delete_atom( atom_index(name) );
	atom_name_to_vd_.erase( atom_name_to_vd_.find(name) );
}

/// @brief flag an atom for deletion by adding its index to the delete_atom_ list
void
ResidueType::delete_atom( Size const index )
{
	finalized_ = false;

	// Delete any atom aliases mentioning this atom.
	// Do this first so we don't attempt to use any deleted atom information
	utility::vector1< std::string > aliases_to_delete;
	for ( auto const elem : atom_aliases_ ) {
		std::string const & rosetta_atom = strip_whitespace( elem.second );
		std::string const & alias        = elem.first;
		if ( elem.second == atom_name( index ) || rosetta_atom == atom_name( index ) ) {
			aliases_to_delete.push_back( alias );
		}
	}
	for ( auto const & alias : aliases_to_delete ) {
		delete_atom_alias( alias );
	}

	VD const vd = ordered_atoms_[index];
	graph_.clear_vertex(vd);
	graph_.remove_vertex(vd);

}

/// @brief Add an alias name for an atom.
void
ResidueType::add_atom_alias( std::string const & rosetta_atom, std::string const & alias ) {
	finalized_ = false;
	if ( ! has( rosetta_atom ) ) {
		utility_exit_with_message( "Unable to add atom alias for non-existent atom "+rosetta_atom );
	}
	std::string stripped_alias( strip_whitespace( alias ) );
	if ( stripped_alias.size() == 0 ) {
		utility_exit_with_message( "Cannot alias atom name to empty or all whitespace string." );
	}
	if ( has( stripped_alias ) || atom_aliases_.count( stripped_alias ) ) {
		utility_exit_with_message( "Cannot add atom alias, residue type already has an atom or alias named "+stripped_alias );
	}

	atom_aliases_[ alias ] = rosetta_atom;
	atom_aliases_[ stripped_alias ] = rosetta_atom;
}


/// @brief Remove a given alias name for an atom.
void
ResidueType::delete_atom_alias( std::string const & alias ) {
	finalized_ = false;
	if ( ! atom_aliases_.count(alias) ) {
		utility_exit_with_message("Cannot remove atom alias "+alias+" as it does not exist as an alias.");
	}
	atom_aliases_.erase( atom_aliases_.find(alias) );
	std::string stripped_alias( strip_whitespace( alias ) );
	if ( atom_aliases_.count(stripped_alias) ) { // Double check, as it might not be there, or alias==stripped_alias
		atom_aliases_.erase( atom_aliases_.find(stripped_alias) );
	}
}


/// @details Set atom type, correctly updating internal state.
/// @note    This method also sets/updates the atom's element and the residue's mass.
void
ResidueType::set_atom_type( VD atom, std::string const & atom_type_name )
{
	debug_assert( has( atom ) );
	Atom & a = graph_[ atom ];

	// If this is a re-setting, update H-bonding counts accordingly....
	core::uint const old_atom_type_index( a.atom_type_index() );
	if ( old_atom_type_index ) {
		if ( ( *atom_types_ )[ old_atom_type_index ].is_acceptor() ) { --n_hbond_acceptors_; }
		if ( ( *atom_types_ )[ old_atom_type_index ].is_donor() ) { --n_hbond_donors_; }
	}

	// ...and update mass accordingly.
	if ( a.element_type() ) {
		mass_ -= a.element_type()->weight();
	}

	// Get the new AtomType and its index.
	// (includes internal check for invalid type name)
	core::uint const atom_type_index = atom_types_->atom_type_index( atom_type_name );
	AtomType const & atype = ( *atom_types_ )[atom_type_index];

	// Set/update AtomType index.
	a.atom_type_index( atom_type_index );

	// Set/update H-bonding counts.
	if ( atype.is_acceptor() ) { ++n_hbond_acceptors_; }
	if ( atype.is_donor() ) { ++n_hbond_donors_; }

	// Set/update element and mass.
	if ( elements_ ) {  // Be robust if elements_ isn't defined.
		std::string const & element_name( ( *atom_types_ )[ atom_type_index ].element() );
		core::uint const element_index = elements_->element_index( element_name );
		//tr.Trace << elements_->element_index( element_name ) << ' ' << element_name << std::endl;
		a.element_type( ( *elements_ )[ element_index ] );
		mass_ += a.element_type()->weight();
	} else {
		tr.Warning << "Elements set undefined." << std::endl;
	}

	// Set/update the atom's properties.
	a.is_polar_hydrogen( atype.is_polar_hydrogen() );
	a.is_hydrogen( atype.is_hydrogen() );
	a.is_haro( atype.is_haro() );
	a.is_acceptor( atype.is_acceptor() );
	a.is_virtual( atype.is_virtual() );
	a.has_orbitals( atype.atom_has_orbital() );
}


void
ResidueType::set_mm_atom_type( std::string const & atom_name, std::string const & mm_atom_type_name )
{
	core::uint const mm_type_index = mm_atom_types_->atom_type_index( mm_atom_type_name );
	Atom & a = graph_[ atom_name_to_vd_[atom_name] ];
	a.mm_atom_type_index( mm_type_index );
	a.mm_name( mm_atom_type_name );
}

/// @brief Get the MM atom_type for this atom by its index number in this residue
MMAtomType const &
ResidueType::mm_atom_type( Size const atomno ) const
{
	MMAtomTypeSetCOP mm_atom_types( mm_atom_types_ );
	return ( *mm_atom_types )[graph_[ ordered_atoms_[atomno] ].mm_atom_type_index() ];
}

/// @brief Manually set the gasteiger typeset - will use the default set otherwise
void ResidueType::set_gasteiger_typeset( gasteiger::GasteigerAtomTypeSetCOP gasteiger_atom_types ) {
	gasteiger_atom_types_ = gasteiger_atom_types;
}

/// @brief set gasteiger atom type
void
ResidueType::set_gasteiger_atom_type(
	std::string const & atom_name,
	std::string const & gasteiger_atom_type_name
)
{
	set_gasteiger_atom_type( atom_name_to_vd_[atom_name] ,gasteiger_atom_type_name);
}

/// @brief set gasteiger atom type
void
ResidueType::set_gasteiger_atom_type(
	VD atom,
	std::string const & gasteiger_atom_type_name
)
{
	gasteiger::GasteigerAtomTypeDataCOP gasteiger_type;
	if ( gasteiger_atom_type_name == "" ) {
		gasteiger_type = nullptr;
	} else {
		if ( ! gasteiger_atom_types_ ) {
			gasteiger_atom_types_ = ChemicalManager::get_instance()->gasteiger_atom_type_set();
		}
		gasteiger_type = gasteiger_atom_types_->atom_type( gasteiger_atom_type_name );
	}
	Atom & a = graph_[ atom ];
	a.gasteiger_atom_type( gasteiger_type );
}

/// @brief Get the Gasteiger atom_type for this atom by its index number in this residue
gasteiger::GasteigerAtomTypeDataCOP
ResidueType::gasteiger_atom_type( Size const atomno ) const
{
	return graph_[ ordered_atoms_[atomno] ].gasteiger_atom_type();
}


gasteiger::GasteigerAtomTypeSetCOP ResidueType::gasteiger_atom_typeset() const {
	return gasteiger_atom_types_;
}

VD
ResidueType::atom_vertex( std::string const & name) const{
	auto atom_name_to_vd_iter( atom_name_to_vd_.find( name ) );
	if ( atom_name_to_vd_iter == atom_name_to_vd_.end() ) {
		tr.Error << "atom name : '" << name << "' not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		utility_exit_with_message("Atom name not found in ResidueType.");
	}

	return atom_name_to_vd_iter->second;
}


void
ResidueType::dump_vd_info() const {
	tr << "Residue " << name() << std::endl;

	for ( core::Size ii(1); ii <= ordered_atoms_.size(); ++ii ) {
		tr << " atom index " << ii << ": " << ordered_atoms_[ii] << " In graph? " << has(ordered_atoms_[ii]) << std::endl;
	}
	tr << "-------------" << std::endl;
	VIterPair allverts( boost::vertices( graph_ ) );
	for ( auto iter(allverts.first); iter != allverts.second; ++iter ) {
		tr << " atom " << atom_name(*iter) << " vd: " << *iter << std::endl;
	}
	tr << "-------------" << std::endl;
}

/// @brief number of orbitals
Size
ResidueType::n_orbitals() const
{
	return orbitals_.size();
}

orbitals::OrbitalType const &
ResidueType::orbital_type(int const orbital_index)const
{
	orbitals::OrbitalTypeSetCOP orbital_types( orbital_types_ );
	return ( *orbital_types )[ orbitals_[ orbital_index ].orbital_type_index() ];
}

// Return a pointer to the object containing the set of ring conformers possible for this residue's nth cycle.
core::chemical::rings::RingConformerSetCOP
ResidueType::ring_conformer_set( core::uint ring_num ) const
{
	return conformer_sets_[ ring_num ];
}

void
ResidueType::clear_orbitals()
{
	orbitals_.clear();
	orbitals_index_.clear();
	for ( Size index=1; index<= natoms(); ++index ) {
		graph_[ordered_atoms_[index]].bonded_orbitals().clear();
	}
}

/// @note this does not set xyz coordinates for the added orbital but sets the index of the orbital and maps
/// it to the type of orbital.
void
ResidueType::add_orbital(
	std::string & orbital_name,
	std::string & orbital_type_name
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	// store the atom type
	// the next call will fail if the orbital type name is unrecognized
	orbitals::OrbitalTypeSetCOP orbital_types( orbital_types_ );
	Size type( orbital_types->orbital_type_index( orbital_type_name ) );

	// store the name
	orbitals_.push_back(Orbital(orbital_name, type, Vector(0.0)));

	orbitals_index_[ orbital_name ] = orbitals_.size();
	orbitals_index_[ strip_whitespace( orbital_name ) ] = orbitals_.size();
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Add an atom to the list of atoms that can potentially form a bond to a metal ion.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ResidueType::add_metalbinding_atom (
	std::string const & atom_name
) {
	if ( !has(atom_name) ) {
		std::string message = "Error in adding metal-binding atom to residue type " + name3() + ". Atom " + atom_name + " was not found.";
		utility_exit_with_message(message);
	}
	metal_binding_atoms_.push_back( atom_name ); //Store names rather than indices, since indices might change.
	return;
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Remove an atom from the list of atoms that can potentially form a bond to a metal ion
/// (used in patching when it kills the valence that is thus used)
/// @author Andrew Watkins (amw579@nyu.edu)
void
ResidueType::delete_metalbinding_atom (
	std::string const & atom_name
) {
	if ( !has( atom_name ) ) {
		std::string message = "Error in removing metal-binding atom from residue type " + name3() + ". Atom " + atom_name + " was not found.";
		utility_exit_with_message(message);
	}
	metal_binding_atoms_.erase( std::remove( metal_binding_atoms_.begin(), metal_binding_atoms_.end(), atom_name ),
		metal_binding_atoms_.end() );

	return;
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Remove an atom from the list of act coord atoms
/// (used in patching when it kills the valence that is thus used)
/// @author Andrew Watkins (amw579@nyu.edu)
void
ResidueType::delete_act_coord_atom (
	std::string const & atom_name
) {
	if ( !has( atom_name ) ) {
		std::string message = "Error in removing act coord atom from residue type " + name3() + ". Atom " + atom_name + " was not found.";
		utility_exit_with_message(message);
	}

	actcoord_atoms_.erase( std::remove( actcoord_atoms_.begin(), actcoord_atoms_.end(), atom_name_to_vd_[ atom_name ] ),
		actcoord_atoms_.end() );
	actcoord_atoms_indices_.erase( std::remove( actcoord_atoms_indices_.begin(), actcoord_atoms_indices_.end(), atom_index( atom_name ) ),
		actcoord_atoms_indices_.end() );

	return;
}

///////////////////////////////////////////////////////////////////////////////

/// @details add a bond between atom1 and atom2 and add a BondType object referencing the bond using the specified bondName
void
ResidueType::add_bond(std::string const & atom_name1, std::string const & atom_name2, BondName bondLabel /*=SingleBond*/)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
		std::string message = "add_bond: atoms " + atom_name1 + " and " + atom_name2 + " don't exist!";
		utility_exit_with_message( message  );
	}

	/////// Standard Version /////////

	Size const i1( atom_index( atom_name1 ) );
	Size const i2( atom_index( atom_name2 ) );

	if ( bonded_neighbor_.size() < Size( std::max( i1, i2) ) ) {
		// ensure enough space for nbrs
		utility_exit_with_message("ResidueType:: shouldnt get here -- resizing in add_bond");
	}

	bonded_neighbor_[ i1 ].push_back( i2 );
	bonded_neighbor_[ i2 ].push_back( i1 );
	bonded_neighbor_type_[i1].push_back(bondLabel);
	bonded_neighbor_type_[i2].push_back(bondLabel);  //bondType_vector_.push_back(BondType(i1,i2,bondLabel));

	NameVDMap::const_iterator source = atom_name_to_vd_.find( atom_name1 );
	NameVDMap::const_iterator target = atom_name_to_vd_.find( atom_name2 );
	debug_assert( source != atom_name_to_vd_.end());
	debug_assert( target != atom_name_to_vd_.end());
	VD const vd_source = source->second;
	VD const vd_target = target->second;

	// check if bond already exists...
	if ( boost::edge(vd_source, vd_target, graph_).second ) {
		utility_exit_with_message( "dont add residue bonds more than once!" );
	}

	ResidueGraph::edge_descriptor e_added;
	bool added;
	boost::tie(e_added, added) = graph_.add_edge( vd_source, vd_target, Bond(-1, bondLabel)); /// -1 means Bond distance not set here. This will be fixed in the future
	debug_assert(added);
}

// @details add a bond between atom1 and atom2 and add a BondType object referencing the bond using the specified bondName
void
ResidueType::add_bond(VD atom1, VD atom2, BondName bondLabel /*=SingleBond*/)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom1 ) || !has( atom2 ) ) {
		utility_exit_with_message( "add_bond: atoms don't exist!" );
	}

	/////// Standard Version /////////

	Size const i1( atom_index( atom1 ) );
	Size const i2( atom_index( atom2 ) );

	if ( bonded_neighbor_.size() < Size( std::max( i1, i2) ) ) {
		// ensure enough space for nbrs
		utility_exit_with_message("ResidueType:: shouldn't get here -- resizing in add_bond");
	}

	bonded_neighbor_[ i1 ].push_back( i2 );
	bonded_neighbor_[ i2 ].push_back( i1 );
	bonded_neighbor_type_[i1].push_back(bondLabel);
	bonded_neighbor_type_[i2].push_back(bondLabel);  //bondType_vector_.push_back(BondType(i1,i2,bondLabel));

	// check if bond already exists...
	if ( boost::edge(atom1, atom2, graph_).second ) {
		utility_exit_with_message( "dont add residue bonds more than once!" );
	}

	ResidueGraph::edge_descriptor e_added;
	bool added;
	boost::tie(e_added, added) = graph_.add_edge( atom1, atom2, Bond(-1, bondLabel)); /// -1 means Bond distance not set here. This will be fixed in the future
	debug_assert(added);
}

// Change the bond type of the given bond from one type to another.
/// @author   Labonte <JWLabonte@jhu.edu>
void
ResidueType::change_bond_type(
	std::string const & atom_name1,
	std::string const & atom_name2,
	BondName const old_bond_label,
	BondName const new_bond_label )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	if ( ! ( has( atom_name1 ) || has( atom_name2 ) ) ) {
		utility_exit_with_message( "change_bond_type: atoms " + atom_name1 + " and " + atom_name2 + " don't exist!" );
	}

	// First, change the bond types stored for each atom across the bond.
	core::uint const atom1( atom_index( atom_name1 ) );
	core::uint const atom2( atom_index( atom_name2 ) );

	utility::vector1< BondName >::iterator it1, it2;
	it1 = std::find( bonded_neighbor_type_[ atom1 ].begin(), bonded_neighbor_type_[ atom1 ].end(), old_bond_label );
	if ( it1 == bonded_neighbor_type_[ atom1 ].end() ) {
		utility_exit_with_message(
			"change_bond_type: atom " + atom_name1 + " does not have the requested bond type!" );
	}
	it2 = std::find( bonded_neighbor_type_[ atom2 ].begin(), bonded_neighbor_type_[ atom2 ].end(), old_bond_label );
	if ( it1 == bonded_neighbor_type_[ atom2 ].end() ) {
		utility_exit_with_message(
			"change_bond_type: atom " + atom_name2 + " does not have the requested bond type!" );
	}

	bonded_neighbor_type_[ atom1 ].erase( it1 );
	bonded_neighbor_type_[ atom2 ].erase( it2 );
	bonded_neighbor_type_[ atom1 ].push_back( new_bond_label );
	bonded_neighbor_type_[ atom2 ].push_back( new_bond_label );

	// Now, change the bond type stored in the edge of the atom graph.
	NameVDMap::const_iterator source = atom_name_to_vd_.find( atom_name1 );
	NameVDMap::const_iterator target = atom_name_to_vd_.find( atom_name2 );
	debug_assert( source != atom_name_to_vd_.end() );
	debug_assert( target != atom_name_to_vd_.end() );
	VD const vd_source = source->second;
	VD const vd_target = target->second;
	graph_.remove_edge( vd_source, vd_target );
	graph_.add_edge( vd_source, vd_target, Bond( -1, new_bond_label ) ); /// -1 means Bond distance not set here.
}

///////////////////////////////////////////////////////////////////////////////
/// @brief add an orbital bond between an atom and an orbital.
/// @note NOTE!!!!! This is indexed based upon atoms, not orbitals. That means that in your params file
/// you must have the atom as the first and orbital as the second.
void
ResidueType::add_orbital_bond(
	std::string const & atom_name1,
	std::string const & orbital_name
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has_orbital( orbital_name ) ) {
		std::string message = "add_bond: atoms " + atom_name1 + " and " + orbital_name + " dont exist!";
		utility_exit_with_message( message  );
	}

	Size const i2( orbitals_index_[ orbital_name ] );

	if ( atom_name_to_vd_.find(atom_name1) == atom_name_to_vd_.end() ) {
		utility_exit_with_message("atom_name: " + atom_name1 +" not found. Improper params file!");

	}

	Size const i1( atom_index( atom_name1 ) );

	graph_[ordered_atoms_[i1]].bonded_orbitals().push_back(i2);

}

orbitals::ICoorOrbitalData const &
ResidueType::orbital_icoor_data(Size const orbital_index) const{
	return orbitals_[orbital_index].icoor();
}


orbitals::ICoorOrbitalData const &
ResidueType::new_orbital_icoor_data(Size const orbital_index) const{
	return orbitals_[orbital_index].new_icoor();
}


/// @details add a cut_bond between atom1 and atom2, which disallows an atom-tree connection,
///            though the atoms are really bonded.
/** update cut_bond_ and resize it as necessary **/
void
ResidueType::add_cut_bond(
	std::string const & atom_name1,
	std::string const & atom_name2
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
		std::string message = "add_cut_bond: atoms " + atom_name1 + " and " + atom_name2 + " don't exist!";
		utility_exit_with_message( message  );
	}

	VD const vd_source = atom_name_to_vd_.find( atom_name1 )->second; //source->second;
	VD const vd_target = atom_name_to_vd_.find( atom_name2 )->second; //target->second;

	std::pair< ED, bool> ed_present_pair( boost::edge(vd_source,vd_target,graph_) );
	if ( ! ed_present_pair.second ) {
		//TODO: Change cut_bond_neighbor_ to derived data status once we make this an error.
		tr.Warning << "Cut bond set for bond which does not exist " << atom_name1 << " -- " << atom_name2 << std::endl;
	} else {
		graph_[ ed_present_pair.first ].cut_bond( true );
	}

	std::map<VD, utility::vector1<VD> >::const_iterator i1_nbrs_iter = cut_bond_neighbor_.find(vd_source);
	if ( i1_nbrs_iter != cut_bond_neighbor_.end() ) {
		utility::vector1<VD> const i1_nbrs(i1_nbrs_iter->second);
		if ( std::find( i1_nbrs.begin(), i1_nbrs.end(), vd_target ) != i1_nbrs.end() ) {
			utility_exit_with_message( "don't add cut bonds more than once!" );
		}
	}
	cut_bond_neighbor_[vd_source].push_back(vd_target);
	cut_bond_neighbor_[vd_target].push_back(vd_source);
}


///////////////////////////////////////////////////////////////////////////////

// Add a chi (side-chain) angle defined by four atoms.
void
ResidueType::add_chi( Size const chino,
	VD atom1,
	VD atom2,
	VD atom3,
	VD atom4)
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	if ( !has( atom1 ) || !has( atom2 ) ||
			!has( atom3 ) || !has( atom4 ) ) {
		utility_exit_with_message("ResidueType::add_chi: atoms don't exist!" );
	}

	utility::vector1<VD> atoms;
	atoms.push_back( atom1 );
	atoms.push_back( atom2 );
	atoms.push_back( atom3 );
	atoms.push_back( atom4 );
	if ( chi_atoms_.size() < chino ) {
		chi_atoms_.resize( chino );
		chi_rotamers_.resize( chino );
		is_proton_chi_.resize( chino );
		chi_2_proton_chi_.resize( chino );
	}
	chi_atoms_[chino] = atoms;
	// Should we be clearing chi_rotamers_ here?
	// No, I don't think so. ~Labonte
	is_proton_chi_[chino] = false;
	chi_2_proton_chi_[ chino ] = 0;
}  // add_chi

void
ResidueType::add_chi(VD atom1,
	VD atom2,
	VD atom3,
	VD atom4)
{
	add_chi(nchi() + 1, atom1, atom2, atom3, atom4);
}

// Add a chi (side-chain) angle defined by four atoms.
void
ResidueType::add_chi(
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
/// @details This method is needed for combinatorial patching of ResidueTypes for which the number of chis is variable.
/// Its primary purpose is to be used with add_chi_rotamer_to_last_chi() that adds rotamer bins to the last chi in the
/// list.  In this way, a new chi can be added by a patch file and its rotamer bins set without needing to designate a
/// chi index.
/// @note    See also add_chi_rotamer_to_last_chi().
/// @author  Labonte
void
ResidueType::add_chi(std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4)
{
	add_chi(nchi() + 1, atom_name1, atom_name2, atom_name3, atom_name4);
}


// Add a nu (internal cyclic) angle defined by four atoms.
void
ResidueType::add_nu( core::uint const nu_index,
	std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4 )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	if ( ! has( atom_name1 ) || ! has( atom_name2 ) || ! has( atom_name3 ) || ! has( atom_name4 ) ) {
		utility_exit_with_message( "ResidueType::add_nu: Requested atoms don't exist in this ResidueType!" );
	}

	utility::vector1< VD > atoms;
	atoms.push_back( ordered_atoms_[ atom_index( atom_name1 ) ] );
	atoms.push_back( ordered_atoms_[ atom_index( atom_name2 ) ] );
	atoms.push_back( ordered_atoms_[ atom_index( atom_name3 ) ] );
	atoms.push_back( ordered_atoms_[ atom_index( atom_name4 ) ] );

	if ( nu_atoms_.size() < nu_index ) {
		nu_atoms_.resize( nu_index );
	}
	nu_atoms_[ nu_index ] = atoms;
}


// Add a ring definition.
void
ResidueType::add_ring( core::uint const ring_num, utility::vector1< std::string > const & ring_atoms )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	Size const ring_size( ring_atoms.size() );
	utility::vector1< VD > atoms( ring_size );
	for ( uint i( 1 ); i <= ring_size; ++i ) {
		if ( ! has( ring_atoms[ i ] ) ) {
			utility_exit_with_message( "ResidueType::add_ring: Requested atoms don't exist in this ResidueType!" );
		}
		atoms[ i ] = ordered_atoms_[ atom_index( ring_atoms[ i ] ) ];
	}

	if ( ring_atoms_.size() < ring_num ) {
		ring_atoms_.resize( ring_num );
	}
	if ( lowest_ring_conformer_.size() < ring_num ) {
		lowest_ring_conformer_.resize( ring_num );
	}
	if ( low_ring_conformers_.size() < ring_num ) {
		low_ring_conformers_.resize( ring_num );
	}
	ring_atoms_[ ring_num ] = atoms;
}

// Set this cyclic residue's lowest-energy ring conformer for the nth ring by IUPAC name.
void
ResidueType::set_lowest_energy_ring_conformer( core::uint const ring_num, std::string const & conformer )
{
	// Signal that we need to update the derived data.
	finalized_ = false;
	lowest_ring_conformer_[ ring_num ] = conformer;
}

// Set this cyclic residue's low-energy ring conformers for the nth ring by IUPAC name.
void
ResidueType::set_low_energy_ring_conformers( core::uint const ring_num, utility::vector1< std::string > const & conformers )
{
	// Signal that we need to update the derived data.
	finalized_ = false;
	low_ring_conformers_[ ring_num ] = conformers;
}


/// @details Describe proton behavior for residue type; where should rotamer samples be considered,
/// and if expanded rotamers are desired, what deviations from the original rotamer samples
/// should be included.
/// E.g. dihedral_samples of 60, -60, and 180 could have an extra_sample of
/// 20 which would produce rotamers at 40 60 & 80, -40 -60 & -80, and -160, 180 & 160.
/// Extra_samples at 10 and 20 would produce 15 different rotamer samples.
void
ResidueType::set_proton_chi(
	Size chino,
	utility::vector1< Real > const & dihedral_samples,
	utility::vector1< Real > const & extra_samples
)
{
	debug_assert( is_proton_chi_.size() >= nchi() );
	debug_assert( chi_2_proton_chi_.size() >= nchi() );
	if ( chino > nchi() ) {
		utility_exit_with_message("Error setting proton chi: Chi to set as proton chi does not exist.");
	}
	if ( is_proton_chi_[ chino ] ) {
		debug_assert( std::find(proton_chis_.begin(),proton_chis_.end(),chino) != proton_chis_.end() );
		core::Size h_chi( chi_2_proton_chi_[ chino ] );
		proton_chi_samples_[ h_chi ] = dihedral_samples;
		proton_chi_extra_samples_[ h_chi ] = extra_samples;
	} else { // New proton chi
		is_proton_chi_[ chino ] = true;
		proton_chis_.push_back( chino );
		proton_chi_samples_.push_back( dihedral_samples );
		proton_chi_extra_samples_.push_back( extra_samples );
		chi_2_proton_chi_[ chino ] = proton_chis_.size();
	}
}

///////////////////////////////////////////////////////////////////////////////

// Add a rotamer bin for a given chi.
/// @details A rotamer bin has the mean and standard deviation.
void
ResidueType::add_chi_rotamer(
	Size const chino,
	Real const mean,
	Real const sdev
)
{
	if ( chi_rotamers_.size() < chino ) { chi_rotamers_.resize( chino ); }
	chi_rotamers_[chino].push_back( std::make_pair( mean, sdev ) );
}

// Adds a chi rotamer bin to the highest-indexed chi in the list of chis for this ResidueType.
/// @details This method is needed for combinatorial patching of ResidueTypes for which the number of chis is variable.
/// Its primary purpose is to be used with the overloaded version of add_chi() that adds a new chi to the end of the
/// list.  In this way, a new chi can be added by a patch file and its rotamer bins set without needing to designate a
/// chi index.
/// @note    See also add_chi().
/// @author  Labonte
void
ResidueType::add_chi_rotamer_to_last_chi(core::Angle const mean, core::Angle const sdev)
{
	add_chi_rotamer(nchi(), mean, sdev);
}


// Delete all of the chi rotamer bins from the specified chi for this ResidueType.
/// @details This method is useful if one has redefined a chi within a patch file such that the old rotamer bins need
/// to be regenerated.
/// @author  Labonte <JWLabonte@jhu.edu>
void
ResidueType::clear_chi_rotamers( core::uint const chi_no )
{
	chi_rotamers_[ chi_no ].clear();
}


/// @brief Regenerate the rotatable chi bonds from the internal graph structure.
/// If the number of proton chi samples would exceed max_proton_chi_samples, don't add extra sampling to proton chis.
/// As a special case, if this is zero don't add any proton chi sampling at all.
///
/// Requires that Icoor and atom base records are up-to-date, and that ring bonds have been annotated.
void
ResidueType::autodetermine_chi_bonds( core::Size max_proton_chi_samples ) {
	// First, clear off the current internal information on rotatable chis;
	chi_atoms_.clear();
	chi_rotamers_.clear();
	is_proton_chi_.clear();
	// TODO: Notate if we have a "non-standard" setting for any of the proton chis
	proton_chis_.clear();
	chi_2_proton_chi_.clear();
	proton_chi_samples_.clear();
	proton_chi_extra_samples_.clear();
	// signal that we need to update the derived data TODO: Is this necessary?
	finalized_ = false;

	// As far as I can tell, there isn't any specified ordering for chis.
	// The canonical residues go from root out, but for ligands
	// there doesn't look to be any ordering guarantee.
	utility::vector1<VDs> found_chis( core::chemical::find_chi_bonds( *this ) );
	utility::vector1< core::Size > proton_chis; // Not the member varaible as set_proton_chi modifies that.

	if ( is_protein() ) {

		utility::vector1<VDs> true_chis; // filtered and ordered from found_chis.
		// Note that this algorithm to get down to the 'real' chis is pretty
		// gross, but when N is < 10 most reasonable big-Os are fine, right?

		// Step 1. Get chi1 (it's the one with N as first or fourth)
		// Other criterion -- atom 3 can't be C (it'll find N CA C O)
		for ( VDs const & chi : found_chis ) {
			tr.Trace << "looking at found chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
			if ( atom_name( chi[ 1 ] ) == "N" && atom_name( chi[ 3 ] ) != "C" ) {
				true_chis.push_back( chi );

				// Third atom of chi1 is first sidechain atom.
				for ( Size ii = 1; ii < atom_index( atom_name( chi[ 3 ] ) ); ++ii ) {
					set_backbone_heavyatom( atom_name( ii ) );
				}
				break;
			}


			// In prior versions of this code, this was necessary because we
			// hadn't re-rooted on N in assign_internal_coordinates (generally
			// called before this function).
			// Now that we've properly re-rooted on N, this probably won't happen.
			/*if ( atom_name( chi[ 4 ] ) == "N" && atom_name( chi[ 2 ] ) != "C" ) {
			VDs reversed_chi;
			reversed_chi.push_back( chi[ 4 ] );
			reversed_chi.push_back( chi[ 3 ] );
			reversed_chi.push_back( chi[ 2 ] );
			reversed_chi.push_back( chi[ 1 ] );
			true_chis.push_back( reversed_chi );
			}*/
		}

		// Step 2. Get remainder of chis by asking each one to start with the
		// second atom of the prior chi[s]. Note that this will potentially
		// confuse branches, but branched sidechains with lots of chis are treated
		// poorly by essentially any chi system.
		std::string target_first_atom = atom_name( true_chis[ 1 ][ 2 ] );
		while ( true ) {

			// this extra loop is to future-proof a bit against branching: multiple
			// chis per pass may start with the target_first_atom and therefore
			// we don't want to update it right away.

			std::string candidate_new_atom = target_first_atom;
			for ( VDs const & chi : found_chis ) {
				tr.Trace << "looking at found chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
				if ( atom_name( chi[ 1 ] ) == target_first_atom ) {
					true_chis.push_back( chi );
					candidate_new_atom = atom_name( chi[ 2 ] );
				}

				// Now that we've properly re-rooted on N, this probably won't happen.
				/*if ( atom_name( chi[ 4 ] ) == target_first_atom ) {
				VDs reversed_chi;
				reversed_chi.push_back( chi[ 4 ] );
				reversed_chi.push_back( chi[ 3 ] );
				reversed_chi.push_back( chi[ 2 ] );
				reversed_chi.push_back( chi[ 1 ] );
				candidate_new_atom = atom_name( reversed_chi[ 2 ] );
				true_chis.push_back( reversed_chi );
				}*/
			}
			if ( candidate_new_atom == target_first_atom ) break;

			// This may have to become a vector -- where we accumulated many
			// candidate_new_atom -- later.
			target_first_atom = candidate_new_atom;
		}

		for ( VDs const & chi : true_chis ) {
			tr.Debug << "looking at true chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
			debug_assert( chi.size() == 4 );
			add_chi( chi[1], chi[2], chi[3], chi[4] );
			if ( atom( chi[4] ).element_type()->element() == core::chemical::element::H ) {
				// proton chi
				proton_chis.push_back( nchi() );
			}
		} // for all found chis
	} else if ( is_RNA() ) {
		utility::vector1<VDs> true_chis; // filtered and ordered from found_chis.
		
		//CHI 1 C2' C1' N9  C4
		//CHI 2 C4' C3' C2' C1'
		//CHI 3 C3' C2' C1' N9 
		//CHI 4 C3' C2' O2' HO2'
		// First base atom is either N1 or N9
		
		VD first_base_atom = atom_vertex( "N9" );
		for ( VDs const & chi : found_chis ) {
			tr.Trace << "looking at found chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
			if ( atom_name( chi[ 1 ] ) == "C2'" && atom_name( chi[ 2 ] ) != "O2'" ) {
				true_chis.push_back( chi );
				first_base_atom = chi[3];
				// Third atom of chi1 is first base atom.
				for ( Size ii = 1; ii < atom_index( atom_name( chi[ 3 ] ) ); ++ii ) {
					set_backbone_heavyatom( atom_name( ii ) );
				}
				break;
			}
		}
		
		// Step 2. Hard-fix three chis: two rings, and proton chi for HO2'.
		VDs chi{atom_vertex("C4'"), atom_vertex("C3'"), atom_vertex("C2'"), atom_vertex("C1'")};
		true_chis.emplace_back( chi );
		chi = VDs{ atom_vertex("C3'"), atom_vertex("C2'"), atom_vertex("C1'"), first_base_atom };
		//true_chis.emplace_back( { atom_vertex("C4'"), atom_vertex("C3'"), atom_vertex("C2'"), atom_vertex("C1'") } );
		true_chis.emplace_back( chi );
		// What to do absent HO2'?
		// answer: whatever else O2' is bonded to that's not C2'
		if ( has( "HO2'" ) ) {
			//chi = ;
			true_chis.emplace_back( VDs{ atom_vertex("C3'"), atom_vertex("C2'"), atom_vertex("O2'"), atom_vertex("HO2'") } );
		} else {
			// implement later
		}
		
		// Theoretical final step (AMW TODO): add all chis that are children of the base
		// or of O2', in a protein-y way, as chis 5+.
		
		std::string target_first_atom = atom_name( true_chis[ 1 ][ 2 ] );
		
		while ( true ) {
			
			// this extra loop is to future-proof a bit against branching: multiple
			// chis per pass may start with the target_first_atom and therefore
			// we don't want to update it right away.
			
			std::string candidate_new_atom = target_first_atom;
			for ( VDs const & chi : found_chis ) {
				tr.Trace << "looking at found chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
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
			tr.Debug << "looking at true chi: " << atom_name( chi[1] ) << " " << atom_name( chi[2] ) << " " << atom_name( chi[3] ) << " " << atom_name( chi[4] ) << std::endl;
			debug_assert( chi.size() == 4 );
			add_chi( chi[1], chi[2], chi[3], chi[4] );
			if ( atom( chi[4] ).element_type()->element() == core::chemical::element::H ) {
				// proton chi
				proton_chis.push_back( nchi() );
			}
		} // for all found chis

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
		tr.Warning << "Number of base proton chi samples (" << num_H_confs << ") for " << name() << " exceeds requested number of samples" << std::endl;
	} else if ( 3* num_H_confs > max_proton_chi_samples ) {
		tr << "Skipping extra samples for proton chis on " << name() << "; would give " << num_H_confs << " conformers." << std::endl;
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


///////////////////////////////////////////////////////////////////////////////

/// @details sets atom_base_[atom1] = atom2
/** resize atom_base_ vector as necessary **/
void
ResidueType::set_atom_base(
	std::string const & atom_name1,
	std::string const & atom_name2
)
{
	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
		if ( ! has( atom_name1 ) ) {
			tr.Error << "Atom '" << atom_name1 << "' does not exist in restype " << name() << std::endl;
		}
		if ( ! has( atom_name2 ) ) {
			tr.Error << "Atom '" << atom_name2 << "' does not exist in restype " << name() << std::endl;
		}
		utility_exit_with_message( "set_atom_base: atom names don't exist!" );
	}

	VD const vd_source = atom_name_to_vd_.find( atom_name1 )->second;
	VD const vd_target = atom_name_to_vd_.find( atom_name2 )->second;

	set_atom_base( vd_source, vd_target );
}

/// @details sets atom_base_[atom1] = atom2
/** resize atom_base_ vector as necessary **/
void
ResidueType::set_atom_base(
	VD const & atom1,
	VD const & atom2
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom1 ) || !has( atom2 ) ) {
		utility_exit_with_message( "set_atom_base: atom vertex descriptors don't exist!" );
	}

	// atom base must be bonded (or the same atom)
	if ( atom1 != atom2 && !boost::edge(atom1, atom2, graph_).second ) {
		tr.Error << "Atoms " << atom_name(atom1) << " and " << atom_name( atom2 ) << " are not bonded. " << std::endl;
		utility_exit_with_message( "in set_atom_base(), atoms must be bonded to set atom base!" );
	}

	//make sure that you do not set an atom base at a cut bond
	utility::vector1<VD> const i1_nbrs(cut_bond_neighbor_[atom1] );
	if ( std::find( i1_nbrs.begin(), i1_nbrs.end(), atom2 ) != i1_nbrs.end() ) {
		utility_exit_with_message( "Don't set atom bases to cut bonds!" );
	}

	atom_base_[atom1] = atom2;
}

/// @brief set indices of all mainchain atoms
void
ResidueType::set_mainchain_atoms( AtomIndices const & mainchain )
{
	utility::vector1<VD> vd_mainchain;
	for ( Size i=1; i<= mainchain.size(); ++i ) {
		vd_mainchain.push_back(ordered_atoms_[mainchain[i]]);
	}
	mainchain_atoms_ = vd_mainchain;
}


///////////////////////////////////////////////////////////////////////////////
// Property-Related Methods
///////////////////////////////////////////////////////////////////////////////

ResidueProperties const &
ResidueType::properties() const
{
	return *properties_;
}

void
ResidueType::set_properties( ResiduePropertiesOP properties ) {
	properties_ = properties;
}

void
ResidueType::add_property( std::string const & property )
{
	// signal that we need to update the derived data.
	finalized_ = false;

	properties_->set_property( property, true );

	// AMW: as of 8/18/15, L_AA and D_AA no longer imply ALPHA_AA.

	// Special "umbrella cases"
	// FIXME: There really shouldn't be as many umbrella cases, IMO. ~Labonte
	if ( property == "PROTEIN" ) {
		properties_->set_property( POLYMER, true );
	} else if ( property == "ALPHA_AA" ) {
		properties_->set_property( PROTEIN, true );
		properties_->set_property( POLYMER, true );
	} else if ( property == "BETA_AA" ) {
		properties_->set_property( PROTEIN, true );
		properties_->set_property( POLYMER, true );
	} else if ( property == "L_AA" ) {
		properties_->set_property( PROTEIN, true );
		properties_->set_property( POLYMER, true );
		//properties_->set_property( ALPHA_AA, true );
	} else if ( property == "D_AA" ) {
		properties_->set_property( PROTEIN, true );
		properties_->set_property( POLYMER, true );
		//properties_->set_property( ALPHA_AA, true );
	} else if ( property == "DNA" ) {
		properties_->set_property( POLYMER, true );
	} else if ( property == "RNA" ) {
		properties_->set_property( POLYMER, true );
	} else if ( property == "PEPTOID" ) {
		properties_->set_property( POLYMER, true );
		// amw: This will no longer be true if we incorporate
		// alanine peptoids (or update our chirality model entirely)
		properties_->set_property( ACHIRAL_BACKBONE, true );
	} else if ( property == "LOWERTERM_TRUNC" ) {
		properties_->set_property( LOWER_TERMINUS, true );
	} else if ( property == "UPPERTERM_TRUNC" ) {
		properties_->set_property( UPPER_TERMINUS, true );
	} else if ( property == "PHOSPHONATE" ) {
		properties_->set_property( POLYMER, true );
	} else if ( property == "PHOSPHONATE_UPPER" ) {
		properties_->set_property( UPPER_TERMINUS, true );
		properties_->set_property( PHOSPHONATE, true );
	} else if ( property == "ACETYLATED_NTERMINUS" ) {
		properties_->set_property( LOWER_TERMINUS, true );
	} else if ( property == "METHYLATED_CTERMINUS" ) {
		properties_->set_property( UPPER_TERMINUS, true );
	}
}

void
ResidueType::set_adduct_flag( bool adduct_in )
{
	properties_->set_property( ADDUCT, adduct_in );
}


void
ResidueType::add_numeric_property(std::string const & tag, core::Real value)
{
	properties_->add_numeric_property( tag, value );
}

void
ResidueType::add_string_property(std::string const & tag, std::string value)
{
	properties_->add_string_property( tag, value );
}

/// @details This is needed for deleting properties, which occurs in certain PTMs.
void
ResidueType::delete_property( std::string const & property )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	properties_->set_property( property, false );
}

/// @brief Is this ResidueType a base type?
/// @details Checks the base_type_cop_ pointer.  If it's null, this is assumed to be a base type.
bool
ResidueType::is_base_type() const {
	return bool( !base_type_cop_ );
}

/// @brief Get a pointer to this ResidueType's base ResidueType.
/// @details Returns the base_type_cop_ pointer if not null, self pointer if null.
ResidueTypeCOP
ResidueType::get_base_type_cop() const {
	if ( base_type_cop_ ) {
		return base_type_cop_;
	}
	return get_self_ptr();
}

/// @brief Reset the base type COP to be null.  This implies that this ResidueType is a base type.
///
void
ResidueType::reset_base_type_cop() {
	base_type_cop_ = ResidueTypeCOP();
}

/// @brief Set the base type COP.  This implies that this ResidueType is NOT a base type.
///
void
ResidueType::set_base_type_cop(
	ResidueTypeCOP new_base_type
) {
	runtime_assert_string_msg( new_base_type, "Error in core::chemical::ResidueType::set_base_type_cop(): A null pointer was passed to this function." );
	base_type_cop_ = new_base_type;
}

bool
ResidueType::is_polymer() const
{
	return properties_->has_property( POLYMER );
}

bool
ResidueType::is_sidechain_thiol() const
{
	return properties_->has_property( SIDECHAIN_THIOL );
}

bool
ResidueType::is_disulfide_bonded() const
{
	return properties_->has_property( DISULFIDE_BONDED );
}

bool
ResidueType::is_sidechain_amine() const
{
	return properties_->has_property( SIDECHAIN_AMINE );
}

bool
ResidueType::is_protein() const
{
	return properties_->has_property( PROTEIN );
}

/// @brief Is this an alpha-amino acid?
///
bool
ResidueType::is_alpha_aa() const {
	return properties_->has_property( ALPHA_AA );
}

/// @brief Is this a beta-amino acid?
///
bool
ResidueType::is_beta_aa() const {
	return properties_->has_property( BETA_AA );
}

/// @brief Is this a gamma-amino acid?
///
bool
ResidueType::is_gamma_aa() const {
	return properties_->has_property( GAMMA_AA );
}


/// @brief Does this type have groups (not just single atoms) that are polymer-bond dependent?
///
bool
ResidueType::has_polymer_dependent_groups() const {
	return is_n_methylated(); //TODO: Update this if other polymer-dependent types are added.
}

/// @brief Is this one of SRI's special heteropolymer building blocks?
///
bool
ResidueType::is_sri() const {
	return properties_->has_property( SRI );
}

/// @brief Is this a triazolemer?
///
bool
ResidueType::is_triazolemer() const
{
	return properties_->has_property( TRIAZOLE_LINKER );
}

bool
ResidueType::is_d_aa() const
{
	return properties_->has_property( D_AA );
}

bool
ResidueType::is_l_aa() const
{
	return properties_->has_property( L_AA );
}

bool
ResidueType::is_d_rna() const
{
	return properties_->has_property( D_RNA );
}

bool
ResidueType::is_l_rna() const
{
	return properties_->has_property( L_RNA );
}

/// @brief Is this residue N-methylated?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::is_n_methylated() const {
	return properties_->has_property( N_METHYLATED );
}

bool
ResidueType::is_achiral_backbone() const
{
	return properties_->has_property( ACHIRAL_BACKBONE );
}

bool
ResidueType::is_DNA() const
{
	return properties_->has_property( DNA );
}

bool
ResidueType::is_RNA() const
{
	return properties_->has_property( RNA );
}

bool
ResidueType::is_coarse() const
{
	return properties_->has_property( COARSE );
}

bool
ResidueType::is_NA() const
{
	return is_DNA() || is_RNA();
}

bool
ResidueType::is_purine() const
{
	return properties_->has_property( PURINE );
}

bool
ResidueType::is_pyrimidine() const
{
	return properties_->has_property( PYRIMIDINE );
}

/// @brief Is this a solvent molecule (SOLVENT property)?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::is_solvent() const
{
	return properties_->has_property( SOLVENT );
}

/// @brief Is this a canonical nucleic acid (CANONICAL_NUCLEIC property)?
/// @details Only the standard nucliec acid types (dA, dC, dG, dT, A, C, G, U) are canonical.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::is_canonical_nucleic() const
{
	return properties_->has_property( CANONICAL_NUCLEIC );
}

/// @brief Is this a canonical amino acid (CANONICAL_AA property)?
/// @details Only the standard amino acid types (ACDEFGHIKLMNPQRSTVWY) are canonical.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
ResidueType::is_canonical_aa() const
{
	return properties_->has_property( CANONICAL_AA );
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
ResidueType::is_peptoid() const
{
	return properties_->has_property( PEPTOID );
}

bool
ResidueType::is_carbohydrate() const
{
	return properties_->has_property( CARBOHYDRATE );
}

bool
ResidueType::is_ligand() const
{
	return properties_->has_property( LIGAND );
}

bool
ResidueType::is_lipid() const
{
	return properties_->has_property( LIPID );
}

/// @details The METAL property is specified in the params file under PROPERTIES.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
ResidueType::is_metal() const
{
	return properties_->has_property( METAL );
}

/// @details The METALBINDING property is specified in the params file under PROPERTIES.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
ResidueType::is_metalbinding() const
{
	return properties_->has_property( METALBINDING );
}

bool
ResidueType::is_membrane() const
{
	return properties_->has_property( MEMBRANE );
}

bool
ResidueType::is_surface() const
{
	return properties_->has_property( SURFACE );
}

bool
ResidueType::has_sc_orbitals() const
{
	return properties_->has_property( SC_ORBITALS );
}

bool
ResidueType::is_polar() const
{
	return properties_->has_property( POLAR );
}

bool
ResidueType::is_charged() const
{
	return properties_->has_property( CHARGED );
}

bool
ResidueType::is_aromatic() const
{
	return properties_->has_property( AROMATIC );
}

bool
ResidueType::is_cyclic() const
{
	return properties_->has_property( CYCLIC );
}

bool
ResidueType::is_terminus() const
{
	return is_upper_terminus() || is_lower_terminus();
}

bool
ResidueType::is_lower_terminus() const
{
	return properties_->has_property( LOWER_TERMINUS );
}

bool
ResidueType::is_upper_terminus() const
{
	return properties_->has_property( UPPER_TERMINUS );
}

bool
ResidueType::is_branch_point() const
{
	return properties_->has_property( BRANCH_POINT );
}

bool
ResidueType::is_branch_lower_terminus() const
{
	return properties_->has_property( BRANCH_LOWER_TERMINUS );
}

bool
ResidueType::is_acetylated_nterminus() const
{
	return properties_->has_property( ACETYLATED_NTERMINUS );
}

bool
ResidueType::is_methylated_cterminus() const
{
	return properties_->has_property( METHYLATED_CTERMINUS );
}

bool
ResidueType::is_virtual_residue() const
{
	return properties_->has_property( VIRTUAL_RESIDUE );
}

/// @brief  Check if residue is 'INVERTING_VIRTUAL_RESIDUE'
/// @details Used by the symmetry machinery for mirror symmetry operations.
bool
ResidueType::is_inverted_virtual_residue() const
{
	return properties_->has_property( INVERTED_VIRTUAL_RESIDUE );
}


bool
ResidueType::is_adduct() const
{
	return properties_->has_property( ADDUCT );
}

/// @brief  Generic property access.
bool
ResidueType::has_property( std::string const & property ) const
{
	return properties_->has_property( property );
}

/// @brief  Generic property access, by ResidueProperty.
///
bool
ResidueType::has_property( ResidueProperty const property ) const
{
	return properties_->has_property( property );
}

core::Real
ResidueType::get_numeric_property(std::string const & tag) const
{
	std::map<std::string, core::Real> const numeric_properties( properties_->numeric_properties() );
	auto property_it( numeric_properties.find( tag ) );
	if ( property_it == numeric_properties.end() ) {
		throw utility::excn::EXCN_KeyError( tag + " does not exist in ResidueType with name " + name3_ );
		return 0.0; //keep compilers happy
	}

	return property_it->second;
}

std::string
ResidueType::get_string_property(std::string const & tag) const
{
	std::map<std::string, std::string> const string_properties( properties_->string_properties() );
	auto property_it(string_properties.find(tag));
	if ( property_it == string_properties.end() ) {
		throw utility::excn::EXCN_KeyError(tag + " does not exist in ResidueType with name " + name3_);
		return "";
	}
	return property_it->second;
}


void
ResidueType::add_variant_type( VariantType const variant_type )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	properties_->set_variant_type( variant_type, true );
}

void
ResidueType::add_variant_type( std::string const & variant_type )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	properties_->set_variant_type( variant_type, true );
}

void
ResidueType::remove_variant_type( VariantType const variant_type )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	properties_->set_variant_type( variant_type, false );
}

void
ResidueType::remove_variant_type( std::string const & variant_type )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	properties_->set_variant_type( variant_type, false );
}

/// @brief  Generic variant access.
bool
ResidueType::has_variant_type( VariantType const variant_type ) const
{
	return properties_->is_variant_type( variant_type );
}

// TODO: Find a way to remove this; it only exists because of how ResidueTypeSelectors are currently written. ~Labonte
/// @brief  Generic variant access by string.
bool
ResidueType::has_variant_type( std::string const & variant_type ) const
{
	return properties_->is_variant_type( variant_type );
}

/// @details "Custom" VariantTypes as strings are permitted for the enzdes and metalloproteins cases.
/// Do not enable unless you have a good reason to, as string look-ups are less efficient and more error-prone.
void
ResidueType::enable_custom_variant_types()
{
	properties_->enable_custom_variant_types();
}

/// @brief Get a vector of VariantType enums for this ResidueType.
/// @details This ONLY includes standard, enum-based variants, not on-the-fly custom variants.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
utility::vector1< VariantType >
ResidueType::variant_type_enums() const {
	return properties_->get_list_of_variant_enums();
}

/// @brief Get a list of custom VariantType strings for this ResidueType (by const reference).
/// @details This ONLY includes custom, on-the-fly variants, not standard variants.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
utility::vector1< std::string > const &
ResidueType::custom_variant_types() const{
	return properties_->get_list_of_custom_variants_by_reference();
}

///////////////////////////////////////////////////////////////////////////////

// delete the terminal chi angle
/// @author Andrew M. Watkins (April 2015)
void
ResidueType::delete_terminal_chi(
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	// if the terminal chi was a proton chi, get rid of it!
	if ( is_proton_chi_[ is_proton_chi_.size() ] ) {
		proton_chis_.resize( proton_chis_.size() - 1 );
	}

	Size new_size = chi_atoms_.size() - 1;
	// resize vectors that include every chi
	chi_atoms_.resize(        new_size );
	chi_rotamers_.resize(     new_size );
	is_proton_chi_.resize(    new_size );
	chi_2_proton_chi_.resize( new_size );


} // delete_terminal_chi

void
ResidueType::delete_child_proton( std::string const & atom ) {
	std::string res_varname( atom + "-PRUNEH" );
	Size count = 0;
	while ( true ) {
		if ( count > 20 ) {
			utility_exit_with_message( "Could not find a new VariantType for ResidueType: " + name() );
		}
		++count;
		if ( count == 1 ) {
			if ( ! has_variant_type( res_varname ) ) break;
		} else {
			res_varname = atom + "-PRUNEH" + utility::to_string( count );
			if ( ! has_variant_type( res_varname ) ) break;
		}
	}
	enable_custom_variant_types();
	add_variant_type( res_varname );

	// AMW: It seems like when we "delete" a proton, or fail to do so and virt
	// instead, it doesn't keep track of it...
	core::Size nhydrogens = number_bonded_hydrogens( atom_index( atom ) );
	if ( nhydrogens == 0 ) {
		tr.Trace << "No bonded hydrogens at " << atom << " in " << name() << std::endl;
	} else {
		// delete last proton
		Size proton_index = attached_H_end( atom_index( atom ) );

		// 1: delete
		tr.Trace << "Removing " << atom_name( proton_index ) << std::endl;
		delete_atom( proton_index );

		// 2: remove any chi containing this H
		for ( Size ii = 1; ii <= nchi(); ++ii ) {
			if ( chi_atoms( ii )[ 4 ] != proton_index )  continue;

			if ( ii == nchi() ) {
				delete_terminal_chi();
			} else {
				// Redefine every chi from ii to nchi - 1 to jj + 1
				// Then delete the terminal one.
				for ( Size jj = ii; jj <= nchi() - 1; ++jj ) {
					redefine_chi( jj,
						atom_name( chi_atoms( jj + 1 )[ 1 ] ),
						atom_name( chi_atoms( jj + 1 )[ 2 ] ),
						atom_name( chi_atoms( jj + 1 )[ 3 ] ),
						atom_name( chi_atoms( jj + 1 )[ 4 ] ) );
				}
				delete_terminal_chi();
			}
		}

		// 3: ensure that the deleted proton is not used to build another atom in the residue
		//    die for now.  If this is a problem this logic could be made smarter.
		for ( Size ii = 1; ii <= natoms(); ++ii ) {
			AtomICoor aicoor = icoor( ii );
			if ( aicoor.stub_atom1().atomno() == proton_index
					|| aicoor.stub_atom2().atomno() == proton_index
					|| aicoor.stub_atom3().atomno() == proton_index
					) {
				utility_exit_with_message( "Deleted proton " + atom_name( proton_index ) + " used to build neighbor atom!" );
			}
		}

		// 4: if there is more than one proton, allow the remain proton to occupy other positions
		if ( nhydrogens > 1 ) {
			Size alt_proton_index = attached_H_begin( atom_index( atom ) );
			AtomICoor aicoor = icoor( alt_proton_index );

			for ( Size ii = 1; ii <= nchi(); ++ii ) {
				if ( chi_atoms( ii )[ 4 ] != proton_index )  continue;
				utility::vector1< Real > dihedral_samples;
				for ( Size jj = 0; jj<nhydrogens; ++jj ) {
					dihedral_samples.push_back( fmod( aicoor.phi() + jj*(360.0/nhydrogens), 360.0) );
				}
				set_proton_chi( ii, dihedral_samples, utility::vector1< Real >() );
			}
		}
	}

	//fd  we need to update the attachedH mappings so call finalize rather than update_derived
	finalize();
}

void
ResidueType::add_metapatch_connect( std::string const & atom ) {
	// Provide unique variant name
	// We have to do this or connections get dropped--not all variants get put
	// back in. This is worse than you think--because they DON'T get dropped by
	// the metal!
	using namespace numeric::conversions;
	std::string res_varname( atom + "-CONNECT" );
	Size count=0;
	while ( true ) {
		if ( count > 20 ) {
			utility_exit_with_message( "Could not find a new VariantType for ResidueType: " + name() );
		}
		++count;
		if ( count == 1 ) {
			if ( ! has_variant_type( res_varname ) ) break;
		} else {
			res_varname = atom + "-CONNECT" + utility::to_string( count );
			if ( ! has_variant_type( res_varname ) ) break;
		}
	}
	enable_custom_variant_types();
	add_variant_type( res_varname );

	if ( number_bonded_hydrogens( atom_index( atom ) ) == 0 ) {
		Size const connid( add_residue_connection( atom ) );
		AtomICoor aicoor = icoor( atom_index( atom ) );

		// These coordinates are generic.
		set_icoor( "CONN"+ObjexxFCL::string_of( connid ), 3.14159, 70.600000*3.14159/180.000000, 1.37, atom, atom_name( aicoor.stub_atom1().atomno() ), atom_name( aicoor.stub_atom2().atomno() ) );
	} else {
		Size proton_index = attached_H_begin( atom_index( atom ) );
		AtomICoor aicoor = icoor( proton_index );

		Size const connid( add_residue_connection( atom ) );
		set_icoor( "CONN"+ObjexxFCL::string_of( connid ),
			aicoor.phi()+radians(180.0),
			aicoor.theta(),
			1.37,
			atom_name( aicoor.stub_atom1().atomno() ),
			atom_name( aicoor.stub_atom2().atomno() ),
			atom_name( aicoor.stub_atom3().atomno() ) );
	}

	update_derived_data();
}

///////////////////////////////////////////////////////////////////////////////

// redefine a chi angle based on four atoms
/// @details This function is almost an exact copy of the add_chi function except that vector resizing does NOT occur.
/// It is needed for certain PTMs that affect proton chis (e.g., phosphorylation and sulfation).
/// @author Andy M. Chen (June 2009)
void
ResidueType::redefine_chi(
	Size const chino,
	std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ||
			!has( atom_name3 ) || !has( atom_name4 ) ) {
		utility_exit_with_message("ResidueType::redefine_chi: atoms dont exist!" );
	}

	utility::vector1<VD> atoms;
	atoms.push_back( ordered_atoms_[ atom_index( atom_name1 ) ]);
	atoms.push_back( ordered_atoms_[ atom_index( atom_name2 ) ]);
	atoms.push_back( ordered_atoms_[ atom_index( atom_name3 ) ]);
	atoms.push_back( ordered_atoms_[ atom_index( atom_name4 ) ]);
	chi_atoms_[chino] = atoms;

	// Assumes that the redefined chi is NOT a proton chi.
	//  (This is adequate in most cases because PTMs tend to replace hydrogens
	//  with functional groups rather than the other way around.)
	is_proton_chi_[ chino ] = false;
	chi_2_proton_chi_[ chino ] = 0;
} // redefine_chi


/////////////////////////////////////////////////////////////////

/// @details add an atom to the list for calculating actcoord center
void
ResidueType::add_actcoord_atom( std::string const & atom )
{
	if ( ! is_protein() ) {
		tr.Warning << "ACT_COORD_ATOM specified for non-protein residue type '" << name() << "' . This doesn't make much sense." << std::endl;
	}
	finalized_ = false;
	tr.Trace << "adding act coord atom: " << name_ << ' ' << atom << std::endl;
	Size atomindex = atom_index( atom );
	actcoord_atoms_.push_back( ordered_atoms_[atomindex] );
}

///////////////////////////////////////////////////////////////////////////////

// set up atom ordering, called by finalize()
/**
@details Because some new heavy atoms are added and removed by patching
we need to reorder the atoms. The atoms are reordered with the backbone
atoms coming first, sidechain atoms second, and finaly hydrogens last.

The graph is iterated over, where bb atoms are defined by the force_bb_
map. If a hydrogen is encountered, it is pushed into the hydrogen vector.
If the backbone is not found in the force_bb_ map, then the atom is a side
chain atom, if not a hydrogen.

During this process, the attached_H_end and attached_H_begin vectors are
filled out.

In order to keep the bb hydrogen atoms first in the hydrogen vector, only the
heavy sidechain atoms are assigned to the sidechain vector while bb hydorgens
atoms are pushed into the hydrogen vectors.

After the bb and sidechain vectors are assigned, the sidechains are iterated
over and their hydrogens are pushed into the hydrogen vector.

Finally, the ordered_atoms_ vector is filled out with the bb atoms, sidechain,
and hydrogen atoms.

Derived data set in this method:
* ordered_atoms_ - (Reset from graph.)
* n_backbone_heavyatoms_
* nheavyatoms_
* attached_H_begin_
* attached_H_end_
**/
void
ResidueType::setup_atom_ordering()
{
	utility::vector1<VD> bb_atoms, sidechain_atoms, hydrogens;
	utility::vector1<Size> h_begin;
	utility::vector1<Size> h_end;
	//update the bonded neighbors here
	for ( VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first ) {
		VD const & vd = *vp.first;
		//we need to iterate through the edges to update
		if ( !graph_[vd].is_hydrogen() ) {
			if ( std::find( force_bb_.begin(), force_bb_.end(), vd ) != force_bb_.end() ) {
				bb_atoms.push_back(vd);
				h_begin.push_back(hydrogens.size() +1 );
				bool h_present(false);
				for ( OutEdgeIterPair ep = boost::out_edges(vd, graph_); ep.first != ep.second; ++ep.first ) { //iterate through the edges for hydrogens
					VD target = boost::target(*ep.first, graph_);
					if ( graph_[target].is_hydrogen() ) {
						hydrogens.push_back(target);
						h_present = true;
					}
				}
				if ( h_present ) { h_end.push_back(hydrogens.size()); } else { h_end.push_back(0); }
			} else {
				sidechain_atoms.push_back(vd);
			}
		}
	}
	//std::cout << "bb_atoms.size(): " << bb_atoms.size() << std::endl;
	//once we set up the sidechains, we have to push back the hydrogens for the sidechains
	for ( Size i= 1; i<= sidechain_atoms.size(); ++i ) {
		h_begin.push_back(hydrogens.size() +1 );
		bool h_present(false);
		for ( OutEdgeIterPair ep = boost::out_edges(sidechain_atoms[i], graph_); ep.first != ep.second; ++ep.first ) { //iterate through the edges for hydrogens
			VD target = boost::target(*ep.first, graph_);
			if ( graph_[target].is_hydrogen() ) {
				hydrogens.push_back(target);
				h_present = true;
			}
		}
		if ( h_present ) {
			h_end.push_back(hydrogens.size());
		} else {
			h_end.push_back(0);
		}
	}

	//create the ordered_atoms_ data by putting bb_atoms, then sidechain atoms, and finally
	//hydrogen atoms into the vector
	ordered_atoms_.clear();
	ordered_atoms_.reserve(bb_atoms.size()+sidechain_atoms.size()+hydrogens.size());
	ordered_atoms_.insert(ordered_atoms_.end(), bb_atoms.begin(), bb_atoms.end());
	ordered_atoms_.insert(ordered_atoms_.end(), sidechain_atoms.begin(), sidechain_atoms.end());
	ordered_atoms_.insert(ordered_atoms_.end(), hydrogens.begin(), hydrogens.end());

	debug_assert( boost::num_vertices( graph_ ) == ordered_atoms_.size() );

	n_backbone_heavyatoms_ = bb_atoms.size();
	nheavyatoms_ = bb_atoms.size() + sidechain_atoms.size();


	//This step is done at the end after the nheavy_atoms_ has been assigned.
	attached_H_begin_.clear();
	attached_H_end_.clear();
	for ( Size i =1; i<= h_begin.size(); ++i ) {
		attached_H_begin_.push_back(nheavyatoms_+h_begin[i]);
		if ( h_end[i] ==0 ) {
			attached_H_end_.push_back(0);
		} else { attached_H_end_.push_back(nheavyatoms_+h_end[i]); }
	}

}

///////////////////////////////////////////////////////////////////////////////

/**
The private data for ResidueType is based on vertex descriptors (VD).
Rosetta relies on atom indices to access the private data. In order to
keep this interface, the VDs must be converted into atom indices. These
are called *_indices, where the primary data is a VD but the indices_data
is the atom indices. Because the atom indices change during the setup
atom ordering, we must regenerate the cached data. This would not have
to occur if Rosetta relied on VDs instead of atom indices.

First, the private data that relies on atom ordering, atom_name_to_vd_ and
vd_to_index_ is created by iterating over the total number of atoms in the
graph.

Second, the bonded neighbors are generated based on the graph structure.

Finally, all the cached data is generated by iterating over the VDs.

-- Assumes that ordered_atoms_ is updated appropriately.
Derived data set in this method:
* atom_name_to_vd_
* vd_to_index_
* atom_base_indices_
* chi_atoms_indices_
* nu_atoms_indices_
* ring_atoms_indices_
* mainchain_atoms_indices_
* nbr_atom_indices_
* actcoord_atoms_indices_
* cut_bond_neighbor_indices_
* atom_shadowed_indices_
* bonded_neighbor_
* bonded_neighbor_type_
* atom_2_residue_connection_map_

**/
void
ResidueType::generate_atom_indices()
{
	atom_2_residue_connection_map_.resize( natoms() );

	atom_name_to_vd_.clear();
	vd_to_index_.clear();
	for ( Size i=1; i<= natoms(); ++i ) {
		atom_name_to_vd_[ graph_[ ordered_atoms_[i]].name() ] = ordered_atoms_[i];
		atom_name_to_vd_[ strip_whitespace( graph_[ ordered_atoms_[i]].name()  ) ] = ordered_atoms_[i];
		vd_to_index_[ordered_atoms_[i]] = i;
	}
	// Add in aliased atoms
	for ( auto iter( atom_aliases_.begin() ), iter_end( atom_aliases_.end() );
			iter != iter_end; ++iter ) {
		debug_assert( !atom_name_to_vd_.count( iter->first ) );
		debug_assert( atom_name_to_vd_.count( iter->second ) );
		VD vd( atom_name_to_vd_[ iter->second ] );
		atom_name_to_vd_[iter->first] = vd;
	}

	bonded_neighbor_.clear();
	bonded_neighbor_type_.clear();
	bonded_neighbor_.resize(natoms());
	bonded_neighbor_type_.resize(natoms());
	//setup bond ordering!
	for ( VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first ) {
		VD const & vd = *vp.first;
		bonded_neighbor_[vd_to_index_[vd]].clear(); //we have to clear the vectors first before pushing back to them
		bonded_neighbor_type_[vd_to_index_[vd]].clear(); //we have to clear the vectors first before pushing back to them
		for ( OutEdgeIterPair ep = boost::out_edges(vd, graph_); ep.first != ep.second; ++ep.first ) { //iterate through the edges for hydrogens
			VD target = boost::target(*ep.first, graph_);
			bonded_neighbor_[vd_to_index_[vd]].push_back(vd_to_index_[target]);
			ED const & edge = *ep.first;
			BondName bond = graph_[edge].bond_name();
			bonded_neighbor_type_[vd_to_index_[vd]].push_back(bond);
		}
	}


	atom_base_indices_.clear();
	atom_shadowed_indices_.clear();
	for ( Size atomno=1; atomno<= natoms(); ++atomno ) {
		{
			if ( atom_base_.find(ordered_atoms_[atomno]) == atom_base_.end() ) {
				atom_base_indices_.push_back(0);
			} else {
				atom_base_indices_.push_back(vd_to_index_.find((atom_base_.find(ordered_atoms_[atomno])->second))->second);
			}
		}

		{
			if ( atom_shadowed_.find(ordered_atoms_[atomno]) == atom_shadowed_.end() ) {
				atom_shadowed_indices_.push_back(0);
			} else {
				atom_shadowed_indices_.push_back(vd_to_index_.find((atom_shadowed_.find(ordered_atoms_[atomno])->second))->second);
			}
		}
	}

	cut_bond_neighbor_indices_.clear();
	cut_bond_neighbor_indices_.resize(natoms());
	AtomIndices atoms;
	for ( Size atomno=1; atomno<= natoms(); ++atomno ) {
		if ( cut_bond_neighbor_.find(ordered_atoms_[atomno]) != cut_bond_neighbor_.end() ) {
			utility::vector1<VD> const & vd_cut_atoms = cut_bond_neighbor_.find( ordered_atoms_[atomno] )->second;
			for ( Size i = 1; i <= vd_cut_atoms.size(); ++i ) { //if you get here, you will have assigned chi atoms
				atoms.push_back(vd_to_index_.find( vd_cut_atoms[i] )->second);
			}
			cut_bond_neighbor_indices_[atomno] =  atoms;
			atoms.clear();
		}
	}

	chi_atoms_indices_.clear();
	atoms.clear();
	utility::vector1<AtomIndices> chi_atoms;
	for ( Size chino=1; chino <= chi_atoms_.size(); ++chino ) {
		for ( Size atom_index=1; atom_index <= chi_atoms_[chino].size(); ++atom_index ) {
			atoms.push_back( vd_to_index_.find( chi_atoms_[chino][atom_index] )->second);
		}
		chi_atoms_indices_.push_back(atoms);
		atoms.clear();
	}
	nu_atoms_indices_.clear();
	atoms.clear();
	for ( Size nu_no = 1; nu_no <= nu_atoms_.size(); ++nu_no ) {
		for ( Size atom_index=1; atom_index <= nu_atoms_[ nu_no ].size(); ++atom_index ) {
			atoms.push_back( vd_to_index_.find( nu_atoms_[ nu_no ][ atom_index ] )->second);
		}
		nu_atoms_indices_.push_back( atoms );
		atoms.clear();
	}
	ring_atoms_indices_.clear();
	atoms.clear();
	for ( uint ring_no( 1 ); ring_no <= ring_atoms_.size(); ++ring_no ) {
		for ( uint atom_index( 1 ); atom_index <= ring_atoms_[ ring_no ].size(); ++atom_index ) {
			atoms.push_back( vd_to_index_.find( ring_atoms_[ ring_no ][ atom_index ] )->second);
		}
		ring_atoms_indices_.push_back( atoms );
		atoms.clear();
	}


	actcoord_atoms_indices_.clear();
	for ( Size i=1; i<= actcoord_atoms_.size(); ++i ) {
		actcoord_atoms_indices_.push_back(vd_to_index_.find(actcoord_atoms_[i])->second);
	}

	mainchain_atoms_indices_.clear();
	for ( Size i=1; i<= mainchain_atoms_.size(); ++i ) {
		mainchain_atoms_indices_.push_back(vd_to_index_.find(mainchain_atoms_[i])->second);
	}

	if ( nbr_atom_ == ResidueGraph::null_vertex() ) {
		nbr_atom_indices_ = 0;
	} else {
		std::map<VD, Size>::const_iterator nbr_translation( vd_to_index_.find(nbr_atom_) );
		debug_assert( nbr_translation != vd_to_index_.end() );
		nbr_atom_indices_ = nbr_translation->second;
	}

	for ( Size index=1; index<= natoms(); ++index ) {
		utility::vector1< core::Size > const orbs(graph_[ordered_atoms_[index]].bonded_orbitals());
		for ( core::Size orb : orbs ) {
			orbitals_[orb].new_icoor().replace_stub1( vd_to_index_[ordered_atoms_[orbitals_[orb].new_icoor().get_stub1()]] );
			orbitals_[orb].new_icoor().replace_stub2( vd_to_index_[ordered_atoms_[orbitals_[orb].new_icoor().get_stub2()]] );
			orbitals_[orb].new_icoor().replace_stub3( vd_to_index_[ordered_atoms_[orbitals_[orb].new_icoor().get_stub3()]] );
		}
	}


	for ( Size index=1; index<= natoms(); ++index ) {
		for ( Size i=1; i<= 3; ++i ) {
			ICoorAtomID & stub_atom( icoor_[ ordered_atoms_[index] ].stub_atom( i )   );
			if ( stub_atom.type() == ICoorAtomID::INTERNAL ) {
				stub_atom.atomno(   vd_to_index_.find(stub_atom.vertex())->second ); //somewhat of a problem. if vertex doesnt exist the map constructor will create a value
				if ( stub_atom.atomno() == 0 ) { // this will trigger if we deleted a stub atom for some other atom
					tr.Error << "===============================================" << std::endl;
					tr.Error << "Internal coordinates for " << name() << std::endl;
					pretty_print_atomicoor(tr.Error, *this);
					tr.Error << "===============================================" << std::endl;
					utility_exit_with_message("ERROR: The internal coordinates of ResidueType '" + name() + "' depend on a missing atom!");
				}
			}
		}
	}

	///TODO fix this problem: deleting atoms can invalidate these residue_connections_
	for ( Size i=1; i<= residue_connections_.size(); ++i ) {
		residue_connections_[i].atomno( vd_to_index_[residue_connections_[i].vertex()]  );
		AtomICoor new_icoor = residue_connections_[i].icoor();
		for ( Size j = 1; j <= 3; ++j ) {
			new_icoor.stub_atom( j ).atomno( vd_to_index_.find(new_icoor.stub_atom(j).vertex())->second );
		}
		residue_connections_[ i ].icoor( new_icoor );
		debug_assert( residue_connections_[i].atomno() ); //this will fail if we deleted an atom involved in an inter-rsd connection
	}
	update_residue_connection_mapping();

}


///////////////////////////////////////////////////////////////////////////////

/// @details update derived data in ResidueType, called by finalize()
/**
after primary data have been reordered, update derived data accordingly.

Derived data updated by this method:
* first_sidechain_hydrogen_
* accpt_pos_
* accpt_pos_sc_
* Haro_index_
* Hpol_index_
* atoms_with_orb_index_
* Hpos_polar_
* Hpos_apolar_
* Hpos_polar_sc_
* all_bb_atoms_
* all_sc_atoms_
* abase2_
* path_distance_
* dihedral_atom_sets_
* dihedrals_for_atom_
* bondangle_atom_sets_
* bondangles_for_atom_
* atoms_within_one_bond_of_a_residue_connection_
* within1bonds_sets_for_atom_
* atoms_within_two_bonds_of_a_residue_connection_
* within2bonds_sets_for_atom_
* last_controlling_chi_
* atoms_last_controlled_by_chi_

* Atom.heavyatom_has_polar_hydrogens_ - for all atoms.
* Some AtomProperties:
** AROMATIC_CARBON_WITH_FREE_VALENCE

* rna_info_  -- Will be reset based on other ResidueType data
* carbohydrate_info_ -- Will be reset based on other ResidueType data
**/
void
ResidueType::update_derived_data()
{
	first_sidechain_hydrogen_ = natoms() + 1;
	for ( Size i= n_backbone_heavyatoms_ + 1; i<= nheavyatoms_; ++i ) {
		if ( attached_H_begin_[i] <= attached_H_end_[i] ) {
			first_sidechain_hydrogen_ = attached_H_begin_[i];
			break;
		}
	}

	// compile atom-index lists of subsets of the atoms
	accpt_pos_.clear();
	accpt_pos_sc_.clear();
	Haro_index_.clear();
	Hpol_index_.clear();
	atoms_with_orb_index_.clear();
	Hpos_polar_.clear();
	Hpos_apolar_.clear();
	Hpos_polar_sc_.clear();
	all_bb_atoms_.clear();
	all_sc_atoms_.clear();

	for ( Size i=1; i<= natoms(); ++i ) {
		Atom const & atom(graph_[ ordered_atoms_[i]]); //get the atom that we are working on

		// info derived from the atom
		if ( atom.has_orbitals() ) atoms_with_orb_index_.push_back(i); //get atoms with orbitals on it
		if ( atom.is_haro() ) Haro_index_.push_back( i ); //get aromatic hydrogens
		if ( atom.is_polar_hydrogen() ) Hpol_index_.push_back( i ); //get polar hydrogens
		if ( atom.is_acceptor() && !atom.is_virtual() ) {
			accpt_pos_.push_back( i );
			if ( i > n_backbone_heavyatoms_ ) {
				accpt_pos_sc_.push_back( i );
			}
		}
		if ( atom.is_polar_hydrogen() && !atom.is_virtual() ) {
			Hpos_polar_.push_back( i );
			if ( i >= first_sidechain_hydrogen_ ) {
				Hpos_polar_sc_.push_back( i );
			}
		}
		if ( atom.is_hydrogen() && !atom.is_polar_hydrogen() ) {
			Hpos_apolar_.push_back( i );
		}

		// Which atoms are backbone and which are sidechain; sometimes nice to just get
		// lists instead of iterating over the subranges.
		if ( atom.is_hydrogen() ) {
			if ( i < first_sidechain_hydrogen_ ) {
				all_bb_atoms_.push_back( i );
			} else {
				all_sc_atoms_.push_back( i );
			}
		} else {
			if ( i <= n_backbone_heavyatoms_ ) {
				all_bb_atoms_.push_back( i );
			} else {
				all_sc_atoms_.push_back( i );
			}
		}
	}

	// setup the hydrogen information
	for ( Size Aindex=1; Aindex<= ordered_atoms_.size(); ++Aindex ) {
		graph_[ordered_atoms_[Aindex]].heavyatom_has_polar_hydrogens(0);
	}

	// donor heavy atoms, acceptor heavy atoms, donor hydrogen atoms setup.
	// Must be executed after Hpos_polar_ and accpt_pos_ have been updated.
	for ( Size ii = 1; ii <= Hpos_polar_.size(); ++ii ) {
		Size hind = Hpos_polar_[ ii ];
		Size base = atom_base(hind);
		graph_[ordered_atoms_[base]].heavyatom_has_polar_hydrogens(1);
	}


	// now setup abase2
	abase2_.clear();
	for ( Size ii=1, ii_end= accpt_pos_.size(); ii<= ii_end; ++ii ) {
		uint const acceptor_position( accpt_pos_[ii] ); //acceptor_position
		uint const acc_base( atom_base(acceptor_position) ); //acceptor_base
		debug_assert(acc_base == atom_base(acceptor_position) );
		debug_assert( acc_base != 0 );
		AtomIndices const & i_nbrs(bonded_neighbor(acceptor_position));
		if ( i_nbrs.size() == 0 ) {
			utility_exit_with_message( "failed to set abase2 for acceptor atom, it has no nbrs!" );
		} else if ( i_nbrs.size() == 1 ) {
			//debug_assert( i_nbrs[1] == acc_base );
			abase2_[ordered_atoms_[acceptor_position]] = ordered_atoms_[atom_base(acc_base)];
			//iwd  The first child of the root is root's atom_base.
			//iwd  But if that child has no children, it ends up as its own abase2.
			//iwd  So instead we use the second child of the parent,
			//iwd  which must exist if there are 3+ atoms in this tree.
			std::map<VD, VD>::const_iterator  acc_abase2( abase2_.find(ordered_atoms_[acceptor_position]) );
			debug_assert( acc_abase2 != abase2_.end() );
			std::map<VD, Size>::const_iterator acc_abase2_index( vd_to_index_.find( acc_abase2->second ) );
			debug_assert( acc_abase2_index != vd_to_index_.end() );
			if ( acc_abase2_index->second == acceptor_position ) {
				AtomIndices const & i_base_nbrs(bonded_neighbor(acc_base) );
				for ( Size jj = 1, jj_end = i_base_nbrs.size(); jj <= jj_end; ++jj ) {
					if ( i_base_nbrs[ jj ] != acceptor_position ) {
						abase2_[ordered_atoms_[acceptor_position] ]  = ordered_atoms_[i_base_nbrs[ jj ]];
						break;
					}
				}
			}
			//debug_assert(abase2(acceptor_position)!=acceptor_position && abase2(acceptor_position) != acc_base && abase2(acceptor_position) != 0 );
		} else if ( i_nbrs[1] == acc_base ) {
			abase2_[ordered_atoms_[acceptor_position]] = ordered_atoms_[i_nbrs[2]];
		} else {
			abase2_[ordered_atoms_[acceptor_position]] = ordered_atoms_[i_nbrs[1] ];
		}
	}

	abase2_indices_.clear();
	for ( Size atomno=1; atomno<= natoms(); ++atomno ) {
		if ( abase2_.find(ordered_atoms_[atomno]) == abase2_.end() ) {
			abase2_indices_.push_back(0);
		} else {
			abase2_indices_.push_back(vd_to_index_.find(  abase2_.find(ordered_atoms_[atomno])->second )->second);
		}
	}


	// bond path distances
	FArray2D_int path_distances( get_residue_path_distances( *this ));
	path_distance_.resize( natoms() );
	for ( Size ii = 1; ii <= natoms(); ++ii ) {
		path_distance_[ ii ].resize( natoms() );
		for ( Size jj = 1; jj <= natoms(); ++jj ) {
			path_distance_[ ii ][ jj ] = path_distances( ii, jj );
		}
	}

	// get dihedral angles
	dihedral_atom_sets_.clear();
	dihedrals_for_atom_.resize( natoms() );
	for ( Size ii = 1; ii <= natoms(); ++ii ) dihedrals_for_atom_[ ii ].clear();

	// get for all pairs of atoms separated by 1 bond
	for ( Size central_atom1 = 1; central_atom1 < natoms(); ++central_atom1 ) {
		for ( Size central_atom2 = central_atom1+1; central_atom2 <= natoms(); ++central_atom2 ) {
			if ( path_distance_[ central_atom1 ][ central_atom2 ] == 1 ) {

				// get all atoms separated from central_atom1/2 by one bond that are not central_atom2/1
				utility::vector1< Size > ca1d1;
				utility::vector1< Size > ca2d1;

				// ca1
				for ( Size i = 1; i <= natoms(); ++i ) {
					if ( ( path_distance_[ central_atom1 ][ i ] == 1 ) && ( i != central_atom2 ) ) {
						ca1d1.push_back( i );
					}
				}
				// ca2
				for ( Size i = 1; i <= natoms(); ++i ) {
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

	// get dihedral angles
	improper_dihedral_atom_sets_.clear();
	improper_dihedrals_for_atom_.resize( natoms() );
	for ( Size ii = 1; ii <= natoms(); ++ii ) improper_dihedrals_for_atom_[ ii ].clear();

	// get for all pairs of atoms separated by 1 bond
	for ( Size central_atom1 = 1; central_atom1 < natoms(); ++central_atom1 ) {
		for ( Size central_atom2 = central_atom1+1; central_atom2 <= natoms(); ++central_atom2 ) {
			if ( path_distance_[ central_atom1 ][ central_atom2 ] == 1 ) {

				// get all atoms separated from central_atom1/2 by one bond that are not central_atom2/1
				utility::vector1< Size > ca1d1;
				utility::vector1< Size > ca2d1;

				// ca1
				for ( Size i = 1; i <= natoms(); ++i ) {
					if ( ( path_distance_[ central_atom1 ][ i ] == 1 ) && ( i != central_atom2 ) ) {
						ca1d1.push_back( i );
					}
				}
				// ca2
				for ( Size i = 1; i <= natoms(); ++i ) {
					if ( ( path_distance_[ central_atom2 ][ i ] == 1 ) && ( i != central_atom1 ) ) {
						ca2d1.push_back( i );
					}
				}

				// for each pair of dihedral angle start or end atoms create a dihedral angle using central atom
				for ( auto terminal_atom1 = ca1d1.begin();
						terminal_atom1 != ca1d1.end(); ++terminal_atom1 ) {
					for ( auto terminal_atom2 = terminal_atom1+1;
							terminal_atom2 != ca1d1.end(); ++terminal_atom2 ) {
						dihedral_atom_set temp( *terminal_atom1, central_atom1, central_atom2, *terminal_atom2 );
						improper_dihedral_atom_sets_.push_back( temp );
						Size const which_dihedral = improper_dihedral_atom_sets_.size();
						improper_dihedrals_for_atom_[ *terminal_atom1 ].push_back( which_dihedral );
						improper_dihedrals_for_atom_[   central_atom1 ].push_back( which_dihedral );
						improper_dihedrals_for_atom_[   central_atom2 ].push_back( which_dihedral );
						improper_dihedrals_for_atom_[ *terminal_atom2 ].push_back( which_dihedral );

					}
				}
				for ( auto terminal_atom1 = ca2d1.begin();
						terminal_atom1 != ca2d1.end(); ++terminal_atom1 ) {
					for ( auto terminal_atom2 = terminal_atom1+1;
							terminal_atom2 != ca2d1.end(); ++terminal_atom2 ) {
						dihedral_atom_set temp( *terminal_atom1, central_atom1, central_atom2, *terminal_atom2 );
						improper_dihedral_atom_sets_.push_back( temp );
						Size const which_dihedral = improper_dihedral_atom_sets_.size();
						improper_dihedrals_for_atom_[ *terminal_atom1 ].push_back( which_dihedral );
						improper_dihedrals_for_atom_[   central_atom1 ].push_back( which_dihedral );
						improper_dihedrals_for_atom_[   central_atom2 ].push_back( which_dihedral );
						improper_dihedrals_for_atom_[ *terminal_atom2 ].push_back( which_dihedral );

					}
				}
			}
		}
	}

	// get bond angles
	bondangle_atom_sets_.clear();
	bondangles_for_atom_.resize( natoms() );
	for ( Size ii = 1; ii <= natoms(); ++ii ) bondangles_for_atom_[ ii ].clear();

	// iterate over all atoms that could be a central atom
	for ( Size central_atom = 1; central_atom <= natoms(); ++central_atom ) {

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

	// Now for inter-residue connections, find the sets of atoms that are within one and within two bonds
	// of a residue connection point.  From these sets, all inter-residue bond angle and bond torsions may
	// be enumerated when evaluating residue pair energies.  Also compute the backwards mapping: a list for
	// each atom of the within-1-bond and within-2-bond sets that the atom is listed as being part of. These
	// lists are needed when evaluating atom derivatives wrt the bond dihedrals and angles.
	atoms_within_one_bond_of_a_residue_connection_.resize( residue_connections_.size() );
	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) atoms_within_one_bond_of_a_residue_connection_[ ii ].clear();

	within1bonds_sets_for_atom_.resize( natoms() );
	for ( Size ii = 1; ii <= natoms(); ++ii ) within1bonds_sets_for_atom_[ ii ].clear();

	atoms_within_two_bonds_of_a_residue_connection_.resize( residue_connections_.size() );
	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) atoms_within_two_bonds_of_a_residue_connection_[ ii ].clear();

	within2bonds_sets_for_atom_.resize( natoms() );
	for ( Size ii = 1; ii <= natoms(); ++ii ) within2bonds_sets_for_atom_[ ii ].clear();

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

	// Assign (a) set(s) of possible ring conformations.
	if ( properties_->has_property( CYCLIC ) ) {
		conformer_sets_.resize( n_rings() );
		for ( uint i( 1 ); i <= n_rings(); ++i ) {
			conformer_sets_[ i ] = rings::RingConformerSetOP( new rings::RingConformerSet(
				ring_atoms_[ i ].size(), lowest_ring_conformer_[ i ], low_ring_conformers_[ i ] ) );
		}
	}

	if ( properties_->has_property( RNA ) ) { //reinitialize and RNA derived data.
		//Reinitialize rna_info_ object! This also make sure rna_info_ didn't inherit anything from the previous update!
		//It appears that the rna_info_ is shared across multiple ResidueType object, if the rna_info_ is not reinitialized here!
		rna_info_ = core::chemical::rna::RNA_InfoOP( new core::chemical::rna::RNA_Info );
		//update_last_controlling_chi is treated separately for RNA case. Parin Sripakdeevong, June 26, 2011
		if ( nchi() >= 4 ) {
			// safety against hypothetical RNA RTs without 4 chi AND vs.
			// premature finalize() calls.
			rna_info_->rna_update_last_controlling_chi( get_self_weak_ptr(), last_controlling_chi_, atoms_last_controlled_by_chi_);
			rna_info_->update_derived_rna_data( get_self_weak_ptr() );
		}
	} else if ( properties_->has_property( CARBOHYDRATE ) ) {
		carbohydrate_info_ =
			carbohydrates::CarbohydrateInfoOP( new carbohydrates::CarbohydrateInfo( get_self_weak_ptr() ) );
		update_last_controlling_chi();
	} else {
		update_last_controlling_chi();
	}

	// Set up some atom properties
	for ( Size ii = 1; ii <= natoms(); ++ii ) {
		// Am I an aromatic carbon atom with a free valence?
		// Can also be specified in params, but add if not already set!
		atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, false );
		if ( atom( ii ).mm_name() == "CA" ) {
			AtomIndices const & ii_bonded_neighbors( bonded_neighbor( ii ) );
			for ( Size jj = 1; jj <= ii_bonded_neighbors.size(); ++jj ) {
				if ( atom( ii_bonded_neighbors[ jj ] ).is_haro() ) {
					atom( ii ).set_property( AROMATIC_CARBON_WITH_FREE_VALENCE, true );
				}
			}
		}
	}

	// Set the RamaPrePro potential name to be this residue's name.
	if ( is_base_type() ) {
		if ( !rama_prepro_map_file_name_.empty() && rama_prepro_mainchain_torsion_potential_name_.empty() ) {
			set_rama_prepro_mainchain_torsion_potential_name( name(), false );
		}
		if ( !rama_prepro_map_file_name_beforeproline_.empty() && rama_prepro_mainchain_torsion_potential_name_beforeproline_.empty() ) {
			set_rama_prepro_mainchain_torsion_potential_name( name(), true );
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
/// @brief Final check of ResidueType data, called by finalize().
/// @details These checks are meant to be quick and low-expense, and are only
/// called on finalize(), so they shouldn't generally add much to Rosetta
/// processing time.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
///////////////////////////////////////////////////////////////////////////////
void
ResidueType::perform_checks()
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

	if ( is_metalbinding() && metal_binding_atoms_.size()==0 ) {
		msg << "A metal-binding residue " << name() << " has no metal binding atoms listed in its params file (PROPERTIES METALBINDING "
			"without METAL_BINDING_ATOMS list)." << std::endl;
		checkspass=false;
	} else if ( !is_metalbinding() && metal_binding_atoms_.size()>0 ) {
		msg << "A residue " << name() << " that has not been declared as a metal-binding residue has metal binding atoms listed in its "
			"params file (METAL_BINDING_ATOMS list without PROPERTIES METALBINDING)." << std::endl;
		checkspass=false;
	}

	if ( properties_->has_property( ALPHA_AA ) && properties_->has_property( BETA_AA ) ) {
		msg << "Error!  A residue type " << name() << " specifies that it is both an alpha and a beta amino acid in its params file." <<
			std::endl;
		checkspass=false;
	}
	if ( properties_->has_property( L_AA ) && properties_->has_property( D_AA ) ) {
		msg << "Error!  A residue type " << name() << " specifies that it is both an L-amino acid and a D-amino acid in its params "
			"file." << std::endl;
		checkspass=false;
	}
	if ( (backbone_aa_ != core::chemical::aa_unk) && ! properties_->has_property( ALPHA_AA ) ) {
		msg << "Error!  A residue type " << name() << " specifies a standard alpha amino acid to use as a template for backbone scoring"
			" (rama and p_aa_pp scoring functions) without specifying that it is itself an alpha amino acid "
			"(PROPERTIES ALPHA_AA)." << std::endl;
		checkspass=false;
	}
	if ( (na_analogue_ != core::chemical::aa_unp) && ( ! properties_->has_property( RNA )  || properties_->has_property( CANONICAL_NUCLEIC ) ) ) {
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

	if ( !checkspass ) {
		utility_exit_with_message(msg.str());
	}

	return;
}


///////////////////////////////////////////////////////////////////////////////

/*
data that we have prior to calling this routine:

name                   type                 setting method
----------------------------------------------------------
ordered_atoms_         v1<AtomAP>           add_atom //from base class
atom_name_             v1<string>           add_atom
atom_name_to_vd_       map<string,VD>       add_atom
atomic_charge          v1<Real>             add_atom
bonded_neighbor_       v1<v1<int>>          add_bond
bonded_neighbor_type   v1<v1<BondName>>     add_bond
atom_base_             v1<int>              set_atom_base
chi_atoms_indices      v1<v1<uint>>         add_chi
nu_atoms_indices       v1<v1<uint>>         add_nu
ring_atoms_indices     v1<v1<uint>>         add_ring
properties_            ResidueProperties    add_property
nbr_atom_              int                  nbr_atom( int )

This routine updates all the derived data.

Atoms_ order will probably change after this call, so if you add a new
property that depends on atom-indices that will have to be updated below.
*/
/// @details recalculate derived data, potentially reordering atom-indices
/**
This routine updates all the derived data.\n
Atom order will probably change after this call, so if you add a new
property that depends on atom-indices that will have to be updated below.
**/
//basic graph structure is first state
//setup_atom_ordering and update_derivied data is second state
//third state
void
ResidueType::finalize()
{

	// We've made substantial changes to this ResidueType.
	// On the off chance any observers has cached data about it,
	// we need to notify them that the ResidueType they were observing is effectively destroyed.
	destruction_obs_hub_( RestypeDestructionEvent( this ) );
	destruction_obs_hub_.clear();

	regenerate_graph_vertex_index( graph_ );

	setup_atom_ordering();

	generate_atom_indices();

	update_derived_data();

	perform_checks();

	// signal that derived data is up to date now
	finalized_ = true;

}

////////////////////////////////////////////////////////////////////
utility::vector1< std::string >
ResidueType::variant_types() const
{
	return properties_->get_list_of_variants();
}

Size
ResidueType::atom_index( std::string const & name ) const
{
	// NOTE: Currently we have to iterate twice because atom_name_to_vd_ stores vertex_descriptors not indices.
	// A substantial change to the interface will fix this, but everyone's code will need to switch too.

	auto graph_iter( atom_name_to_vd_.find( name ) );
	if ( graph_iter == atom_name_to_vd_.end() ) {
#if defined BOINC
		// chu temporary graphic fix for boinc
		if ( name == "CA" && !is_protein() ) return 1;
#endif
		if ( name == "CA" && is_membrane() ) return 2;
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		utility_exit_with_message("unknown atom_name: '" + this->name3() + "'  in residue " + name );
	}
	VD const & vd = graph_iter->second;

	Size ordered_index = 0;
	for ( Size i = 1; i <= ordered_atoms_.size(); ++i ) {
		if ( &graph_[ordered_atoms_[i]] == &graph_[vd] ) {
			ordered_index = i;
			break;
		}
	}

	if ( ordered_index == 0 ) {
#if defined BOINC
		// chu temporary graphic fix for boinc
		if ( name == "CA" && !is_protein() ) return 1;
#endif
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		tr.Error << "printing ordered_atomns names" << std::endl;
		for ( Size index=1; index <= ordered_atoms_.size(); ++index ) {
			tr.Error << graph_[ordered_atoms_[index]].name() << " " << &graph_[ordered_atoms_[index]] << std::endl;
		}
		tr.Error << "vd memory address: " << &graph_[vd] << std::endl;
		utility_exit_with_message("unknown atom_name: " + this->name() + "  " + name );
	}

	return ordered_index;
}

Size
ResidueType::atom_index( VD const & vd ) const
{
	if ( vd == ResidueGraph::null_vertex() ) return 0;

	if ( ! vd_to_index_.count(vd) ) {
		utility_exit_with_message("Attempted to get atom index for atom that doesn't appear to be in ResidueType.");
	}
	return vd_to_index_.find(vd)->second;
}

VD ResidueType::atom_vertex(Size const & atomno) const{

	if ( ! atomno ) return ResidueGraph::null_vertex();
	return ordered_atoms_[atomno];
}

core::Size
ResidueType::orbital_index( std::string const & name ) const
{
	auto iter( orbitals_index_.find( name ) );
	if ( iter == orbitals_index_.end() ) {
		utility_exit_with_message("unknown orbital_name: " + name3() + "  " + name );
	}
	return iter->second;
}


void
ResidueType::set_backbone_heavyatom( std::string const & name )
{
	finalized_ = false;
	if ( !has(name) ) {
		utility_exit_with_message("Trying to set bb atom that does not exist in residuetype");
	}
	force_bb_.push_back( atom_name_to_vd_[name]);
}

/// @brief AtomICoord of an atom
AtomICoor const &
ResidueType::icoor( Size const atm ) const
{
	debug_assert( 1 <= atm && atm <= ordered_atoms_.size() );
	auto itr( icoor_.find(ordered_atoms_[atm]) );
	debug_assert( itr != icoor_.end() );
	return itr->second;
}

/// @brief AtomICoord of an atom
AtomICoor const &
ResidueType::icoor( VD const atm ) const
{
	debug_assert( has(atm) );
	auto itr( icoor_.find(atm) );
	debug_assert( itr != icoor_.end() );
	return itr->second;
}

Size
ResidueType::add_residue_connection( std::string const & atom_name )
{
	finalized_ = false;

	++n_non_polymeric_residue_connections_;
	residue_connections_.push_back(
		ResidueConnection( atom_index( atom_name ), ordered_atoms_[atom_index( atom_name )] ) );
	update_residue_connection_mapping();
	return residue_connections_.size();
}


bool
ResidueType::requires_actcoord() const
{
	return properties_->has_property( PROTEIN ) &&
		( properties_->has_property( POLAR ) || properties_->has_property( AROMATIC ) ) &&
		actcoord_atoms_.size() != 0;
}

/// @details average geometrical center of the set of actcoord_atoms_
void
ResidueType::update_actcoord( conformation::Residue & rot ) const
{
	rot.actcoord().zero();
	core::Size const n_actcoord_atoms( actcoord_atoms_.size() );
	if ( n_actcoord_atoms > 0 ) {
		for ( Size ii = 1; ii <= n_actcoord_atoms; ++ii ) {
			rot.actcoord() += rot.atoms()[ vd_to_index_.find(actcoord_atoms_[ ii ])->second ].xyz();
		}
		rot.actcoord() /= n_actcoord_atoms;
	}
}

/// @details set AtomICoor for an atom
///
/// will update the xyz coords as well if desired, useful inside a patching operation where new
/// atoms are being added.
void
ResidueType::set_icoor(
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
	ICoorAtomID id( atm, *this );
	AtomICoor const ic( atm, phi, theta, d, stub_atom1, stub_atom2, stub_atom3, *this );
	set_icoor_private( atm, id, ic, update_xyz );
}

void
ResidueType::set_icoor(
	VD const & atm,
	Real const phi,
	Real const theta,
	Real const d,
	VD const & stub_atom1,
	VD const & stub_atom2,
	VD const & stub_atom3,
	bool const update_xyz // = false
)
{
	debug_assert( has(atm) && has(stub_atom1) && has(stub_atom2) && has(stub_atom3) );

	ICoorAtomID id( atm, *this );
	AtomICoor const ic( atm, phi, theta, d, stub_atom1, stub_atom2, stub_atom3, *this );

	debug_assert( id.type() == ICoorAtomID::INTERNAL ); //It should be, as we're using vertex descriptors

	// previously, a huge chunk of code was just copy and pasted by lazy ass developers.
	set_icoor_private( graph_[atm].name(), id, ic, update_xyz );
}

void
ResidueType::set_icoor(
	std::string const & atm,
	Real const phi,
	Real const theta,
	Real const d,
	ICoorAtomID const & stub_atom1,
	ICoorAtomID const & stub_atom2,
	ICoorAtomID const & stub_atom3,
	bool const update_xyz // = false
)
{
	ICoorAtomID id( atm, *this );
	AtomICoor const ic( atm, phi, theta, d, stub_atom1, stub_atom2, stub_atom3, *this );
	set_icoor_private( atm, id, ic, update_xyz );
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
ResidueType::reset_bond_distance_to_atom( std::string const & atm, core::Distance const d )
{
	AtomICoor const & atm_ic( icoor( atom_index( atm ) ) );
	if ( atm_ic.is_internal() ) {
		set_icoor( atm_ic.built_atom_vertex(), atm_ic.phi(), atm_ic.theta(), d,
			atm_ic.stub_atom1().vertex(), atm_ic.stub_atom2().vertex(), atm_ic.stub_atom3().vertex(), true );
	} else {
		utility_exit_with_message( "ResidueType::reset_bond_distance_to_atom( std::string const & atm, "
			"core:Distance const d ): This function cannot reset ResidueConnections." );
	}
}


//set the orbital icoor data.
void
ResidueType::set_orbital_icoor_id(
	std::string const & orbital,
	Real const phi,
	Real const theta,
	Real const d,
	std::string const & stub_atom1,
	std::string const & stub_atom2,
	std::string const & stub_atom3
)
{

	Size orb_indx(orbital_index(orbital));
	core::Size s1(atom_index( stub_atom1 ));
	core::Size s2(atom_index( stub_atom2 ));
	core::Size s3(atom_index( stub_atom3 ));
	orbitals::ICoorOrbitalData new_icoor(phi, theta, d, s1, s2, s3, ordered_atoms_[s1], ordered_atoms_[s2], ordered_atoms_[s3]);

	orbitals_[ orb_indx ].new_icoor( new_icoor );

}

void
ResidueType::rotamer_library_specification( rotamers::RotamerLibrarySpecificationOP rotlibspec) {
	rotamer_library_specification_ = rotlibspec;
}

rotamers::RotamerLibrarySpecificationCOP
ResidueType::rotamer_library_specification() const {
	return rotamer_library_specification_;
}

/// @brief Nonconst access to the RotamerLibrarySpecification.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
rotamers::RotamerLibrarySpecificationOP
ResidueType::rotamer_library_specification_nonconst() {
	return rotamer_library_specification_;
}


/// @brief Remove any rotamer library specifications attached to this ResidueType.
/// @details After this operation, the rotamer_library_specification() method returns a NULL pointer.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
ResidueType::strip_rotamer_library_specification() {
	rotamer_library_specification_.reset();
}

void ResidueType::assign_neighbor_atom()
{
	//calculate the geometric center of all atoms in the residue

	Vector total(0.0,0.0,0.0);
	for ( core::Size index = 1; index <= ordered_atoms_.size(); ++index ) {
		total += graph_[ordered_atoms_[index]].ideal_xyz();
	}

	Vector center = total/ordered_atoms_.size();

	//locate the atom which is closest to the center
	core::Size min_index = 0;
	core::Real min_distance = 50000.0;

	for ( core::Size index = 1; index <= ordered_atoms_.size(); ++index ) {
		core::Real distance = center.distance(graph_[ordered_atoms_[index]].ideal_xyz());
		if ( (distance < min_distance) && (!atom_is_hydrogen(index)) ) {
			min_distance = distance;
			min_index = index;
		}
	}
	debug_assert(min_index != 0);
	//set neighbor atom
	nbr_atom(graph_[ordered_atoms_[min_index]].name());
}

void ResidueType::assign_internal_coordinates()
{
	// Reuse the existing root, or failing that, the neighbor atom
	// As a last resort, just use atom #1
	VD new_root( root_atom_ );
	if ( new_root == ResidueType::null_vertex ) {
		new_root = nbr_atom_;
	}
	if ( new_root == ResidueType::null_vertex ) {
		new_root = ordered_atoms_[1];
	}
	debug_assert( new_root != ResidueType::null_vertex );
	assign_internal_coordinates( new_root );
}

void ResidueType::assign_internal_coordinates(core::chemical::VD new_root)
{
	//%TODO: right now we're ignoring M FRAG lines and M SPLT lines in molfiles
	if ( n_possible_residue_connections() != 0 ) {
		tr.Error << "Residue " << name() << " has connections - can't assign internal coordinates.";
		utility_exit_with_message("Cannot currently assign internal coordinates for polymeric residue.");
	}
	debug_assert( new_root != ResidueType::null_vertex );
	// Reset the root atom so we can re-root the tree
	root_atom_ = ResidueType::null_vertex;
	reroot_restype(*this, graph_, new_root);
}

void
ResidueType::set_ideal_xyz(
	std::string const & atm,
	Vector const & xyz_in
)
{
	Size const index( atom_index(atm) );
	set_ideal_xyz(index,xyz_in);
}

void
ResidueType::set_ideal_xyz(
	Size index,
	Vector const & xyz_in
)
{
	if ( index > ordered_atoms_.size() ) {
		utility_exit_with_message("Cannot set ideal coordinates for non-existent atom.");
	}
	Atom & a = graph_[ ordered_atoms_[index] ];
	a.ideal_xyz( xyz_in );
}

void
ResidueType::set_ideal_xyz(
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

void
ResidueType::fill_ideal_xyz_from_icoor() {
	debug_assert( natoms() == icoor_.size() );
	core::chemical::fill_ideal_xyz_from_icoor(*this, graph_);
}

void
ResidueType::set_shadowing_atom(
	std::string const & atom,
	std::string const & atom_being_shadowed
)
{

	VD const index_shadower( ordered_atoms_[atom_index(atom)] );
	VD const index_shadowee( ordered_atoms_[ atom_index(atom_being_shadowed) ]);
	atom_shadowed_[ index_shadower ] = index_shadowee;
}

// Return the CarbohydrateInfo object containing sugar-specific properties for this residue.
core::chemical::carbohydrates::CarbohydrateInfoCOP
ResidueType::carbohydrate_info() const
{
	return carbohydrate_info_;
}

void
ResidueType::print_dihedrals() const
{
	core::Size ndihe( dihedral_atom_sets_.size() );
	tr.Debug << "START DIHEDRAL ANGLES ATOM NAMES" << std::endl;
	tr.Debug << "Number of dihe: " << ndihe << std::endl;
	for ( Size i = 1; i <= ndihe; ++i ) {

		AtomType at1 = atom_type( dihedral_atom_sets_[ i ].key1() );
		AtomType at2 = atom_type( dihedral_atom_sets_[ i ].key2() );
		AtomType at3 = atom_type( dihedral_atom_sets_[ i ].key3() );
		AtomType at4 = atom_type( dihedral_atom_sets_[ i ].key4() );
		MMAtomType at5 = mm_atom_type( dihedral_atom_sets_[ i ].key1() );
		MMAtomType at6 = mm_atom_type( dihedral_atom_sets_[ i ].key2() );
		MMAtomType at7 = mm_atom_type( dihedral_atom_sets_[ i ].key3() );
		MMAtomType at8 = mm_atom_type( dihedral_atom_sets_[ i ].key4() );

		tr.Debug << "PDB:" << "\t"
			<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key1() ]].name() << "\t"
			<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key2() ]].name() << "\t"
			<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key3() ]].name() << "\t"
			<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key4() ]].name() << "\t"
			<< "MM:" << "\t"
			<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key1() ]].mm_name() << "\t"
			<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key2() ]].mm_name() << "\t"
			<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key3() ]].mm_name() << "\t"
			<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key4() ]].mm_name() << "\t"
			<< "MM2:" << "\t"
			<< at5.name() << "\t"
			<< at6.name() << "\t"
			<< at7.name() << "\t"
			<< at8.name() << "\t"
			<< "ROS:" << "\t"
			<< at1.name() << "\t"
			<< at2.name() << "\t"
			<< at3.name() << "\t"
			<< at4.name() << "\t"
			<< std::endl;
	}
	tr.Debug << "END DIHEDRAL ANGLES ATOM NAMES" << std::endl;
}

void
ResidueType::print_bondangles() const
{
	tr.Debug << "START BOND ANGLES ATOM NAMES" << std::endl;
	tr.Debug << "Number of bond angles: " << bondangle_atom_sets_.size() << std::endl;
	for ( Size i = 1; i <= bondangle_atom_sets_.size(); ++i ) {
		AtomType at1 = atom_type( bondangle_atom_sets_[ i ].key1() );
		AtomType at2 = atom_type( bondangle_atom_sets_[ i ].key2() );
		AtomType at3 = atom_type( bondangle_atom_sets_[ i ].key3() );
		MMAtomType at5 = mm_atom_type( bondangle_atom_sets_[ i ].key1() );
		MMAtomType at6 = mm_atom_type( bondangle_atom_sets_[ i ].key2() );
		MMAtomType at7 = mm_atom_type( bondangle_atom_sets_[ i ].key3() );

		tr.Debug << "PDB:" << "\t"
			<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key1() ]].name() << "\t"
			<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key2() ]].name() << "\t"
			<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key3() ]].name() << "\t"
			<< "MM:" << "\t"
			<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key1() ]].mm_name() << "\t"
			<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key2() ]].mm_name() << "\t"
			<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key3() ]].mm_name() << "\t"
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
	tr.Debug << "END BOND ANGLES ATOM NAMES" << std::endl;
}

void
ResidueType::print_pretty_path_distances() const
{
	tr.Debug << "START PATH DISTANCES" << std::endl;
	// print header line
	for ( Size i = 1; i <= natoms(); ++i ) {
		tr.Debug << "\t" << graph_[ordered_atoms_[i]].name();
	}
	tr.Debug << std::endl;

	for ( Size j = 1; j <= natoms(); ++j ) {
		tr.Debug << graph_[ordered_atoms_[j]].name() << "\t";
		for ( Size k = 1; k <= natoms(); ++k ) {
			tr.Debug << path_distance_[j][k] << "\t";
		}
		tr.Debug << std::endl;
	}
	tr.Debug << "END PATH DISTANCES" << std::endl;
}

void
ResidueType::update_residue_connection_mapping()
{
	//std::fill( atom_2_residue_connection_map_.begin(), atom_2_residue_connection_map_.end(), 0 );
	for ( Size ii = 1; ii <= natoms(); ++ii ) { atom_2_residue_connection_map_[ ii ].clear(); }

	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) {
		atom_2_residue_connection_map_[ vd_to_index_.find(residue_connections_[ ii ].vertex())->second ].push_back( ii );
		residue_connections_[ ii ].index( ii );
	}
}

void
ResidueType::update_last_controlling_chi() {
	last_controlling_chi_.resize( natoms() );
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
	for ( Size ii = 1; ii <= natoms(); ++ii ) {
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
	for ( Size ii = 1; ii <= natoms(); ++ii ) {
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

void
ResidueType::set_icoor_private(
	std::string const & atm,
	ICoorAtomID const & id,
	AtomICoor const & ic,
	bool update_xyz
)
{
	VD atom_vd( id.vertex() );
	switch ( id.type() ) {
	case ICoorAtomID::INTERNAL :
		//debug_assert( atom_vd == atom_vertex( atm ) );
		if ( atom_vd == ic.stub_atom1().vertex() ) {
			//Root atom case
			if ( root_atom_ != ResidueType::null_vertex ) {
				tr.Error << "Can't reset root. Was " << atom_name( root_atom_ ) << " Resetting to " << atm << std::endl;
				utility_exit_with_message( "Attempted to inappropriately reset ICOOR root atom." );
			}
			root_atom_ = atom_vd;
			// Now we can continue as normal.
		}
		icoor_[ atom_vd ] = ic;

		// update atom_base if it isn't set yet (or is set to itself)
		if ( ! atom_base_.count(atom_vd) || atom_base_[atom_vd] == atom_vd ) {
			if ( atom_vd == ic.stub_atom1().vertex() ) {
				//root of tree
				if ( natoms() == 1 ) {
					set_atom_base( atm, atm );
				} else {
					set_atom_base( atom_vd, ic.stub_atom2().vertex() );
				}
			} else {
				set_atom_base( atom_vd, ic.stub_atom1().vertex() );
			}
		}
		if ( update_xyz ) {
			set_ideal_xyz( atm, ic.build( *this ) );
			//std::cout << "building coords for atm " << name_ << ' ' << atm << ' ' <<
			//  ic.build(*this)(1) << ' ' <<
			//  ic.build(*this)(2) << ' ' <<
			//  ic.build(*this)(3) << std::endl;
		}
		break;
	case ICoorAtomID::CONNECT :
		// For CONNECT, the atomno is repurposed as the connection number
		residue_connections_[ id.atomno() ].icoor( ic );
		break;
	case ICoorAtomID::POLYMER_LOWER :
		debug_assert( lower_connect_id_ != 0 );
		residue_connections_[ lower_connect_id_ ].icoor( ic );
		break;
	case ICoorAtomID::POLYMER_UPPER :
		debug_assert( upper_connect_id_ != 0 );
		residue_connections_[ upper_connect_id_ ].icoor( ic );
		break;
	default :
		utility_exit_with_message( "unrecognized stub atom id type!" );
		break; //to silence warning
	}
}

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
			utility_exit_with_message("Cannot superimpose residues of type "+name());
		}
		//std::cout << "Superimposing on " << atom_name(center) << " " << atom_name(nbr1) << " " << atom_name(nbr2) << "\n";

	} else {
		// look for a backbone atom, one of whose neighbors is a sidechain atom
		// center will be this atom
		// nbr1 and nbr2 will be the backbone heavyatom nbrs of this atom
		// eg center = CA, nbr1 = N. nbr2 = C in the protein case
		for ( Size atom_index(1); atom_index <= natoms(); ++atom_index ) {
			if ( atom_is_backbone( atom_index ) ) {
				AtomIndices const & nbrs( bonded_neighbor( atom_index ) );
				center = 0; nbr1 = 0; nbr2 = 0;
				for ( Size nbr_index(1); nbr_index <= nbrs.size(); ++nbr_index ) {
					Size const nbr( nbrs[ nbr_index ] );
					if ( !atom_is_backbone( nbr ) && atom_base( nbr ) == atom_index ) {
						// nbr is a sidechain atom that branches from the atom at atom_index
						center = atom_index;
					} else if ( atom_is_backbone( nbr ) && nbr <= nheavyatoms() ) {
						// nbr is a backbone heavy atom neighbor of the atom at atom_index
						if ( nbr1 ) nbr2 = nbr;
						else nbr1 = nbr;
					}
				}
			} // atom_index is backbone
			if ( center && nbr1 && nbr2 ) break;
		} // atom_index
	}
}

// A graph-based function to determine the size of the smallest ring that involves a given atom.
core::Size
ResidueType::smallest_ring_size( VD const & atom, core::Size const & max_size /*= 999999*/ ) const
{
	return utility::graph::smallest_ring_size(atom, graph_, max_size);
}

void
ResidueType::report_adducts()
{
	if ( defined_adducts_.size() == 0 ) return;

	for ( Size ii = 1 ; ii <= defined_adducts_.size() ; ++ii ) {
		Adduct & add( defined_adducts_[ii] );
		tr.Debug << "Residue: " << name3() << " Adduct: " << add.adduct_name() <<
			" Atom name: " << add.atom_name() << std::endl;
	}
}

void
ResidueType::debug_dump_icoor() const
{

	tr.Debug << "ICoor for " << name3() << std::endl;
	for ( Size ii = 1 ; ii <= natoms() ; ++ii ) {
		tr.Debug << " Atom name: " << atom_name( ii ) << "vertex: " << atom_vertex( ii ) << " ideal xyz " << atom(ii).ideal_xyz()[0] << "  " << atom(ii).ideal_xyz()[1] << "  " << atom(ii).ideal_xyz()[2] << std::endl;
	}
	pretty_print_atomicoor(tr.Debug, *this);
}

/// @brief Get the key name for the mainchain torsion potential.
/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
/// residue types (in which case this function returns the string stored in the base ResidueType), though this can be overridden.
std::string const &
ResidueType::get_rama_prepro_mainchain_torsion_potential_name( bool const pre_proline_position ) const {

	if ( pre_proline_position ) {
		if ( !rama_prepro_mainchain_torsion_potential_name_beforeproline_.empty() ) return rama_prepro_mainchain_torsion_potential_name_beforeproline_;
		//If the string is empty...:
		if ( !is_base_type() ) {
			std::string const & basestring( get_base_type_cop()->get_rama_prepro_mainchain_torsion_potential_name(pre_proline_position) );
			if ( !basestring.empty() ) return basestring;
		}
		return rama_prepro_mainchain_torsion_potential_name_beforeproline_; //Returns an empty string if this is empty AND the base type is empty.
	}

	//Otherwise...

	if ( !rama_prepro_mainchain_torsion_potential_name_.empty() ) return rama_prepro_mainchain_torsion_potential_name_;
	//If the string is empty...:
	if ( !is_base_type() ) {
		std::string const & basestring( get_base_type_cop()->get_rama_prepro_mainchain_torsion_potential_name(pre_proline_position) );
		if ( !basestring.empty() ) return basestring;
	}
	return rama_prepro_mainchain_torsion_potential_name_; //Returns an empty string if this is empty AND the base type is empty.
}

/// @brief Do the rama_prepro mainchain torsion potentials of this residue match another?
///
bool
ResidueType::mainchain_potentials_match(
	ResidueType const &other
) const {
	return
		rama_prepro_mainchain_torsion_potential_name_.compare( other.rama_prepro_mainchain_torsion_potential_name_ ) == 0 &&
		rama_prepro_mainchain_torsion_potential_name_beforeproline_.compare( other.rama_prepro_mainchain_torsion_potential_name_beforeproline_ ) == 0 &&
		rama_prepro_map_file_name_.compare( other.rama_prepro_map_file_name_ ) == 0 &&
		rama_prepro_map_file_name_beforeproline_.compare( other.rama_prepro_map_file_name_beforeproline_ ) == 0
		;
}


/// @brief Set the key name for the mainchain torsion potential.
/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
/// residue types (pointing the function to the base type), though this can be overridden using this function.
void
ResidueType::set_rama_prepro_mainchain_torsion_potential_name(
	std::string const &name_in,
	bool const pre_proline_position
) {
	if ( pre_proline_position ) {
		rama_prepro_mainchain_torsion_potential_name_beforeproline_ = name_in;
	} else {
		rama_prepro_mainchain_torsion_potential_name_ = name_in;
	}
}

/// @brief Get the file name for the mainchain torsion potential used by the RamaPrePro score term.
/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
/// residue types (in which case this function returns the string stored in the base ResidueType), though this can be overridden.
std::string const &
ResidueType::get_rama_prepro_map_file_name( bool const pre_proline_position ) const {
	if ( pre_proline_position ) {
		if ( !rama_prepro_map_file_name_beforeproline_.empty() ) return rama_prepro_map_file_name_beforeproline_;
		//If the string is empty...:
		if ( !is_base_type() ) {
			std::string const & basestring( get_base_type_cop()->get_rama_prepro_map_file_name(pre_proline_position) );
			if ( !basestring.empty() ) return basestring;
		}
		return rama_prepro_map_file_name_beforeproline_; //Returns an empty string if this is empty AND the base type is empty.
	}

	//Otherwise...

	if ( !rama_prepro_map_file_name_.empty() ) return rama_prepro_map_file_name_;
	//If the string is empty...:
	if ( !is_base_type() ) {
		std::string const & basestring( get_base_type_cop()->get_rama_prepro_map_file_name(pre_proline_position) );
		if ( !basestring.empty() ) return basestring;
	}
	return rama_prepro_map_file_name_; //Returns an empty string if this is empty AND the base type is empty.
}

/// @brief Set the file name for the mainchain torsion potential used by the RamaPrePro score term.
/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
/// residue types (pointing the function to the base type), though this can be overridden using this function.
void
ResidueType::set_rama_prepro_map_file_name(
	std::string const &filename_in,
	bool const pre_proline_position
) {
	if ( pre_proline_position ) {
		rama_prepro_map_file_name_beforeproline_ = filename_in;
	} else {
		rama_prepro_map_file_name_ = filename_in;
	}
}

/// @brief Returns true if and only if (a) this is not a base type, AND (b) there is a rama_prepro_mainchain_torsion_map_file_name_
/// defined for this ResidueType (which is presumably different from that of the base type).
/// @details If pre_proline_position is true, checks rama_prepro_mainchain_torsion_map_file_name_beforeproline_ instead of
/// rama_prepro_mainchain_torsion_potential_name_.
bool
ResidueType::defines_custom_rama_prepro_map( bool const pre_proline_position ) const {
	if ( is_base_type() ) return false;
	if ( pre_proline_position ) {
		return ( !rama_prepro_map_file_name_beforeproline_.empty() );
	}
	return !rama_prepro_map_file_name_.empty();
}

/// @brief Set the names of the mainchain torsion potential maps to use to "".
///
void
ResidueType::reset_mainchain_torsion_potential_names() {
	rama_prepro_mainchain_torsion_potential_name_.clear();
	rama_prepro_map_file_name_.clear();
	rama_prepro_mainchain_torsion_potential_name_beforeproline_.clear();
	rama_prepro_map_file_name_beforeproline_.clear();
}


void
ResidueType::show_all_atom_names( std::ostream & out ) const {

	for ( VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first ) {
		auto v_iter= vp.first;
		VD vd = *v_iter;
		Atom a = graph_[vd];
		out << "'" << a.name() << "' " << &graph_[vd] << std::endl;
	}

}

/// @brief  Check if atom is virtual.
bool
ResidueType::is_virtual( Size const & atomno ) const
{
	return ( atom_type( atomno ).is_virtual() );
}

/// @brief  Check if atom is repulsive.
bool
ResidueType::is_repulsive( Size const & atomno ) const
{
	return ( atom_type( atomno ).is_repulsive() );
}


///////////////////////////////////////////////////////////////
core::chemical::rna::RNA_Info const &
ResidueType::RNA_info() const{
	return ( *rna_info_ );
}

/// @author Labonte <JWLabonte@jhu.edu>
void
ResidueType::show( std::ostream & output, bool output_atomic_details ) const
{
	using namespace std;
	using namespace utility;

	output << name_ << " (" << name3_ << ", " << name1_ << "):" << endl;

	output << "Base: " << base_name_ << std::endl;

	properties_->show( output );

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
			Size const n_ring_atoms( ring_atoms_indices_[ i ].size() );
			for ( uint j = 1; j <= n_ring_atoms; ++j ) {
				output << ' ' << atom_name( ring_atoms_indices_[ i ][ j ] );
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

	if ( properties_->has_property( CARBOHYDRATE ) ) {
		carbohydrate_info_->show( output );
	}

	if ( output_atomic_details ) {
		output << " Atomic Details:" << endl;
		Size const n_atoms( natoms() );
		for ( uint i( 1 ); i <= n_atoms; ++i ) {
			output << "  Atom " << i << ": ";
			atom( i ).show( output );
		}
	}
}

void
ResidueType::set_atom_type_set( AtomTypeSetCOP setting ) {
	debug_assert( setting );
	atom_types_ = setting;
	mode_ = atom_types_->mode();
}

//////////////////////////////////////////////////////
/// Make all atoms virtual
/// @author Sebastian Rämisch <raemisch@scripps.edu>
void
ResidueType::real_to_virtual() {
	std::string VIRT = "VIRT";
	for ( Size i=1; i<=this->natoms(); ++i ) {
		this->set_atom_type( (this->atom_name(i) ), VIRT);
		this->atom(i).charge(0.0);
		this->atom(i).is_virtual( true );
	}
	this->add_property("VIRTUAL_RESIDUE");
	this->finalize();

}



// Helper methods //////////////////////////////////////////////////////////////
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

	arc( CEREAL_NVP( atom_types_ ) ); // AtomTypeSetCOP
	arc( CEREAL_NVP( elements_ ) ); // ElementSetCOP
	arc( CEREAL_NVP( mm_atom_types_ ) ); // MMAtomTypeSetCOP
	arc( CEREAL_NVP( gasteiger_atom_types_ ) ); // gasteiger::GasteigerAtomTypeSetCOP
	arc( CEREAL_NVP( orbital_types_ ) ); // orbitals::OrbitalTypeSetCOP

	runtime_assert( finalized_ ); // can't serialize if it's not finalized
	// EXEMPT graph_ ordered_atoms_ vd_to_index_ icoor_
	// We can't serialize the Boost graph directly, and, besides, the VDs will change
	// when we reconstruct the graph. To circumvent this, we serialize atoms directly
	// in index order (as atom names might not be unique ... theoretically).
	arc( CEREAL_NVP_( "natoms", ordered_atoms_.size() ) );
	for( VD vd : ordered_atoms_ ) {
		SERIALIZE_VD( arc, vd );
		arc( graph_[ vd ] ); // EXEMPT graph_ ordered_atoms_ vd_to_index_ atom_name_to_vd_
	}

	//std::cout << "DONE ATOMS" << std::endl;

	arc( CEREAL_NVP_( "bonds", boost::num_edges(graph_) ) );
	for( EIterPair biter_pair( boost::edges(graph_) ); biter_pair.first != biter_pair.second; ++biter_pair.first ) {
		ED ed( *biter_pair.first );
		VD source( boost::source( ed, graph_ ) ), target( boost::target( ed, graph_ ) );
		SERIALIZE_VD( arc, source );
		SERIALIZE_VD( arc, target );
		arc( graph_[ ed ] );
	}

	for( VD vd : ordered_atoms_ ) {
		arc( icoor_.at( vd ) ); // EXEMPT icoor_
	}

	// Need to be VD-adjusted on the back end.
	arc( CEREAL_NVP( residue_connections_ ) ); // utility::vector1<ResidueConnection>
	arc( CEREAL_NVP( orbitals_ ) ); // utility::vector1<Orbital>

	SERIALIZE_VD( arc, root_atom_, "root_atom_" ); // EXEMPT root_atom_
	SERIALIZE_VD( arc, nbr_atom_, "nbr_atom_" ); // EXEMPT nbr_atom_

	SERIALIZE_VD_VD_MAP( arc, atom_base_ ); // EXEMPT atom_base_
	SERIALIZE_VD_VD_MAP( arc, abase2_ ); // EXEMPT abase2_
	SERIALIZE_VD_VD_MAP( arc, atom_shadowed_ ); // EXEMPT atom_shadowed_

	SERIALIZE_VD_VD_VECTOR_MAP( arc, cut_bond_neighbor_ ); // EXEMPT cut_bond_neighbor_

	SERIALIZE_VD_VECTOR( arc, mainchain_atoms_ ); // EXEMPT mainchain_atoms_
	SERIALIZE_VD_VECTOR( arc, actcoord_atoms_ ); // EXEMPT actcoord_atoms_
	SERIALIZE_VD_VECTOR( arc, force_bb_ ); // EXEMPT force_bb_

	SERIALIZE_NESTED_VD_VECTOR( arc, chi_atoms_ ); // EXEMPT chi_atoms_
	SERIALIZE_NESTED_VD_VECTOR( arc, nu_atoms_ ); // EXEMPT nu_atoms_
	SERIALIZE_NESTED_VD_VECTOR( arc, ring_atoms_ ); // EXEMPT ring_atoms_

	arc( rings_and_their_edges_.size() ); // EXEMPT rings_and_their_edges_
	for ( utility::vector1<ED> const & innervec : rings_and_their_edges_ ) {
		arc( innervec.size() );
		for ( ED ed : innervec ) {
			SERIALIZE_VD( arc, boost::source( ed, graph_ ) );
			SERIALIZE_VD( arc, boost::target( ed, graph_ ) );
		}
	}

	// Many of the following will be reset in finalize(), but I serialize them here anyway because it's just as easy as not doing so

	// RingConformerSets are reset in finalize()
	// EXEMPT conformer_sets_
	arc( CEREAL_NVP( mode_ ) ); // enum core::chemical::TypeSetMode
	arc( CEREAL_NVP( nheavyatoms_ ) ); // Size
	arc( CEREAL_NVP( n_hbond_acceptors_ ) ); // Size
	arc( CEREAL_NVP( n_hbond_donors_ ) ); // Size
	arc( CEREAL_NVP( n_backbone_heavyatoms_ ) ); // Size
	arc( CEREAL_NVP( first_sidechain_hydrogen_ ) ); // Size
	arc( CEREAL_NVP( bonded_neighbor_ ) ); // utility::vector1<AtomIndices>
	arc( CEREAL_NVP( bonded_neighbor_type_ ) ); // utility::vector1<utility::vector1<BondName> >
	arc( CEREAL_NVP( attached_H_begin_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( attached_H_end_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( dihedral_atom_sets_ ) ); // utility::vector1<dihedral_atom_set>
	arc( CEREAL_NVP( dihedrals_for_atom_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( improper_dihedral_atom_sets_ ) ); // utility::vector1<dihedral_atom_set>
	arc( CEREAL_NVP( improper_dihedrals_for_atom_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( bondangle_atom_sets_ ) ); // utility::vector1<bondangle_atom_set>
	arc( CEREAL_NVP( bondangles_for_atom_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( last_controlling_chi_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( atoms_last_controlled_by_chi_ ) ); // utility::vector1<AtomIndices>
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
	arc( CEREAL_NVP( metal_binding_atoms_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( disulfide_atom_name_ ) ); // std::string
	arc( CEREAL_NVP( is_proton_chi_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( proton_chis_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( chi_2_proton_chi_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( proton_chi_samples_ ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( proton_chi_extra_samples_ ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( path_distance_ ) ); // utility::vector1<utility::vector1<int> >
	arc( CEREAL_NVP( atom_aliases_ ) ); // std::map<std::string, std::string>
	arc( CEREAL_NVP( orbitals_index_ ) ); // std::map<std::string, int>
	arc( CEREAL_NVP( chi_rotamers_ ) ); // utility::vector1<utility::vector1<std::pair<Real, Real> > >
	arc( CEREAL_NVP( rotamer_library_specification_ ) ); // rotamers::RotamerLibrarySpecificationOP
	arc( CEREAL_NVP( ring_sizes_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( lowest_ring_conformer_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( low_ring_conformers_ ) ); // utility::vector1<utility::vector1<std::string> >
	arc( CEREAL_NVP( properties_ ) ); // ResiduePropertiesOP
	arc( CEREAL_NVP( aa_ ) ); // enum core::chemical::AA
	arc( CEREAL_NVP( rotamer_aa_ ) ); // enum core::chemical::AA
	arc( CEREAL_NVP( backbone_aa_ ) ); // enum core::chemical::AA
	arc( CEREAL_NVP( na_analogue_ ) ); // enum core::chemical::AA
	arc( CEREAL_NVP( base_name_ ) ); // std::string
	arc( CEREAL_NVP( base_type_cop_ ) ); // ResidueTypeCOP
	arc( CEREAL_NVP( name_ ) ); // std::string
	arc( CEREAL_NVP( name3_ ) ); // std::string
	arc( CEREAL_NVP( name1_ ) ); // char
	arc( CEREAL_NVP( interchangeability_group_ ) ); // std::string
	arc( CEREAL_NVP( nbr_radius_ ) ); // Real
	arc( CEREAL_NVP( force_nbr_atom_orient_ ) ); // _Bool
	arc( CEREAL_NVP( remap_pdb_atom_names_ ) ); // _Bool
	arc( CEREAL_NVP( mass_ ) ); // Real
	arc( CEREAL_NVP( atom_2_residue_connection_map_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( atoms_within_one_bond_of_a_residue_connection_ ) ); // utility::vector1<utility::vector1<two_atom_set> >
	arc( CEREAL_NVP( within1bonds_sets_for_atom_ ) ); // utility::vector1<utility::vector1<std::pair<Size, Size> > >
	arc( CEREAL_NVP( atoms_within_two_bonds_of_a_residue_connection_ ) ); // utility::vector1<utility::vector1<three_atom_set> >
	arc( CEREAL_NVP( within2bonds_sets_for_atom_ ) ); // utility::vector1<utility::vector1<std::pair<Size, Size> > >
	arc( CEREAL_NVP( lower_connect_id_ ) ); // Size
	arc( CEREAL_NVP( upper_connect_id_ ) ); // Size
	arc( CEREAL_NVP( n_non_polymeric_residue_connections_ ) ); // Size
	arc( CEREAL_NVP( n_polymeric_residue_connections_ ) ); // Size
	// RNA_ResidueType and CarbohydrateInfo are reset in finalize()
	// EXEMPT rna_info_ carbohydrate_info_
	arc( CEREAL_NVP( atom_base_indices_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( abase2_indices_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( chi_atoms_indices_ ) ); // utility::vector1<AtomIndices>
	arc( CEREAL_NVP( nu_atoms_indices_ ) ); // utility::vector1<AtomIndices>
	arc( CEREAL_NVP( ring_atoms_indices_ ) ); // utility::vector1<AtomIndices>
	arc( CEREAL_NVP( mainchain_atoms_indices_ ) ); // AtomIndices
	arc( CEREAL_NVP( nbr_atom_indices_ ) ); // Size
	arc( CEREAL_NVP( actcoord_atoms_indices_ ) ); // AtomIndices
	arc( CEREAL_NVP( cut_bond_neighbor_indices_ ) ); // utility::vector1<AtomIndices>
	arc( CEREAL_NVP( atom_shadowed_indices_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( rama_prepro_mainchain_torsion_potential_name_ ) ); // std::string
	arc( CEREAL_NVP( rama_prepro_mainchain_torsion_potential_name_beforeproline_ ) ); // std::string
	arc( CEREAL_NVP( rama_prepro_map_file_name_ ) );// std::string
	arc( CEREAL_NVP( rama_prepro_map_file_name_beforeproline_ ) );// std::string
	// EXEMPT finalized_
	// ( will call finalization function on load )
	arc( CEREAL_NVP( defined_adducts_ ) ); // utility::vector1<Adduct>
	arc( CEREAL_NVP( nondefault_ ) ); // _Bool

	// Observers aren't being serialized - any observer on the remote side will have to reattach itself.
	// EXEMPT destruction_obs_hub_
	// EXEMPT destruction_obs_mutex_;
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ResidueType::load( Archive & arc ) {
	using namespace core::chemical;

	arc( atom_types_ ); // AtomTypeSetCOP
	arc( elements_ ); // ElementSetCOP
	arc( mm_atom_types_ ); // MMAtomTypeSetCOP
	arc( gasteiger_atom_types_ ); // gasteiger::GasteigerAtomTypeSetCOP
	arc( orbital_types_ ); // orbitals::OrbitalTypeSetCOP

	std::map< VD, VD > old_to_new;
	old_to_new[ ResidueType::null_vertex ] = ResidueType::null_vertex; // Null vertex self-converts.

	// EXEMPT graph_
	core::Size natoms; arc( natoms );
	for( core::Size ii(1); ii <= natoms; ++ii ) {
		VD old_vd; DESERIALIZE_VD( arc, old_vd );
		Atom atom; arc( atom );
		atom.update_typesets( *this );
		VD new_vd = graph_.add_vertex(atom);
		old_to_new[ old_vd ] = new_vd;
		//std::cout << "Old to new " << old_vd << " : " << new_vd << std::endl;
		ordered_atoms_.push_back(new_vd); // EXEMPT ordered_atoms_
		vd_to_index_[ new_vd ] = ii; // EXEMPT vd_to_index_
		atom_name_to_vd_[ atom.name() ] = new_vd; // EXEMPT atom_name_to_vd_
		atom_name_to_vd_[ strip_whitespace( atom.name() ) ] = new_vd;
	}

	//std::cout << "DONE ATOMS" << std::endl;

	//std::map< ED, ED > old_to_new_bonds;

	core::Size nbonds; arc( nbonds );
	for( core::Size ii(1); ii <= nbonds; ++ii ) {
		VD source, target;
		DESERIALIZE_VD( arc, source );
		source = old_to_new.at(source);
		DESERIALIZE_VD( arc, target );
		target = old_to_new.at(target);
		Bond bond; arc( bond );
		bool added; ED new_ed;
		boost::tie( new_ed, added ) = graph_.add_edge( source, target, bond );
		runtime_assert( added );
		//old_to_new_bonds[ old_ed ] = new_ed;
	}

	// Now we have a completely set old_to_nes
	for( core::Size ii(1); ii <= natoms; ++ii ) {
		VD vd( ordered_atoms_[ii] );

		AtomICoor icoor; arc( icoor );
		icoor.remap_atom_vds( old_to_new );
		icoor_[ vd ] = icoor; // EXEMPT icoor_

	}

	arc( residue_connections_ ); // utility::vector1<ResidueConnection>
	for( ResidueConnection & rescon : residue_connections_ ) {
		rescon.remap_atom_vds( old_to_new );
	}
	arc( orbitals_ ); // utility::vector1<Orbital>
	for( Orbital & orbital : orbitals_ ) {
		orbital.remap_atom_vds( old_to_new );
	}

	DESERIALIZE_VD( arc, root_atom_, old_to_new ); // EXEMPT root_atom_
	DESERIALIZE_VD( arc, nbr_atom_, old_to_new ); // EXEMPT nbr_atom_

	DESERIALIZE_VD_VD_MAP( arc, atom_base_, old_to_new ); // EXEMPT atom_base_
	DESERIALIZE_VD_VD_MAP( arc, abase2_, old_to_new ); // EXEMPT abase2_
	DESERIALIZE_VD_VD_MAP( arc, atom_shadowed_, old_to_new ); // EXEMPT atom_shadowed_

	DESERIALIZE_VD_VD_VECTOR_MAP( arc, cut_bond_neighbor_, old_to_new ); // EXEMPT cut_bond_neighbor_

	DESERIALIZE_VD_VECTOR( arc, mainchain_atoms_, old_to_new ); // EXEMPT mainchain_atoms_
	DESERIALIZE_VD_VECTOR( arc, actcoord_atoms_, old_to_new ); // EXEMPT actcoord_atoms_
	DESERIALIZE_VD_VECTOR( arc, force_bb_, old_to_new ); // EXEMPT force_bb_

	DESERIALIZE_NESTED_VD_VECTOR( arc, chi_atoms_, old_to_new ); // EXEMPT chi_atoms_
	DESERIALIZE_NESTED_VD_VECTOR( arc, nu_atoms_, old_to_new ); // EXEMPT nu_atoms_
	DESERIALIZE_NESTED_VD_VECTOR( arc, ring_atoms_, old_to_new ); // EXEMPT ring_atoms_

	core::Size r_a_t_e_size; arc( r_a_t_e_size ); // EXEMPT rings_and_their_edges_
	rings_and_their_edges_.clear();
	for ( core::Size ii(1); ii <= r_a_t_e_size; ++ii ) {
		utility::vector1<ED> innervec;
		core::Size innervec_size; arc( innervec_size );
		for ( core::Size jj(1); jj <= innervec_size; ++jj ) {
			VD source, target;
			DESERIALIZE_VD( arc, source, old_to_new );
			DESERIALIZE_VD( arc, target, old_to_new );
			bool has_edge; ED edge;
			boost::tie( edge, has_edge ) = boost::edge( source, target, graph_ );
			debug_assert( has_edge );
			innervec.push_back( edge );
		}
		rings_and_their_edges_.push_back( innervec );
	}

	// Many of the following will be reset in finalize(), but I deserialize them here anyway because it's just as easy as not doing so

	// RingConformerSets are reset in finalize()
	// EXEMPT conformer_sets_
	arc( mode_ ); // enum core::chemical::TypeSetMode
	arc( nheavyatoms_ ); // Size
	arc( n_hbond_acceptors_ ); // Size
	arc( n_hbond_donors_ ); // Size
	arc( n_backbone_heavyatoms_ ); // Size
	arc( first_sidechain_hydrogen_ ); // Size
	arc( bonded_neighbor_ ); // utility::vector1<AtomIndices>
	arc( bonded_neighbor_type_ ); // utility::vector1<utility::vector1<BondName> >
	arc( attached_H_begin_ ); // utility::vector1<Size>
	arc( attached_H_end_ ); // utility::vector1<Size>
	arc( dihedral_atom_sets_ ); // utility::vector1<dihedral_atom_set>
	arc( dihedrals_for_atom_ ); // utility::vector1<utility::vector1<Size> >
	arc( improper_dihedral_atom_sets_ ); // utility::vector1<dihedral_atom_set>
	arc( improper_dihedrals_for_atom_ ); // utility::vector1<utility::vector1<Size> >
	arc( bondangle_atom_sets_ ); // utility::vector1<bondangle_atom_set>
	arc( bondangles_for_atom_ ); // utility::vector1<utility::vector1<Size> >
	arc( last_controlling_chi_ ); // utility::vector1<Size>
	arc( atoms_last_controlled_by_chi_ ); // utility::vector1<AtomIndices>
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
	arc( metal_binding_atoms_ ); // utility::vector1<std::string>
	arc( disulfide_atom_name_ ); // std::string
	arc( is_proton_chi_ ); // utility::vector1<_Bool>
	arc( proton_chis_ ); // utility::vector1<Size>
	arc( chi_2_proton_chi_ ); // utility::vector1<Size>
	arc( proton_chi_samples_ ); // utility::vector1<utility::vector1<Real> >
	arc( proton_chi_extra_samples_ ); // utility::vector1<utility::vector1<Real> >
	arc( path_distance_ ); // utility::vector1<utility::vector1<int> >
	arc( atom_aliases_ ); // std::map<std::string, std::string>
	arc( orbitals_index_ ); // std::map<std::string, int>
	arc( chi_rotamers_ ); // utility::vector1<utility::vector1<std::pair<Real, Real> > >
	arc( rotamer_library_specification_ ); // rotamers::RotamerLibrarySpecificationOP
	arc( ring_sizes_ ); // utility::vector1<Size>
	arc( lowest_ring_conformer_ ); // utility::vector1<std::string>
	arc( low_ring_conformers_ ); // utility::vector1<utility::vector1<std::string> >
	arc( properties_ ); // ResiduePropertiesOP
	arc( aa_ ); // enum core::chemical::AA
	arc( rotamer_aa_ ); // enum core::chemical::AA
	arc( backbone_aa_ ); // enum core::chemical::AA
	arc( na_analogue_ ); // enum core::chemical::AA
	arc( base_name_ ); // std::string
	arc( base_type_cop_ ); // ResidueTypeCOP
	arc( name_ ); // std::string
	arc( name3_ ); // std::string
	arc( name1_ ); // char
	arc( interchangeability_group_ ); // std::string
	arc( nbr_radius_ ); // Real
	arc( force_nbr_atom_orient_ ); // _Bool
	arc( remap_pdb_atom_names_ ); // _Bool
	arc( mass_ ); // Real
	arc( atom_2_residue_connection_map_ ); // utility::vector1<utility::vector1<Size> >
	arc( atoms_within_one_bond_of_a_residue_connection_ ); // utility::vector1<utility::vector1<two_atom_set> >
	arc( within1bonds_sets_for_atom_ ); // utility::vector1<utility::vector1<std::pair<Size, Size> > >
	arc( atoms_within_two_bonds_of_a_residue_connection_ ); // utility::vector1<utility::vector1<three_atom_set> >
	arc( within2bonds_sets_for_atom_ ); // utility::vector1<utility::vector1<std::pair<Size, Size> > >
	arc( lower_connect_id_ ); // Size
	arc( upper_connect_id_ ); // Size
	arc( n_non_polymeric_residue_connections_ ); // Size
	arc( n_polymeric_residue_connections_ ); // Size
	// RNA_ResidueType and CarbohydrateInfo are reset in finalize()
	// EXEMPT rna_info_ carbohydrate_info_
	arc( atom_base_indices_ ); // utility::vector1<Size>
	arc( abase2_indices_ ); // utility::vector1<Size>
	arc( chi_atoms_indices_ ); // utility::vector1<AtomIndices>
	arc( nu_atoms_indices_ ); // utility::vector1<AtomIndices>
	arc( ring_atoms_indices_ ); // utility::vector1<AtomIndices>
	arc( mainchain_atoms_indices_ ); // AtomIndices
	arc( nbr_atom_indices_ ); // Size
	arc( actcoord_atoms_indices_ ); // AtomIndices
	arc( cut_bond_neighbor_indices_ ); // utility::vector1<AtomIndices>
	arc( atom_shadowed_indices_ ); // utility::vector1<Size>
	arc( rama_prepro_mainchain_torsion_potential_name_ ); // std::string
	arc( rama_prepro_mainchain_torsion_potential_name_beforeproline_ ); // std::string
	arc( rama_prepro_map_file_name_ ); // std::string
	arc( rama_prepro_map_file_name_beforeproline_ ); // std::string
	arc( defined_adducts_ ); // utility::vector1<Adduct>
	arc( nondefault_ ); // _Bool

	// Observers aren't being serialized - any observer on the remote side will have to reattach itself.
	// EXEMPT destruction_obs_hub_
	// EXEMPT destruction_obs_mutex_;

	//finalize(); // Make sure all the derived data is up-to-date
	// EXEMPT finalized_
}
SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ResidueType );
CEREAL_REGISTER_TYPE( core::chemical::ResidueType )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_ResidueType )
#endif // SERIALIZATION
