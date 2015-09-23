// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ResidueType.cc
/// @brief Method definitions for ResidueType
/// @author
/// Phil Bradley
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to properties, lipids, carbohydrates, and other non-AAs)

// Unit headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <boost/graph/graph_utility.hpp>
#include <core/chemical/util.hh>

// Package Headers
#include <core/chemical/ResidueProperties.hh>
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/chemical/residue_support.hh>
#include <core/chemical/icoor_support.hh>
#include <core/chemical/AtomType.hh>
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
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/bond_support.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/NumericTraits.hh>

// Basic headers
#include <basic/Tracer.hh>
// Options and Option key includes (needed for protonated versions of the residues - pH mode)
#include <basic/options/option.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>

// Utility headers
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>
#include <utility/graph/ring_detection.hh>

// External headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <boost/graph/graph_utility.hpp>

// C++ headers
#include <algorithm>


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

namespace core {
namespace chemical {

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

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
	orbital_types_( orbital_types),
	conformer_set_(/* NULL */),
	residue_type_set_( /* 0 */ ),
	graph_(),
	orbitals_(),
	nheavyatoms_(0),
	n_hbond_acceptors_(0),
	n_hbond_donors_(0),
	n_backbone_heavyatoms_(0),
	first_sidechain_hydrogen_( 0 ),
	rotamer_library_specification_( 0 ),
	lowest_ring_conformer_( "" ),
	low_ring_conformers_(),
	properties_( ResiduePropertiesOP( new ResidueProperties( this ) ) ),
	aa_( aa_unk ),
	rotamer_aa_( aa_unk ),
	backbone_aa_( aa_unk ),
	name_(),
	name3_(),
	name1_(),
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
	finalized_(false),
	nondefault_(false),
	base_restype_name_(""),
	serialized_(false)
{}

ResidueType::~ResidueType()
{
	tr.Trace << "Residue dstor" << std::endl;
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
	atom_types_ = residue_type.atom_types_;
	elements_ = residue_type.elements_;
	mm_atom_types_ = residue_type.mm_atom_types_;
	gasteiger_atom_types_ = residue_type.gasteiger_atom_types_;
	orbital_types_ = residue_type.orbital_types_;
	conformer_set_ = residue_type.conformer_set_;
	residue_type_set_ = residue_type.residue_type_set_;
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
	ring_atoms_ = residue_type.ring_atoms_;
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
	rna_residue_type_ = residue_type.rna_residue_type_;
	// CarbohydrateInfo has a back-pointer to ResidueType and must be reset during finalize
	carbohydrate_info_ = 0; /* NULL */
	//rings_and_their_edges_=residue_type.rings_and_their_edges_; //Apparently updated below
	atom_base_indices_ = residue_type.atom_base_indices_;
	abase2_indices_ = residue_type.abase2_indices_;
	chi_atoms_indices_ = residue_type.chi_atoms_indices_;
	nu_atoms_indices_ = residue_type.nu_atoms_indices_;
	mainchain_atoms_indices_ = residue_type.mainchain_atoms_indices_;
	nbr_atom_indices_ = residue_type.nbr_atom_indices_;
	actcoord_atoms_indices_ = residue_type.actcoord_atoms_indices_;
	cut_bond_neighbor_indices_ = residue_type.cut_bond_neighbor_indices_;
	atom_shadowed_indices_ = residue_type.atom_shadowed_indices_;
	finalized_ = residue_type.finalized_;
	defined_adducts_ = residue_type.defined_adducts_;
	nondefault_ = residue_type.nondefault_;
	base_restype_name_ = residue_type.base_restype_name_;
	serialized_ = residue_type.serialized_;

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
		VIter v_iter= vp.first;
		VIter old_v_iter= old_vp.first;
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
	for ( std::map<std::string,std::string>::iterator iter( atom_aliases_.begin() ), iter_end( atom_aliases_.end() );
			iter != iter_end; ++iter ) {
		debug_assert( !atom_name_to_vd_.count( iter->first ) );
		debug_assert( atom_name_to_vd_.count( iter->second ) );
		VD vd( atom_name_to_vd_[ iter->second ] );
		atom_name_to_vd_.insert( NameVDPair( iter->first, vd) );
	}
	// Setup the temporary ordered_atoms_ vector for refactor
	// VKM, 5 Sept 2015: I'm not sure the code below is working as desired.  Valgrind issues might be caused here.
	VDs::const_iterator begin = residue_type.ordered_atoms_.begin();
	VDs::const_iterator const end = residue_type.ordered_atoms_.end();
	for ( ; begin != end; ++begin ) {
		VD old_vd = *begin;
		VD vd = old_to_new[old_vd];
		ordered_atoms_.push_back(vd);
		vd_to_index_[vd] = ordered_atoms_.size();
	}
	std::map<VD, VD> old_atom_base(atom_base_);
	atom_base_.clear();
	for ( std::map<VD, VD>::iterator it = old_atom_base.begin(); it != old_atom_base.end(); ++it ) {
		VD old_key = it->first;
		VD old_value = it->second;
		VD new_key = old_to_new[old_key];
		VD new_value = old_to_new[old_value];
		atom_base_[new_key] = new_value;
	}

	std::map<VD, VD> old_atom_shadowed(atom_shadowed_);
	atom_shadowed_.clear();
	for ( std::map<VD, VD>::iterator it = old_atom_shadowed.begin(); it != old_atom_shadowed.end(); ++it ) {
		VD old_key = it->first;
		VD old_value = it->second;
		VD new_key = old_to_new[old_key];
		VD new_value = old_to_new[old_value];
		atom_shadowed_[new_key] = new_value;
	}

	utility::vector1<utility::vector1<VD> >  old_chi_atoms(chi_atoms_);
	chi_atoms_.clear();
	//chi_atoms_.resize(old_chi_atoms.size());
	for ( utility::vector1<utility::vector1<VD> >::const_iterator it= old_chi_atoms.begin(); it != old_chi_atoms.end(); ++it ) {
		utility::vector1<VD> old_vector = *it;
		debug_assert(old_vector.size() == 4);
		utility::vector1<VD> new_vector;
		for ( Size i= 1; i<= old_vector.size(); ++i ) {
			new_vector.push_back(old_to_new[old_vector[i]]);
		}
		chi_atoms_.push_back(new_vector);
	}

	utility::vector1<utility::vector1<VD> >  old_nu_atoms(nu_atoms_);
	nu_atoms_.clear();
	//chi_atoms_.resize(old_chi_atoms.size());
	for ( utility::vector1<utility::vector1<VD> >::const_iterator it= old_nu_atoms.begin(); it != old_nu_atoms.end(); ++it ) {
		utility::vector1<VD> old_vector = *it;
		debug_assert(old_vector.size() == 4);
		utility::vector1<VD> new_vector;
		for ( Size i= 1; i<= old_vector.size(); ++i ) {
			new_vector.push_back(old_to_new[old_vector[i]]);
		}
		nu_atoms_.push_back(new_vector);
	}

	utility::vector1<VD> old_mainchain(mainchain_atoms_);
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
	for ( std::map<VD, utility::vector1<VD> >::const_iterator it = old_cut_bonds.begin(); it != old_cut_bonds.end(); ++it ) {
		VD old_key = it->first;
		utility::vector1<VD> old_value = it->second;
		VD new_key = old_to_new[old_key];
		utility::vector1<VD> new_value;
		for ( Size i=1; i<= old_value.size(); ++i ) {
			new_value.push_back(old_to_new[ old_value[i] ] );
		}
		cut_bond_neighbor_[ new_key ] = new_value;
	}

	std::map< VD, AtomICoor > old_icoor(icoor_);
	icoor_.clear();
	for ( std::map<VD, AtomICoor >::const_iterator it= old_icoor.begin(); it != old_icoor.end(); ++it ) {
		VD old_key = it->first;
		VD new_key = old_to_new[old_key];
		AtomICoor old_icoor = it->second; //now we have to change the vertex descriptors within icoor. They are pointing to an old vertex descriptor
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


	for ( Size i=1; i<= residue_connections_.size(); ++i ) {
		VD old_vertex = residue_connections_[i].vertex();
		residue_connections_[i].vertex(old_to_new[old_vertex]);
		residue_connections_[i].atomno(vd_to_index_.find(residue_connections_[i].vertex())->second);

		AtomICoor new_icoor = residue_connections_[i].icoor();
		for ( Size j = 1; j <= 3; ++j ) {
			VD old_vd = new_icoor.stub_atom( j ).vertex();
			new_icoor.stub_atom( j ).vertex( old_to_new[old_vd] );
			new_icoor.stub_atom( j ).atomno( vd_to_index_.find(new_icoor.stub_atom(j).vertex())->second );
		}
		residue_connections_[ i ].icoor( new_icoor );

	}

	for ( Size index=1; index<= natoms(); ++index ) {
		utility::vector1< core::Size > const orbs(graph_[ordered_atoms_[index]].bonded_orbitals());
		for ( utility::vector1< core::Size >::const_iterator orb = orbs.begin(); orb != orbs.end(); ++orb ) {
			orbitals_[*orb].new_icoor().vertex1( old_to_new[ orbitals_[*orb].new_icoor().vertex1()  ] );
			orbitals_[*orb].new_icoor().vertex2(old_to_new[ orbitals_[*orb].new_icoor().vertex2()  ] );
			orbitals_[*orb].new_icoor().vertex3( old_to_new[ orbitals_[*orb].new_icoor().vertex3()  ] );
		}
	}

	utility::vector1<VD> old_bb(force_bb_);
	force_bb_.clear();
	for ( Size i=1; i<= old_bb.size(); ++i ) {
		force_bb_.push_back(old_to_new[ old_bb[i] ]);
	}

	return *this;
}


ResidueTypeSet const &
ResidueType::residue_type_set() const
{
	ResidueTypeSetCOP residue_type_set = residue_type_set_.lock();
	if ( !residue_type_set ) {
		utility_exit_with_message( "ResidueType::residue_type_set: pointer is not set!");
	}
	return *residue_type_set;
}

bool
ResidueType::in_residue_type_set() const
{
	return ! residue_type_set_.expired();
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
	ResidueTypeOP rsd( new ResidueType( 0, 0, 0, 0 ) );
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

void
ResidueType::residue_type_set( ResidueTypeSetCAP set_in )
{
	residue_type_set_ = set_in;
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
	NameVDMap::const_iterator found = atom_name_to_vd_.find( atom_name );
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

Bond & ResidueType::bond(ED const ed){
	return graph_[ ed ];
}
Bond const & ResidueType::bond(ED const ed) const{
	return graph_[ ed ];
}

Bond &
ResidueType::bond(std::string const & atom1, std::string const & atom2) {
	VD vd1( atom_vertex( atom1 ) ), vd2( atom_vertex( atom2 ) );
	ED ed;
	bool found;
	boost::tie( ed, found ) = boost::edge( vd1, vd2 , graph_ );
	if ( ! found ) {
		utility_exit_with_message( "Cannot find bond between " + atom1 + " and " + atom2 + " in residue " + name() );
	}
	return graph_[ ed ];
}

Bond const &
ResidueType::bond(std::string const & atom1, std::string const & atom2) const {
	VD vd1( atom_vertex( atom1 ) ), vd2( atom_vertex( atom2 ) );
	ED ed;
	bool found;
	boost::tie( ed, found ) = boost::edge( vd1, vd2 , graph_ );
	if ( ! found ) {
		utility_exit_with_message( "Cannot find bond between " + atom1 + " and " + atom2 + " in residue " + name() );
	}
	return graph_[ ed ];
}


/// @brief set the atom which connects to the lower connection
void
ResidueType::set_lower_connect_atom( std::string const & atm_name )
{
	finalized_ = false;
	if ( atm_name == "NONE" ) {
		if ( lower_connect_id_ != 0 ) {
			tr.Debug << "ERASING LOWER_CONNECT: " << lower_connect_id_ << " lcid: " << upper_connect_id_ << std::endl;
			utility::vector1< ResidueConnection >::iterator to_erase( residue_connections_.begin() );
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

/// @brief set the atom which connects to the upper connection
void
ResidueType::set_upper_connect_atom( std::string const & atm_name )
{
	finalized_ = false;
	if ( atm_name == "NONE" ) {
		if ( upper_connect_id_ != 0 ) {
			tr.Debug << "ERASING UPPER_CONNECT: " << upper_connect_id_ << " lcid: " << lower_connect_id_  << std::endl;
			utility::vector1< ResidueConnection >::iterator to_erase( residue_connections_.begin() );
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

ResidueConnection const &
ResidueType::upper_connect() const
{
	debug_assert( properties_->has_property( POLYMER ) );
	debug_assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ];
}

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

/// @brief index number of the atom which connects to the upper connection
Size
ResidueType::upper_connect_atom() const
{
	debug_assert( properties_->has_property( POLYMER ) );
	debug_assert( upper_connect_id_ != 0 );
	return vd_to_index_.find(residue_connections_[ upper_connect_id_ ].vertex())->second;
}

/// @brief number of ResidueConnections, counting polymeric residue connections
Size
ResidueType::n_residue_connections() const
{
	return residue_connections_.size();
}

Size
ResidueType::residue_connect_atom_index( Size const resconn_id ) const {
	return vd_to_index_.find(residue_connections_[ resconn_id ].vertex())->second;
}


/// @brief get a ResidueConection
ResidueConnection const &
ResidueType::residue_connection( Size const i ) const
{
	return residue_connections_[i];
}

ResidueConnection &
ResidueType::residue_connection( Size const i )
{
	return residue_connections_[ i ];
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
		assert(atom_name_to_vd_.find(atom_name) == atom_name_to_vd_.end());
		assert(atom_name_to_vd_.find( strip_whitespace(atom_name)) == atom_name_to_vd_.end());
		atom_name_to_vd_[ atom_name ] = v;
		atom_name_to_vd_[ strip_whitespace( atom_name ) ] = v;
	}

	ordered_atoms_.push_back(v);

	// allocate space for the new atom !!!!!!!!!!!!!!!!!!!!!!
	// eg, in the atom/resconn map
	assert( atom_2_residue_connection_map_.size() == ordered_atoms_.size()-1 );

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
}

/// @brief flag an atom for deletion by adding its index to the delete_atom_ list
void
ResidueType::delete_atom( Size const index )
{
	finalized_ = false;
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
		tr.Warning << "WARNING Elements set undefined." << std::endl;
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
		gasteiger_type = 0;
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
	NameVDMap::const_iterator atom_name_to_vd_iter( atom_name_to_vd_.find( name ) );
	if ( atom_name_to_vd_iter == atom_name_to_vd_.end() ) {
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		debug_assert(false);
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
	for ( VIter iter(allverts.first); iter != allverts.second; ++iter ) {
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

// Return a pointer to the object containing the set of ring conformers possible for this saccharide.
core::chemical::rings::RingConformerSetCOP
ResidueType::ring_conformer_set() const
{
	return conformer_set_;
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

/// @note this does not set xyz coordiates for the added orbital but sets the index of the orbital and maps
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
	std::string const atom_name
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
		std::string message = "add_bond: atoms " + atom_name1 + " and " + atom_name2 + " dont exist!";
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
ResidueType::add_nu(core::uint const nu_index,
	std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4)
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	if ( !has(atom_name1) || !has(atom_name2) || !has(atom_name3) || !has(atom_name4) ) {
		utility_exit_with_message("ResidueType::add_nu: Requested atoms don't exist in this ResidueType!");
	}

	utility::vector1<VD> atoms;
	atoms.push_back(ordered_atoms_[atom_index(atom_name1)]);
	atoms.push_back(ordered_atoms_[atom_index(atom_name2)]);
	atoms.push_back(ordered_atoms_[atom_index(atom_name3)]);
	atoms.push_back(ordered_atoms_[atom_index(atom_name4)]);

	if ( nu_atoms_.size() < nu_index ) {
		nu_atoms_.resize(nu_index);
	}
	nu_atoms_[nu_index] = atoms;
}


// Set this cyclic residue's lowest-energy ring conformer by IUPAC name.
void
ResidueType::set_lowest_energy_ring_conformer( std::string const & conformer )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	lowest_ring_conformer_ = conformer;
}

// Set this cyclic residue's low-energy ring conformers by IUPAC name.
void
ResidueType::set_low_energy_ring_conformers( utility::vector1< std::string > const & conformers )
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	low_ring_conformers_ = conformers;
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
	for ( utility::vector1<VDs>::const_iterator iter( found_chis.begin() ); iter != found_chis.end(); ++iter ) {
		VDs const & chi( *iter );
		debug_assert( chi.size() == 4 );
		add_chi( chi[1], chi[2], chi[3], chi[4] );
		if ( atom( chi[4] ).element_type()->element() == core::chemical::element::H ) {
			// proton chi
			proton_chis.push_back( nchi() );
		}
	} // for all found chis

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
		tr.Warning << "Warning: Number of base proton chi samples (" << num_H_confs << ") for " << name() << " exceeds requested number of samples" << std::endl;
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

bool
ResidueType::is_adduct() const
{
	return properties_->has_property( ADDUCT );
}


core::Real
ResidueType::get_numeric_property(std::string const & tag) const
{
	std::map<std::string, core::Real> const numeric_properties( properties_->numeric_properties() );
	std::map<std::string, core::Real>::const_iterator property_it( numeric_properties.find( tag ) );
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
	std::map<std::string, std::string>::const_iterator property_it(string_properties.find(tag));
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


/// @details "Custom" VariantTypes as strings are permitted for the enzdes and metalloproteins cases.
/// Do not enable unless you have a good reason to, as string look-ups are less efficient and more error-prone.
void
ResidueType::enable_custom_variant_types()
{
	properties_->enable_custom_variant_types();
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
		tr.Warning << "WARNING: ACT_COORD_ATOM specified for non-protein residue type '" << name() << "' . This doesn't make much sense." << std::endl;
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
	for ( std::map<std::string,std::string>::iterator iter( atom_aliases_.begin() ), iter_end( atom_aliases_.end() );
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
			atoms.push_back( vd_to_index_.find( chi_atoms_[chino][atom_index] )->second) ;
		}
		chi_atoms_indices_.push_back(atoms);
		atoms.clear();
	}
	nu_atoms_indices_.clear();
	atoms.clear();
	for ( Size nu_no=1; nu_no <= nu_atoms_.size(); ++nu_no ) {
		for ( Size atom_index=1; atom_index <= nu_atoms_[nu_no].size(); ++atom_index ) {
			atoms.push_back( vd_to_index_.find( nu_atoms_[nu_no][atom_index] )->second) ;
		}
		nu_atoms_indices_.push_back(atoms);
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
		for ( utility::vector1< core::Size >::const_iterator orb = orbs.begin(); orb != orbs.end(); ++orb ) {
			orbitals_[*orb].new_icoor().replace_stub1( vd_to_index_[ordered_atoms_[orbitals_[*orb].new_icoor().get_stub1()]] );
			orbitals_[*orb].new_icoor().replace_stub2( vd_to_index_[ordered_atoms_[orbitals_[*orb].new_icoor().get_stub2()]] );
			orbitals_[*orb].new_icoor().replace_stub3( vd_to_index_[ordered_atoms_[orbitals_[*orb].new_icoor().get_stub3()]] );
		}
	}


	for ( Size index=1; index<= natoms(); ++index ) {
		for ( Size i=1; i<= 3; ++i ) {
			ICoorAtomID & stub_atom( icoor_[ ordered_atoms_[index] ].stub_atom( i )   );
			if ( stub_atom.type() == ICoorAtomID::INTERNAL ) {
				stub_atom.atomno(   vd_to_index_.find(stub_atom.vertex())->second ); //somewhat of a problem. if vertex doesnt exist the map constructor will create a value
				if( stub_atom.atomno() == 0 ) { // this will trigger if we deleted a stub atom for some other atom
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
* ring_atoms_
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

* rna_residue_type_  -- Will be reset based on other ResidueType data
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
	ring_atoms_.clear();

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

	// Set the ring atoms.
	// The logic here is tricky. The nu torsion definitions contain all the ring atoms, of course.
	// However, the first nu definition will include a virtual atom as its first atom, if defined properly.
	// The last will include a virtual atom as its final atom and will not contain any atoms already included in
	// earlier definitions.
	// Hence, we can take the last three atoms indices from the first nu definition, and then the last index from the
	// rest of the definitions except for the last, which we can completely ignore.
	if ( properties_->has_property( CYCLIC ) ) {
		Size const n_nus( nu_atoms_indices_.size() );
		if ( ! graph_[ ordered_atoms_[ nu_atoms_indices_[ 1 ][ 1 ] ] ].is_virtual() ||
				! graph_[ ordered_atoms_[ nu_atoms_indices_[ n_nus ][ 4 ] ] ].is_virtual() ) {
			utility_exit_with_message( "The nu angles for this ResidueType are not properly defined.  "
				"The first atom of the first nu and the last atom of the last nu must be virtual atoms." );
		}
		for ( uint j( 2 ); j <= 4; ++j ) {
			ring_atoms_.push_back( nu_atoms_indices_[ 1 ][ j ] );
		}
		for ( uint i( 2 ); i < n_nus; ++i ) {
			ring_atoms_.push_back( nu_atoms_indices_[ i ][ 4 ] );
		}
		// You always need 1 fewer nu angles to define a ring than the number of atoms in that ring.
		debug_assert( ring_atoms_.size() == nu_atoms_indices_.size() + 1 );
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
				for ( utility::vector1< Size >::iterator terminal_atom1 = ca1d1.begin();
						terminal_atom1 != ca1d1.end(); ++terminal_atom1 ) {
					for ( utility::vector1< Size >::iterator terminal_atom2 = ca2d1.begin();
							terminal_atom2 != ca2d1.end(); ++terminal_atom2 ) {
						dihedral_atom_set temp( *terminal_atom1, central_atom1, central_atom2, *terminal_atom2 );
						dihedral_atom_sets_.push_back( temp );
						Size const which_dihedral = dihedral_atom_sets_.size();
						dihedrals_for_atom_[ *terminal_atom1 ].push_back( which_dihedral );
						dihedrals_for_atom_[   central_atom1 ].push_back( which_dihedral );
						dihedrals_for_atom_[   central_atom2 ].push_back( which_dihedral );
						dihedrals_for_atom_[ *terminal_atom2 ].push_back( which_dihedral );
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
				for ( utility::vector1< Size >::iterator terminal_atom1 = ca1d1.begin();
						terminal_atom1 != ca1d1.end(); ++terminal_atom1 ) {
					for ( utility::vector1< Size >::iterator terminal_atom2 = terminal_atom1+1;
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
				for ( utility::vector1< Size >::iterator terminal_atom1 = ca2d1.begin();
						terminal_atom1 != ca2d1.end(); ++terminal_atom1 ) {
					for ( utility::vector1< Size >::iterator terminal_atom2 = terminal_atom1+1;
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

	// Assign a set of possible ring conformations.
	// Ring size is determined by the number of NU angles listed in the .params file, which should always be 1 less
	// than the size of the ring.
	if ( properties_->has_property( CYCLIC ) ) {
		conformer_set_ = rings::RingConformerSetOP( new rings::RingConformerSet(
			ring_atoms_.size(), lowest_ring_conformer_, low_ring_conformers_ ) );
	}

	if ( properties_->has_property( RNA ) ) { //reinitialize and RNA derived data.
		//Reinitialize rna_residue_type_ object! This also make sure rna_residue_type_ didn't inherit anything from the previous update!
		//It appears that the rna_residue_type_ is shared across multiple ResidueType object, if the rna_residue_type_ is not reinitialized here!
		rna_residue_type_ = core::chemical::rna::RNA_ResidueTypeOP( new core::chemical::rna::RNA_ResidueType );
		//update_last_controlling_chi is treated separately for RNA case. Parin Sripakdeevong, June 26, 2011
		rna_residue_type_->rna_update_last_controlling_chi( get_self_weak_ptr(), last_controlling_chi_, atoms_last_controlled_by_chi_);
		rna_residue_type_->update_derived_rna_data( get_self_weak_ptr() );
	} else if ( properties_->has_property( CARBOHYDRATE ) ) {
		carbohydrate_info_ =
			carbohydrates::CarbohydrateInfoOP( new carbohydrates::CarbohydrateInfo( get_self_weak_ptr() ) );
		update_last_controlling_chi();
	} else {
		update_last_controlling_chi();
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
	msg << "One or more internal errors have occurred in residue type setup:" << std::endl;

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
chi_atoms_             v1<v1<uint>>         add_chi
nu_atoms_              v1<v1<uint>>         add_nu
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

	NameVDMap::const_iterator graph_iter( atom_name_to_vd_.find( name ) );
	if ( graph_iter == atom_name_to_vd_.end() ) {
#if defined BOINC
		// chu temporary graphic fix for boinc
		if ( name == "CA" && !is_protein() ) return 1;
#endif
		if ( name == "CA" && is_membrane() ) return 2;
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		// AMW: this really should be the residue's name, so that people trying to debug can see what the issue is
		utility_exit_with_message("unknown atom_name: " + this->name() + "  " + name );
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
		// AMW: this really should be the residue's name, so that people trying to debug can see what the issue is
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

	std::map< std::string, int >::const_iterator iter
		( orbitals_index_.find( name ) );
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
	std::map< VD, AtomICoor >::const_iterator itr( icoor_.find(ordered_atoms_[atm]) );
	debug_assert( itr != icoor_.end() );
	return itr->second;
}

/// @brief AtomICoord of an atom
AtomICoor const &
ResidueType::icoor( VD const atm ) const
{
	debug_assert( has(atm) );
	std::map< VD, AtomICoor >::const_iterator itr( icoor_.find(atm) );
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

	Size atomno;
	VD atom_vd( id.vertex() );
	switch ( id.type() ) {
	case ICoorAtomID::INTERNAL :
		debug_assert( atom_vd == atom_vertex( atm ) );
		if ( atm == stub_atom1 ) {
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
			if ( atm == stub_atom1 ) {
				//root of tree
				if ( natoms() == 1 ) {
					set_atom_base( atm, atm );
				} else {
					set_atom_base( atm, stub_atom2 );
				}
			} else {
				set_atom_base( atm, stub_atom1 );
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
		atomno = id.atomno(); // For CONNECT, the atomno is repurposed as the connection number
		residue_connections_[ atomno ].icoor( ic );
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
	if ( atm == stub_atom1 ) {
		//Root atom case
		if ( root_atom_ != ResidueType::null_vertex ) {
			tr.Error << "Can't reset root. Was " << atom_name( root_atom_ ) << " Resetting to " << atom_name( atm ) << std::endl;
			utility_exit_with_message( "Attempted to inappropriately reset ICOOR root atom." );
		}
		root_atom_ = atm;
		// Now we can continue as normal.
	}
	icoor_[ atm ] = ic;
	// update atom_base if it isn't set yet (or is set to itself)
	if ( ! atom_base_.count(atm) || atom_base_[atm] == atm ) {
		if ( atm == stub_atom1 ) {
			//root of tree
			if ( natoms() == 1 ) {
				set_atom_base( atm, atm );
			} else {
				set_atom_base( atm, stub_atom2 );
			}
		} else {
			set_atom_base( atm, stub_atom1 );
		}
	}
	if ( update_xyz ) {
		set_ideal_xyz( atm, ic.build( *this ) );
	}
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
	if ( n_residue_connections() != 0 ) {
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


void
ResidueType::show_all_atom_names( std::ostream & out ) const {

	for ( VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first ) {
		VIter v_iter= vp.first;
		VD vd = *v_iter;
		Atom a = graph_[vd];
		out << a.name() << " " << &graph_[vd] << std::endl;
	}

}

/// @brief  Check if atom is virtual.
bool
ResidueType::is_virtual( Size const & atomno ) const
{
	return ( atom_type( atomno ).is_virtual() );
}


///////////////////////////////////////////////////////////////
core::chemical::rna::RNA_ResidueType const &
ResidueType::RNA_type() const{
	return ( *rna_residue_type_ );
}

/// @author Labonte <JWLabonte@jhu.edu>
void
ResidueType::show( std::ostream & output, bool output_atomic_details ) const
{
	using namespace std;
	using namespace utility;

	output << name_ << " (" << name3_ << ", " << name1_ << "):" << endl;

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

	output << " Ring atoms:  ";
	Size const n_ring_atoms( ring_atoms_.size() );
	for ( uint i = 1; i <= n_ring_atoms; ++i ) {
		output << ' ' << atom_name( ring_atoms_[ i ] );
	}
	output << endl;

	output << " Side-chain atoms:";
	Size const n_sc_atoms( all_sc_atoms_.size() );
	for ( uint i = 1; i <= n_sc_atoms; ++i ) {
		output << ' ' << atom_name( all_sc_atoms_[ i ] );
	}
	output << endl;

	if ( properties_->has_property( CARBOHYDRATE ) ) {
		carbohydrate_info_->show( output );
	}

	if ( output_atomic_details ) {
		output << " Atomic Details:" << endl;
		Size n_atoms = natoms();
		for ( core::uint i = 1; i <= n_atoms; ++i ) {
			output << "  Atom " << i << ": ";
			atom( i ).show( output );
		}
	}
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
