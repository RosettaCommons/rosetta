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
/// Jason W. Labonte (code related to lipids, carbohydrates, and other non-AAs)

// Unit headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>

// Package Headers
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/chemical/ResidueSupport.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/RingConformerSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/NumericTraits.hh>

// Basic headers
#include <basic/Tracer.hh>
// Options and Option key includes (needed for protonated versions of the residues - pH mode)
#include <basic/options/option.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>

// Utility headers
#include <utility/PyAssert.hh>
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

static basic::Tracer tr("core.chemical.ResidueType");

// must be a better place for this, probably already exists!
inline
std::string
strip_whitespace( std::string const & name )
{
	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to dothis?
	return trimmed_name;
}

ResidueType::ResidueType(
		AtomTypeSetCAP atom_types,
		ElementSetCAP elements,
		MMAtomTypeSetCAP mm_atom_types,
		orbitals::OrbitalTypeSetCAP orbital_types//, CSDAtomTypeSetCAP csd_atom_types kwk commenting out csd atom types until they are fully functional
) : utility::pointer::ReferenceCount(),
		atom_types_( atom_types ),
		elements_( elements ),
		mm_atom_types_( mm_atom_types ),
		orbital_types_( orbital_types),
		conformer_set_(NULL),
		residue_type_set_( 0 ),
		graph_(),
		orbitals_(),
		nheavyatoms_(0),
		n_hbond_acceptors_(0),
		n_hbond_donors_(0),
		n_orbitals_(0),
		n_backbone_heavyatoms_(0),
		first_sidechain_hydrogen_( 0 ),
		ndihe_( 0 ),
		nbonds_(0),
		rotamer_library_name_( "" ),
		use_ncaa_rotlib_( false ),
		ncaa_rotlib_n_rots_( 0 ),
		is_polymer_( false ),
		is_protein_( false ),
		is_alpha_aa_(false),
		is_beta_aa_(false),
		is_l_aa_(false),
		is_d_aa_(false),
		is_charged_( false ),
		is_polar_( false ),
		has_sc_orbitals_(false),
		is_aromatic_( false ),
		is_cyclic_( false ),
		is_DNA_( false ),
		is_RNA_( false ),
		is_NA_( false ),
		is_carbohydrate_( false ),
		is_lipid_( false ),
		is_ligand_( false ),
		is_metal_( false ), //Is this residue type a metal ion?
		is_metalbinding_( false ), //Is this residue type a type that has the potential to bind to metal ions?
		is_surface_( false ),
		is_terminus_( false ),
		is_lower_terminus_( false ),
		is_upper_terminus_( false ),
		is_branch_lower_terminus_( false ),
		is_phosphonate_( false ),
		is_phosphonate_upper_( false ),
		is_acetylated_nterminus_( false ),
		is_methylated_cterminus_( false ),
		is_coarse_( false ), //currently for coarse_RNA only
		is_adduct_( false ),
		aa_( aa_unk ),
		rotamer_aa_( aa_unk ),
		backbone_aa_( aa_unk ),
		name_(),
		name3_(),
		name1_(),
		interchangeability_group_(),
		nbr_atom_( ResidueGraph::null_vertex() ),
		nbr_radius_( 0 ),
		force_nbr_atom_orient_(false),
		mass_(0),
		n_actcoord_atoms_( 0 ),
		lower_connect_id_( 0 ),
		upper_connect_id_( 0 ),
		n_non_polymeric_residue_connections_( 0 ),
		n_polymeric_residue_connections_( 0 ),
		carbohydrate_info_(NULL),
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

ResidueType::ResidueType(ResidueType const & residue_type):
		utility::pointer::ReferenceCount(),
		atom_types_( residue_type.atom_types_ ),
		elements_( residue_type.elements_ ),
		mm_atom_types_( residue_type.mm_atom_types_ ),
		orbital_types_( residue_type.orbital_types_ ),
		conformer_set_(residue_type.conformer_set_),
		residue_type_set_( residue_type.residue_type_set_ ),
		graph_(residue_type.graph_),
		vd_to_index_(),
		atom_base_(residue_type.atom_base_),
		abase2_(residue_type.abase2_),
		orbitals_(residue_type.orbitals_),
		nheavyatoms_(residue_type.nheavyatoms_),
		n_hbond_acceptors_(residue_type.n_hbond_acceptors_),
		n_hbond_donors_(residue_type.n_hbond_donors_),
		n_orbitals_(residue_type.n_orbitals_),
		n_backbone_heavyatoms_(residue_type.n_backbone_heavyatoms_),
		first_sidechain_hydrogen_( residue_type.first_sidechain_hydrogen_ ),
		ndihe_( residue_type.ndihe_ ),
		nbonds_(residue_type.nbonds_),
		//orbital_bonded_neighbor_(residue_type.orbital_bonded_neighbor_),
		bonded_neighbor_(residue_type.bonded_neighbor_),
		bonded_neighbor_type_(residue_type.bonded_neighbor_type_),
		cut_bond_neighbor_(residue_type.cut_bond_neighbor_),
		attached_H_begin_(residue_type.attached_H_begin_),
		attached_H_end_(residue_type.attached_H_end_),
		parents_(residue_type.parents_),
		icoor_(residue_type.icoor_),
		dihedral_atom_sets_(residue_type.dihedral_atom_sets_),
		dihedrals_for_atom_(residue_type.dihedrals_for_atom_),
		bondangle_atom_sets_(residue_type.bondangle_atom_sets_),
		bondangles_for_atom_(residue_type.bondangles_for_atom_),
		atom_shadowed_(residue_type.atom_shadowed_),
		last_controlling_chi_(residue_type.last_controlling_chi_),
		atoms_last_controlled_by_chi_(residue_type.atoms_last_controlled_by_chi_),
		atoms_with_orb_index_(residue_type.atoms_with_orb_index_),
		Haro_index_(residue_type.Haro_index_),
		Hpol_index_(residue_type.Hpol_index_),
		accpt_pos_(residue_type.accpt_pos_),
		Hpos_polar_(residue_type.Hpos_polar_),
		Hpos_apolar_(residue_type.Hpos_apolar_),
		accpt_pos_sc_(residue_type.accpt_pos_sc_),
		Hpos_polar_sc_(residue_type.Hpos_polar_sc_),
		all_bb_atoms_(residue_type.all_bb_atoms_),
		all_sc_atoms_(residue_type.all_sc_atoms_),
		metal_binding_atoms_(residue_type.metal_binding_atoms_),
		mainchain_atoms_(residue_type.mainchain_atoms_),
		actcoord_atoms_(residue_type.actcoord_atoms_),
		chi_atoms_(residue_type.chi_atoms_),
		is_proton_chi_(residue_type.is_proton_chi_),
		proton_chis_(residue_type.proton_chis_),
		chi_2_proton_chi_(residue_type.chi_2_proton_chi_),
		proton_chi_samples_(residue_type.proton_chi_samples_),
		proton_chi_extra_samples_(residue_type.proton_chi_extra_samples_),
		nu_atoms_(residue_type.nu_atoms_),
		path_distance_(residue_type.path_distance_),
		atom_name_to_vd_(), // This must be regenerated below to hold the new new vertex_descriptors
		ordered_atoms_(), // This must be regenerated to hold the new vertex_descriptors
		orbitals_index_(residue_type.orbitals_index_),
		chi_rotamers_(residue_type.chi_rotamers_),
		rotamer_library_name_( residue_type.rotamer_library_name_ ),
		use_ncaa_rotlib_( residue_type.use_ncaa_rotlib_ ),
		ncaa_rotlib_path_( residue_type.ncaa_rotlib_path_),
		ncaa_rotlib_n_rots_( residue_type.ncaa_rotlib_n_rots_ ),
		ncaa_rotlib_n_bins_per_rot_(residue_type.ncaa_rotlib_n_bins_per_rot_),
		properties_(residue_type.properties_),
		is_polymer_( residue_type.is_polymer_ ),
		is_protein_( residue_type.is_protein_ ),
		is_alpha_aa_( residue_type.is_alpha_aa_),
		is_beta_aa_( residue_type.is_beta_aa_),
		is_l_aa_( residue_type.is_l_aa_),
		is_d_aa_( residue_type.is_d_aa_),
		is_charged_( residue_type.is_charged_ ),
		is_polar_( residue_type.is_polar_ ),
		has_sc_orbitals_(residue_type.has_sc_orbitals_),
		is_aromatic_( residue_type.is_aromatic_ ),
		is_cyclic_( residue_type.is_cyclic_ ),
		is_DNA_( residue_type.is_DNA_ ),
		is_RNA_( residue_type.is_RNA_ ),
		is_NA_( residue_type.is_NA_ ),
		is_carbohydrate_( residue_type.is_carbohydrate_ ),
		is_lipid_( residue_type.is_lipid_ ),
		is_ligand_( residue_type.is_ligand_ ),
		is_metal_( residue_type.is_metal_ ), //Is this residue type a metal ion?
		is_metalbinding_( residue_type.is_metalbinding_ ), //Is this residue type capable of binding to a metal ion?
		is_surface_( residue_type.is_surface_ ),
		is_terminus_( residue_type.is_terminus_ ),
		is_lower_terminus_( residue_type.is_lower_terminus_ ),
		is_upper_terminus_( residue_type.is_upper_terminus_ ),
		is_branch_lower_terminus_( residue_type.is_branch_lower_terminus_ ),
		is_acetylated_nterminus_( residue_type.is_acetylated_nterminus_ ),
		is_methylated_cterminus_( residue_type.is_methylated_cterminus_ ),
		is_coarse_( residue_type.is_coarse_ ), //currently for coarse_RNA only
		is_adduct_( residue_type.is_adduct_ ),
		variant_types_( residue_type.variant_types_ ),
		numeric_properties_(residue_type.numeric_properties_),
		string_properties_(residue_type.string_properties_),
		aa_( residue_type.aa_ ),
		rotamer_aa_( residue_type.rotamer_aa_ ),
		backbone_aa_( residue_type.backbone_aa_ ),
		name_( residue_type.name_),
		name3_( residue_type.name3_),
		name1_(residue_type.name1_),
		interchangeability_group_( residue_type.interchangeability_group_ ),
		nbr_atom_(residue_type.nbr_atom_),
		nbr_radius_( residue_type.nbr_radius_ ),
		force_nbr_atom_orient_(residue_type.force_nbr_atom_orient_),
		mass_(residue_type.mass_),
		n_actcoord_atoms_( residue_type.n_actcoord_atoms_ ),
		mol_data_(residue_type.mol_data_),
		residue_connections_(residue_type.residue_connections_),
		atom_2_residue_connection_map_(residue_type.atom_2_residue_connection_map_),
		atoms_within_one_bond_of_a_residue_connection_(residue_type.atoms_within_one_bond_of_a_residue_connection_),
		within1bonds_sets_for_atom_(residue_type.within1bonds_sets_for_atom_),
		atoms_within_two_bonds_of_a_residue_connection_(residue_type.atoms_within_two_bonds_of_a_residue_connection_),
		within2bonds_sets_for_atom_(residue_type.within2bonds_sets_for_atom_),
		lower_connect_id_( residue_type.lower_connect_id_ ),
		upper_connect_id_( residue_type.upper_connect_id_ ),
		n_non_polymeric_residue_connections_( residue_type.n_non_polymeric_residue_connections_ ),
		n_polymeric_residue_connections_( residue_type.n_polymeric_residue_connections_ ),
		force_bb_(residue_type.force_bb_),
		rna_residue_type_(residue_type.rna_residue_type_),
		carbohydrate_info_(residue_type.carbohydrate_info_),
		atom_base_indices_(residue_type.atom_base_indices_),
		abase2_indices_(residue_type.abase2_indices_),
		chi_atoms_indices_(residue_type.chi_atoms_indices_),
		nu_atoms_indices_(residue_type.nu_atoms_indices_),
		mainchain_atoms_indices_(residue_type.mainchain_atoms_indices_),
		nbr_atom_indices_(residue_type.nbr_atom_indices_),
		actcoord_atoms_indices_(residue_type.actcoord_atoms_indices_),
		cut_bond_neighbor_indices_(residue_type.cut_bond_neighbor_indices_),
		atom_shadowed_indices_(residue_type.atom_shadowed_indices_),
		finalized_(residue_type.finalized_),
		defined_adducts_(residue_type.defined_adducts_),
		nondefault_(residue_type.nondefault_),
		base_restype_name_(residue_type.base_restype_name_),
		serialized_(residue_type.serialized_)
{
	// When you copy vertex descriptors from cached data, the vertex descriptors are pointing to the old copied graph.
	// New vertices are assigned.  You have to map the old vertex to the new vertex.

	std::map<ED, ED> old_edges_to_new_edges;
	for(
		EIterPair ep = boost::edges(graph_), old_ep = boost::edges(residue_type.graph_);
		ep.first != ep.second;
		++ep.first, ++old_ep.first
	){
		EIter e_iter = ep.first;
		EIter old_e_iter = old_ep.first;
		ED  ed = *e_iter;
		ED old_ed = *old_e_iter;
		old_edges_to_new_edges[old_ed] = ed;
	}

	std::list<utility::vector1<ED> >  old_rings(rings_and_their_edges_);
	rings_and_their_edges_.clear();
	for(std::list<utility::vector1<ED> >::const_iterator it = old_rings.begin(); it != old_rings.end(); ++it){
		utility::vector1<ED> old_edges = *it;
		utility::vector1<ED> new_edges;
		for(Size i = 1; i <= old_edges.size(); ++i){
			new_edges.push_back(old_edges_to_new_edges[old_edges[i]]);
		}
		rings_and_their_edges_.push_back(new_edges);
	}


	// Setup the atom_name_to_vd_ and map old to new...
	std::map<VD, VD> old_to_new;
	for(
			VIterPair vp = boost::vertices(graph_), old_vp= boost::vertices(residue_type.graph_);
			vp.first != vp.second;
			++vp.first, ++old_vp.first
	){
		VIter v_iter= vp.first;
		VIter old_v_iter= old_vp.first;
		VD vd = *v_iter;
		VD old_vd = *old_v_iter;
		old_to_new[old_vd] = vd; //Assuming the boost::graph copy preserves ordering within the vertices list
		Atom & a = graph_[vd];
		assert( a == residue_type.graph_[old_vd]);
#ifndef NDEBUG
		NameVDInserted const name_vd_inserted =
#endif
				atom_name_to_vd_.insert(NameVDPair(a.name(), vd));
		assert(name_vd_inserted.second); // Don't add atoms with the same name
		atom_name_to_vd_.insert( NameVDPair( strip_whitespace( a.name() ), vd) );
		//assert(strip_name_vd_inserted.second); // If this is 4 chars, than it will be the same as before.
	}
	// Setup the temporary ordered_atoms_ vector for refactor
	VDs::const_iterator begin = residue_type.ordered_atoms_.begin();
	VDs::const_iterator const end = residue_type.ordered_atoms_.end();
	for(; begin != end; ++begin){
		VD old_vd = *begin;
		VD vd = old_to_new[old_vd];
		ordered_atoms_.push_back(vd);
		vd_to_index_[vd] = ordered_atoms_.size();
	}
	std::map<VD, VD> old_atom_base(atom_base_);
	atom_base_.clear();
	for(std::map<VD, VD>::iterator it = old_atom_base.begin(); it != old_atom_base.end(); ++it){
		VD old_key = it->first;
		VD old_value = it->second;
		VD new_key = old_to_new[old_key];
		VD new_value = old_to_new[old_value];
		atom_base_[new_key] = new_value;
	}

	std::map<VD, VD> old_atom_shadowed(atom_shadowed_);
	atom_shadowed_.clear();
	for(std::map<VD, VD>::iterator it = old_atom_shadowed.begin(); it != old_atom_shadowed.end(); ++it){
		VD old_key = it->first;
		VD old_value = it->second;
		VD new_key = old_to_new[old_key];
		VD new_value = old_to_new[old_value];
		atom_shadowed_[new_key] = new_value;
	}

	utility::vector1<utility::vector1<VD> >  old_chi_atoms(chi_atoms_);
	chi_atoms_.clear();
	//chi_atoms_.resize(old_chi_atoms.size());
	for(utility::vector1<utility::vector1<VD> >::const_iterator it= old_chi_atoms.begin(); it != old_chi_atoms.end(); ++it){
		utility::vector1<VD> old_vector = *it;
		assert(old_vector.size() == 4);
		utility::vector1<VD> new_vector;
		for(Size i= 1; i<= old_vector.size(); ++i){
			new_vector.push_back(old_to_new[old_vector[i]]);
		}
		chi_atoms_.push_back(new_vector);
	}

	utility::vector1<utility::vector1<VD> >  old_nu_atoms(nu_atoms_);
	nu_atoms_.clear();
	//chi_atoms_.resize(old_chi_atoms.size());
	for(utility::vector1<utility::vector1<VD> >::const_iterator it= old_nu_atoms.begin(); it != old_nu_atoms.end(); ++it){
		utility::vector1<VD> old_vector = *it;
		assert(old_vector.size() == 4);
		utility::vector1<VD> new_vector;
		for(Size i= 1; i<= old_vector.size(); ++i){
			new_vector.push_back(old_to_new[old_vector[i]]);
		}
		nu_atoms_.push_back(new_vector);
	}

	utility::vector1<VD> old_mainchain(mainchain_atoms_);
	mainchain_atoms_.clear();
	for(Size i= 1; i <= old_mainchain.size(); ++i){
		mainchain_atoms_.push_back( old_to_new[ old_mainchain[i] ]);
	}

	if( nbr_atom_ != ResidueGraph::null_vertex() ) {
		VD old_nbr = nbr_atom_;
		nbr_atom_ = old_to_new[old_nbr];
	}

	utility::vector1<VD> old_act(actcoord_atoms_);
	actcoord_atoms_.clear();
	for(Size i=1; i<= old_act.size(); ++i){
		actcoord_atoms_.push_back(old_to_new[old_act[i]]);
	}


	std::map<VD, utility::vector1<VD> > old_cut_bonds(cut_bond_neighbor_);
	cut_bond_neighbor_.clear();
	for(std::map<VD, utility::vector1<VD> >::const_iterator it = old_cut_bonds.begin(); it != old_cut_bonds.end(); ++it){
		VD old_key = it->first;
		utility::vector1<VD> old_value = it->second;
		VD new_key = old_to_new[old_key];
		utility::vector1<VD> new_value;
		for(Size i=1; i<= old_value.size(); ++i){
			new_value.push_back(old_to_new[ old_value[i] ] );
		}
		cut_bond_neighbor_[ new_key ] = new_value;
	}

	std::map< VD, AtomICoor > old_icoor(icoor_);
	icoor_.clear();
	for(std::map<VD, AtomICoor >::const_iterator it= old_icoor.begin(); it != old_icoor.end(); ++it){
		VD old_key = it->first;
		VD new_key = old_to_new[old_key];
		AtomICoor old_icoor = it->second; //now we have to change the vertex descriptors within icoor. They are pointing to an old vertex descriptor
		for(Size i=1; i<= 3; ++i){
			ICoorAtomID & stub_atom( old_icoor.stub_atom(i) );
			if ( stub_atom.type() == ICoorAtomID::INTERNAL ) {
				VD old_vertex = stub_atom.vertex();
				VD new_vertex = old_to_new[old_vertex];
				stub_atom.vertex(new_vertex);
			}
		}
		icoor_[new_key] = old_icoor;
	}


	for(Size i=1; i<= residue_connections_.size(); ++i){
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

	for(Size index=1; index<= natoms(); ++index){
		utility::vector1< core::Size > const orbs(graph_[ordered_atoms_[index]].bonded_orbitals());
		for(utility::vector1< core::Size >::const_iterator orb = orbs.begin(); orb != orbs.end(); ++orb)
		{
			orbitals_[*orb].new_icoor().vertex1( old_to_new[ orbitals_[*orb].new_icoor().vertex1()  ] );
			orbitals_[*orb].new_icoor().vertex2(old_to_new[ orbitals_[*orb].new_icoor().vertex2()  ] );
			orbitals_[*orb].new_icoor().vertex3( old_to_new[ orbitals_[*orb].new_icoor().vertex3()  ] );
		}
	}

	utility::vector1<VD> old_bb(force_bb_);
	force_bb_.clear();
	for(Size i=1; i<= old_bb.size(); ++i){
		force_bb_.push_back(old_to_new[ old_bb[i] ]);
	}

}


///
ResidueTypeSet const &
ResidueType::residue_type_set() const
{
	if ( residue_type_set_ == 0 ) {
		utility_exit_with_message( "ResidueType::residue_type_set: pointer is not set!");
	}
	return *residue_type_set_;
}

//////////////////////////////////////////////////////////////////////////////

/// @brief make a copy
ResidueTypeOP
ResidueType::clone() const
{
	ResidueTypeOP rsd_ptr( new ResidueType( *this ) );
	return rsd_ptr;
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
Atom & ResidueType::atom(VD const vd){
	return graph_[ vd ];
}
Atom const & ResidueType::atom(VD const vd) const{
	return graph_[ vd ];
}
Atom & ResidueType::atom(std::string const & atom_name){
	return graph_[ atom_name_to_vd_[atom_name] ];
}
Atom const & ResidueType::atom(std::string const & atom_name) const{
	NameVDMap::const_iterator found = atom_name_to_vd_.find( atom_name );
	assert( found != atom_name_to_vd_.end());
	return graph_[  found->second ];
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
			assert( n_polymeric_residue_connections_ != 0 );
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
			assert( n_polymeric_residue_connections_ != 0 );
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
	assert( is_polymer_ );
	assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ];
}

ResidueConnection const &
ResidueType::lower_connect() const
{
	assert( is_polymer_ );
	assert( lower_connect_id_ != 0 );
	return residue_connections_[ lower_connect_id_ ];
}

Size
ResidueType::lower_connect_atom() const {
	assert( is_polymer_ );
	assert( lower_connect_id_ != 0 );
	return vd_to_index_.find(residue_connections_[ lower_connect_id_ ].vertex())->second;
}

/// @brief index number of the atom which connects to the upper connection
Size
ResidueType::upper_connect_atom() const
{
	assert( is_polymer_ );
	assert( upper_connect_id_ != 0 );
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

///////////////////////////////////////////////////////////////////////////////
///
/// @brief Get the chemical atom_type for this atom by it index number in this residue
///
/// @details If we want the atom_type index (integer), we get this from
/// the conformation::Atom itself, as seen in the code below
AtomType const &
ResidueType::atom_type( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <= ordered_atoms_.size()), "ResidueType::atom_type( Size const atomno ): atomno is not in this ResidueType!");
	assert( (atomno > 0) && (atomno <= ordered_atoms_.size()) );
	return ( *atom_types_ )[ graph_[ ordered_atoms_[atomno] ].atom_type_index() ];
}

AtomType const &
ResidueType::atom_type( VD const vd ) const
{
	//TODO: Is there a better way of validating the VD instead of checking all the ordered atoms?
	PyAssert( (vd != ResidueGraph::null_vertex() ) &&
			( std::find(ordered_atoms_.begin(), ordered_atoms_.end(), vd) != ordered_atoms_.end()),
			"ResidueType::atom_type( VD const vd ): vd is not in this ResidueType!");
	assert( (vd != ResidueGraph::null_vertex()) &&
			( std::find(ordered_atoms_.begin(), ordered_atoms_.end(), vd) != ordered_atoms_.end()) );
	return ( *atom_types_ )[ graph_[ vd ].atom_type_index() ];
}

/// @brief Get the atom name by index
std::string const &
ResidueType::atom_name( Size const index ) const
{
	PyAssert((index > 0) && (index <= ordered_atoms_.size()), "ResidueType::atom_name( Size const index ): index is not in this ResidueType!");
	assert((index > 0) && (index <= ordered_atoms_.size()));
	return graph_[ ordered_atoms_[index] ].name();
}

/// @brief get index of an atom's base atom
Size
ResidueType::atom_base( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <= ordered_atoms_.size()), "ResidueType::atom_base( Size const atomno ): atomno is not in this ResidueType!");
	return atom_base_indices_[atomno];
}

/// @brief get index of an atom's second base atom
Size
ResidueType::abase2( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <=  ordered_atoms_.size()), "ResidueType::abase2( Size const atomno ): atomno is not in this ResidueType!");
	return abase2_indices_[atomno];
}

///@brief Counts the number of virtual atoms and returns the count.
///@details The virtual count is not stored in the resiude type.  This count is performed on the fly, and
///can hurt performance if reapeatedly carried out.  Not intended for use in large loops -- instead, call
///once and store the value.
///@author Vikram K. Mulligan (vmullig@uw.edu)
Size
ResidueType::n_virtual_atoms () const
{
	core::Size virtcount = 0;
	for(core::Size ia=1, iamax=natoms(); ia<=iamax; ++ia) {
		if(is_virtual(ia)) ++virtcount;
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
void
ResidueType::add_atom(
		std::string const & atom_name,
		std::string const & atom_type_name,
		std::string const & mm_atom_type_name,
		Real const charge
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	assert(atom_name_to_vd_.find(atom_name) == atom_name_to_vd_.end());
	assert(atom_name_to_vd_.find( strip_whitespace(atom_name)) == atom_name_to_vd_.end());

	// index lookup by name
	// store the atom types
	// the next calls will fail if the atom type name is unrecognized
	Size const type( atom_types_->atom_type_index( atom_type_name ) );
	Size const mm_type( mm_atom_types_->atom_type_index( mm_atom_type_name ) );

	if ( (*atom_types_)[type].is_acceptor() ) ++n_hbond_acceptors_;
	if ( (*atom_types_)[type].is_donor() ) ++n_hbond_donors_;


	//get the element information
	ElementCOP element;
	if( elements_ )	{ // Be robust if elements_ isn't defined.
		std::string const & element_name= (*atom_types_)[type].element();
		int const element_index = elements_->element_index(element_name);
		//std::cout << elements_->element_index(element_name) << " " << element_name << std::endl;
		element = (*elements_)[element_index];
		mass_ += element->weight();
		//mass_ += (*elements_)[element_index]->weight();
	} else {
		tr.Warning << "WARNING Elements set undefined." << std::endl;
	}



	Atom atom(
			atom_name,
			mm_atom_type_name,
			type, mm_type,
			element, charge,
			Vector(0.0)
	);

	VD v = graph_.add_vertex( atom );
	atom_name_to_vd_[ atom_name ] = v;
	atom_name_to_vd_[ strip_whitespace( atom_name ) ] = v;
	ordered_atoms_.push_back(v);
	AtomAP graph_atom = &graph_[v];




	parents_.push_back(0);

	//setup data for atoms
	AtomType const & atype = (*atom_types_)[type];
	graph_atom->is_polar_hydrogen(atype.is_polar_hydrogen());
	graph_atom->is_hydrogen(atype.is_hydrogen());
	graph_atom->is_haro(atype.is_haro());
	graph_atom->is_acceptor(atype.is_acceptor());
	graph_atom->is_virtual(atype.is_virtual());
	graph_atom->has_orbitals(atype.atom_has_orbital());


	// allocate space for the new atom !!!!!!!!!!!!!!!!!!!!!!
	// eg, in the atom/resconn map
	assert( atom_2_residue_connection_map_.size() == ordered_atoms_.size()-1 );

	atom_2_residue_connection_map_.resize( ordered_atoms_.size() );
	bonded_neighbor_.resize(ordered_atoms_.size());
	bonded_neighbor_type_.resize(ordered_atoms_.size());
	vd_to_index_[v] = ordered_atoms_.size();
	atom_base_[v] =v;  //graph_atom->atom_base(ordered_atoms_.size()); // base defaults to self
	icoor_[v] =  AtomICoor( 0.0, 0.0, 0.0, atom_name, atom_name, atom_name, *this );


}

/// @brief flag an atom for deletion by adding its index to the delete_atom_ list
void
ResidueType::delete_atom( std::string const & name )
{
	assert( has( name ) );
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

/// @brief set atom type
void
ResidueType::set_atom_type(
		std::string const & atom_name,
		std::string const & atom_type_name
)
{
	int const atom_type_index = atom_types_->atom_type_index( atom_type_name );
	Atom & a = graph_[ atom_name_to_vd_[atom_name] ];
	AtomType const & atype = (*atom_types_)[atom_type_index];
	a.atom_type_index( atom_type_index );
	a.is_polar_hydrogen(atype.is_polar_hydrogen());
	a.is_hydrogen(atype.is_hydrogen());
	a.is_haro(atype.is_haro());
	a.is_acceptor(atype.is_acceptor());
	a.is_virtual(atype.is_virtual());
	a.has_orbitals(atype.atom_has_orbital());
}


/// @brief set mm atom type
void
ResidueType::set_mm_atom_type(
		std::string const & atom_name,
		std::string const & mm_atom_type_name
)
{
	int const mm_type_index = mm_atom_types_->atom_type_index( mm_atom_type_name );
	Atom & a = graph_[ atom_name_to_vd_[atom_name] ];
	a.mm_atom_type_index( mm_type_index );
}

/// @brief Get the MM atom_type for this atom by its index number in this residue
MMAtomType const &
ResidueType::mm_atom_type( Size const atomno ) const
{
	return ( *mm_atom_types_ )[graph_[ ordered_atoms_[atomno] ].mm_atom_type_index() ];
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
	if ( ! gasteiger_atom_types_ ) {
		gasteiger_atom_types_ = ChemicalManager::get_instance()->gasteiger_atom_type_set();
	}
	gasteiger::GasteigerAtomTypeDataCOP gasteiger_type = gasteiger_atom_types_->atom_type( gasteiger_atom_type_name );
	Atom & a = graph_[ atom_name_to_vd_[atom_name] ];
	a.gasteiger_atom_type( gasteiger_type );
}

/// @brief Get the Gasteiger atom_type for this atom by its index number in this residue
gasteiger::GasteigerAtomTypeDataCOP
ResidueType::gasteiger_atom_type( Size const atomno ) const
{
	return graph_[ ordered_atoms_[atomno] ].gasteiger_atom_type();
}
VD
ResidueType::vd_from_name( std::string const & name) const{
	NameVDMap::const_iterator atom_name_to_vd_iter( atom_name_to_vd_.find( name ) );
	if ( atom_name_to_vd_iter == atom_name_to_vd_.end() ) {
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		assert(false);
	}

	return atom_name_to_vd_iter->second;
}


orbitals::OrbitalType const &
ResidueType::orbital_type(int const orbital_index)const
{
	return ( *orbital_types_ )[ orbitals_[ orbital_index ].orbital_type_index() ];
}

// Return a pointer to the object containing the set of ring conformers possible for this saccharide.
core::chemical::RingConformerSetCOP
ResidueType::ring_conformer_set() const
{
	return conformer_set_;
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

	// increment orbital count
	++n_orbitals_;

	// store the atom type
	// the next call will fail if the orbital type name is unrecognized
	Size type( orbital_types_->orbital_type_index( orbital_type_name ) );

	// store the name
	orbitals_.push_back(Orbital(orbital_name, type, Vector(0.0)));
	assert( orbitals_.size() == n_orbitals_ );

	orbitals_index_[ orbital_name ] = n_orbitals_;
	orbitals_index_[ strip_whitespace( orbital_name ) ] = n_orbitals_;
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Add an atom to the list of atoms that can potentially form a bond to a metal ion.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ResidueType::add_metalbinding_atom (
	std::string const atom_name
) {
	if(!has(atom_name)) {
		std::string message = "Error in adding metal-binding atom to residue type " + name3() + ". Atom " + atom_name + " was not found.";
		utility_exit_with_message(message);
	}
	metal_binding_atoms_.push_back( atom_name ); //Store names rather than indices, since indices might change.
	return;
}

///////////////////////////////////////////////////////////////////////////////

/// @details add a bond between atom1 and atom2 and add a BondType object referencing the bond (default bond type of SingleBond)
/** update bonded_neighbor_ and resize it as necessary **/
void
ResidueType::add_bond(std::string const & atom_name1, std::string const & atom_name2)
{
	add_bond(atom_name1, atom_name2, SingleBond);
}

/// @details add a bond between atom1 and atom2 and add a BondType object referencing the bond using the specified bondName
void
ResidueType::add_bond(std::string const & atom_name1, std::string const & atom_name2, BondName bondLabel)
{
	++nbonds_;
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
		utility_exit_with_message("ResidueType:: shouldnt get here -- resizing in add_atom");
	}

	bonded_neighbor_[ i1 ].push_back( i2 );
	bonded_neighbor_[ i2 ].push_back( i1 );
	bonded_neighbor_type_[i1].push_back(bondLabel);
	bonded_neighbor_type_[i2].push_back(bondLabel); 	//bondType_vector_.push_back(BondType(i1,i2,bondLabel));

	NameVDMap::const_iterator source = atom_name_to_vd_.find( atom_name1 );
	NameVDMap::const_iterator target = atom_name_to_vd_.find( atom_name2 );
	assert( source != atom_name_to_vd_.end());
	assert( target != atom_name_to_vd_.end());
	VD const vd_source = source->second;
	VD const vd_target = target->second;

	// check if bond already exists...
	if( boost::edge(vd_source, vd_target, graph_).second ){
		utility_exit_with_message( "dont add residue bonds more than once!" );
	}

	ResidueGraph::edge_descriptor e_added;
	bool added;
	boost::tie(e_added, added) = graph_.add_edge( vd_source, vd_target, Bond(-1, bondLabel)); /// -1 means Bond distance not set here. This will be fixed in the future
	assert(added);
}

///////////////////////////////////////////////////////////////////////////////
///@brief add an orbital bond between an atom and an orbital.
///@note NOTE!!!!! This is indexed based upon atoms, not orbitals. That means that in your params file
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

	if(atom_name_to_vd_.find(atom_name1) == atom_name_to_vd_.end()){
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

	utility::vector1<VD> const i1_nbrs(cut_bond_neighbor_.find(vd_source)->second );
	if ( std::find( i1_nbrs.begin(), i1_nbrs.end(), vd_target ) != i1_nbrs.end() ) {
		utility_exit_with_message( "don't add residue bonds more than once!" );
	}

	cut_bond_neighbor_[vd_source].push_back(vd_target);
	cut_bond_neighbor_[vd_target].push_back(vd_source);
}


///////////////////////////////////////////////////////////////////////////////

// Add a chi (side-chain) angle defined by four atoms.
void
ResidueType::add_chi(
		Size const chino,
		std::string const & atom_name1,
		std::string const & atom_name2,
		std::string const & atom_name3,
		std::string const & atom_name4
)
{
	// Signal that we need to update the derived data.
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ||
			!has( atom_name3 ) || !has( atom_name4 ) ) {
		utility_exit_with_message("ResidueType::add_chi: atoms don't exist!" );
	}

	utility::vector1<VD> atoms;
	atoms.push_back( ordered_atoms_[ atom_index( atom_name1 ) ]);
	atoms.push_back( ordered_atoms_[ atom_index( atom_name2 ) ]);
	atoms.push_back( ordered_atoms_[ atom_index( atom_name3 ) ]);
	atoms.push_back( ordered_atoms_[ atom_index( atom_name4 ) ]);
	if ( chi_atoms_.size() < chino ) {
		chi_atoms_.resize( chino );
		chi_rotamers_.resize( chino );
		chi_2_proton_chi_.resize( chino );
	}
	chi_atoms_[chino] = atoms;
	is_proton_chi_.push_back( false );
	chi_2_proton_chi_[ chino ] = 0;
}  // add_chi

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

	if (!has(atom_name1) || !has(atom_name2) || !has(atom_name3) || !has(atom_name4)) {
		utility_exit_with_message("ResidueType::add_nu: Requested atoms don't exist in this ResidueType!");
	}

	utility::vector1<VD> atoms;
	atoms.push_back(ordered_atoms_[atom_index(atom_name1)]);
	atoms.push_back(ordered_atoms_[atom_index(atom_name2)]);
	atoms.push_back(ordered_atoms_[atom_index(atom_name3)]);
	atoms.push_back(ordered_atoms_[atom_index(atom_name4)]);

	if (nu_atoms_.size() < nu_index) {
		nu_atoms_.resize(nu_index);
	}
	nu_atoms_[nu_index] = atoms;
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
		utility::vector1< Real > dihedral_samples,
		utility::vector1< Real > extra_samples
)
{
	assert( is_proton_chi_.size() >= nchi() );
	assert( chi_2_proton_chi_.size() >= nchi() );
	if( chino > nchi() ) {
		utility_exit_with_message("Error setting proton chi: Chi to set as proton chi does not exist.");
	}
	is_proton_chi_[ chino ] = true;
	proton_chis_.push_back( chino );
	proton_chi_samples_.push_back( dihedral_samples );
	proton_chi_extra_samples_.push_back( extra_samples );
	chi_2_proton_chi_[ chino ] = proton_chis_.size();
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
	if ( chi_rotamers_.size() < chino ) chi_rotamers_.resize( chino );
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


///////////////////////////////////////////////////////////////////////////////

/// @details sets atom_base_[atom1] = atom2
/** resize atom_base_ vector as necessary **/
void
ResidueType::set_atom_base(
		std::string const & atom_name1,
		std::string const & atom_name2
)
{

	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
		utility_exit_with_message( "set_atom_base: atoms dont exist!" );
	}

	VD const vd_source = atom_name_to_vd_.find( atom_name1 )->second; //source->second;
	VD const vd_target = atom_name_to_vd_.find( atom_name2 )->second; //target->second;

	// atom base must be bonded.
	if( !boost::edge(vd_source, vd_target, graph_).second ){
		utility_exit_with_message( "in set_atom_base(), atoms must be bonded to set atom base!" );
	}

	//make sure that you do not set an atom base at a cut bond
	utility::vector1<VD> const i1_nbrs(cut_bond_neighbor_.find(vd_source)->second );
	if ( std::find( i1_nbrs.begin(), i1_nbrs.end(), vd_target ) != i1_nbrs.end() ) {
		utility_exit_with_message( "Don't set atom bases to cut bonds!" );
	}

	atom_base_[vd_source] = vd_target;

}

/// @brief set indices of all mainchain atoms
void
ResidueType::set_mainchain_atoms( AtomIndices const & mainchain )
{
	utility::vector1<VD> vd_mainchain;
	for(Size i=1; i<= mainchain.size(); ++i){
		vd_mainchain.push_back(ordered_atoms_[mainchain[i]]);
	}
	mainchain_atoms_ = vd_mainchain;
}


///////////////////////////////////////////////////////////////////////////////

/// @details get all specified properties for this residue type
utility::vector1< std::string > const &
ResidueType::properties() const
{
	return properties_;
}

///////////////////////////////////////////////////////////////////////////////

/// @details add a property to this residue
/** update boolean property member data accordingly **/
void
ResidueType::add_property( std::string const & property )
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( property == "POLYMER" ) {
		is_polymer_ = true;
	} else if ( property == "PROTEIN" ) {
		is_protein_ = true;
		is_polymer_ = true;
	} else if ( property == "ALPHA_AA" ) {
		is_protein_ = true;
		is_polymer_ = true;
		is_alpha_aa_ = true;
	} else if ( property == "BETA_AA" ) {
		is_protein_ = true;
		is_polymer_ = true;
		is_beta_aa_ = true;
	} else if ( property == "L_AA" ) {
		is_protein_ = true;
		is_polymer_ = true;
		is_alpha_aa_ = true;
		is_l_aa_ = true;
	} else if ( property == "D_AA" ) {
		is_protein_ = true;
		is_polymer_ = true;
		is_alpha_aa_ = true;
		is_d_aa_ = true;
	} else if ( property == "POLAR" ) {
		is_polar_ = true;
	} else if( property == "SC_ORBITALS"){
		has_sc_orbitals_ = true;
	} else if ( property == "CHARGED" ) {
		is_charged_ = true;
	} else if ( property == "AROMATIC" ) {
		is_aromatic_ = true;
	} else if ( property == "CYCLIC" ) {
		is_cyclic_ = true;
	} else if ( property == "COARSE" ) {
		is_coarse_ = true; //currently only for RNA
	} else if ( property == "DNA" ) {
		is_DNA_ = true;
		is_NA_ = true;
		is_polymer_ = true;
	} else if ( property == "RNA" ) {
		is_RNA_ = true;
		is_NA_ = true;
		is_polymer_ = true;
	} else if ( property == "CARBOHYDRATE") {
		is_carbohydrate_ = true;
	} else if ( property == "LIPID" ) {
		is_lipid_ = true;
	} else if ( property == "LIGAND" ) {
		is_ligand_ = true;
	} else if ( property == "METAL" ) { //Is this a metal ion?
		is_metal_ = true;
	} else if ( property == "METALBINDING" ) { //Can this amino acid residue bind metals?
		is_metalbinding_ = true;
	} else if ( property == "SURFACE" ) {
		is_surface_ = true;
	} else if ( property == "LOWER_TERMINUS" ) {
		is_terminus_ = true;
		is_lower_terminus_ = true;
	} else if ( property == "UPPER_TERMINUS" ) {
		is_terminus_ = true;
		is_upper_terminus_ = true;
	} else if ( property == "BRANCH_LOWER_TERMINUS" ) {
		is_branch_lower_terminus_ = true;
	} else if ( property == "LOWERTERM_TRUNC" ) {
		is_terminus_ = true;
		is_lower_terminus_ = true;
	} else if ( property == "UPPERTERM_TRUNC" ) {
		is_terminus_ = true;
		is_upper_terminus_ = true;
	} else if ( property == "PHOSPHONATE" ) {
		is_polymer_ = true;
		is_phosphonate_ = true;
	} else if ( property == "PHOSPHONATE_UPPER" ) {
		is_terminus_ = true;
		is_upper_terminus_ = true;
		is_phosphonate_ = true;
		is_phosphonate_upper_ = true;
	} else if ( property == "TERMINUS" ) {
		is_terminus_ = true;
	} else if ( property == "ACETYLATED_NTERMINUS" ) {
		is_terminus_ = true;
		is_lower_terminus_ = true;
		is_acetylated_nterminus_ = true;
	} else if ( property == "METHYLATED_CTERMINUS" ) {
		is_terminus_ = true;
		is_upper_terminus_ = true;
		is_methylated_cterminus_ = true;
	} else if (property == "TAUTOMER") {
		; // this is for HIS_D, following someone's suggestion in Patch applications. -- rhiju
	} else if (property == "BRANCH_POINT") {
		;  // Null statement for now.... ~ Labonte
	} else if (carbohydrates::CarbohydrateInfo::sugar_properties().contains(property)) {
		;  // Null statement -- these properties will be added to carbohydrate_info_ by update_derived_data().
	} else {
		tr.Warning << "WARNING:: unrecognized residue type property: " << property << std::endl;
	}

	properties_.push_back( property );
}

void
ResidueType::add_numeric_property(std::string const & tag, core::Real value)
{
	numeric_properties_.insert(std::make_pair(tag,value));
}


void
ResidueType::add_string_property(std::string const & tag, std::string value)
{
	string_properties_.insert(std::make_pair(tag,value));
}
///////////////////////////////////////////////////////////////////////////////

/// @details delete a property to this residue
/** update boolean property member data accordingly **/
//    Added by Andy M. Chen in June 2009
//    This is needed for deleting properties, which occurs in certain PTM's
void
ResidueType::delete_property( std::string const & property )
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( property == "POLYMER" ) {
		is_polymer_ = false;
	} else if ( property == "PROTEIN" ) {
		is_protein_ = false;
	} else if ( property == "ALPHA_AA" ) {
		is_alpha_aa_ = false;
	} else if ( property == "BETA_AA" ) {
		is_beta_aa_ = false;
	} else if ( property == "L_AA" ) {
		is_l_aa_ = false;
	} else if ( property == "D_AA" ) {
		is_d_aa_ = false;
	} else if ( property == "POLAR" ) {
		is_polar_ = false;
	}else if(property == "SC_ORBITALS"){
		has_sc_orbitals_ = false;
	}else if ( property == "CHARGED" ) {
		is_charged_ = false;
	} else if ( property == "AROMATIC" ) {
		is_aromatic_ = false;
	} else if ( property == "CYCLIC" ) {
		is_cyclic_ = false;
	} else if ( property == "COARSE" ) {
		is_coarse_ = false;
	} else if ( property == "DNA" ) {
		is_DNA_ = false;
	} else if ( property == "RNA" ) {
		is_RNA_ = false;
	} else if ( property == "CARBOHYDRATE") {
		is_carbohydrate_ = false;
	} else if ( property == "LIPID") {
		is_lipid_ = false;
	} else if ( property == "LIGAND" ) {
		is_ligand_ = false;
	} else if ( property == "METAL" ) {
		is_metal_ = false;
	} else if ( property == "METALBINDING" ) {
		is_metalbinding_ = false;
	} else if ( property == "SURFACE" ) {
		is_surface_ = false;
	} else if ( property == "LOWER_TERMINUS" ) {
		// could add an is_lower_terminus_ bool if needed?
		is_lower_terminus_ = false;
	} else if ( property == "UPPER_TERMINUS" ) {
		is_upper_terminus_ = false;
	} else if ( property == "BRANCH_LOWER_TERMINUS" ) {
		is_branch_lower_terminus_ = false;
	} else if ( property == "TERMINUS" ) {
		is_terminus_ = false;
	} else if ( property == "PHOSPHONATE" ) {
		is_phosphonate_ = false;
	} else if ( property == "PHOSPHONATE_UPPER" ) {
		is_phosphonate_upper_ = false;
	} else if ( property == "ACETYLATED_NTERMINUS" ) {
		is_acetylated_nterminus_ = false;
	} else if ( property == "METHYLATED_CTERMINUS" ) {
		is_methylated_cterminus_ = false;
	} else {
		tr.Warning << "WARNING:: unrecognized residue type property: " << property << std::endl;
	}

	utility::vector1<std::string>::iterator i = std::find(properties_.begin(), properties_.end(), property);
	properties_.erase(i);
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
	if( ! is_protein() ) {
		tr.Warning << "WARNING: ACT_COORD_ATOM specified for non-protein residue type '" << name() << "' . This doesn't make much sense." << std::endl;
	}
	finalized_ = false;
	tr.Trace << "adding act coord atom: " << name_ << ' ' << atom << std::endl;
	Size atomindex = atom_index( atom );
	actcoord_atoms_.push_back( ordered_atoms_[atomindex] );
	++n_actcoord_atoms_;
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

 **/
void
ResidueType::setup_atom_ordering()
{
	utility::vector1<VD> bb_atoms, sidechain_atoms, hydrogens;
	utility::vector1<Size> h_begin;
	utility::vector1<Size> h_end;
	//update the bonded neighbors here
	for(VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first){
		VD const & vd = *vp.first;
		//we need to iterate through the edges to update
		if(!graph_[vd].is_hydrogen()){
			if(std::find( force_bb_.begin(), force_bb_.end(), vd ) != force_bb_.end() ){
				bb_atoms.push_back(vd);
				h_begin.push_back(hydrogens.size() +1 );
				bool h_present(false);
				for(OutEdgeIterPair ep = boost::out_edges(vd, graph_); ep.first != ep.second; ++ep.first){ //iterate through the edges for hydrogens
					VD target = boost::target(*ep.first, graph_);
					if(graph_[target].is_hydrogen()){
						hydrogens.push_back(target);
						h_present = true;
					}
				}
				if(h_present){ h_end.push_back(hydrogens.size()); } else{ h_end.push_back(0); }
			} else {
				sidechain_atoms.push_back(vd);
			}
		}
	}
	//std::cout << "bb_atoms.size(): " << bb_atoms.size() << std::endl;
	//once we set up the sidechains, we have to push back the hydrogens for the sidechains
	for(Size i= 1; i<= sidechain_atoms.size(); ++i){
		h_begin.push_back(hydrogens.size() +1 );
		bool h_present(false);
		for(OutEdgeIterPair ep = boost::out_edges(sidechain_atoms[i], graph_); ep.first != ep.second; ++ep.first){ //iterate through the edges for hydrogens
			VD target = boost::target(*ep.first, graph_);
			if(graph_[target].is_hydrogen()){
				hydrogens.push_back(target);
				h_present = true;
			}
		}
		if(h_present){
			h_end.push_back(hydrogens.size());
		} else{
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


	n_backbone_heavyatoms_ = bb_atoms.size();
	nheavyatoms_ = bb_atoms.size() + sidechain_atoms.size();


	//This step is done at the end after the nheavy_atoms_ has been assigned.
	attached_H_begin_.clear();
	attached_H_end_.clear();
	for(Size i =1; i<= h_begin.size(); ++i){
		attached_H_begin_.push_back(nheavyatoms_+h_begin[i]);
		if(h_end[i] ==0){
			attached_H_end_.push_back(0);
		}else{ attached_H_end_.push_back(nheavyatoms_+h_end[i]); }
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

	bonded_neighbor_.clear();
	bonded_neighbor_type_.clear();
	bonded_neighbor_.resize(natoms());
	bonded_neighbor_type_.resize(natoms());
	//setup bond ordering!
	for(VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first){
		VD const & vd = *vp.first;
		bonded_neighbor_[vd_to_index_[vd]].clear(); //we have to clear the vectors first before pushing back to them
		bonded_neighbor_type_[vd_to_index_[vd]].clear(); //we have to clear the vectors first before pushing back to them
		for(OutEdgeIterPair ep = boost::out_edges(vd, graph_); ep.first != ep.second; ++ep.first){ //iterate through the edges for hydrogens
			VD target = boost::target(*ep.first, graph_);
			bonded_neighbor_[vd_to_index_[vd]].push_back(vd_to_index_[target]);
			ED const & edge = *ep.first;
			BondName bond = graph_[edge].bond_name();
			bonded_neighbor_type_[vd_to_index_[vd]].push_back(bond);
		}
	}



	atom_base_indices_.clear();
	for ( Size atomno=1; atomno<= natoms(); ++atomno ) {
		{
			if(atom_base_.find(ordered_atoms_[atomno]) == atom_base_.end() ){
				atom_base_indices_.push_back(0);
			} else {
				atom_base_indices_.push_back(vd_to_index_.find((atom_base_.find(ordered_atoms_[atomno])->second))->second);
			}
		}

		{
			if(atom_shadowed_.find(ordered_atoms_[atomno]) == atom_shadowed_.end()) {
				atom_shadowed_indices_.push_back(0);
			} else {
				atom_shadowed_indices_.push_back(vd_to_index_.find((atom_shadowed_.find(ordered_atoms_[atomno])->second))->second);
			}
		}
	}

	cut_bond_neighbor_indices_.clear();
	cut_bond_neighbor_indices_.resize(natoms());
	AtomIndices atoms;
	for ( Size atomno=1; atomno<= natoms(); ++atomno ){
		if(cut_bond_neighbor_.find(ordered_atoms_[atomno]) != cut_bond_neighbor_.end()){
			utility::vector1<VD> const & vd_cut_atoms = cut_bond_neighbor_.find( ordered_atoms_[atomno] )->second;
			for(Size i = 1; i <= vd_cut_atoms.size(); ++i){ //if you get here, you will have assigned chi atoms
				atoms.push_back(vd_to_index_.find( vd_cut_atoms[i] )->second);
			}
			cut_bond_neighbor_indices_[atomno] =  atoms;
			atoms.clear();
		}
	}

	chi_atoms_indices_.clear();
	atoms.clear();
	utility::vector1<AtomIndices> chi_atoms;
	for(Size chino=1; chino <= chi_atoms_.size(); ++chino){
		for(Size atom_index=1; atom_index <= chi_atoms_[chino].size(); ++atom_index){
			atoms.push_back( vd_to_index_.find( chi_atoms_[chino][atom_index] )->second) ;
		}
		chi_atoms_indices_.push_back(atoms);
		atoms.clear();
	}
	nu_atoms_indices_.clear();
	atoms.clear();
	for(Size nu_no=1; nu_no <= nu_atoms_.size(); ++nu_no){
		for(Size atom_index=1; atom_index <= nu_atoms_[nu_no].size(); ++atom_index){
			atoms.push_back( vd_to_index_.find( nu_atoms_[nu_no][atom_index] )->second) ;
		}
		nu_atoms_indices_.push_back(atoms);
		atoms.clear();
	}


	actcoord_atoms_indices_.clear();
	for(Size i=1; i<= actcoord_atoms_.size(); ++i){
		actcoord_atoms_indices_.push_back(vd_to_index_.find(actcoord_atoms_[i])->second);
	}

	mainchain_atoms_indices_.clear();
	for(Size i=1; i<= mainchain_atoms_.size(); ++i){
		mainchain_atoms_indices_.push_back(vd_to_index_.find(mainchain_atoms_[i])->second);
	}

	if( nbr_atom_ == ResidueGraph::null_vertex() ) {
		nbr_atom_indices_ = 0;
	} else {
		std::map<VD, Size>::const_iterator nbr_translation( vd_to_index_.find(nbr_atom_) );
		assert( nbr_translation != vd_to_index_.end() );
		nbr_atom_indices_ = nbr_translation->second;
	}

	for(Size index=1; index<= natoms(); ++index){
		utility::vector1< core::Size > const orbs(graph_[ordered_atoms_[index]].bonded_orbitals());
		for(utility::vector1< core::Size >::const_iterator orb = orbs.begin(); orb != orbs.end(); ++orb)
		{
			orbitals_[*orb].new_icoor().replace_stub1( vd_to_index_[ordered_atoms_[orbitals_[*orb].new_icoor().get_stub1()]] );
			orbitals_[*orb].new_icoor().replace_stub2( vd_to_index_[ordered_atoms_[orbitals_[*orb].new_icoor().get_stub2()]] );
			orbitals_[*orb].new_icoor().replace_stub3( vd_to_index_[ordered_atoms_[orbitals_[*orb].new_icoor().get_stub3()]] );
		}
	}


	for(Size index=1; index<= natoms(); ++index){
		for(Size i=1; i<= 3; ++i){
			ICoorAtomID & stub_atom( icoor_[ ordered_atoms_[index] ].stub_atom( i )   );
			if ( stub_atom.type() == ICoorAtomID::INTERNAL ) {
				stub_atom.atomno(   vd_to_index_.find(stub_atom.vertex())->second ); //somewhat of a problem. if vertex doesnt exist the map constructor will create a value
				assert( stub_atom.atomno() ); // this will fail if we deleted a stub atom for some other atom
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
		assert( residue_connections_[i].atomno() ); //this will fail if we deleted an atom involved in an inter-rsd connection
	}
	update_residue_connection_mapping();

}


///////////////////////////////////////////////////////////////////////////////

/// @details update derived data in ResidueType, called by finalize()
/// after primary data have been reordered, update derived data accordingly,
/// including\n, Hbond donor and acceptors, path_distance etc.
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
		if(atom.has_orbitals()) atoms_with_orb_index_.push_back(i); //get atoms with orbitals on it
		if(atom.is_haro()) Haro_index_.push_back( i ); //get aromatic hydrogens
		if(atom.is_polar_hydrogen()) Hpol_index_.push_back( i ); //get polar hydrogens
		if ( atom.is_acceptor() && !atom.is_virtual() ){
			accpt_pos_.push_back( i );
			if ( i > n_backbone_heavyatoms_ ) {
				accpt_pos_sc_.push_back( i );
			}
		}
		if ( atom.is_polar_hydrogen() && !atom.is_virtual() ){
			Hpos_polar_.push_back( i );
			if ( i >= first_sidechain_hydrogen_ ) {
				Hpos_polar_sc_.push_back( i );
			}
		}
		if ( atom.is_hydrogen() && !atom.is_polar_hydrogen() ){
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
	for(Size Aindex=1; Aindex<= ordered_atoms_.size(); ++Aindex){
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
		assert(acc_base == atom_base(acceptor_position) );
		assert( acc_base != 0 );
		AtomIndices const & i_nbrs(bonded_neighbor(acceptor_position));
		if ( i_nbrs.size() == 0 ) {
			utility_exit_with_message( "failed to set abase2 for acceptor atom, it has no nbrs!" );
		} else if ( i_nbrs.size() == 1 ) {
			//assert( i_nbrs[1] == acc_base );
			abase2_[ordered_atoms_[acceptor_position]] = ordered_atoms_[atom_base(acc_base)];
			//iwd  The first child of the root is root's atom_base.
			//iwd  But if that child has no children, it ends up as its own abase2.
			//iwd  So instead we use the second child of the parent,
			//iwd  which must exist if there are 3+ atoms in this tree.
			if(vd_to_index_.find( abase2_.find(ordered_atoms_[acceptor_position])->second )->second == acceptor_position ) {
				AtomIndices const & i_base_nbrs(bonded_neighbor(acc_base) );
				for(Size jj = 1, jj_end = i_base_nbrs.size(); jj <= jj_end; ++jj) {
					if(i_base_nbrs[ jj ] != acceptor_position) {
						abase2_[ordered_atoms_[acceptor_position] ]  = ordered_atoms_[i_base_nbrs[ jj ]];
						break;
					}
				}
			}
			// assert(abase2(acceptor_position)!=acceptor_position && abase2(acceptor_position) != acc_base && abase2(acceptor_position) != 0 );
		} else if ( i_nbrs[1] == acc_base ) {
			abase2_[ordered_atoms_[acceptor_position]] = ordered_atoms_[i_nbrs[2]];
		} else {
			abase2_[ordered_atoms_[acceptor_position]] = ordered_atoms_[i_nbrs[1] ];
		}
	}

	abase2_indices_.clear();
	for ( Size atomno=1; atomno<= natoms(); ++atomno ) {
		if(abase2_.find(ordered_atoms_[atomno]) == abase2_.end() ){
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
		for (Size jj = 1; jj <= natoms(); ++jj ) {
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
	ndihe_ = dihedral_atom_sets_.size();


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
	// Ring size is determined by the number of NU angles listed in the .params file, which should always be 2 less
	// than the size of the ring.
	if ( is_cyclic_ ) {
		// ring_size could be made a private datum, but it only really makes sense for monocyclics.  Since its only use
		// for the time being is to set the proper RingConformerSet, I'll just leave it as a local variable here.
		// ~Labonte
		Size ring_size = nu_atoms_indices_.size() + 2;
		conformer_set_ = new RingConformerSet(ring_size);
	}

	if(is_RNA_){ //reinitialize and RNA derived data.
		//Reinitialize rna_residue_type_ object! This also make sure rna_residue_type_ didn't inherit anything from the previous update!
		//It appears that the rna_residue_type_ is shared across multiple ResidueType object, if the rna_residue_type_ is not reinitialized here!
		rna_residue_type_ = new core::chemical::rna::RNA_ResidueType;
		//update_last_controlling_chi is treated separately for RNA case. Parin Sripakdeevong, June 26, 2011
		rna_residue_type_->rna_update_last_controlling_chi( this, last_controlling_chi_, atoms_last_controlled_by_chi_);
		rna_residue_type_->update_derived_rna_data( this );
	} else if (is_carbohydrate_) {
		carbohydrate_info_ = new chemical::carbohydrates::CarbohydrateInfo(this);
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

	if(is_metal() && (1 > nheavyatoms_ || is_virtual(1) )) {
		msg << "A metal residue type has a non-metal atom as atom 1." << std::endl;
		checkspass=false;
	}

	if(is_metalbinding() && metal_binding_atoms_.size()==0) {
		msg << "A metal-binding residue has no metal binding atoms listed in its params file (PROPERTIES METALBINDING without METAL_BINDING_ATOMS list)." << std::endl;
		checkspass=false;
	} else if (!is_metalbinding() && metal_binding_atoms_.size()>0) {
		msg << "A residue that has not been declared as a metal-binding residue has metal binding atoms listed in its params file (METAL_BINDING_ATOMS list without PROPERTIES METALBINDING)." << std::endl;
		checkspass=false;
	}

	if(is_alpha_aa_ && is_beta_aa_) {
		msg << "Error!  A residue type specifies that it is both an alpha and a beta amino acid in its params file." << std::endl;
		checkspass=false;
	}
	if(is_l_aa_ && is_d_aa_) {
		msg << "Error!  A residue type specifies that it is both an L-amino acid and a D-amino acid in its params file." << std::endl;
		checkspass=false;
	}
	if( (backbone_aa_ != core::chemical::aa_unk) && !is_alpha_aa_) {
		msg << "Error!  A residue type specifies a standard alpha amino acid to use as a template for backbone scoring (rama and p_aa_pp scoring functions) without specifying that it is itself an alpha amino acid (PROPERTIES ALPHA_AA)." << std::endl;
		checkspass=false;
	}

	if(!checkspass) {
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
         atom_name_to_vd_      map<string,VD>       add_atom
         atomic_charge          v1<Real>             add_atom
         bonded_neighbor_       v1<v1<int>>          add_bond
         bonded_neighbor_type   v1<v1<BondName>>     add_bond
         atom_base_             v1<int>              set_atom_base
         chi_atoms_             v1<v1<uint>>         add_chi
         nu_atoms_              v1<v1<uint>>         add_nu
         properties             bools                add_property
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


bool
ResidueType::variants_match( ResidueType const & other ) const
{
	if ( ! basic::options::option[ basic::options::OptionKeys::pH::pH_mode ].user() ) {
		for ( Size ii = 1; ii <= variant_types_.size(); ++ii ) {
			if ( ! other.has_variant_type( variant_types_[ ii ] ) ) {
				return false;
			}
		}
		return (variant_types_.size() == other.variant_types_.size());
	}

	//needed for protonated versions of the residues
	else {
		int this_variant_count_offset( 0 );
		for ( Size ii = 1; ii <= variant_types_.size(); ++ii ) {
			if ( variant_types_[ii] == PROTONATED || variant_types_[ii] == DEPROTONATED ) {
				this_variant_count_offset = 1;
				continue;
			}
			if ( ! other.has_variant_type( variant_types_[ ii ] ) ) {
				return false;
			}
		}

		int other_variant_count_offset( 0 );
		if( other.has_variant_type( PROTONATED ) || other.has_variant_type( DEPROTONATED ) ) {
			other_variant_count_offset = 1;
		}

		return ( ( variant_types_.size() - this_variant_count_offset ) ==
				( other.variant_types_.size() - other_variant_count_offset ) );
	}
}

bool
ResidueType::nonadduct_variants_match( ResidueType const & other ) const
{
	int this_variant_count_offset( 0 );
	for ( Size ii = 1; ii <= variant_types_.size(); ++ii ) {
		if ( variant_types_[ii] == ADDUCT ) {
			this_variant_count_offset = 1;
			continue;
		}
		if ( ! other.has_variant_type( variant_types_[ ii ] ) ) {
			return false;
		}
	}

	int other_variant_count_offset( 0 );
	if( other.has_variant_type( ADDUCT ) ) {
		other_variant_count_offset = 1;
	}

	return ( ( variant_types_.size() - this_variant_count_offset ) ==
			( other.variant_types_.size() - other_variant_count_offset ) );
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
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		utility_exit_with_message("unknown atom_name: " + name3() + "  " + name );
	}
	VD const & vd = graph_iter->second;

	Size ordered_index = 0;
	for( Size i = 1; i <= ordered_atoms_.size(); ++i ) {
		if( &graph_[ordered_atoms_[i]] == &graph_[vd]) {
			ordered_index = i;
			break;
		}
	}

	if( ordered_index == 0 ) {
		#if defined BOINC
		// chu temporary graphic fix for boinc
		if ( name == "CA" && !is_protein() ) return 1;
		#endif
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		tr.Error << "printing ordered_atomns names" << std::endl;
		for(Size index=1; index <= ordered_atoms_.size(); ++index){
			tr.Error << graph_[ordered_atoms_[index]].name() << " " << &graph_[ordered_atoms_[index]] << std::endl;
		}
		tr.Error << "vd memory address: " << &graph_[vd] << std::endl;
		utility_exit_with_message("unknown atom_name: " + name3() + "  " + name );
	}

	return ordered_index;
}

/// @brief get atom index by VD
Size
ResidueType::atom_index( VD const & vd) const{
	/// NOTE: Currently we have to iterate twice because atom_name_to_vd_ stores vertex_descriptors not indices.
	/// A substantial change to the interface will fix this but everyone's code will need to switch to
	if( vd == ResidueGraph::null_vertex()) return 0;

	for( core::Size ii(1); ii <= ordered_atoms_.size(); ++ii ) {
		if( ordered_atoms_[ii] == vd ) { return ii; }
	}
	tr.Error << "VD" << vd << " not available in residue " << name3() << std::endl;
	show_all_atom_names( tr.Error );
	tr.Error << std::endl;
	utility_exit_with_message("unknown vertex descriptor in " + name3() );
	return 0;
}

VD ResidueType::vd_from_index(Size const & atomno) const{

	if( ! atomno ) return ResidueGraph::null_vertex();
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
	if(!has(name)){
		utility_exit_with_message("Trying to set bb atom that does not exist in residuetype");
	}
	force_bb_.push_back( atom_name_to_vd_[name]);
}

/// @brief AtomICoord of an atom
AtomICoor const &
ResidueType::icoor( Size const atm ) const
{
	return icoor_.find(ordered_atoms_[atm])->second;

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

///@details update actcoord
/** average geometrical center of the set of actcoord_atoms_ */
void
ResidueType::update_actcoord( conformation::Residue & rot ) const
{
	rot.actcoord().zero();
	if ( n_actcoord_atoms_ > 0 ) {
		for ( Size ii = 1; ii <= n_actcoord_atoms_; ++ii )
		{
			rot.actcoord() += rot.atoms()[ vd_to_index_.find(actcoord_atoms_[ ii ])->second ].xyz();
		}
		rot.actcoord() /= n_actcoord_atoms_;
	}
}

/// @details set AtomICoor for an atom
///
/// will update the xyz coords as well if desired, useful inside a patching operation where new
/// atoms are being added.
void
ResidueType::set_icoor(
		Size const & index,
		std::string const & atm,
		Real const phi,
		Real const theta,
		Real const d,
		std::string const & stub_atom1,
		std::string const & stub_atom2,
		std::string const & stub_atom3,
		bool const update_xyz /* = false*/)
{
	ICoorAtomID id( atm, *this );
	AtomICoor const ic( index, phi, theta, d, stub_atom1, stub_atom2, stub_atom3, *this );
	Size atomno;

	switch ( id.type() ) {
		case ICoorAtomID::INTERNAL:
			atomno = vd_to_index_.find(id.vertex())->second;
			if ( graph_.num_vertices() < atomno ) utility_exit_with_message("ResidueType:: shoudnt get here!");//icoor_.resize(atomno);
			if ( ordered_atoms_.size() < atomno ) utility_exit_with_message("ResidueType:: shoudnt get here!");//icoor_.resize(atomno);
			//graph_[ordered_atoms_[ atomno ]].icoor( ic );
			icoor_[ ordered_atoms_[atomno] ] = ic;

			// update atom_base?
			if ( ( stub_atom1 != atm ) && has( stub_atom1 ) &&
					( atom_base_.find(ordered_atoms_[atomno]) == atom_base_.end() ||
							vd_to_index_.find((atom_base_.find(ordered_atoms_[atomno])->second))->second == atomno ) ) {
				set_atom_base( atm, stub_atom1 );
			}
			if ( update_xyz ) {
				set_ideal_xyz( atm, ic.build( *this ) );
				//std::cout << "building coords for atm " << name_ << ' ' << atm << ' ' <<
				//		ic.build(*this)(1) << ' ' <<
				//		ic.build(*this)(2) << ' ' <<
				//		ic.build(*this)(3) << std::endl;
			}
			break;
		case ICoorAtomID::CONNECT:
			atomno = vd_to_index_.find(id.vertex())->second;
			residue_connections_[ atomno ].icoor( ic );
			break;
		case ICoorAtomID::POLYMER_LOWER:
			assert( lower_connect_id_ != 0 );
			residue_connections_[ lower_connect_id_ ].icoor( ic );
			break;
		case ICoorAtomID::POLYMER_UPPER:
			assert( upper_connect_id_ != 0 );
			residue_connections_[ upper_connect_id_ ].icoor( ic );
			break;
		default:
			utility_exit_with_message( "unrecognized stub atom id type!" );
	}
}

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
	AtomICoor const ic( phi, theta, d, stub_atom1, stub_atom2, stub_atom3, *this );
	Size atomno;

	switch ( id.type() ) {
		case ICoorAtomID::INTERNAL:
			atomno = vd_to_index_.find(id.vertex())->second;
			if ( ordered_atoms_.size() < atomno ) utility_exit_with_message("ResidueType:: shoudnt get here!");//icoor_.resize(atomno);
			icoor_[ ordered_atoms_[atomno] ] = ic;
			//graph_[ordered_atoms_[ atomno ]].icoor( ic );
			// update atom_base?
			if ( ( stub_atom1 != atm ) && has( stub_atom1 ) &&
					(  atom_base_.find(ordered_atoms_[atomno]) == atom_base_.end() ||
							vd_to_index_.find((atom_base_.find(ordered_atoms_[atomno])->second))->second == atomno) ) {
				set_atom_base( atm, stub_atom1 );
			}
			if ( update_xyz ) {
				set_ideal_xyz( atm, ic.build( *this ) );
				//std::cout << "building coords for atm " << name_ << ' ' << atm << ' ' <<
				//		ic.build(*this)(1) << ' ' <<
				//		ic.build(*this)(2) << ' ' <<
				//		ic.build(*this)(3) << std::endl;
			}
			break;
		case ICoorAtomID::CONNECT:
			atomno = vd_to_index_.find(id.vertex())->second;
			residue_connections_[ atomno ].icoor( ic );
			break;
		case ICoorAtomID::POLYMER_LOWER:
			assert( lower_connect_id_ != 0 );
			residue_connections_[ lower_connect_id_ ].icoor( ic );
			break;
		case ICoorAtomID::POLYMER_UPPER:
			assert( upper_connect_id_ != 0 );
			residue_connections_[ upper_connect_id_ ].icoor( ic );
			break;
		default:
			utility_exit_with_message( "unrecognized stub atom id type!" );
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


void ResidueType::set_RotamerLibraryName( std::string const & filename )
{
	rotamer_library_name_ = filename;
}

/// @brief A residue parameter file can refer to a set of "pdb rotamers" that can be
/// superimposed onto a starting position for use in the packer.  These rotamers
/// are loaded into the pack::dunbrack::RotamerLibrary at the time of their first use.
std::string ResidueType::get_RotamerLibraryName() const
{
	return rotamer_library_name_;
}


void ResidueType::assign_neighbor_atom()
{
	//calculate the geometric center of all atoms in the residue

	Vector total(0.0,0.0,0.0);
	for(core::Size index = 1; index <= ordered_atoms_.size();++index)
	{
		total += graph_[ordered_atoms_[index]].ideal_xyz();
	}

	Vector center = total/ordered_atoms_.size();

	//locate the atom which is closest to the center
	core::Size min_index = 0;
	core::Real min_distance = 50000.0;

	for(core::Size index = 1; index <= ordered_atoms_.size();++index)
	{
		core::Real distance = center.distance(graph_[ordered_atoms_[index]].ideal_xyz());
		if( (distance < min_distance) && (!atom_is_hydrogen(index)) )
		{
			min_distance = distance;
			min_index = index;
		}
	}
	assert(min_index != 0);
	//set neighbor atom
	nbr_atom(graph_[ordered_atoms_[min_index]].name());
}

void ResidueType::assign_internal_coordinates()
{
	assign_internal_coordinates( atom_name(nbr_atom()) );
	AtomICoor nbr_icoor = icoor(nbr_atom());
	set_icoor(atom_name(nbr_atom()),0.0,0.0,0.0,atom_name(nbr_icoor.stub_atom1().atomno()),
			atom_name(nbr_icoor.stub_atom2().atomno()),atom_name(nbr_icoor.stub_atom3().atomno()));

}

void ResidueType::assign_internal_coordinates(std::string const & current_atom)
{
	//%TODO: right now i'm ignoring M FRAG lines and M SPLT lines in molfiles
	core::Size current_atom_index = atom_index(current_atom);
	std::string parent_stub1;
	std::string parent_stub2;
	std::string parent_stub3;

	//the root atom has dummy stubs and icoords of 0
	if(current_atom_index == nbr_atom())
	{
		core::Size first_child = bonded_neighbor(current_atom_index)[1];
		parent_stub1 = atom_name(current_atom_index);
		parent_stub2 = atom_name(first_child);

		if(bonded_neighbor(first_child).size() >	0)
		{
			parent_stub3 = atom_name(bonded_neighbor(first_child)[1]);
		}else
		{
			parent_stub3 = atom_name(bonded_neighbor(current_atom_index)[2]);
		}
	}else
	{
		core::Size parent_index = vd_to_index_.find( parents_[ current_atom_index ] )->second;
		AtomICoor parent_icoor = icoor(parent_index);
		parent_stub1 = atom_name(current_atom_index);
		parent_stub2 = atom_name(parent_index);
		parent_stub3 = atom_name(parent_icoor.stub_atom2().atomno());
	}

	std::string previous_sibling = parent_stub3;
	AtomIndices children = bonded_neighbor(current_atom_index);
	for(core::Size index = 1; index <children.size();++index)
	{
		core::Size child_index = children[index];

		if(( vd_to_index_.find(parents_[child_index])->second != 0) && (child_index != nbr_atom()))
		{
			continue;
		}

		std::string child_stub1 = parent_stub1;
		std::string child_stub2 = parent_stub2;
		std::string child_stub3 = previous_sibling;
		parents_[child_index] = ordered_atoms_[current_atom_index];
		if((current_atom_index == nbr_atom()) && (previous_sibling == parent_stub2))
		{
			child_stub3 = parent_stub3;
		}
		calculate_icoor(atom_name(child_index),child_stub1,child_stub2,child_stub3);
		//set_atom_base(atom_name(child_index),)
                		assign_internal_coordinates(atom_name(child_index) );
		previous_sibling = atom_name(child_index);
	}
}

void ResidueType::calculate_icoor(std::string const & child,
		std::string const & stub_atom1,
		std::string const & stub_atom2,
		std::string const & stub_atom3)
{
	//std::cout <<child << " \""<<stub_atom1 << "\" \""<<stub_atom2<< "\" \""<<stub_atom3 << std::endl;
	// This is basically a direct port of calc_internal_coords()
	// found in /python/apps/public/molfile_to_params.py
	Vector const child_xyz = atom(atom_index(child)).ideal_xyz();
	Vector const stub1_xyz = atom(atom_index(stub_atom1)).ideal_xyz();
	Vector const stub2_xyz = atom(atom_index(stub_atom2)).ideal_xyz();
	Vector const stub3_xyz = atom(atom_index(stub_atom3)).ideal_xyz();

	core::Real distance = child_xyz.distance(stub1_xyz);
	core::Real theta = 0.0;
	core::Real phi = 0.0;
	if(distance <1e-2)
	{
		tr << "WARNING: extremely small distance=" << distance << " for " <<
				child << " ,using 0.0 for theta and phi."<<
				" If you were not expecting this warning, something is very wrong" <<std::endl;
	}else
	{
		theta = numeric::angle_radians<core::Real>(child_xyz,stub1_xyz,stub2_xyz);
		if( (theta < 1e-2) || (theta > numeric::NumericTraits<Real>::pi()-1e-2) )
		{
			phi = 0.0;
		}else
		{
			phi = numeric::dihedral_radians<core::Real>(child_xyz,stub1_xyz,stub2_xyz,stub3_xyz);
		}

	}
	//tr << child << " " << stub_atom1 << " "<< stub_atom2 << " " <<stub_atom3 << " " <<distance << " " << phi << " " << theta <<std::endl;
	set_icoor(child,phi,theta,distance,stub_atom1,stub_atom2,stub_atom3);
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
ResidueType::set_shadowing_atom(
		std::string const & atom,
		std::string const & atom_being_shadowed
)
{

	VD const index_shadower( ordered_atoms_[atom_index(atom)] );
	VD const index_shadowee( ordered_atoms_[ atom_index(atom_being_shadowed) ]);
	atom_shadowed_[ index_shadower ] = index_shadowee;
}

void
ResidueType::set_ideal_xyz(
		Size index,
		Vector const & xyz_in
)
{
	if ( index > ordered_atoms_.size() ) ordered_atoms_.resize(index); ///TODO REMOVE THIS!!! This would create NULL vertex_descriptors
	Atom & a = graph_[ ordered_atoms_[index] ];
	a.ideal_xyz( xyz_in );
}

////////// Utility functions for retype_atoms

/// @brief Should the element be considered to be a virtual atom?
bool
retype_is_virtual( std::string const & element ) {
	return (element == "*" || element == "X" || element == "V"); // TODO: Permit Vandium to be an actual element
}

std::string
retype_get_element(VD const & vd, Atom const & a, ElementMap const & emap, AtomTypeSet const & atom_type_set ) {
	ElementMap::const_iterator emiter( emap.find(vd) );
	//std::string element( emap[ vd ]  );
	if( emiter == emap.end() ) {
		if( a.atom_type_index() != 0 ) {
			// Assume we're keeping the same element.
			return atom_type_set[ a.atom_type_index() ].element();
		} else {
			utility_exit_with_message("Cannot retype atoms - element unknown.");
		}
	} else {
		return emiter->second;
	}
}

/// @brief An atom is aromatic if it has any aromatic bonds to a non-virtual atom.
/// TODO: We need better aromatic ring detection.
bool
retype_is_aromatic(VD const & atom, ResidueGraph const & graph, ElementMap const & emap, AtomTypeSet const & atom_type_set ) {
	OutEdgeIter bonds, bonds_end;
	for( boost::tie(bonds, bonds_end) = boost::out_edges(atom,graph); bonds != bonds_end; ++bonds ) {
		if( graph[ *bonds ].bond_name() == AromaticBond ) {
			VD const & tvd( boost::target( *bonds, graph) );
			Atom const & t( graph[tvd] );
			std::string t_element( retype_get_element(tvd,t,emap,atom_type_set) );
			if( ! retype_is_virtual( t_element ) ) {
				return true;
			}
		}
	}
	return false;
}

/// @brief Reassign Rosetta atom types based on the current heuristics.
/// emap is a map of VD->element strings. If an atom is not present in the element map,
/// attempt to get the element string from the current type (it's an error if it doesn't have one.)
/// If preserve is true, only retype those atoms which have an atom_type_index of zero.
/// @details The logic here comes from molfile_to_params.py
/// Which is itself based on Rosetta++ ligand_ns.cc set_rosetta_atom_types(),
/// and has been validated against the Meiler and Baker 2006 cross docking test set
/// assignments
///
/// I'm not saying the logic is good, but it's the logic we're using.
///
/// This function assumes that:
///   * All bonds and atoms exist.
///   * Bond types (bond_name) are correctly set
///   * The element symbols are either provided in emap, or are available through the currently set types.
void
ResidueType::retype_atoms(ElementMap const & emap, bool preserve) {
	// For each atom, analyze bonding pattern to determine type
	VDs aroCs; // Atoms assigned as aroC - need to change all attached hydrogens to Haro.
	VIter itr, itr_end;
	for( boost::tie(itr, itr_end) = vertices(graph_); itr != itr_end; ++itr) {
		Atom & a( graph_[*itr] );
		if( preserve && a.atom_type_index() != 0 ) {
			continue;
		}
		std::string element( retype_get_element(*itr,a,emap,*atom_types_)  );
		// H, C, O, N have complicated rules.
		// Everything else maps to a single atom type.
		if ( retype_is_virtual( element ) ) {
			a.atom_type_index( atom_types_->atom_type_index("VIRT") );
		} else if( element == "H" ) {
			OutEdgeIter bonds, bonds_end;
			core::Size num_aro_C(0), num_NOS(0);
			for( boost::tie(bonds, bonds_end) = boost::out_edges(*itr,graph_); bonds != bonds_end; ++bonds ) {
				VD const & tvd( boost::target( *bonds, graph_) );
				Atom const & t( graph_[tvd] );
				std::string t_element( retype_get_element(tvd,t,emap,*atom_types_) );

				if( t_element == "N" || t_element == "O" || t_element == "S" ) { ++num_NOS; }
				// Instead of also counting number of aroC's here (which may depend on atom iteration ordering,
				// we annotate the ones we've assigned, and then adjust them afterwards.
				// We still include the following test here, though, as it may catch hydrogens on carbons which
				// don't get the aroC label.
				if( t_element == "C" && retype_is_aromatic(tvd,graph_,emap,*atom_types_) ) {
					++num_aro_C;
				}
			}
			if( num_NOS >=1 ) {
				a.atom_type_index( atom_types_->atom_type_index("Hpol") );
			} else if ( num_aro_C >= 1 ) {
				a.atom_type_index( atom_types_->atom_type_index("Haro") );
			} else {
				a.atom_type_index( atom_types_->atom_type_index("Hapo") );
			}
		} else if( element == "C") {
			OutEdgeIter bonds, bonds_end;
			bool saturated(true);
			core::Size num_H(0), num_dbl_nonO(0), num_aro_nonO(0), num_aro_N(0);
			for( boost::tie(bonds, bonds_end) = boost::out_edges(*itr,graph_); bonds != bonds_end; ++bonds ) {
				VD const & tvd( boost::target( *bonds, graph_) );
				Atom const & t( graph_[tvd] );
				std::string t_element( retype_get_element(tvd,t,emap,*atom_types_) );

				if( retype_is_virtual(t_element) ) { continue; }
				switch( graph_[*bonds].bond_name() ) {
					case SingleBond:
						if( t_element == "H" ) { ++num_H; }
						break;
					case DoubleBond:
						saturated = false;
						if ( t_element != "O" ) { ++num_dbl_nonO; }
						break;
					case TripleBond:
						saturated = false;
						break;
					case AromaticBond:
						saturated = false;
						if ( t_element != "O" ) { ++num_aro_nonO; }
						if ( t_element == "N" ) { ++num_aro_N; } // really if, not else if
						break;
					default:
						break;
				}
			}
			if( saturated ) {
				if( num_H >= 3 ) {
					a.atom_type_index( atom_types_->atom_type_index("CH3 ") );
				} else if( num_H == 2 ) {
					a.atom_type_index( atom_types_->atom_type_index("CH2 ") );
				} else {
					a.atom_type_index( atom_types_->atom_type_index("CH1 ") );
				}
			} else { // unsaturated
				if( num_aro_nonO >= 2 ) {
					a.atom_type_index( atom_types_->atom_type_index("aroC") );
					aroCs.push_back( *itr ); // for later attached H annotation
				} else if( num_dbl_nonO >= 1 ) {
					a.atom_type_index( atom_types_->atom_type_index("aroC") );
					aroCs.push_back( *itr ); // for later attached H annotation
				} else if( num_aro_N >= 1 ) {
					a.atom_type_index( atom_types_->atom_type_index("CNH2") );
				} else {
					a.atom_type_index( atom_types_->atom_type_index("COO ") );
				}
			}
		} else if( element == "N" ) {
			OutEdgeIter bonds, bonds_end;
			bool saturated(true);
			core::Size num_H(0), heavy_nbrs(0);
			for( boost::tie(bonds, bonds_end) = boost::out_edges(*itr,graph_); bonds != bonds_end; ++bonds ) {
				VD const & tvd( boost::target( *bonds, graph_) );
				Atom const & t( graph_[tvd] );
				std::string t_element( retype_get_element(tvd,t,emap,*atom_types_) );

				if( retype_is_virtual(t_element) ) { continue; }
				if( t_element == "H" ) { ++num_H; }
				else { ++heavy_nbrs; } // We've already ignored all the virtual atoms.
				if( graph_[*bonds].bond_name() != SingleBond ) { saturated = false; }
			}

			if( num_H >= 3 ) {
				a.atom_type_index( atom_types_->atom_type_index("Nlys") ); // carries a VERY high desolvation penalty
			} else if( num_H == 2 ) {
				// Not totally sure about this one, may want Ntrp instead if more than one heavy neighbor:
				a.atom_type_index( atom_types_->atom_type_index("NH2O") ); // Narg would also be a possibility, but they're fairly similar
			} else if( num_H == 1 ) {
				if( heavy_nbrs <= 2 ) {
					a.atom_type_index( atom_types_->atom_type_index("Ntrp") ); // should always be 2 neighbors, not less
				} else {
					a.atom_type_index( atom_types_->atom_type_index("Ntrp") ); // Npro? protonated tertiary amine
				} // I know they're the same -- I'm just copying molfile_to_params, which splits the case.
			} else {
				if( heavy_nbrs <= 2 ) {
					a.atom_type_index( atom_types_->atom_type_index("Nhis") );
				} else if ( heavy_nbrs == 3 ) {
					if( saturated ) {
						a.atom_type_index( atom_types_->atom_type_index("Nhis") ); // deprotonated tertiary amine; need an sp3 hybrid H-bond acceptor type...
					} else { // This also catches nitro groups -- is that what we want here?
						a.atom_type_index( atom_types_->atom_type_index("Npro") ); // X=[N+](X)X, including nitro groups
					}
				} else {
					a.atom_type_index( atom_types_->atom_type_index("Npro") ); // quaternary amine
				}
			}
		} else if( element == "O" ) {
			OutEdgeIter bonds, bonds_end;
			bool saturated(true);
			core::Size num_H(0), num_bonds(0), bonded_to_N(0), bonded_to_C_to_N(0), unsat_nbrs(0);
			for( boost::tie(bonds, bonds_end) = boost::out_edges(*itr,graph_); bonds != bonds_end; ++bonds ) {
				VD const & tvd( boost::target( *bonds, graph_) );
				Atom const & t( graph_[tvd] );
				std::string t_element( retype_get_element(tvd,t,emap,*atom_types_) );

				if( retype_is_virtual(t_element) ) { continue; }
				++num_bonds; // Bonds to non-virtual atoms.
				if( graph_[*bonds].bond_name() != SingleBond ) { saturated = false; }
				if( t_element == "H" ) { ++num_H; }
				else if( t_element == "N" ) { ++bonded_to_N; }
				OutEdgeIter bonds2, bonds_end2; // second degree bonds.
				bool sat_neighbor = true;
				for( boost::tie(bonds2, bonds_end2) = boost::out_edges(tvd,graph_); bonds2 != bonds_end2; ++bonds2 ) {
					// Ignore the bond back to the atom we're typing.
					//if( boost::target( *bonds2, graph_) == *itr ) { continue; }
					VD const & tvd2( boost::target( *bonds2, graph_) );
					Atom const & t2( graph_[tvd2] );
					std::string t2_element( retype_get_element(tvd2,t2,emap,*atom_types_) );

					if( retype_is_virtual(t2_element) ) { continue; }
					if( t_element == "C" && t2_element == "N") { ++bonded_to_C_to_N; }
					if( graph_[*bonds2].bond_name() != SingleBond ) { sat_neighbor = false; }
				}
				if( ! sat_neighbor ) { ++unsat_nbrs; }
			}
			if( saturated ) {
				if( num_bonds < 2 ) {
					a.atom_type_index( atom_types_->atom_type_index("OOC ") ); // catches C(=O)[O-] (Kekule form) -- new rule by IWD
				} else {
					core::Size ring_size( smallest_ring_size( *itr ) );
					if( num_H > 0 ) {
						a.atom_type_index( atom_types_->atom_type_index("OH  ") ); // catches C(=O)OH (Kekule form)
					} else if ( ring_size < 5 ) {
						a.atom_type_index( atom_types_->atom_type_index("OH  ") ); // small, strained rings leave the O more exposed? (IWD, see 1p8d)
					} else if ( ring_size < 999999 && unsat_nbrs > 0 ) {
						a.atom_type_index( atom_types_->atom_type_index("Oaro") ); // catches aromatic O in furan-like rings, though I rarely see these H-bond (IWD)
					} else {
						a.atom_type_index( atom_types_->atom_type_index("OH  ") ); // catches ethers, ROR (IWD, see comment)
						// The lone pairs on ethers are capable of H-bonding in the same way that alcohols are.
						// While alkyl ethers are quite non-polar, many others seem to make Hbonds,
						// such as those attached to phosphates (R-O-PO3), methyls (R-O-CH3), and aromatic rings (R-O-Ph).
						// It is unclear from the literature how strong these are, and is probably very situation dependent.
					}
				}
			} else if ( num_H > 0 ) {
				a.atom_type_index( atom_types_->atom_type_index("OH  ") ); // catches c(o)oH (aromatic bonds to both O)
			} else if ( bonded_to_N ) {
				a.atom_type_index( atom_types_->atom_type_index("ONH2") );
			} else if ( bonded_to_C_to_N ) { // This is a non-standard rule introduced by IWD, agreed to by KWK:
				a.atom_type_index( atom_types_->atom_type_index("ONH2") );
			} else {
				a.atom_type_index( atom_types_->atom_type_index("OOC ") );
			}
		} else if( "S"  == element ) {
			a.atom_type_index( atom_types_->atom_type_index("S   ") );
		} else if( "P"  == element ) {
			a.atom_type_index( atom_types_->atom_type_index("Phos") );
		} else if( "F"  == element ) {
			a.atom_type_index( atom_types_->atom_type_index("F   ") );
		} else if( "CL" == element ) {
			a.atom_type_index( atom_types_->atom_type_index("Cl  ") );
		} else if( "BR" == element ) {
			a.atom_type_index( atom_types_->atom_type_index("Br  ") );
		} else if( "I"  == element ) {
			a.atom_type_index( atom_types_->atom_type_index("I   ") );
		} else if( "NA" == element ) {
			a.atom_type_index( atom_types_->atom_type_index("Na1p") );
		} else if( "K"  == element ) {
			a.atom_type_index( atom_types_->atom_type_index("K1p ") );
		} else if( "MG" == element ) {
			a.atom_type_index( atom_types_->atom_type_index("Mg2p") );
		} else if( "FE" == element ) {
			a.atom_type_index( atom_types_->atom_type_index("Fe3p") );
		} else if( "CA" == element ) {
			a.atom_type_index( atom_types_->atom_type_index("Ca2p") );
		} else if( "ZN" == element ) {
			a.atom_type_index( atom_types_->atom_type_index("Zn2p") );
		} else {
			utility_exit_with_message("Cannot type atom with element '"+element+"'");
		}
	} // For vertices in graph

	// Hydrogens attached to aroCs == Haro.
	// Technically doesn't match molfile_to_params, as a hydrogen simultaneously bonded to an aroC and an N/O/S
	// would be typed Hpol there, but is typed Haro here. Though if you're silly enough to make two bonds to a hydrogen,
	// you really can't complain when things come out mucked up.
	for( VDs::const_iterator aroit(aroCs.begin()), aroend(aroCs.end()); aroit != aroend; ++aroit ) {
		OutEdgeIter bonds, bonds_end;
		for( boost::tie(bonds, bonds_end) = boost::out_edges(*aroit,graph_); bonds != bonds_end; ++bonds ) {
			VD const & tvd( boost::target( *bonds, graph_) );
			Atom & t( graph_[tvd] );
			std::string t_element( retype_get_element(tvd,t,emap,*atom_types_) );
			if( t_element == "H" && graph_[ *bonds ].bond_name() != UnknownBond  ) {
				if( preserve && t.atom_type_index() != 0 ) { continue; }
				t.atom_type_index( atom_types_->atom_type_index("Haro") );
			}
		}
	}
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
	tr.Debug << "START DIHEDRAL ANGLES ATOM NAMES" << std::endl;
	tr.Debug << "Number of dihe: " << ndihe_ << " " << dihedral_atom_sets_.size() << std::endl;
	for ( Size i = 1; i <= ndihe_; ++i )
	{

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
	for ( Size i = 1; i <= bondangle_atom_sets_.size(); ++i )
	{

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
	for ( Size i = 1; i <= natoms(); ++i )
	{
		tr.Debug << "\t" << graph_[ordered_atoms_[i]].name();
	}
	tr.Debug << std::endl;

	for ( Size j = 1; j <= natoms(); ++j )
	{
		tr.Debug << graph_[ordered_atoms_[j]].name() << "\t";
		for ( Size k = 1; k <= natoms(); ++k )
		{
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
			if (atom_base(jj_atom) == iiat3 && iiat3base != jj_atom ) {
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
	assert(  atom_base(atomno) != atomno );

	/// End the recursion: this atom already has had it's last chi identified, and it's not
	/// the chi we're currently labeling atoms with.
	if ( last_controlling_chi_[ atomno ] != 0 && last_controlling_chi_[ atomno ] != chi ) return;

	last_controlling_chi_[ atomno ] = chi;

	AtomIndices const & nbrs(bonded_neighbor(atomno) );
	for ( Size ii = 1; ii <= nbrs.size(); ++ii ) {
		/// descend into atoms who list atomno as their parent;
		/// atom_base_ defines a tree except at the root, where
		/// atom_base_[ atom_base_[ ii ]] == ii
		if (atom_base(nbrs[ii]) == atomno ) {
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
		while( center > nheavyatoms() || bonded_neighbor(center).size() < 2 ) {
			center = atom_base(center);
		}
		AtomIndices const & nbrs( bonded_neighbor(center) );
		// First try to find two neighbors that are heavyatoms
		for( Size j=1; j<= nbrs.size(); ++j ) {
			Size const nbr( nbrs[j] );
			if( nbr <= nheavyatoms() ) {
				if ( nbr1 ) nbr2 = nbr;
				else nbr1 = nbr;
			}
		}
		// Failing that, just try for two neighbors!
		if( !( center && nbr1 && nbr2 ) ) {
			for( Size j=1; j<= nbrs.size(); ++j ) {
				Size const nbr( nbrs[j] );
				if ( nbr1 ) nbr2 = nbr;
				else nbr1 = nbr;
			}
		}
		if( !( center && nbr1 && nbr2 ) ) {
			// assert() isn't enough for these cases b/c they're typically ligands
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
	if( defined_adducts_.size() == 0 ) return;

	for( Size ii = 1 ; ii <= defined_adducts_.size() ; ++ii) {
		Adduct & add( defined_adducts_[ii] );
		tr.Debug << "Residue: " << name3() << " Adduct: " << add.adduct_name() <<
				" Atom name: " << add.atom_name() << std::endl;
	}
}

void
ResidueType::debug_dump_icoor()
{

	tr.Debug << "ICoor for " << name3() << std::endl;
	for( Size ii = 1 ; ii <= natoms() ; ++ii) {
		tr.Debug << " Atom name: " << atom_name( ii ) << " ideal xyz " << atom(ii).ideal_xyz()[0] << "  " << atom(ii).ideal_xyz()[1] << "  " << atom(ii).ideal_xyz()[2] << std::endl;
	}

}


void
ResidueType::show_all_atom_names( std::ostream & out ) const {

	for( VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first){
		VIter v_iter= vp.first;
		VD vd = *v_iter;
		Atom a = graph_[vd];
		out << a.name() << " " << &graph_[vd] << std::endl;
	}

}

void
ResidueType::set_ncaa_rotlib_n_bin_per_rot( utility::vector1<Size> n_bins_per_rot )
{
	assert( ncaa_rotlib_n_rots_ == n_bins_per_rot.size() );
	ncaa_rotlib_n_bins_per_rot_.resize( ncaa_rotlib_n_rots_ );
	for( Size i = 1; i <= ncaa_rotlib_n_rots_; ++i ) {
		ncaa_rotlib_n_bins_per_rot_[i] = n_bins_per_rot[i];
	}
}

/// @brief  Check if atom is virtual.
bool
ResidueType::is_virtual( Size const & atomno ) const
{
	return ( atom_type( atomno ).is_virtual() );
}

/// @brief  Check if residue is 'VIRTUAL_RESIDUE'
bool
ResidueType::is_virtual_residue() const{
	return ( has_variant_type( "VIRTUAL_RESIDUE" ) );
}

///////////////////////////////////////////////////////////////
core::chemical::rna::RNA_ResidueType const &
ResidueType::RNA_type() const{
			return ( *rna_residue_type_ );
		}


/// @author Labonte
void
ResidueType::show( std::ostream & output ) const
{
	using namespace std;

	output << name_ << " (" << name3_ << ", " << name1_ << "):" << endl;

	output << " Properties:";
	Size n_properties = properties_.size();
	for ( uint i = 1; i <= n_properties; ++i ) {
		output << ' ' << properties_[ i ];
	}
	output << endl;

	output << " Main-chain atoms:";
	Size n_mainchain_atoms = mainchain_atoms_indices_.size();
	for ( uint i = 1; i <= n_mainchain_atoms; ++i ) {
		output << ' ' << atom_name( mainchain_atoms_indices_[ i ] );
	}
	output << endl;

	output << " Backbone atoms:  ";
	Size n_bb_atoms = all_bb_atoms_.size();
	for ( uint i = 1; i <= n_bb_atoms; ++i ) {
		output << ' ' << atom_name( all_bb_atoms_[ i ] );
	}
	output << endl;

	output << " Side-chain atoms:";
	Size n_sc_atoms = all_sc_atoms_.size();
	for ( uint i = 1; i <= n_sc_atoms; ++i ) {
		output << ' ' << atom_name( all_sc_atoms_[ i ] );
	}
	output << endl;

	if ( is_carbohydrate_ ) {
		carbohydrate_info_->show(output);
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
