// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/RDMolToRestype.cc
///
/// @brief A class for creating a residuetype based on a RDKit fragment
//
/// @details This class takes a RDKit molecule and converts it into a ResidueType.
///
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/rdkit/RDKit.fwd.hh>

#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <numeric/xyzVector.hh>
//support type classes
#include <core/chemical/bond_support.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/atomtype_support.hh>

#include <utility/numbers.hh>
#include <basic/Tracer.hh>

#include <rdkit/GraphMol/AtomIterators.h>
#include <rdkit/GraphMol/BondIterators.h>
#include <rdkit/GraphMol/MolOps.h>

namespace core {
namespace chemical {
namespace rdkit {

static basic::Tracer TR("core.chemical.rdkit.RDMolToRestype");

RDMolToRestype::RDMolToRestype( ::RDKit::ROMol const & rdmol) :
	nbr_(utility::get_undefined_size()),
	index_to_vd_(utility::get_undefined_size(), ResidueGraph::null_vertex()),
	rdmol_(rdmol)
{}

MutableResidueTypeOP
RDMolToRestype::generate_restype(VDIndexMapping const & mapping /*= {}*/) {
	ChemicalManager * cm(core::chemical::ChemicalManager::get_instance());
	std::string const tag("fa_standard");
	AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
	ElementSetCOP element_types = cm->element_set("default");
	MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
	orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);

	MutableResidueType restype(atom_types, element_types, mm_atom_types, orbital_types);
	restype.name( "UNKNOWN" );
	restype.name3( "UNK" );
	restype.name1( 'X' );
	restype.interchangeability_group( "UNK" );

	return generate_restype(restype, mapping);
}

MutableResidueTypeOP
RDMolToRestype::generate_restype(MutableResidueType const & orig_restype, VDIndexMapping const & mapping) {
	IndexNameMapping name_mapping;
	for ( VD atm: orig_restype.all_atoms() ) {
		name_mapping[ mapping[ atm ] ] = orig_restype.atom_name( atm );
	}
	return generate_restype(orig_restype, name_mapping);
}

MutableResidueTypeOP
RDMolToRestype::generate_restype(IndexNameMapping const & mapping) {
	ChemicalManager * cm(core::chemical::ChemicalManager::get_instance());
	std::string const tag("fa_standard");
	AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
	ElementSetCOP element_types = cm->element_set("default");
	MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
	orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);

	MutableResidueType restype(atom_types, element_types, mm_atom_types, orbital_types);
	restype.name( "UNKNOWN" );
	restype.name3( "UNK" );
	restype.name1( 'X' );
	restype.interchangeability_group( "UNK" );

	return generate_restype(restype, mapping);
}

MutableResidueTypeOP
RDMolToRestype::generate_restype(MutableResidueType const & orig_restype, IndexNameMapping const & mapping) { // INTERACTIVE RM
	MutableResidueTypeOP residue( new core::chemical::MutableResidueType(
		orig_restype.atom_type_set_ptr(),
		orig_restype.element_set_ptr(),
		orig_restype.mm_atom_types_ptr(),
		orig_restype.orbital_types_ptr()) );

	residue->set_gasteiger_atom_typeset(orig_restype.gasteiger_atom_typeset());

	residue->name(  orig_restype.name() );
	residue->name3( orig_restype.name3() );
	residue->name1( orig_restype.name1() );
	residue->interchangeability_group( orig_restype.interchangeability_group() );

	index_to_vd_.clear();

	debug_assert( rdmol_.getNumConformers() > 0 );

	// Add hydrogens if there aren't any already
	// We assume that all valid indexes in the input structure correspond to the same index in the added H structure
	::RDKit::ROMolOP rdmol( ::RDKit::MolOps::addHs(rdmol_,false,true) ); // Add physical hydrogens with coordinates - both for the "explicit" and "implicit" hydrogens on the atoms

	::RDKit::Conformer const & conf( rdmol->getConformer() ); // Just the first is fine.

	for ( ::RDKit::ROMol::AtomIterator aitr( rdmol->beginAtoms() ), aitr_end( rdmol->endAtoms() ); aitr != aitr_end; ++aitr ) {
		::RDKit::Atom const & atom( *(*aitr) );
		std::string element_name( atom.getSymbol() );
		int charge( atom.getFormalCharge() );

		std::string name = "";
		// Try copying the atom name
		if ( mapping.count( atom.getIdx() ) ) {
			name = mapping[ atom.getIdx() ];
		}

		VD atom_vd = residue->add_atom(name);
		Atom & restype_atom( residue->atom( atom_vd ) );
		numeric::xyzVector<core::Real> xyz_coords;
		::RDGeom::Point3D const & pos( conf.getAtomPos( atom.getIdx() ) );
		xyz_coords.x() = pos.x;
		xyz_coords.y() = pos.y;
		xyz_coords.z() = pos.z;

		restype_atom.element_type( residue->element_set().element(element_name) );
		restype_atom.formal_charge( charge );
		restype_atom.ideal_xyz( xyz_coords );
		restype_atom.mm_name( "VIRT" ); // We need to do better on this typing.

		index_to_vd_[ atom.getIdx() ] = atom_vd;
	}

	for ( ::RDKit::ROMol::BondIterator bitr( rdmol->beginBonds() ), bitr_end( rdmol->endBonds() ); bitr != bitr_end; ++bitr ) {
		::RDKit::Bond const & bond( *(*bitr) );
		core::chemical::BondName bond_name( convert_from_rdkit_bondtype(bond.getBondType()) );

		residue->add_bond(
			index_to_vd_[ bond.getBeginAtom()->getIdx() ],
			index_to_vd_[ bond.getEndAtom()->getIdx() ],
			bond_name
		);
	}

	core::chemical::rename_atoms(*residue, /*preserve=*/true);
	rosetta_retype_fullatom(*residue, false); //need to do this, fails currently
	rosetta_recharge_fullatom(*residue);
	core::chemical::find_bonds_in_rings( *residue );
	VD nbr_atom = MutableResidueType::null_vertex;
	if ( ! utility::is_undefined( nbr_ ) ) {
		nbr_atom = index_to_vd_[ nbr_ ];
	}
	residue->nbr_radius( find_nbr_dist( *residue, nbr_atom ) );
	residue->nbr_atom( nbr_atom );
	residue->assign_internal_coordinates(); // Also sets atom base. Needs nbr atom assignment
	residue->autodetermine_chi_bonds(); // Must be after internal coordinate setup

	// We want to remove indexes which correspond to the added hydrogens, if any.
	for ( core::Size ii(rdmol_.getNumAtoms()); ii <= rdmol->getNumAtoms(); ++ii ) {
		index_to_vd_.erase(ii);
	}

	return residue;
}



}
}
}
