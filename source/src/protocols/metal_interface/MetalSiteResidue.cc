// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/metal_interface/MetalSiteResidue.cc
/// @brief   Stores data that describes a metal-coordinating residue.
/// @details The intended use is for a utility::vector1 of MetalSiteResidues to describe a multiple-residue metal site.  The atom id's make it convenient generalize the process of adding metalsite constraints.
/// @author Bryan Der

// Headers
#include <protocols/metal_interface/MetalSiteResidue.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


typedef numeric::xyzVector<core::Real> point;

namespace protocols {
namespace metal_interface {

MetalSiteResidue::MetalSiteResidue()
{
}

MetalSiteResidue::~MetalSiteResidue()
{
}


core::Size
MetalSiteResidue::get_seqpos() {
	return seqpos_;
}
void
MetalSiteResidue::set_seqpos( core::Size seqpos ) {
	seqpos_ = seqpos;
}

point
MetalSiteResidue::get_ligand_atom_xyz() { // the atom that forms a coordinating bond to metal
	return ligand_atom_xyz_;
}
void
MetalSiteResidue::set_ligand_atom_xyz( point ligand_atom_xyz ) {
	ligand_atom_xyz_ = ligand_atom_xyz;
}

std::string
MetalSiteResidue::get_ligand_atom_name() {
	return ligand_atom_name_;
}
void
MetalSiteResidue::set_ligand_atom_name( std::string ligand_atom_name ) {
	ligand_atom_name_ = ligand_atom_name;
}

core::id::AtomID
MetalSiteResidue::get_ligand_atom_id() { //atom_id( atomno, residue )
	return ligand_atom_id_;
}
void
MetalSiteResidue::set_ligand_atom_id( core::id::AtomID ligand_atom_id ) {
	ligand_atom_id_ = ligand_atom_id;
}

core::id::AtomID
MetalSiteResidue::get_pre_ligand_atom_id() { //one atom nearer to the CA than the coordinating atom, required for angle and dihedral constraints
	return pre_ligand_atom_id_;
}
void
MetalSiteResidue::set_pre_ligand_atom_id( core::id::AtomID pre_ligand_atom_id ) {
	pre_ligand_atom_id_ = pre_ligand_atom_id;
}

core::id::AtomID
MetalSiteResidue::get_pre_pre_ligand_atom_id() { //two atoms nearer to the CA than the coordinating atom, required for dihedral constraints.  Cysteine does not need a pre_pre_ligan_atom_id because the dihedral is free
	return pre_pre_ligand_atom_id_;
}
void
MetalSiteResidue::set_pre_pre_ligand_atom_id( core::id::AtomID pre_pre_ligand_atom_id ) {
	pre_pre_ligand_atom_id_ = pre_pre_ligand_atom_id;
}

std::string
MetalSiteResidue::get_resname() {
	return resname_;
}
void
MetalSiteResidue::set_resname( std::string resname ) {
	resname_ = resname;
}

}//namespace metal_interface
}//namespace protocols
