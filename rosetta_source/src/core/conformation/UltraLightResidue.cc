// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/UltraLightResidue.cc
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)
/// @detail Look at protocols/ligand_docking/Transform.cc for a usage example.
/// Basically, this is just a bunch of atoms and a reference to the Residue they originated from.

#include <core/conformation/UltraLightResidue.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/xyz.functions.hh>

namespace core {
namespace conformation {

UltraLightResidue::UltraLightResidue(ResidueCOP residue)
{
	residue_ = residue;
	core::Size resnum = residue->seqpos();
	for(core::Size atom_index = 1; atom_index <= residue->natoms();++atom_index)
	{
		coords_.push_back(residue->xyz(atom_index));
		id::AtomID new_atom_id(atom_index,resnum);
		atom_ids_.push_back(new_atom_id);
	}
	center_ =numeric::center_of_mass(coords_);
}

UltraLightResidue::UltraLightResidue(UltraLightResidue const & src) : atom_ids_(src.atom_ids_), coords_(src.coords_),residue_(src.residue_)
{

}

void UltraLightResidue::update_conformation(Conformation & conformation) const
{
	conformation.batch_set_xyz(atom_ids_,coords_);
}

void UltraLightResidue::transform(numeric::xyzMatrix<core::Real> const & rotation_matrix, core::Vector const & translation_vector)
{
	center_ = numeric::center_of_mass(coords_)+translation_vector;
	numeric::xyzTransform<core::Real> transformer(numeric::xyzTransform<core::Real>::rot(rotation_matrix,center_));
	for(utility::vector1<PointPosition>::iterator it = coords_.begin(); it != coords_.end(); ++it)
	{
		*it = transformer*(*it);
	}
}

void UltraLightResidue::slide(core::Vector const & translation_vector)
{
	numeric::xyzMatrix<core::Real> identity(numeric::xyzMatrix<core::Real>::identity());
	numeric::xyzTransform<core::Real> transformer(numeric::xyzTransform<core::Real>::rot(identity,translation_vector));
	for(utility::vector1<PointPosition>::iterator it = coords_.begin(); it != coords_.end(); ++it)
	{
		*it = transformer*(*it);
	}
}

}
}
