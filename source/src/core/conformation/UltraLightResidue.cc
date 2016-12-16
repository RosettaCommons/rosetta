// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/UltraLightResidue.cc
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)
/// @detail Look at protocols/ligand_docking/Transform.cc for a usage example.
/// Basically, this is just a bunch of atoms and a reference to the Residue they originated from.

#include <core/conformation/UltraLightResidue.hh>
#include <core/conformation/Conformation.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>

namespace core {
namespace conformation {

UltraLightResidue::UltraLightResidue(ResidueCOP residue)
{
	residue_ = residue;
	core::Size resnum = residue->seqpos();
	for ( core::Size atom_index = 1; atom_index <= residue->natoms(); ++atom_index ) {
		coords_.push_back(residue->xyz(atom_index));
		id::AtomID new_atom_id(atom_index,resnum);
		atom_ids_.push_back(new_atom_id);
	}
	center_ =numeric::center_of_mass(coords_);
}

UltraLightResidue::UltraLightResidue(UltraLightResidue const & src) : ReferenceCount(),
	atom_ids_(src.atom_ids_), coords_(src.coords_), center_(src.center_), residue_(src.residue_)
{}

void UltraLightResidue::update_conformation(Conformation & conformation) const
{
	conformation.batch_set_xyz(atom_ids_,coords_);
}

void UltraLightResidue::transform(numeric::xyzMatrix<core::Real> const & rotation_matrix, core::Vector const & translation_vector)
{
	PointPosition old_center = numeric::center_of_mass(coords_);
	center_ = old_center+translation_vector;
	numeric::xyzTransform<core::Real> transformer(numeric::xyzTransform<core::Real>::rot(rotation_matrix,old_center,center_));
	for ( auto & coord : coords_ ) {
		coord = transformer*coord;
	}
}

void UltraLightResidue::align_to_residue(UltraLightResidue const & other_residue)
{

	utility::vector1<PointPosition> target_coords(other_residue.coords_vector());
	utility::vector1<core::Real> weights(target_coords.size(),1.0); //weight every atom equally
	numeric::xyzMatrix<core::Real> rot_matrix;
	core::Real sigma3 = 0.0; //unused but findUU() needs it

	numeric::model_quality::findUU(coords_,target_coords,weights,target_coords.size(),rot_matrix,sigma3);

	//Coords are shifted to be centered at origin after findUU

	//original code
	//numeric::xyzTransform<core::Real> transformer(rot_matrix,center_);

	//transform to other residue center

	numeric::xyzTransform<core::Real> transformer(rot_matrix,other_residue.center());

	center_ = numeric::center_of_mass(coords_);

	//coords_ and target_coords get recentered to 0,0,0.  rot_matrix gets set to the correct rotation matrix. sigma3 is set but nobody cares
	//setup the transform, use the last move center to recenter away from zero
	for ( auto & coord : coords_ ) {
		coord = transformer*coord;
	}

	center_ = numeric::center_of_mass(coords_);

}

void UltraLightResidue::slide(core::Vector const & translation_vector)
{
	numeric::xyzMatrix<core::Real> identity(numeric::xyzMatrix<core::Real>::identity());
	numeric::xyzTransform<core::Real> transformer(numeric::xyzTransform<core::Real>::rot(identity,translation_vector));
	for ( auto & coord : coords_ ) {
		coord = translation_vector+coord;
	}
	center_ = numeric::center_of_mass(coords_);
}

}
}
