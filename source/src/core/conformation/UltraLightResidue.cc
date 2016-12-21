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
#include <ObjexxFCL/FArray2D.hh>

namespace core {
namespace conformation {

UltraLightResidue::UltraLightResidue(ResidueCOP residue)
{
	residue_ = residue;
	core::Size resnum = residue->seqpos();
	for ( core::Size atom_index = 1; atom_index <= residue->natoms(); ++atom_index ) {
		coords_.push_back(residue->xyz(atom_index));
		if(!residue->type().atom_is_hydrogen(atom_index))
		{heavy_coords_.push_back(residue->xyz(atom_index));}
		id::AtomID new_atom_id(atom_index,resnum);
		atom_ids_.push_back(new_atom_id);
	}
	center_ =numeric::center_of_mass(coords_);
}

UltraLightResidue::UltraLightResidue(UltraLightResidue const & src) : ReferenceCount(),
		atom_ids_(src.atom_ids_), coords_(src.coords_),heavy_coords_(src.heavy_coords_),residue_(src.residue_),center_(src.center_)

{}

void UltraLightResidue::update_conformation(Conformation & conformation) const
{
	conformation.batch_set_xyz(atom_ids_,coords_);
}

//Transform for Single Ligand
void UltraLightResidue::transform(numeric::xyzMatrix<core::Real> const & rotation_matrix, core::Vector const & translation_vector)
{

	PointPosition new_center = center_ + translation_vector;

	numeric::xyzTransform<core::Real> transformer(numeric::xyzTransform<core::Real>::rot(rotation_matrix,center_,new_center));
	for ( auto & coord : coords_ ) {
		coord = transformer*coord;
	}
	for ( auto & heavy_coord : heavy_coords_ ) {
		heavy_coord = transformer*heavy_coord;
	}
	center_ = numeric::center_of_mass(coords_);

}

//Transform for Ensemble of Ligands
void UltraLightResidue::transform(numeric::xyzMatrix<core::Real> const & rotation_matrix, core::Vector const & translation_vector, core::Vector group_center)
{
	PointPosition new_center = group_center + translation_vector;

	numeric::xyzTransform<core::Real> transformer(numeric::xyzTransform<core::Real>::rot(rotation_matrix,group_center,new_center));

	for ( auto & coord : coords_ ) {
		coord = transformer*coord;
	}
	for ( auto & heavy_coord : heavy_coords_ ) {
		heavy_coord = transformer*heavy_coord;
	}

	center_ = numeric::center_of_mass(coords_);
}

void UltraLightResidue::align_to_residue(UltraLightResidue const & other_residue)
{

	center_ = numeric::center_of_mass(coords_);

	utility::vector1<PointPosition> reference_coords(other_residue.coords_vector()); //FOR RMSD COMPARISON AT THE END
	//utility::vector1<PointPosition> reference_heavy(other_residue.heavy_coords_vector());

	//OLD CODE
	//utility::vector1<core::Real> weights(target_coords.size(),1.0); //weight every atom equally
	//numeric::xyzMatrix<core::Real> rot_matrix;
	//core::Real sigma3 = 0.0; //unused but findUU() needs it
	// numeric::model_quality::findUU(coords_,target_coords,weights,target_coords.size(),rot_matrix,sigma3);
	//OLD CODEE


	utility::vector1<core::PointPosition > target_coords = other_residue.coords_vector();
	utility::vector1<core::PointPosition > coords = coords_;

	//Convert to arrays
	ObjexxFCL::FArray2D< numeric::Real > copy_Farray(numeric::vector_of_xyzvectors_to_FArray<numeric::Real>(coords));
	ObjexxFCL::FArray2D< numeric::Real > target_Farray(numeric::vector_of_xyzvectors_to_FArray<numeric::Real>(target_coords));
	ObjexxFCL::FArray1D< numeric::Real > ww(coords.size(),1.0);
	core::Real bogus = 0;

	//rms_Fit(size, to - unchanged, from - changed, weights, size, bogus)
	numeric::model_quality::rms_fit(coords.size(), target_Farray, copy_Farray, ww, coords.size(), bogus);
	coords = numeric::FArray_to_vector_of_xyzvectors(copy_Farray);
	coords_ = coords;
	center_ = numeric::center_of_mass(coords_);

	numeric::xyzMatrix<core::Real> identity(numeric::xyzMatrix<core::Real>::identity());
	numeric::xyzTransform<core::Real> transformer(identity,other_residue.center());

	for ( auto & coord : coords_ ) {
		coord = transformer*coord;
	}

	for ( auto & heavy_coord : heavy_coords_ ) {
		heavy_coord = transformer*heavy_coord;
	}

	center_ = numeric::center_of_mass(coords_);

	//Calculate RMSD and compare

	core::Real deviation = 0;
	//core::Real heavy = 0;

	for(core::Size i=1; i <= coords_.size(); ++i)
	{
		core::Real deviation_x = ((coords_[i][0]-reference_coords[i][0]) * (coords_[i][0]-reference_coords[i][0]));
		core::Real deviation_y = ((coords_[i][1]-reference_coords[i][1]) * (coords_[i][1]-reference_coords[i][1]));
		core::Real deviation_z = ((coords_[i][2]-reference_coords[i][2]) * (coords_[i][2]-reference_coords[i][2]));

		core::Real total_dev = deviation_x + deviation_y + deviation_z;
		deviation += total_dev;
	}

	deviation /= (core::Real)coords_.size();
	deviation = sqrt(deviation);

	//	std::cout<<"Deviation after conformer change is: " << deviation << std::endl;

	//	for(core::Size i=1; i <= coords_.size(); ++i)
	//	{
	//		heavy += ((coords_[i][0]-reference_heavy[i][0]) * (coords_[i][0]-reference_heavy[i][0]));
	//		heavy += ((coords_[i][1]-reference_heavy[i][1]) * (coords_[i][1]-reference_heavy[i][1]));
	//		heavy += ((coords_[i][2]-reference_heavy[i][2]) * (coords_[i][2]-reference_heavy[i][2]));
	//	}

	//	heavy /= (core::Real)heavy_coords_.size();
	//	heavy = sqrt(heavy);

//		std::cout<<"heavy atom Deviation after conformer change is: " << heavy << std::endl;


}

void UltraLightResidue::slide(core::Vector const & translation_vector)
{
	numeric::xyzMatrix<core::Real> identity(numeric::xyzMatrix<core::Real>::identity());
	numeric::xyzTransform<core::Real> transformer(numeric::xyzTransform<core::Real>::rot(identity,translation_vector));
	for ( auto & coord : coords_ ) {
		coord = translation_vector+coord;
	}
	for(utility::vector1<PointPosition>::iterator it = heavy_coords_.begin(); it != heavy_coords_.end(); ++it)
		{
			*it = translation_vector+(*it);
		}
	center_ = numeric::center_of_mass(coords_);
}


}
}
