// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/grid/ProteinGrid.hh
/// @author Ari Ginsparg

#ifndef INCLUDED_core_grid_ProteinGrid_fwd_hh
#define INCLUDED_core_grid_ProteinGrid_fwd_hh


#include <core/grid/ProteinGrid.fwd.hh>
#include <core/types.hh>

#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/VirtualBase.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/Binary_Util.hh>
#include <core/pose/Pose.hh>

#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

#include <algorithm>
#include <sstream>

namespace core {
namespace grid {

class ProteinGrid : public utility::VirtualBase
{
public:
	/// @brief condensing the data type name for the main 3D vector (or matrix) that represents the protein
	// currently laid out in that the matrix internal data can only be of Size. If flexibility is needed, this could be expanded upon later. The main reason for not doing this now is that operations that use the matrix data rely on the data being a discrete set of positive integers
	typedef utility::vector1<utility::vector1<utility::vector1<core::Size>>> ProteinMatrix;

	//tracer for class prints
	static basic::Tracer ms_tr( "core.grid.ProteinGrid" );

	//parameterized constructors

	//simple constructor that only takes in pose
	ProteinGrid(core::pose::PoseOP in_pose)
	{
		protein_matrix_ = in_pose;

		//wrap matrix around pose
		wrap_matrix_around_pose();
	}

	//constructor that takes in pose and new resolution value
	ProteinGrid(core::pose::PoseOP in_pose, core::Real resolution)
	{
		protein_matrix_ = in_pose;
		resolution_ = resolution;

		//wrap matrix around pose
		wrap_matrix_around_pose();
	}

	//constructor that takes in pose and sub_region_min/max vectors
	ProteinGrid(core::pose::PoseOP in_pose, utility::vector1<core::Size> sub_region_max, utility::vector1<core::Size> sub_region_min)
	{
		protein_matrix_ = in_pose;

		sub_region_max_[1] = sub_region_max[1];
		sub_region_max_[2] = sub_region_max[2];
		sub_region_max_[3] = sub_region_max[3];

		sub_region_min_[1] = sub_region_min[1];
		sub_region_min_[2] = sub_region_min[2];
		sub_region_min_[3] = sub_region_min[3];

		//wrap matrix around pose
		wrap_matrix_around_pose();
	}


	//constructor that takes in pose, resolution values, and subregion vectors
	ProteinGrid(core::pose::PoseOP in_pose, core::Real resolution, utility::vector1<core::Size> sub_region_max, utility::vector1<core::Size> sub_region_min)
	{
		protein_matrix_ = in_pose;
		resolution_ = resolution;

		sub_region_max_[1] = sub_region_max[1];
		sub_region_max_[2] = sub_region_max[2];
		sub_region_max_[3] = sub_region_max[3];

		sub_region_min_[1] = sub_region_min[1];
		sub_region_min_[2] = sub_region_min[2];
		sub_region_min_[3] = sub_region_min[3];

		//wrap matrix around pose
		wrap_matrix_around_pose();
	}

	//copy constructor
	ProteinGrid( ProteinGrid const & other ) 
	{
		other.clone( *this );
	}

	//@brief = operator overload
	ProteinGrid & operator=( ProteinGrid const & other ) {
		other.clone( *this );
		return *this;
	}

	// destructor
	~ProteinGrid() override = default;

	// @brief function to clone the current ProteinGrid into inputted ProteinGrid
	void clone(ProteinGrid & copy) const {
		// copy over gross information
		copy.protein_matrix_ = this->protein_matrix_;
		copy.working_pose_ = this->working_pose_;
		copy.xyz_shift_ = this->xyz_shift_;
		copy.xyz_bound_ = this->xyz_bound_;
		copy.resolution_ = this->resolution_;
		copy.sub_region_max_ = this->sub_region_max_;
		copy.sub_region_min_ = this->sub_region_min_;
	}

	// @brief function to elaborate upon the protein_matrix_, and will review the pose and update occupied cells by projecting atom lennard jobes radii and marking cells within the radius as occupied
	// if a sub area boundary is defined, will define that area with different values

	//(function for determining space fill difference by also including a ligand residuetype)



private:

	// @brief default constructor
	//will need to use class functions to seed values for input pose and other potential input data
	//should only use parameterized or copy constructor
	ProteinGrid();

	// @brief critical function that builds/rebuilds (overwrites existing) the proteingrid protein_matrix_ around the pose in working_pose_
	//called on creation of the object, and also called in cases like changing the pose or resolution
	//this function will simply fill out cells with either a 0 if there is no atom present in the corresponding cell, or a 1 if there is
	//for a more in-depth space fill (project LJ radii), the space fill function will need to be called
	void wrap_matrix_around_pose()
	{
		//run through all atoms to derive a range of dimensions to contain the protein in a 3D  space
		//since we can't have negative indices, we need to normalize the coordinate values so that everything is positive
		//derive constant values based  on the most negative values in each dimension, and then add that constant to all coordinates
		//to be safe, values need to be seeded with an initial value, or errors could be thrown when deriving shift values

		int smallest_x = 1;
		int smallest_y = 1;
		int smallest_z = 1;

		int largest_x = 1;
		int largest_y = 1;
		int largest_z = 1;

		//create a list of coordinates of each atom to hold and work with to fill the protein_representation_matrix
		//can't seem to make a vector of xyzVector objects, so will need to just make a custome 2D vector  to  hold the data
		utility::vector1<numeric::xyzVector<int>> atom_coordinates;

		//determine largest and smallest x,y,z  values to determine dimensions of matrix
		for ( core::Size res_num = 1; res_num <= working_pose_->size(); ++res_num ) {
			for ( core::Size atom_num = 1; atom_num <= working_pose_->residue(res_num).natoms(); ++atom_num ) {
				//get the x,y,z data of the atom, rounded to the closest value
				numeric::xyzVector<int> atom_xyz;
				//floor the coordinates down for a constant negative directional shift
				atom_xyz.x() = std::floor(working_pose_->residue(res_num).xyz(atom_num).x());
				atom_xyz.y() = std::floor(working_pose_->residue(res_num).xyz(atom_num).y());
				atom_xyz.z() = std::floor(working_pose_->residue(res_num).xyz(atom_num).z());

				//safe handling for the first atom encountered to be set as the smallest and largest values
				if ( res_num == 1 && atom_num == 1)
				{
					smallest_x = atom_xyz.x();
					smallest_y = atom_xyz.y();
					smallest_z = atom_xyz.z();
					largest_x = atom_xyz.x();
					largest_y = atom_xyz.y();
					largest_z = atom_xyz.z();
					continue;
				}

				//determine if any of the values  are the smallest
				if ( smallest_x > atom_xyz.x() ) {
					smallest_x = atom_xyz.x();
				}
				if ( smallest_y > atom_xyz.y() ) {
					smallest_y = atom_xyz.y();
				}
				if ( smallest_z > atom_xyz.z() ) {
					smallest_z = atom_xyz.z();
				}

				//determine if any  are the largest
				if ( largest_x < atom_xyz.x() ) {
					largest_x = atom_xyz.x();
				}
				if ( largest_y < atom_xyz.y() ) {
					largest_y = atom_xyz.y();
				}
				if ( largest_z < atom_xyz.z() ) {
					largest_z = atom_xyz.z();
				}

				atom_coordinates.push_back(atom_xyz);

			}
		}

		//take negative values of the smallest values and then add 1 to derive the constants
		//the logic here should apply, whether the smallest value is positive or negative
		//for the smallest value in the system to be indexed to 1, you add the negative of itself + 1; this shift would be applied to all other atom coordinates
		xyz_shift_[1] = std::floor(((smallest_x * -1) + 1) * resolution_);
		xyz_shift_[2] = std::floor(((smallest_y * -1) + 1) * resolution_);
		xyz_shift_[3] = std::floor(((smallest_z * -1) + 1) * resolution_);

		//apply shift values to largest to get boundaries
		xyz_bound_[1] = std::floor((x_shift + largest_x) * resolution_);
		xyz_bound_[2] = std::floor((y_shift + largest_y) * resolution_);
		xyz_bound_[3] = std::floor((z_shift + largest_z) * resolution_);

		//create 3D matrix to roughly represent 3D coordinate space of protein
		ms_tr.Debug << "Creating protein clash coordinate matrix. Dimensions of matrix are " << xyz_bound_[1] << "," << xyz_bound_[2] << "," << xyz_bound_[3] << std::endl;
		ms_tr.Debug << "Shift from from original coordinates, and multiplied by current resolution factor (" << resolution_ << ") are: " << xyz_shift_[1] << "," << xyz_shift_[2] << "," << xyz_shift_[3] << std::endl;

		for ( core::Size x = 1; x <= xyz_bound_[1]; ++x ) {
			//make a 2D  matrix
			utility::vector1<utility::vector1<core::Size>> sub_matrix;

			for ( core::Size y = 1; y <= xyz_bound_[2]; ++y ) {

				//make a 1D matrix, seed with false values
				utility::vector1<core::Size> sub_sub_matrix(xyz_bound_[3],  0);
				//push 1D  matrix into 2D
				sub_matrix.push_back(sub_sub_matrix);

			}
			//push a 2D  matrix into the 3D matrix
			protein_matrix_.push_back(sub_matrix);
		}

		//seed the matrix with approximate coordinates of each atom
		//apply the shift to the coordinates
		//approximated by flooring coordinates down
		for ( core::Size xyzVec = 1; xyzVec <= atom_coordinates.size(); ++xyzVec ) {
			protein_matrix_[atom_coordinates[xyzVec].x() + xyz_shift_[1]][atom_coordinates[xyzVec].y() + xyz_shift_[2]][atom_coordinates[xyzVec].z() + xyz_shift_[3]] = 1;
		}
	}

	// @brief 3D matrix of voxelized representation of atoms in pose; individual indices contain data values that correspond to whether the voxel is occupied by pose atoms
	// the coordinates of pose atoms are used to correspond to voxels in this matrix. Atom coordinates are all shifted so that the minimum coordinate value of all pose atoms are shifted to 1
	ProteinMatrix protein_matrix_;

	// @brief the corresponding pose that the protein_matrix_ is wrapped around
	core::pose::PoseOP working_pose_;

	// @brief vector to hold the values that atom coordinates shift by in the x,y,z directions
	utility::vector1<int> xyz_shift_(3,0);

	// @brief vector to hold the xyz boundaries of the matrix (not the pose coordinates, but the boundaries that are shifted and potentially scaled within the matrix)
	utility::vector1<core::Size> xyz_bound_(3,0);

	// @brief floating value that can be used to scale the resolution of the protein grid, if desired
	// default resolution of voxels are at 1 cubic angstrom (generally recommended)
	//values <1 decrease the resolution of the matrix, increasing the likelihood of multiple atoms appearing in the same voxel
	//values >1 increase the resolution, which can help better define the sphere shape of lj radii from atoms; this comes at a substantial time/memory code
	//resolution value must be positive
	core::Real resolution_ = 1;

	// @brief vectors of maximum and minimum xyz values for a sub-area to potentially investigate in methods like space filling
	// the region is a rectangular prism that is defined by coordinates outlined in sub_refoun_min_ and sub_region_max_
	//seed values to be zero to indicate that we should not be using the subregion unless these values get set to be >=1
	utility::vector1<core::Size> sub_region_max_(3,0);
	utility::vector1<core::Size> sub_region_min_(3,0);
};

}
}


#endif 
