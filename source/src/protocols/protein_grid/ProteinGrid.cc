// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/protein_grid/ProteinGrid.hh
/// @author Ari Ginsparg

#include <protocols/protein_grid/ProteinGrid.fwd.hh>
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

namespace protocols {
namespace protein_grid {


	/// @brief condensing the data type name for the main 3D vector (or matrix) that represents the protein
	// currently laid out in that the matrix internal data can only be of Size. If flexibility is needed, this could be expanded upon later. The main reason for not doing this now is that operations that use the matrix data rely on the data being a discrete set of positive integers
	typedef utility::vector1<utility::vector1<utility::vector1<core::Size>>> ProteinMatrix;

	//tracer for class prints
	static basic::Tracer ms_tr( "core.grid.ProteinGrid" );

	//parameterized constructors

	//simple constructor that only takes in pose
	ProteinGrid::ProteinGrid(core::pose::PoseOP in_pose)
	{
		protein_matrix_ = in_pose;

		//reset the subregion vectors
		//shift and bound vectors will be set/reset in the wrap function
		reset_sub_region_vectors();

		//wrap matrix around pose
		wrap_matrix_around_pose();
	}

	//constructor that takes in pose and new resolution value
	ProteinGrid::ProteinGrid(core::pose::PoseOP in_pose, core::Real resolution)
	{
		protein_matrix_ = in_pose;
		resolution_ = resolution;
		//validate the resolution
		validate_resolution();

		//reset the subregion vectors
		//shift and bound vectors will be set/reset in the wrap function
		reset_sub_region_vectors();

		//wrap matrix around pose
		wrap_matrix_around_pose();
	}

	//constructor that takes in pose and sub_region_min/max vectors
	ProteinGrid::ProteinGrid(core::pose::PoseOP in_pose, utility::vector1<core::Size> sub_region_max, utility::vector1<core::Size> sub_region_min)
	{
		protein_matrix_ = in_pose;

		//reset the subregion vectors
		//shift and bound vectors will be set/reset in the wrap function
		reset_sub_region_vectors();

		//set sub regions
		set_sub_regions(sub_region_max,sub_region_min);

		//wrap matrix around pose
		wrap_matrix_around_pose();
	}

	//constructor that takes in pose, resolution values, and subregion vectors
	ProteinGrid::ProteinGrid(core::pose::PoseOP in_pose, core::Real resolution, utility::vector1<core::Size> sub_region_max, utility::vector1<core::Size> sub_region_min)
	{
		protein_matrix_ = in_pose;
		resolution_ = resolution;
		//validate the resolution
		validate_resolution();

		//reset the subregion vectors
		//shift and bound vectors will be set/reset in the wrap function
		reset_sub_region_vectors();

		//set sub regions
		set_sub_regions(sub_region_max,sub_region_min);

		//wrap matrix around pose
		wrap_matrix_around_pose();
	}

	//copy constructor
	ProteinGrid::ProteinGrid( ProteinGrid const & other ) 
	{
		other.clone( *this );
	}

	//@brief = operator overload
	ProteinGrid::ProteinGrid & operator=( ProteinGrid const & other ) {
		other.clone( *this );
		return *this;
	}

	// destructor
	~ProteinGrid::ProteinGrid() override = default;

	// @brief function to clone the current ProteinGrid into inputted ProteinGrid
	void ProteinGrid::clone(ProteinGrid & copy) const {
		// copy over gross information
		copy.protein_matrix_ = this->protein_matrix_;
		copy.working_pose_ = this->working_pose_;
		copy.xyz_shift_ = this->xyz_shift_;
		copy.xyz_bound_ = this->xyz_bound_;
		copy.resolution_ = this->resolution_;
		copy.sub_region_max_ = this->sub_region_max_;
		copy.sub_region_min_ = this->sub_region_min_;
		copy.matrix_fullness_ = this->matrix_fullness_;
		copy.fullness_ratio_ = this->fullness_ratio_;
		copy.matrix_volume_ = this->matrix_volume_;
	}

	// @brief simple function to derive the volume of the matrix
	core::Size ProteinGrid::get_grid_volume()
	{
		return (xyz_bound_[0] * xyz_bound_[1] * xyz_bound_[2]);
	}

	// @brief function to set the xyz coordinates of the sub_region_min
	// reminder, these values should directly relate to coordinates in the pose, and not be directly aimed at the matrix indices
	void ProteinGrid::set_sub_region_min( utility::vector1<core::Size> region_in )
	{
		sub_region_min_[1] = region_in[1];
		sub_region_min_[2] = region_in[2];
		sub_region_min_[3] = region_in[3];
	}

	// @brief function to set the xyz coordinates of the sub_region_max
	// reminder, these values should directly relate to coordinates in the pose, and not be directly aimed at the matrix indices
	void ProteinGrid::set_sub_region_max( utility::vector1<core::Size> region_in )
	{
		sub_region_max_[1] = region_in[1];
		sub_region_max_[2] = region_in[2];
		sub_region_max_[3] = region_in[3];
	}

	// @brief function to set the xyz coordinates of sub_region_max and sun_region_min
	// reminder, these values should directly relate to coordinates in the pose, and not be directly aimed at the matrix indices
	void ProteinGrid::set_sub_regions( utility::vector1<core::Size> region_max, utility::vector1<core::Size> region_min )
	{
		set_sub_region_max(region_max);
		set_sub_region_min(region_min);
	}

	// @brief function to elaborate upon the protein_matrix_, and will review the pose and update occupied cells by projecting atom lennard jobes radii and marking cells within the radius as occupied
	// if a sub area boundary is defined, will define that area with different values
	void ProteinGrid::project_lj_radii(){
		//iterate over each atom in working_pose
		for ( core::Size res_num = 1; res_num <= working_pose_->size(); ++res_num ) {
			for ( core::Size atom_num = 1; atom_num <= working_pose_->residue(res_num).natoms(); ++atom_num ) {
			


			}
		}
		//
	}

	//(function for determining space fill difference by also including a ligand residuetype)

	// @brief default constructor
	//will need to use class functions to seed values for input pose and other potential input data
	//should only use parameterized or copy constructor
	ProteinGrid::ProteinGrid();

	// @brief critical function that builds/rebuilds (overwrites existing) the proteingrid protein_matrix_ around the pose in working_pose_
	//called on creation of the object, and also called in cases like changing the pose or resolution
	//this function will simply fill out cells with either a 0 if there is no atom present in the corresponding cell, or a 1 if there is
	//for a more in-depth space fill (project LJ radii), the space fill function will need to be called
	ProteinGrid::void wrap_matrix_around_pose()
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

		//set (or reset) the shift and bound vectors to be 3 values of 0 before setting them
		//this is necessary at least for the initial setting when the object is constructed
		reset_xyz_vectors();

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

		//wipe the current contents of the protein_matrix_ and reset fullness values
		//create an empty dummy vector, and then assign protein_matrix_ with it
		ProteinMatrix dummy_matrix;
		protein_matrix_ = dummy_matrix;

		matrix_volume_ = get_grid_volume();

		matrix_fullness_ = 0;
		fullness_ratio_ = 0;

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

			//increment occupied cell count by 1
			++matrix_fullness_;
		}

		//safety check to ensure that we actually had atoms in the pose, and that the matrix has a nonzero volume
		if(matrix_volume_ == 0)
		{
			ms_tr.Warning << "The matrix has no volume, likely because the inputted pose has no atoms!" << std::endl;
		}

		//derive the fullness ratio as the number of occupied cells over the total cell number
		fullness_ratio_ = static_cast<core::Real>(matrix_fullness_) / matrix_volume_;

	}

	// @brief function to check and ensure that the resolution value is valid (> 0)
	// if invalid, will throw a warning and set resolution to 1
	void ProteinGrid::validate_resolution()
	{
		if(resolution_ <= 0)
		{
			ms_tr.Warning << "Invalid resolution scale. Resolution value must be >0. Setting resolution factor to default of 1." << std::endl;
			resolution_ = 1;
		}
	}

	// @brief reset (or set) the xyz shift and bound matrices to be 3 values of zeroes
	void ProteinGrid::reset_xyz_vectors()
	{
		xyz_shift_ = utility::vector1<int>(3, 0);
		xyz_bound_ = utility::vector1<int>(3, 0);
	}

	// @brief reset (or set) the sub-region matrices to be 3 values of zeroes
	void ProteinGrid::reset_sub_region_vectors()
	{
		sub_region_max_ = utility::vector1<int>(3, 0);
		sub_region_min_ = utility::vector1<int>(3, 0);
	}
}
}