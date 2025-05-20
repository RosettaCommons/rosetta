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

#include <protocols/protein_grid/ProteinGrid.hh>
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
#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <basic/options/keys/protein_grid.OptionKeys.gen.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <utility/Binary_Util.hh>

#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

#include <algorithm>
#include <sstream>

namespace protocols {
namespace protein_grid {


//key to values that can exist in the ProteinMatrix:
//if the matrix is written to a pdb file, the different states will be translated to different atoms
//notably, even numbers are relegated to define empty space, and odd numbers are relegated to define occupied space
//0 = empty and out of sub area, carbon, black
//1 = pose and out of sub area, fluorine, icy blue
//2 = empty and in sub area, oxygen, red
//3 = pose and in sub area, nitrogen, blue
//4 = do not use, keep even for unoccupied space
//5 = secondary ligand and out of sub area, sulphur, yellow
//6 = do not use, keep even for unoccupied space
//7 = secondary ligand and in sub area, chlorine, green
//8 = do not use, keep even for unoccupied space
//9 = secondary ligand and pose out of sub area, phosphorous, orange
//10 = do not use, keep even for unoccupied space
//11 = secondary ligand and pose in sub area, iodine, purple

//tracer for class prints
static basic::Tracer ms_tr( "core.grid.ProteinGrid" );

//parameterized constructors

//simple constructor that only takes in pose
ProteinGrid::ProteinGrid(core::pose::PoseOP in_pose)
{
	//seed initial values on member variables
	initialize_from_options();

	working_pose_ = in_pose;

	//reset the subregion vectors
	//shift and bound vectors will be set/reset in the wrap function
	reset_sub_region_vectors();

	//wrap matrix around pose
	wrap_matrix_around_pose();
}

//constructor that takes in pose and new resolution value
ProteinGrid::ProteinGrid(core::pose::PoseOP in_pose, core::Real resolution)
{
	//seed initial values on member variables
	initialize_from_options();

	working_pose_ = in_pose;
	resolution_ = resolution;
	//validate the resolution
	validate_resolution();

	//reset the subregion vectors
	//shift and bound vectors will be set/reset in the wrap function
	reset_sub_region_vectors();

	//wrap matrix around pose
	wrap_matrix_around_pose();
}

//constructor that takes in pose and sub region vectors
ProteinGrid::ProteinGrid(core::pose::PoseOP in_pose, numeric::xyzVector<int> sub_area_center, utility::vector1<core::Size> sub_region_dimensions)
{
	//seed initial values on member variables
	initialize_from_options();

	//set the working_pose_ to be the input pose
	working_pose_ = in_pose;

	//reset the subregion vectors
	//shift and bound vectors will be set/reset in the wrap function
	reset_sub_region_vectors();

	//set dimensions
	//throw a warning if the inputted subregion is not 3 entries long, will try to move forward anyway
	if ( sub_region_dimensions.size() != 3 ) {
		ms_tr.Warning << "Sub region dimensions vector that was provided does not have 3 expected values. Contents of vector: [";
		for ( auto val : sub_region_dimensions ) {
			ms_tr.Warning << val << ",";
		}
		ms_tr.Warning << "]" << std::endl;
	}

	//set values of true/absolute sub areas (center and dimensions before shifting)
	true_sub_area_center_ = sub_area_center;
	true_sub_region_dimensions_ = sub_region_dimensions;

	using_sub_area_ = true;

	//wrap matrix around pose
	wrap_matrix_around_pose();
}

//constructor that takes in pose, resolution values, and subregion vectors
ProteinGrid::ProteinGrid(core::pose::PoseOP in_pose, core::Real resolution, numeric::xyzVector<int> sub_area_center, utility::vector1<core::Size> sub_region_dimensions)
{
	//seed initial values on member variables
	initialize_from_options();

	working_pose_ = in_pose;
	resolution_ = resolution;
	//validate the resolution
	validate_resolution();

	//reset the subregion vectors
	//shift and bound vectors will be set/reset in the wrap function
	reset_sub_region_vectors();

	//set dimensions
	//throw a warning if the inputted subregion is not 3 entries long, will try to move forward anyway
	if ( sub_region_dimensions.size() != 3 ) {
		ms_tr.Warning << "Sub region dimensions vector that was provided does not have 3 expected values. Contents of vector: [";
		for ( auto val : sub_region_dimensions ) {
			ms_tr.Warning << val << ",";
		}
		ms_tr.Warning << "]" << std::endl;
	}

	//set values of true/absolute sub areas (center and dimensions before shifting)
	true_sub_area_center_ = sub_area_center;
	true_sub_region_dimensions_ = sub_region_dimensions;

	using_sub_area_ = true;

	//wrap matrix around pose
	wrap_matrix_around_pose();
}

// destructor
ProteinGrid::~ProteinGrid() = default;

// @brief function to clone the current ProteinGrid into inputted ProteinGrid
void ProteinGrid::clone(const ProteinGrid & copy) {
	// copy over all member variables from copy to "this" object
	this->protein_matrix_ = copy.protein_matrix_;
	this->working_pose_ = copy.working_pose_;
	this->xyz_shift_ = copy.xyz_shift_;
	this->xyz_bound_ = copy.xyz_bound_;
	this->matrix_volume_ = copy.matrix_volume_;
	this->sub_matrix_volume_ = copy.sub_matrix_volume_;
	this->resolution_ = copy.resolution_;
	this->using_sub_area_ = copy.using_sub_area_;
	this->using_lj_radii_ = copy.using_lj_radii_;
	this->true_sub_area_center_ = copy.true_sub_area_center_;
	this->adjusted_sub_area_center_ = copy.adjusted_sub_area_center_;
	this->true_sub_region_dimensions_ = copy.true_sub_region_dimensions_;
	this->adjusted_sub_region_dimensions_ = copy.adjusted_sub_region_dimensions_;
	this->sub_region_max_ = copy.sub_region_max_;
	this->sub_region_min_ = copy.sub_region_min_;
	this->matrix_fullness_ = copy.matrix_fullness_;
	this->fullness_ratio_ = copy.fullness_ratio_;
	this->sub_matrix_fullness_ = copy.sub_matrix_fullness_;
	this->sub_fullness_ratio_ = copy.sub_fullness_ratio_;
	this->print_whole_matrix_ = copy.print_whole_matrix_;
	this->print_empty_space_ = copy.print_empty_space_;
}

// @brief function to initialize member variables at constructor calls
//this sets initial and default values for member variables that we can build the matrix upon
void ProteinGrid::initialize_from_options()
{
	//seed variables that are relevant for how to output the matrix to a pdb for visualization
	//not planning to have devoted constructors to define these variables
	print_whole_matrix_ = basic::options::option[ basic::options::OptionKeys::protein_grid::output_whole_matrix ]();
	print_empty_space_ = basic::options::option[ basic::options::OptionKeys::protein_grid::output_empty_space ]();
}

// @brief simple function to derive the volume of the matrix
core::Size ProteinGrid::get_grid_volume()
{
	return (xyz_bound_[1] * xyz_bound_[2] * xyz_bound_[3]);
}

// @brief simple function to derive the volume of the sub area matrix
core::Size ProteinGrid::get_sub_area_grid_volume()
{
	return ((sub_region_max_[1] - sub_region_min_[1] + 1) * (sub_region_max_[2] - sub_region_min_[2] + 1) * (sub_region_max_[3] - sub_region_min_[3] + 1));
}

// @brief simple function to derive the pre-calculated occpuied ratio of the full matrix
core::Real ProteinGrid::get_grid_occupied_cell_ratio()
{
	return fullness_ratio_;
}

// @brief simple function to derive the pre-calculated occpuied ratio of the sub area matrix
core::Real ProteinGrid::get_sub_grid_occupied_cell_ratio()
{
	return sub_fullness_ratio_;
}

// @brief simple function to derive the pre-calculated occpuied cell count of the full matrix
core::Real ProteinGrid::get_grid_occupied_cell_count()
{
	return matrix_fullness_;
}

// @brief simple function to derive the pre-calculated occpuied cell count of the sub area matrix
core::Real ProteinGrid::get_sub_grid_occupied_cell_count()
{
	return sub_matrix_fullness_;
}

// @brief overwrite the true sub area center and dimensions
//makes a call to reset the wrap matrix around pose to account for the change in the sub area
void ProteinGrid::set_sub_regions( numeric::xyzVector<int> sub_area_center, utility::vector1<core::Size> sub_region_dimensions )
{
	using_sub_area_ = true;

	//set center
	true_sub_area_center_ = sub_area_center;

	//set dimensions
	//throw a warning if the inputted subregion is not 3 entries long, will try to move forward anyway
	if ( sub_region_dimensions.size() != 3 ) {
		ms_tr.Warning << "Sub region dimensions vector that was provided does not have 3 expected values. Contents of vector: [";
		for ( auto val : sub_region_dimensions ) {
			ms_tr.Warning << val << ",";
		}
		ms_tr.Warning << "]" << std::endl;
	}


	true_sub_region_dimensions_ = sub_region_dimensions;

	wrap_matrix_around_pose();
}

// @brief function to set up the sub region vectors
// this sets up the following:
//xyz center vector -> fills out the adjusted xyz center vector (shift and resolution)
//xyz dimension vector -> adjusted xyz dimensions vector (resolution)
//xyz min and max values
//leaving this function as public
void ProteinGrid::define_sub_regions()
{


	//first, set the center and then adjust based on the resolution and shift values
	adjusted_sub_area_center_.x() = std::floor((xyz_shift_[1] + true_sub_area_center_.x()) * resolution_);
	adjusted_sub_area_center_.y() = std::floor((xyz_shift_[2] + true_sub_area_center_.y()) * resolution_);
	adjusted_sub_area_center_.z() = std::floor((xyz_shift_[3] + true_sub_area_center_.z()) * resolution_);

	//throw a warning if the adjusted center falls outside the matrix
	if ( adjusted_sub_area_center_.x() > xyz_bound_[1] ) {
		ms_tr.Warning << "Sub region center x outside of matrix. True (inputted): " << true_sub_area_center_.x() << " Adjusted: " << adjusted_sub_area_center_.x() << std::endl;
	}
	//throw a warning if the adjusted center falls outside the matrix
	if ( adjusted_sub_area_center_.y() > xyz_bound_[2] ) {
		ms_tr.Warning << "Sub region center y outside of matrix. True (inputted): " << true_sub_area_center_.y() << " Adjusted: " << adjusted_sub_area_center_.y() << std::endl;
	}
	//throw a warning if the adjusted center falls outside the matrix
	if ( adjusted_sub_area_center_.z() > xyz_bound_[3] ) {
		ms_tr.Warning << "Sub region center z outside of matrix. True (inputted): " << true_sub_area_center_.z() << " Adjusted: " << adjusted_sub_area_center_.z() << std::endl;
	}




	//now, set the adjusted dimensions
	adjusted_sub_region_dimensions_[1] = std::floor(true_sub_region_dimensions_[1] * resolution_);
	adjusted_sub_region_dimensions_[2] = std::floor(true_sub_region_dimensions_[2] * resolution_);
	adjusted_sub_region_dimensions_[3] = std::floor(true_sub_region_dimensions_[3] * resolution_);



	//derive the min and max
	//if the min or max are below 1 or beyond the dimension boundary, set the value to be that boundary

	//max
	//x
	if ( adjusted_sub_area_center_.x() + adjusted_sub_region_dimensions_[1] > xyz_bound_[1] ) {
		//greater than boundary, set to boundary
		sub_region_max_[1] = xyz_bound_[1];
	} else if ( adjusted_sub_area_center_.x() + adjusted_sub_region_dimensions_[1] < 1 ) {
		//smaller than minumum (1), set to 1
		sub_region_max_[1] = 1;
	} else {
		//set normally, as the center + the dimension
		sub_region_max_[1] = adjusted_sub_area_center_.x() + adjusted_sub_region_dimensions_[1];
	}
	//y
	if ( adjusted_sub_area_center_.y() + adjusted_sub_region_dimensions_[2] > xyz_bound_[2] ) {
		sub_region_max_[2] = xyz_bound_[2];
	} else if ( adjusted_sub_area_center_.y() + adjusted_sub_region_dimensions_[2] < 1 ) {
		sub_region_max_[2] = 1;
	} else {
		sub_region_max_[2] = adjusted_sub_area_center_.y() + adjusted_sub_region_dimensions_[2];
	}
	//z
	if ( adjusted_sub_area_center_.z() + adjusted_sub_region_dimensions_[3] > xyz_bound_[3] ) {
		sub_region_max_[3] = xyz_bound_[3];
	} else if ( adjusted_sub_area_center_.z() + adjusted_sub_region_dimensions_[3] < 1 ) {
		sub_region_max_[3] = 1;
	} else {
		sub_region_max_[3] = adjusted_sub_area_center_.z() + adjusted_sub_region_dimensions_[3];
	}

	//min
	//x
	if ( adjusted_sub_area_center_.x() - adjusted_sub_region_dimensions_[1] > xyz_bound_[1] ) {
		//greater than boundary, set to boundary
		sub_region_min_[1] = xyz_bound_[1];
	} else if ( adjusted_sub_area_center_.x() - adjusted_sub_region_dimensions_[1] < 1 ) {
		//smaller than minumum (1), set to 1
		sub_region_min_[1] = 1;
	} else {
		//set normally, as the center - the dimension
		sub_region_min_[1] = adjusted_sub_area_center_.x() - adjusted_sub_region_dimensions_[1];
	}
	//y
	if ( adjusted_sub_area_center_.y() - adjusted_sub_region_dimensions_[2] > xyz_bound_[2] ) {
		sub_region_min_[2] = xyz_bound_[2];
	} else if ( adjusted_sub_area_center_.y() - adjusted_sub_region_dimensions_[2] < 1 ) {
		sub_region_min_[2] = 1;
	} else {
		sub_region_min_[2] = adjusted_sub_area_center_.y() - adjusted_sub_region_dimensions_[2];
	}
	//z
	if ( adjusted_sub_area_center_.z() - adjusted_sub_region_dimensions_[3] > xyz_bound_[3] ) {
		sub_region_min_[3] = xyz_bound_[3];
	} else if ( adjusted_sub_area_center_.z() - adjusted_sub_region_dimensions_[3] < 1 ) {
		sub_region_min_[3] = 1;
	} else {
		sub_region_min_[3] = adjusted_sub_area_center_.z() - adjusted_sub_region_dimensions_[3];
	}



	//with maxima and minima defined, set values within the sub area to 2 instead of 0
	//2 is defined as empty and within the sub area
	for ( core::Size i = sub_region_min_[1]; i <= sub_region_max_[1]; ++i ) {
		for ( core::Size j = sub_region_min_[2]; j <= sub_region_max_[2]; ++j ) {
			for ( core::Size k = sub_region_min_[3]; k <= sub_region_max_[3]; ++k ) {
				protein_matrix_[i][j][k] = 2;
			}
		}
	}


	//set the sub_matrix_volume_
	sub_matrix_volume_ = get_sub_area_grid_volume();
}

// @brief function to elaborate upon the protein_matrix_, and will review the pose and update occupied cells by projecting atom lennard jobes radii and marking cells within the radius as occupied
// if a sub area boundary is defined, will define that area with different values
// note: we are not going to re-size the proteinmatrix of measure the volume that exists outside the ProteinMatrix; this could be updated later, but not doing this simplifies the operation
// a main argument for retaining this as-is is that there is less reason to look at the hwole pose, and usually working with the sub-area only is better, which is less likely to be cut off (only if it is at any edges)
void ProteinGrid::project_lj_radii(){

	//set using lj radii to true
	using_lj_radii_ = true;

	//reset the fullness count
	matrix_fullness_ = 0;
	sub_matrix_fullness_ = 0;

	//iterate over each atom in working_pose
	//translate the atom coordinates to where the center would be in the matrix
	//extract the atom lj radius, and then project it from the center
	for ( core::Size res_num = 1; res_num <= working_pose_->size(); ++res_num ) {
		for ( core::Size atom_num = 1; atom_num <= working_pose_->residue(res_num).natoms(); ++atom_num ) {

			//get the x,y,z data of the atom, rounded to the closest value
			numeric::xyzVector<int> atom_xyz;
			//floor the coordinates down for a constant negative directional shift
			atom_xyz.x() = std::floor(working_pose_->residue(res_num).xyz(atom_num).x());
			atom_xyz.y() = std::floor(working_pose_->residue(res_num).xyz(atom_num).y());
			atom_xyz.z() = std::floor(working_pose_->residue(res_num).xyz(atom_num).z());

			//extract the lj radius as a floating point value, which we will use to project the area about teh atom
			core::Real atom_lj_radius = working_pose_->residue(res_num).atom_type(atom_num).lj_radius();

			//apply the xyz shift and resolution to the atom coordinates
			atom_xyz.x() = std::floor((atom_xyz.x() + xyz_shift_[1]) * resolution_);
			atom_xyz.y() = std::floor((atom_xyz.y() + xyz_shift_[2]) * resolution_);
			atom_xyz.z() = std::floor((atom_xyz.z() + xyz_shift_[3]) * resolution_);

			//apply the resolution to the lj radius and then floor it
			atom_lj_radius *= resolution_;
			atom_lj_radius = std::floor(atom_lj_radius);

			//derive maxima and minima for the radius about the atom
			//if a value would go out of bounds (<1 or > dimension bound), set the the appropriate boundary
			core::Size atom_x_min = 0;
			core::Size atom_x_max = 0;
			core::Size atom_y_min = 0;
			core::Size atom_y_max = 0;
			core::Size atom_z_min = 0;
			core::Size atom_z_max = 0;

			//xmin
			//if below the matrix minimum, set to 1
			if ( atom_xyz.x() - atom_lj_radius < 1 ) {
				atom_x_min = 1;
			} else if ( atom_xyz.x() - atom_lj_radius > xyz_bound_[1] ) {
				//if above the matrix maximum, set to the maximum
				atom_x_min = xyz_bound_[1];
			} else {
				//if within boundaries, keep as is
				atom_x_min = atom_xyz.x() - atom_lj_radius;
			}

			//xmax
			//if below the matrix minimum, set to 1
			if ( atom_xyz.x() + atom_lj_radius < 1 ) {
				atom_x_max = 1;
			} else if ( atom_xyz.x() + atom_lj_radius > xyz_bound_[1] ) {
				//if above the matrix maximum, set to the maximum
				atom_x_max = xyz_bound_[1];
			} else {
				//if within boundaries, keep as is
				atom_x_max = atom_xyz.x() + atom_lj_radius;
			}

			//ymin
			//if below the matrix minimum, set to 1
			if ( atom_xyz.y() - atom_lj_radius < 1 ) {
				atom_y_min = 1;
			} else if ( atom_xyz.y() - atom_lj_radius > xyz_bound_[2] ) {
				//if above the matrix maximum, set to the maximum
				atom_y_min = xyz_bound_[2];
			} else {
				//if within boundaries, keep as is
				atom_y_min = atom_xyz.y() - atom_lj_radius;
			}

			//ymax
			//if below the matrix minimum, set to 1
			if ( atom_xyz.y() + atom_lj_radius < 1 ) {
				atom_y_max = 1;
			} else if ( atom_xyz.y() + atom_lj_radius > xyz_bound_[2] ) {
				//if above the matrix maximum, set to the maximum
				atom_y_max = xyz_bound_[2];
			} else {
				//if within boundaries, keep as is
				atom_y_max = atom_xyz.y() + atom_lj_radius;
			}

			//zmin
			//if below the matrix minimum, set to 1
			if ( atom_xyz.z() - atom_lj_radius < 1 ) {
				atom_z_min = 1;
			} else if ( atom_xyz.z() - atom_lj_radius > xyz_bound_[3] ) {
				//if above the matrix maximum, set to the maximum
				atom_z_min = xyz_bound_[3];
			} else {
				//if within boundaries, keep as is
				atom_z_min = atom_xyz.z() - atom_lj_radius;
			}

			//zmax
			//if below the matrix minimum, set to 1
			if ( atom_xyz.z() + atom_lj_radius < 1 ) {
				atom_z_max = 1;
			} else if ( atom_xyz.z() + atom_lj_radius > xyz_bound_[3] ) {
				//if above the matrix maximum, set to the maximum
				atom_z_max = xyz_bound_[3];
			} else {
				//if within boundaries, keep as is
				atom_z_max = atom_xyz.z() + atom_lj_radius;
			}

			//iterate over a cube around the atom, whose side length is 2x the lj radius (with the atom xyz coordinate in the center)
			//iterate voxel by voxel within the cube to determine which cubes are filled
			//effectively inscribe a sphere within the cube, and the voxels that compose the sphere will be appropriately marked as being occupied
			for ( core::Size i = atom_x_min; i <= atom_x_max; ++i ) {
				for ( core::Size j = atom_y_min; j <= atom_y_max; ++j ) {
					for ( core::Size k = atom_z_min; k <= atom_z_max; ++k ) {
						//check if the coordinate is within the sphere, and if so, provide the proper assignment
						//get the distance of the current point from the atom center, and determine if it is less than the radius
						core::Real atom_cell_distance = sqrt( ((i - atom_xyz.x()) * (i - atom_xyz.x())) + ((j - atom_xyz.y()) * (j - atom_xyz.y())) + ((k - atom_xyz.z()) * (k - atom_xyz.z())) );

						//if the distance from the atom center to the current voxel is within the lj radius, then we have an occupied voxel that needs to be appropriately labeled
						if ( atom_cell_distance <= atom_lj_radius ) {
							//check if the coordinate is within the sub-area, if so, then adjust the value if the point is within the sub area
							if ( using_sub_area_ && is_coordinate_in_sub_area(i,j,k) ) {
								//increment the both fullness if this cell was previously empty (and in the sub area)
								if ( protein_matrix_[i][j][k] == 2 ) {
									++sub_matrix_fullness_;
									++matrix_fullness_;
								}

								protein_matrix_[i][j][k] = 3;

							} else {
								//increment the matrix fullness if this cell was previously empty
								if ( protein_matrix_[i][j][k] == 0 ) {
									++matrix_fullness_;
								}

								protein_matrix_[i][j][k] = 1;
							}


						}
					}
				}
			}

		}
	}




	//safety check to ensure that we actually had atoms in the pose, and that the matrix has a nonzero volume
	if ( matrix_volume_ == 0 ) {
		ms_tr.Warning << "The matrix has no volume, likely because the inputted pose has no atoms!" << std::endl;
	}

	//derive the fullness ratio as the number of occupied cells over the total cell number
	fullness_ratio_ = static_cast<core::Real>(matrix_fullness_) / matrix_volume_;


	//calculate the sub area fullness if using it
	if ( using_sub_area_ ) {
		sub_fullness_ratio_ = static_cast<core::Real>(sub_matrix_fullness_) / sub_matrix_volume_;
	}


}

//@brief function where a residue object (i.e. ligand) outside of a pose can be imposed upon the matrix, and the space filling volume of the system with the ligand can be analyzed
//this forces activation of the space fill data on the class
//the imposed ligand space fill data is retained in object until a function is called that wipes data (i.e. wrap_matrix_around_pose)
//this does allow the placement of multiple ligands if called multiple times, past ligand data will be retained until an operation that called wrap_protein_matrix is performed
void ProteinGrid::placed_ligand_space_fill_analysis(core::conformation::ResidueOP ligresOP)
{
	//begin by forcing the working pose to have the space fill projected about pose atoms, if it is not already
	if ( using_lj_radii_ == false ) {
		using_lj_radii_ = true;
		project_lj_radii();
	}

	//read over each atom in the inputted residue
	//project a sphere around the atom and determine the matrix value to correspond to space occupied by the ligand atoms

	//reminder of what the different values ProteinMatrix values are, since we can use them all here:
	//if the matrix is written to a pdb file, the different states will be translated to different atoms
	//0 = empty and out of sub area, carbon, black
	//1 = pose and out of sub area, fluorine, icy blue
	//2 = empty and in sub area, oxygen, red
	//3 = pose and in sub area, nitrogen, blue
	//4 = do not use, keep even for unoccupied space
	//5 = secondary ligand and out of sub area, sulphur, yellow
	//6 = do not use, keep even for unoccupied space
	//7 = secondary ligand and in sub area, chlorine, green
	//8 = do not use, keep even for unoccupied space
	//9 = secondary ligand and pose out of sub area, phosphorous, orange
	//10 = do not use, keep even for unoccupied space
	//11 = secondary ligand and pose in sub area, iodine, purple

	for ( core::Size residue_atom_iterator = 1; residue_atom_iterator <= ligresOP->natoms(); ++residue_atom_iterator ) {
		//convert coordinates of current atom into format that can be read into space fill matrix matrix
		numeric::xyzVector<int> atom_xyz = ligresOP->xyz(residue_atom_iterator);

		//apply the xyz shift to each value (and do before applying the resolution increase factor)
		atom_xyz.x() = std::floor((atom_xyz.x() + xyz_shift_[1]) * resolution_);
		atom_xyz.y() = std::floor((atom_xyz.y() + xyz_shift_[2]) * resolution_);
		atom_xyz.z() = std::floor((atom_xyz.z() + xyz_shift_[3]) * resolution_);

		//derive the atom lj radius and apply the resolution to it, then floor
		core::Real atom_lj_radius = ligresOP->atom_type(residue_atom_iterator).lj_radius();
		atom_lj_radius *= resolution_;
		atom_lj_radius = std::floor(atom_lj_radius);

		//derive maxima and minima for the radius about the atom
		//if a value would go out of bounds (<1 or > dimension bound), set the the appropriate boundary
		core::Size atom_x_min = 0;
		core::Size atom_x_max = 0;
		core::Size atom_y_min = 0;
		core::Size atom_y_max = 0;
		core::Size atom_z_min = 0;
		core::Size atom_z_max = 0;

		//xmin
		//if below the matrix minimum, set to 1
		if ( atom_xyz.x() - atom_lj_radius < 1 ) {
			atom_x_min = 1;
		} else if ( atom_xyz.x() - atom_lj_radius > xyz_bound_[1] ) {
			//if above the matrix maximum, set to the maximum
			atom_x_min = xyz_bound_[1];
		} else {
			//if within boundaries, keep as is
			atom_x_min = atom_xyz.x() - atom_lj_radius;
		}

		//xmax
		//if below the matrix minimum, set to 1
		if ( atom_xyz.x() + atom_lj_radius < 1 ) {
			atom_x_max = 1;
		} else if ( atom_xyz.x() + atom_lj_radius > xyz_bound_[1] ) {
			//if above the matrix maximum, set to the maximum
			atom_x_max = xyz_bound_[1];
		} else {
			//if within boundaries, keep as is
			atom_x_max = atom_xyz.x() + atom_lj_radius;
		}

		//ymin
		//if below the matrix minimum, set to 1
		if ( atom_xyz.y() - atom_lj_radius < 1 ) {
			atom_y_min = 1;
		} else if ( atom_xyz.y() - atom_lj_radius > xyz_bound_[2] ) {
			//if above the matrix maximum, set to the maximum
			atom_y_min = xyz_bound_[2];
		} else {
			//if within boundaries, keep as is
			atom_y_min = atom_xyz.y() - atom_lj_radius;
		}

		//ymax
		//if below the matrix minimum, set to 1
		if ( atom_xyz.y() + atom_lj_radius < 1 ) {
			atom_y_max = 1;
		} else if ( atom_xyz.y() + atom_lj_radius > xyz_bound_[2] ) {
			//if above the matrix maximum, set to the maximum
			atom_y_max = xyz_bound_[2];
		} else {
			//if within boundaries, keep as is
			atom_y_max = atom_xyz.y() + atom_lj_radius;
		}

		//zmin
		//if below the matrix minimum, set to 1
		if ( atom_xyz.z() - atom_lj_radius < 1 ) {
			atom_z_min = 1;
		} else if ( atom_xyz.z() - atom_lj_radius > xyz_bound_[3] ) {
			//if above the matrix maximum, set to the maximum
			atom_z_min = xyz_bound_[3];
		} else {
			//if within boundaries, keep as is
			atom_z_min = atom_xyz.z() - atom_lj_radius;
		}

		//zmax
		//if below the matrix minimum, set to 1
		if ( atom_xyz.z() + atom_lj_radius < 1 ) {
			atom_z_max = 1;
		} else if ( atom_xyz.z() + atom_lj_radius > xyz_bound_[3] ) {
			//if above the matrix maximum, set to the maximum
			atom_z_max = xyz_bound_[3];
		} else {
			//if within boundaries, keep as is
			atom_z_max = atom_xyz.z() + atom_lj_radius;
		}

		//iterate over a cube around the atom, whose side length is 2x the lj radius (with the atom xyz coordinate in the center)
		//iterate voxel by voxel within the cube to determine which cubes are filled
		//effectively inscribe a sphere within the cube, and the voxels that compose the sphere will be appropriately marked as being occupied
		for ( core::Size i = atom_x_min; i <= atom_x_max; ++i ) {
			for ( core::Size j = atom_y_min; j <= atom_y_max; ++j ) {
				for ( core::Size k = atom_z_min; k <= atom_z_max; ++k ) {
					//check if the coordinate is within the sphere, and if so, provide the proper assignment
					//get the distance of the current point from the atom center, and determine if it is less than the radius
					core::Real atom_cell_distance = sqrt( ((i - atom_xyz.x()) * (i - atom_xyz.x())) + ((j - atom_xyz.y()) * (j - atom_xyz.y())) + ((k - atom_xyz.z()) * (k - atom_xyz.z())) );

					//if the distance from the atom center to the current voxel is within the lj radius, then we have an occupied voxel that needs to be appropriately labeled
					if ( atom_cell_distance <= atom_lj_radius ) {
						//cases to test for appropriately filling the cell

						//outside sub area and not occupied by pose, should be initially 0
						//set to 5
						if ( is_coordinate_in_sub_area(i,j,k) == false && protein_matrix_[i][j][k] == 0 ) {
							protein_matrix_[i][j][k] = 5;
							//increment occupied cell count by 1
							++matrix_fullness_;
						}

						//inside sub area and not occupied by pose, should be initially 2
						//set to 7
						if ( using_sub_area_ && is_coordinate_in_sub_area(i,j,k) && protein_matrix_[i][j][k] == 2 ) {
							protein_matrix_[i][j][k] = 7;
							//increment occupied cell count by 1
							++matrix_fullness_;
							++sub_matrix_fullness_;
						}

						//outside sub area and occupied by pose, should be initially 1
						//set to 9
						//do not increment fullness because it is already occupied by the pose
						if ( is_coordinate_in_sub_area(i,j,k) == false && protein_matrix_[i][j][k] == 1 ) {
							protein_matrix_[i][j][k] = 9;
						}

						//inside sub area and occupied by pose, should be initially 3
						//set to 11
						//do not increment fullness because it is already occupied by the pose
						if ( using_sub_area_ && is_coordinate_in_sub_area(i,j,k) && protein_matrix_[i][j][k] == 3 ) {
							protein_matrix_[i][j][k] = 11;

						}
					}
				}
			}
		}
	}


	//derive the fullness ratio as the number of occupied cells over the total cell number
	fullness_ratio_ = static_cast<core::Real>(matrix_fullness_) / matrix_volume_;

	//calculate the sub area fullness if using it
	if ( using_sub_area_ ) {
		sub_fullness_ratio_ = static_cast<core::Real>(sub_matrix_fullness_) / sub_matrix_volume_;
	}


}

//@brief function that takes in a residue object (i.e. ligand) and determines if the ligand clashes with the pose in this class
//returns true if there is a clash, and false if there is no clash
//this can work with with and without space fill, and is unaffected by a sub area
//this function does not modify the ProteinMatrix, and instead either iterates over all atoms in the ligresOP and returns true if there is no clash, or returns false at the first clashing hit
//to be fast, this function does not utilize lj radii on the ligand atoms
bool ProteinGrid::placed_ligand_clash_analysis(core::conformation::ResidueOP ligresOP)
{
	//iterate over all atoms in the ligresOP to see if any clash with the pose
	for ( core::Size residue_atom_iterator = 1; residue_atom_iterator <= ligresOP->natoms(); ++residue_atom_iterator ) {
		//convert coordinates of current atom into format that can be read into space fill matrix matrix
		numeric::xyzVector<int> atom_xyz = ligresOP->xyz(residue_atom_iterator);

		//apply the xyz shift to each value (and do before applying the resolution increase factor)
		atom_xyz.x() = std::floor((atom_xyz.x() + xyz_shift_[1]) * resolution_);
		atom_xyz.y() = std::floor((atom_xyz.y() + xyz_shift_[2]) * resolution_);
		atom_xyz.z() = std::floor((atom_xyz.z() + xyz_shift_[3]) * resolution_);

		//plug the atom xyz coordinates into the proteinmatrix and see what the value is
		//if the value is odd, then it is occupied by the pose (or technically could be a placed ligand, but that would indicate a clash as well)
		if ( protein_matrix_[atom_xyz.x()][atom_xyz.y()][atom_xyz.z()] % 2 == 1 ) {
			return true;
		}
	}

	//if we iterated over all atoms and made it to this point, return true
	return false;
}

// @ brief function that prints out the current state of the ProteinMatrix as a Pose, so the user can do things like write the pose to a pdb
core::pose::Pose ProteinGrid::export_protein_matrix_to_pose()
{
	//create a dummy mutable residue type that we can use to add atoms to
	//it's a huge pain to build a new residue from scratch, so we are going to assume that working_pose_ actually has contents
	//we'll pluck the first residue from the pose, and use it to make our dummy mrt

	//warning and handling in case working_pose_ is actually empty, just return the working pose
	if ( working_pose_->size() == 0 ) {
		ms_tr.Warning << "There is nothing in the working_pose_! Load something in so we can export the protein matrix. Returning the empty working_pose_." << std::endl;
		return *working_pose_;
	}

	//pluck the first residue in working_pose_
	core::conformation::Residue dummyligres = working_pose_->residue(1);

	//derive a mutable residue type from the residue
	core::chemical::MutableResidueTypeOP dummylig_mrt( new core::chemical::MutableResidueType( dummyligres.type() ) );

	//run through each atom in dummylig_mrt and delete it so it is effectively clear
	//this definitely goes against what the data structure is designed for (considering that the default constructor is blocked from use), but hopefully this just does what I need it to do
	//see if this works: keep deleting the atom at index 1 until the number of atoms in the mrt is 0
	while ( dummylig_mrt->natoms() > 0 )
			{
		//debug print of the first atom in the list of verticies and number of atoms

		ms_tr.Trace << "Number of atoms in dummy mrt is now: " <<  dummylig_mrt->natoms() << std::endl;

		//pop the first atom in the list of vertices
		dummylig_mrt->delete_atom(dummylig_mrt->all_atoms()[1]);
	}


	ms_tr.Trace << "Number of atoms in dummy mrt is now: " <<  dummylig_mrt->natoms() << std::endl;

	//delete all chis in the dummy ligand
	//probably easiest to use delete_terminal_chi
	while ( dummylig_mrt->nchi() > 0 )
	{
		ms_tr.Trace << "Number of chi in residue is currently: " << dummylig_mrt->nchi() << std::endl;

		//delete the terminal chi
		dummylig_mrt->delete_terminal_chi();
	}

	ms_tr.Trace << "Number of chi in residue is currently: " << dummylig_mrt->nchi() << std::endl;

	//we now have a blank dummy residue type that we can rebuild as a way to represent the ProteinMatrix, and could spit out to a pdb
	//declare an additional mrt object that we will fill with atoms. every 100 atoms, we'll add the build mrt to the pose, and then wipe it using the dummy so we can build again
	core::chemical::MutableResidueType my_mrt( *dummylig_mrt );

	//declare strings to use for the atom name, type, and mm_type
	std::string atom_name = "";
	std::string atom_type_name = "";
	std::string mm_atom_type_name = "";

	//declare vector to hold atom coordinates
	utility::vector1<core::Real> atom_xyz(3,0);

	//declare a vector to hold the names of first 3 atoms so we can apply icoor data (especialy for the first 3 atoms)
	utility::vector1<std::string> first_three_atoms;

	//counter to determine the first 3 atoms so we can assign icoor_data to the first 3 atoms after we reach them (can't do beforehand, as we do not know atoms before we reach them)
	core::Size atom_counter = 0;

	//counter to count the total atom number for unique name tracking
	core::Size total_atom_counter = 0;

	//create pose for matrix
	core::pose::Pose matrix_pose;

	//reminder of what the different values ProteinMatrix values are, since we can use them all here:
	//if the matrix is written to a pdb file, the different states will be translated to different atoms
	//0 = empty and out of sub area, carbon, black
	//1 = pose and out of sub area, fluorine, icy blue
	//2 = empty and in sub area, oxygen, red
	//3 = pose and in sub area, nitrogen, blue
	//4 = do not use, keep even for unoccupied space
	//5 = secondary ligand and out of sub area, sulphur, yellow
	//6 = do not use, keep even for unoccupied space
	//7 = secondary ligand and in sub area, chlorine, green
	//8 = do not use, keep even for unoccupied space
	//9 = secondary ligand and pose out of sub area, phosphorous, orange
	//10 = do not use, keep even for unoccupied space
	//11 = secondary ligand and pose in sub area, iodine, purple

	//iterate over each cell in the matrix to create atoms to add to the mrt
	for ( core::Size x = 1; x <= xyz_bound_[1]; ++x ) {
		for ( core::Size y = 1; y <= xyz_bound_[2]; ++y ) {
			for ( core::Size z = 1; z <= xyz_bound_[3]; ++z ) {
				//check whether to print the atom if it is outside the sub area and if we want to print outside the sub area
				//if not using a sub area, we will assume that we want to print the whole pose
				if ( using_sub_area_ && print_whole_matrix_ == false && (protein_matrix_[x][y][z] == 0 || protein_matrix_[x][y][z] == 1 || protein_matrix_[x][y][z] == 5 || protein_matrix_[x][y][z] == 9) ) {
					//move to next voxel
					continue;
				}

				//check whether to print the atom if it is empty space (printing with empty space is slower and makes a more memory intensive file)
				//values that correspond to emptiness are even (0 and 2)
				if ( print_empty_space_ == false && protein_matrix_[x][y][z] % 2 == 0 ) {
					//move to next voxel
					continue;
				}

				//derive unique atom name based on the atom count
				atom_name = base_10_to_base_62(total_atom_counter);

				//carbon
				if ( protein_matrix_[x][y][z] == 0 ) {
					atom_type_name = "aroC";
					mm_atom_type_name = "C";
				}
				//hydrogen
				//hydrogen may be a problem, try fluorine instead (icy blue)
				if ( protein_matrix_[x][y][z] == 1 ) {
					//atom_type_name = "Haro";
					//mm_atom_type_name = "H";
					atom_type_name = "F";
					mm_atom_type_name = "F1";
				}
				//oxygen
				if ( protein_matrix_[x][y][z] == 2 ) {
					atom_type_name = "OH";
					mm_atom_type_name = "O";
				}
				//nitrogen
				if ( protein_matrix_[x][y][z] == 3 ) {
					atom_type_name = "NH2O";
					mm_atom_type_name = "N";
				}
				//sulphur
				if ( protein_matrix_[x][y][z] == 5 ) {
					atom_type_name = "S";
					mm_atom_type_name = "S";
				}
				//chlorine
				if ( protein_matrix_[x][y][z] == 7 ) {
					atom_type_name = "Cl";
					mm_atom_type_name = "CL";
				}
				//phosphorous
				if ( protein_matrix_[x][y][z] == 9 ) {
					atom_type_name = "Pha";
					mm_atom_type_name = "P";
				}
				//iodine
				if ( protein_matrix_[x][y][z] == 11 ) {
					atom_type_name = "I";
					mm_atom_type_name = "I";
				}

				//add the atom to the mrt
				//default charge to 0
				my_mrt.add_atom(atom_name, atom_type_name, mm_atom_type_name, 0);

				//derive the coordinates of the atom
				//divide matrix coordinates by resolution scalar and then subtract corresponding shift (un-apply the resolution and shift so that this atom aligns with a position in the pose)
				atom_xyz[1] = (static_cast<core::Real>(x)/resolution_) - xyz_shift_[1];
				atom_xyz[2] = (static_cast<core::Real>(y)/resolution_) - xyz_shift_[2];
				atom_xyz[3] = (static_cast<core::Real>(z)/resolution_) - xyz_shift_[3];

				//set the coordinates of the atom in the mrt
				my_mrt.set_ideal_xyz(atom_name, numeric::xyzVector<core::Real>(atom_xyz[1],atom_xyz[2],atom_xyz[3]) );

				//handle adding the icoor data for the atom
				//skip if the atom name is "" (which will happen if you only look at the sub-matrix)
				if ( atom_name != "" ) {
					++atom_counter;
					++total_atom_counter;

					//if atom counter is <= 3, add the atom name to the first_three_atoms vector
					//atom count must be >3 so we know the first 3 atoms and can retroactively add icoor for those first 3 atoms
					if ( atom_counter <= 3 ) {
						first_three_atoms.push_back(atom_name);

						//if the count is 3, we can retroactively add icoor data for the first 3 atoms now
						if ( atom_counter == 3 ) {

							my_mrt.set_icoor(first_three_atoms[1],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);

							my_mrt.set_icoor(first_three_atoms[2],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
							//flip order of stub2 and stub3 for 3rd atom (looking at some params files, it seems like they do this)
							my_mrt.set_icoor(first_three_atoms[3],0,0,0,first_three_atoms[1],first_three_atoms[3],first_three_atoms[2],false);

							//add bonds to connect these 3
							my_mrt.add_bond(first_three_atoms[1],first_three_atoms[2]);
							my_mrt.add_bond(first_three_atoms[2],first_three_atoms[3]);
						}
					} else {
						//add icoor data for later atoms, using the first 3 atoms as the same stub atoms each time
						my_mrt.set_icoor(first_three_atoms[1],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);

						//bond everything else to the 3rd atom and hope things dont get wonky
						//if they do, I can try another means to attach atoms to each other looking at previous atoms
						my_mrt.add_bond(first_three_atoms[3],atom_name);

					}

					//cut off when my_mrt becomes 100 atoms and start a new residue so that the conversion for an rt doesn't take too long
					//there might need to be some handling at the end if my_mrt only has 1-2 atoms
					//handling would probably just be to add fake atoms in the place of the first atom to fill out to 3 for the icoor data
					if ( atom_counter == 100 ) {
						//convert my_mrt to a data type that can be added to a pose

						//assign internal coordinates
						my_mrt.assign_internal_coordinates();

						//convert to residuetype
						core::chemical::ResidueTypeCOP my_rt(core::chemical::ResidueType::make(my_mrt));

						//convert to residue
						core::conformation::ResidueOP my_res( core::conformation::ResidueFactory::create_residue(*my_rt));

						//append residue to pose
						matrix_pose.append_residue_by_jump( *my_res, 1 );

						//reset atom counter
						atom_counter = 0;

						//pop back the entries to first_three_atoms so we can start again
						first_three_atoms.pop_back();
						first_three_atoms.pop_back();
						first_three_atoms.pop_back();

						//reset my_mrt
						my_mrt = *dummylig_mrt;
					}

					//set name back to "" so that we don't keep going after we finish the sub area
					atom_name = "";

				}


			}
		}
	}

	ms_tr.Trace << "Made all atoms for this matrix" << std::endl;

	//handle rare case where the my_mrt that reaches this point only has 0-2 atoms and would have no icoor data
	//if no atoms, just return the pose
	if ( my_mrt.natoms() == 0 ) {
		ms_tr.Trace << "Returning system pose." << std::endl;
		return matrix_pose;
	} else if ( my_mrt.natoms() == 1 ) {
		//if only 1 atom, make 2 "copies" (atom in the same position as original) of the atom and assign icoor data
		//make name of atom
		atom_name = base_10_to_base_62(total_atom_counter + 1);
		//add atom to mrt
		//use same types as what was used last
		my_mrt.add_atom(atom_name, atom_type_name, mm_atom_type_name, 0);
		//add atom to the first three atoms list for tracking
		first_three_atoms.push_back(atom_name);

		//make name of atom
		atom_name = base_10_to_base_62(total_atom_counter + 2);
		//add atom to mrt
		//use same types as what was used last
		my_mrt.add_atom(atom_name, atom_type_name, mm_atom_type_name, 0);
		//add atom to the first three atoms list for tracking
		first_three_atoms.push_back(atom_name);

		//set icoor data for these atoms
		my_mrt.set_icoor(first_three_atoms[1],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
		my_mrt.set_icoor(first_three_atoms[2],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
		my_mrt.set_icoor(first_three_atoms[3],0,0,0,first_three_atoms[1],first_three_atoms[3],first_three_atoms[2],false);

		//add bonds to connect these 3
		my_mrt.add_bond(first_three_atoms[1],first_three_atoms[2]);
		my_mrt.add_bond(first_three_atoms[2],first_three_atoms[3]);

	} else if ( my_mrt.natoms() == 2 ) {
		//if 2 atoms, make a 3rd atom that is a "copy" of the first atom
		//make name of atom
		atom_name = base_10_to_base_62(total_atom_counter + 1);
		//add atom to mrt
		//use same types as what was used last
		my_mrt.add_atom(atom_name, atom_type_name, mm_atom_type_name, 0);
		//add atom to the first three atoms list for tracking
		first_three_atoms.push_back(atom_name);

		//set icoor data for these atoms
		my_mrt.set_icoor(first_three_atoms[1],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
		my_mrt.set_icoor(first_three_atoms[2],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
		my_mrt.set_icoor(first_three_atoms[3],0,0,0,first_three_atoms[1],first_three_atoms[3],first_three_atoms[2],false);

		//add bonds to connect the 3 atoms
		my_mrt.add_bond(first_three_atoms[1],first_three_atoms[2]);
		my_mrt.add_bond(first_three_atoms[2],first_three_atoms[3]);
	}

	//now need to finalize the final mrt and add it to the pose

	//assign internal coordinates
	my_mrt.assign_internal_coordinates();

	//convert my_mrt to a data type that can be added to a pose
	//convert to residuetype
	core::chemical::ResidueTypeCOP my_rt(core::chemical::ResidueType::make(my_mrt));

	//convert to residue
	core::conformation::ResidueOP my_res( core::conformation::ResidueFactory::create_residue(*my_rt));

	//append residue to pose
	matrix_pose.append_residue_by_jump( *my_res, 1 );

	ms_tr.Trace << "Returning system pose." << std::endl;

	//return the matrix_pose that is full of atoms in a grid that represents teh ProteinMatrix
	return matrix_pose;
}

// @ brief function that prints out the current state of the ProteinMatrix as a pdb; calls export_protein_matrix_to_pose and goes the extra step to print out the pose to a pdb without the user having to do more
//takes in a string to use to help assign a name to the created pose and pdb
void ProteinGrid::export_protein_matrix_to_pdb(std::string pdb_name_prefix)
{
	//create a file name to output the pose to
	std::string matrix_pdb_name = pdb_name_prefix + "_WholeRatio_" + std::to_string(fullness_ratio_);



	//if using a sub area, tack that onto the name
	if ( using_sub_area_ ) {
		matrix_pdb_name = matrix_pdb_name + "_SubRatio_" + std::to_string(sub_fullness_ratio_);


	}

	//tack on end
	matrix_pdb_name = matrix_pdb_name + ".pdb";

	ms_tr.Trace << "Preparing to make visualization pose for " << matrix_pdb_name << std::endl;

	//create a pose to print to a pdb
	core::pose::Pose matrix_pose = export_protein_matrix_to_pose();

	//dump the pose to a pdb
	matrix_pose.dump_pdb(matrix_pdb_name);

	//add a comment to the working_pose_ that notes the matrix file that was made
	core::pose::add_comment(*working_pose_, "Corresponding space fill matrix file:", matrix_pdb_name);
}


// @brief function to be used to convert a base 10 number to base 62 (as a string with characters represented by upper+lower case letters and digits)
//used in export_protein_matrix_to_pose to assign a unique name to an atom
//Due to limitations in atom icoor data, an atom name can be no longer than 4 characters, so this provides 62^4 (~14.7M) unique atom names
std::string ProteinGrid::base_10_to_base_62(core::Size starting_num)
{
	//bool to indicate whether to keep processing the number
	//stop processing once we fully build the string, which occurs when dividing the current number by 62 is <1 (integer would be 0)
	bool keep_processing = true;

	//string to hold the base 62 representation of the number
	//characters represented by base_62_cipher_ vector in the .hh file
	std::string base_62_number = "";

	core::Size curr_num = starting_num;

	while ( keep_processing )
			{
		//take the quotient of the current number by 62
		core::Size quotient = curr_num / 62;

		//derive the modulus of the current number by 62
		core::Size mod = curr_num % 62;

		//use the mod value to get the corresponding character from the cipher and append to the string
		//character appends to the front of the string
		base_62_number.insert(base_62_number.begin(),utility::code_to_6bit(mod + 1));

		//if the quotient is under 62, get the number from the cipher, otherwise we have to repeat the processing operation
		if ( quotient < 62 ) {
			// use the quotient to add to the number, unless the quotient is 0 (no need to add a placeholder 0)
			if ( quotient != 0 ) {
				base_62_number.insert(base_62_number.begin(),utility::code_to_6bit(quotient + 1));
			}

			//we have now fully derived the base 62 number and can stop processing
			keep_processing = false;
		} else {
			//set the quotient to the current number and we continue the operation off of it
			curr_num = quotient;
		}
	}

	return base_62_number;
}

// @brief critical function that builds/rebuilds (overwrites existing) the proteingrid protein_matrix_ around the pose in working_pose_
//called on creation of the object, and also called in cases like changing the pose or resolution
//this function will simply fill out cells with either a 0 if there is no atom present in the corresponding cell, or a 1 if there is
//for a more in-depth space fill (project LJ radii), the space fill function will need to be called
void ProteinGrid::wrap_matrix_around_pose()
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
			if ( res_num == 1 && atom_num == 1 ) {
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
	xyz_bound_[1] = std::floor((xyz_shift_[1] + largest_x) * resolution_);
	xyz_bound_[2] = std::floor((xyz_shift_[2] + largest_y) * resolution_);
	xyz_bound_[3] = std::floor((xyz_shift_[3] + largest_z) * resolution_);

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
	sub_matrix_fullness_ = 0;
	sub_fullness_ratio_ = 0;

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

	//this needs to be defined after filling out the protein_matrix_
	//if the true_sub_area_dimension_ values are not 0, then we have a sub area that we want to work with and consider for investigation
	if ( using_sub_area_ ) {
		//apply the shift and resolution to the true center and dimension, and derive the sub area max and min
		define_sub_regions();
	}

	//seed the matrix with approximate coordinates of each atom
	//apply the shift to the coordinates
	//approximated by flooring coordinates down
	for ( core::Size xyzVec = 1; xyzVec <= atom_coordinates.size(); ++xyzVec ) {

		//declaring variables to more easily access coordinates in the protein_matrix_
		core::Size x_idx = floor((atom_coordinates[xyzVec].x() + xyz_shift_[1]) * resolution_);
		core::Size y_idx = floor((atom_coordinates[xyzVec].y() + xyz_shift_[2]) * resolution_);
		core::Size z_idx = floor((atom_coordinates[xyzVec].z() + xyz_shift_[3]) * resolution_);

		//throw in a warning/check in case any of these values floor to 0
		if ( x_idx == 0 ) {
			ms_tr.Warning << "Warning, we floored an atom coordinate to 0. Going to automatically nudge it to 1." << std::endl;
			++x_idx;
		}
		if ( y_idx == 0 ) {
			ms_tr.Warning << "Warning, we floored an atom coordinate to 0. Going to automatically nudge it to 1." << std::endl;
			++y_idx;
		}
		if ( z_idx == 0 ) {
			ms_tr.Warning << "Warning, we floored an atom coordinate to 0. Going to automatically nudge it to 1." << std::endl;
			++z_idx;
		}

		//check if the coordinate is within the sub-area, if so, then adjust the value if the point is within the sub area
		if ( using_sub_area_ && is_coordinate_in_sub_area(x_idx,y_idx,z_idx) ) {
			//increment the both fullness if this cell was previously empty (and in the sub area)
			if ( protein_matrix_[x_idx][y_idx][z_idx] == 2 ) {
				++sub_matrix_fullness_;
				++matrix_fullness_;
			}

			protein_matrix_[x_idx][y_idx][z_idx] = 3;
		} else {
			//increment the main fullness if this cell was previously empty
			if ( protein_matrix_[x_idx][y_idx][z_idx] == 2 ) {
				++matrix_fullness_;
			}

			//set to 1
			protein_matrix_[x_idx][y_idx][z_idx] = 1;
		}
	}

	//safety check to ensure that we actually had atoms in the pose, and that the matrix has a nonzero volume
	if ( matrix_volume_ == 0 ) {
		ms_tr.Warning << "The matrix has no volume, likely because the inputted pose has no atoms!" << std::endl;
	}

	//derive the fullness ratio as the number of occupied cells over the total cell number
	fullness_ratio_ = static_cast<core::Real>(matrix_fullness_) / matrix_volume_;

	//calculate the sub area fullness if using it
	if ( using_sub_area_ ) {
		sub_fullness_ratio_ = static_cast<core::Real>(sub_matrix_fullness_) / sub_matrix_volume_;
	}
}

// @brief function to check and ensure that the resolution value is valid (> 0)
// if invalid, will throw a warning and set resolution to 1
void ProteinGrid::validate_resolution()
{
	if ( resolution_ <= 0 ) {
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
	true_sub_area_center_.x() = 0;
	true_sub_area_center_.y() = 0;
	true_sub_area_center_.z() = 0;

	adjusted_sub_area_center_.x() = 0;
	adjusted_sub_area_center_.y() = 0;
	adjusted_sub_area_center_.z() = 0;

	true_sub_region_dimensions_ = utility::vector1<int>(3, 0);
	adjusted_sub_region_dimensions_ = utility::vector1<int>(3, 0);

	sub_region_max_ = utility::vector1<int>(3, 0);
	sub_region_min_ = utility::vector1<int>(3, 0);
}

// @brief run a check to see if a coordinate within the protein_matrix is within the sub area
//returns true if the coordinate is, false otherwise
//I don't think we want this to be a public function, since this uses coordinates relative to the protein matrix, and not pose coordinates (which are shifted and potentially stretched)
bool ProteinGrid::is_coordinate_in_sub_area(core::Size x, core::Size y, core::Size z)
{
	//greedily return false if we are not even using a sub area, since we are not even using a sub area
	if ( using_sub_area_ == false ) {
		return false;
	}

	if ( x >= sub_region_min_[1] && x <= sub_region_max_[1] && y >= sub_region_min_[2] && y <= sub_region_max_[2] && z >= sub_region_min_[3] && z <= sub_region_max_[3] ) {
		return true;
	} else {
		return false;
	}
}

// @brief function that turns off using a sub area (until a function is called that turns it back on, like passing in new sub area dimensions)
// calls a rewrap on the pose that will now ignore the sub area
void ProteinGrid::ignore_sub_area()
{
	//throw warning and return if using_sub_area_ was already false; return
	if ( using_sub_area_ == false ) {
		ms_tr.Warning << "Warning! We were already not using a sub area." << std::endl;
		return;
	}

	using_sub_area_ = false;
	wrap_matrix_around_pose();

	//if using the lj radius to define volume, get the lj data again, because it is wiped by effectively reseting with the wrap_matrix function
	if ( using_lj_radii_ ) {
		project_lj_radii();
	}
}

// @ brief function that turns off lj radii by recalling the wrap_matrix function
void ProteinGrid::ignore_lj_radii()
{
	//throw warning and return if using_lj_radii_ was already false; return
	if ( using_lj_radii_ == false ) {
		ms_tr.Warning << "Warning! We were already not lj radii." << std::endl;
		return;
	}

	//set to false and reset wrap around matrix
	using_lj_radii_ = false;
	wrap_matrix_around_pose();
}

// @brief function to set the value of print_whole_matrix_
void ProteinGrid::set_print_whole_matrix(bool setter)
{
	print_whole_matrix_ = setter;
}

// @brief function to set the value of print_empty_space_
void ProteinGrid::set_print_empty_space(bool setter)
{
	print_empty_space_ = setter;
}
}
}
