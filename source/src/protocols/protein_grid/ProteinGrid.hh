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

#ifndef INCLUDED_protocols_protein_grid_ProteinGrid_hh
#define INCLUDED_protocols_protein_grid_ProteinGrid_hh

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
#include <core/pose/Pose.fwd.hh>

#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

#include <algorithm>
#include <sstream>

namespace protocols {
namespace protein_grid {

class ProteinGrid
{
public:
	/// @brief condensing the data type name for the main 3D vector (or matrix) that represents the protein
	// currently laid out in that the matrix internal data can only be of Size. If flexibility is needed, this could be expanded upon later. The main reason for not doing this now is that operations that use the matrix data rely on the data being a discrete set of positive integers
	typedef utility::vector1<utility::vector1<utility::vector1<core::Size>>> ProteinMatrix;

	//parameterized constructors

	//simple constructor that only takes in pose
	ProteinGrid(core::pose::PoseOP in_pose);

	//constructor that takes in pose and new resolution value
	ProteinGrid(core::pose::PoseOP in_pose, core::Real resolution);

	//constructor that takes in pose and sub_region_min/max vectors
	ProteinGrid(core::pose::PoseOP in_pose, utility::vector1<core::Size> sub_region_max, utility::vector1<core::Size> sub_region_min);

	//constructor that takes in pose, resolution values, and subregion vectors
	ProteinGrid(core::pose::PoseOP in_pose, core::Real resolution, utility::vector1<core::Size> sub_region_max, utility::vector1<core::Size> sub_region_min);

	//copy constructor
	ProteinGrid( ProteinGrid const & other );

	//@brief = operator overload
	ProteinGrid & operator=( ProteinGrid const & other );

	// destructor
	~ProteinGrid();

	// @brief function to clone the current ProteinGrid into inputted ProteinGrid
	void clone(ProteinGrid & copy) const ;

	// @brief simple function to derive the volume of the matrix
	core::Size get_grid_volume();

	// @brief function to set the xyz coordinates of the sub_region_min
	// reminder, these values should directly relate to coordinates in the pose, and not be directly aimed at the matrix indices
	void set_sub_region_min( utility::vector1<core::Size> region_in );

	// @brief function to set the xyz coordinates of the sub_region_max
	// reminder, these values should directly relate to coordinates in the pose, and not be directly aimed at the matrix indices
	void set_sub_region_max( utility::vector1<core::Size> region_in );

	// @brief function to set the xyz coordinates of sub_region_max and sun_region_min
	// reminder, these values should directly relate to coordinates in the pose, and not be directly aimed at the matrix indices
	void set_sub_regions( utility::vector1<core::Size> region_max, utility::vector1<core::Size> region_min );	

	// @brief function to elaborate upon the protein_matrix_, and will review the pose and update occupied cells by projecting atom lennard jobes radii and marking cells within the radius as occupied
	// if a sub area boundary is defined, will define that area with different values
	void project_lj_radii();

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
	void wrap_matrix_around_pose();

	// @brief function to check and ensure that the resolution value is valid (> 0)
	// if invalid, will throw a warning and set resolution to 1
	void validate_resolution();

	// @brief reset (or set) the xyz shift and bound matrices to be 3 values of zeroes
	void reset_xyz_vectors();

	// @brief reset (or set) the sub-region matrices to be 3 values of zeroes
	void reset_sub_region_vectors();	

	// @brief 3D matrix of voxelized representation of atoms in pose; individual indices contain data values that correspond to whether the voxel is occupied by pose atoms
	// the coordinates of pose atoms are used to correspond to voxels in this matrix. Atom coordinates are all shifted so that the minimum coordinate value of all pose atoms are shifted to 1
	ProteinMatrix protein_matrix_;

	// @brief the corresponding pose that the protein_matrix_ is wrapped around
	core::pose::PoseOP working_pose_;

	// @brief vector to hold the values that atom coordinates shift by in the x,y,z directions
	utility::vector1<int> xyz_shift_;

	// @brief vector to hold the xyz boundaries of the matrix (not the pose coordinates, but the boundaries that are shifted and potentially scaled within the matrix)
	utility::vector1<core::Size> xyz_bound_;

	// @brief value to hold the number of cells/voxels within the matrix, derived by the product of the xyz_bound_ dimensions
	core::Size matrix_volume_ = 0;

	// @brief floating value that can be used to scale the resolution of the protein grid, if desired
	// default resolution of voxels are at 1 cubic angstrom (generally recommended)
	//values <1 decrease the resolution of the matrix, increasing the likelihood of multiple atoms appearing in the same voxel
	//values >1 increase the resolution, which can help better define the sphere shape of lj radii from atoms; this comes at a substantial time/memory code
	//resolution value must be positive
	core::Real resolution_ = 1;

	// @brief vectors of maximum and minimum xyz values for a sub-area to potentially investigate in methods like space filling
	// the region is a rectangular prism that is defined by coordinates outlined in sub_refoun_min_ and sub_region_max_
	//seed values to be zero to indicate that we should not be using the subregion unless these values get set to be >=1
	//these are meant to be coordinats in the pose (rounded to an int), not cells in the ProteinMatrix
	//the values will be scaled and shifted when used
	utility::vector1<core::Size> sub_region_max_;
	utility::vector1<core::Size> sub_region_min_;

	// @brief values to track how full the matrix is with atoms, and a ratio to track the fullness value compared to the size of the matrix
	core::Size matrix_fullness_ = 0;
	core::Real fullness_ratio_ = 0;
};

}
}


#endif 
