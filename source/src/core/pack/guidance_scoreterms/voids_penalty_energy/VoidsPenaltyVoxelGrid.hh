// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyVoxelGrid.hh
/// @brief A 3D boolean array used for identifying core voxels in the VoidsPenaltyEnergy.
/// @details The order of operations for using this are:
/// 1. Create an instance.
/// 2. Set voxel size and padding (set_voxel_size_and_padding()), plus cone parameters (set_cone_parameters()).
/// 3. Initialize from a pose (set_up_voxel_grid_and_compute_burial()).
/// 4. Optionally, prune voxels for positions that are not packable by calling prune_voxels_for_fixed_residues().
/// 5. Pass in a residue set and compute volumes of buried rotamers (compute_volumes_of_buried_rotamers()).
/// 6. The total_buried_volume() and reachable_buried_volume() functions have been computed and can be called at this point.
/// 7. It should be safe to call step 4 again without calling reset() or repeating step 3.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_voids_penalty_energy_VoidsPenaltyVoxelGrid_hh
#define INCLUDED_core_pack_guidance_scoreterms_voids_penalty_energy_VoidsPenaltyVoxelGrid_hh

#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyVoxelGrid.fwd.hh>
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyVoxelGridTests.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/MathNTensor.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace voids_penalty_energy {

///@brief A 3D boolean array used for identifying core voxels in the VoidsPenaltyEnergy.
class VoidsPenaltyVoxelGrid : public utility::pointer::ReferenceCount {

	friend class ::VoidsPenaltyVoxelGridTests;
	friend class ::VoidsPenaltyVoxelGridTests_rotamer_setup;
	friend class ::VoidsPenaltyVoxelGridTests_rotamer_setup_2;

public:

	/// @brief Default constructor.
	VoidsPenaltyVoxelGrid();

	/// @brief Default destructor.
	virtual ~VoidsPenaltyVoxelGrid();

	/// @brief Copy constructor.
	VoidsPenaltyVoxelGrid(VoidsPenaltyVoxelGrid const & src);

	/// @brief Clone operator: copy this object and return an owning pointer to the copy.
	VoidsPenaltyVoxelGridOP clone() const;

public: //Public functions

	/// @brief Sets the voxel size and padding.
	/// @details Calls reset(), necessitating subsequent initialization from pose.  Computes half_voxel_size_ and voxel_volume_.
	void set_voxel_size_and_padding( core::Real const &size_in, core::Real const &padding_in );

	/// @brief Set the cone dotproduct cutoff and the cone distance cutoff.
	/// @detals The cone dotproduct cutoff is the cutoff value for the dot product of a cone vector and a cone base-test point vector below which
	/// we declare the test point not to be within the cone.  Effectively, this is the cone width.  Lower values make broader cones.  Default 0.1.
	/// Can range from 1.0 (infinitely thin cone) to -1.0 (full spherical volume), with 0.0 represeting all points on one side of the plane perpendicular
	/// to the cone vector.  The cone distance cutoff is the distance from the cone base that a test point must lie within to be in the cone.  This is
	/// automatically converted to (and stored as) the square of the distance by this function for efficient comparison operations.  the containing cones
	/// cutoff is the minimum number of cones that must contain a given voxel in order for that voxel to be considered to be "buried".
	/// @note Calls reset(), necessitating subsequent initialization from pose.
	void set_cone_parameters( core::Real const &cone_dotproduct_cutoff_in, core::Real const cone_distance_cutoff_in, core::Size const containing_cones_cutoff_in );

	/// @brief Reset this object.
	/// @details Clears stored volumetric data.  Does not clear settings, such as voxel size or cone dotproduct cutoff.
	void reset();

	/// @brief Set up voxel grid from a pose and determine which voxels are buried.
	/// @details This function must be called before compute_volumes_of_buried_rotamers().
	void set_up_voxel_grid_and_compute_burial( core::pose::Pose const &pose );

	/// @brief Remove voxels from voxel grid for positions that cannot repack.
	/// @details set_up_voxel_grid_and_compute_burial() must be called first.
	/// @note If symminfo == nullptr, then we assume this is the asymmetric case.  Otherwise, we
	/// use the symmetry information to aid the pruning.
	void prune_voxels_for_fixed_residues( core::pose::Pose const &pose, core::pack::rotamer_set::RotamerSets const &rotsets, core::conformation::symmetry::SymmetryInfoCOP symminfo, core::conformation::symmetry::SymmetricConformationCOP symmconf );

	/// @brief Given a rotamer set and a place to store rotamer volume information, computes volumes of buried parts of rotamers.
	/// @details Requires that set_up_voxel_grid_and_compute_burial() was called first. Input is a RotamersSet object.
	/// Output is a vector of maps of ResidueOP to floats.  The float values are the total volume of each rotamer
	/// that lies within the buried volume.
	/// @note If symminfo == nullptr, then we assume this is the asymmetric case (in which case symmconf can also be nullptr).  Otherwise, we
	/// use the symmetry information to simplify the calculation, returning rotamer volumes multiplied by the number of symmetry copies ONLY
	/// for the independent positions.
	void compute_volumes_of_buried_rotamers( core::pose::Pose const & pose, core::pack::rotamer_set::RotamerSets const &rotsets, utility::vector1< std::map< core::conformation::ResidueCOP, core::Real > > &rotamer_volume_maps, core::conformation::symmetry::SymmetryInfoCOP symminfo, core::conformation::symmetry::SymmetricConformationCOP symmconf );

	/// @brief Given a pose, compute the total volume of its buried residues.
	/// @details Requires that set_up_voxel_grid_and_compute_burial() was called first.
	core::Real compute_total_volume_of_current_residues( core::pose::Pose const &pose ) const;

	/// @brief Given a rotamer, compute the volume of the buried atoms.
	/// @details This is a bit approximate.  The algorithm is as follows:
	/// - Deterine the bounding box of the residue.
	/// - Loop through every voxel in the bounding box.
	/// - Check whether any atom in the residue overlaps the voxel.  If it does,
	/// add 1 voxel volume to the total volume.
	/// - Return the total volume.
	core::Real compute_volume_of_this_buried_rotamer( core::conformation::ResidueCOP rot, numeric::MathNTensor< bool, 3 > &untouched_voxels, core::Size &untouched_voxel_count ) const;

	/// @brief Given a rotamer, determine the bounding box (in voxel grid indices) of that rotamer.
	/// @details Requires that set_up_voxel_grid_and_compute_burial() was called first.
	/// @note If a rotamer exceeds the bounding box of the overall voxel grid, the edge indices are returned.  Outputs are bounding_box_start and bounding_box_end.
	void get_rotamer_bounding_box( core::conformation::ResidueCOP rot, utility::fixedsizearray1< core::Size, 3 > &bounding_box_start, utility::fixedsizearray1< core::Size, 3 > &bounding_box_end ) const;

	/// @brief Given a voxel index, the coordinates of an atom centre, and an atom radius, determine
	/// whether the atom overlaps the voxel.  Returns true if it does and false otherwise.
	/// @details Requires that set_up_voxel_grid_and_compute_burial() was called first.
	/// @note This works in three steps.  First, we determine whether the bonuding box of the atom
	/// and the voxel overlap, and return false if they don't.  Next, we determine whether the atom
	/// centre is in the voxel, and return true if it is.  Finally, we determine whether the voxel
	/// centre to atom centre distance is less than the atomic radius, and return true if and only if
	/// it is.
	bool atom_overlaps_voxel( utility::fixedsizearray1< core::Size, 3 > const &voxel_index, numeric::xyzVector< core::Real > const &atom_xyz, core::Real const atomic_radius ) const;

	/*
	// @brief Given two atomic radii and their positions, calculate and return the overlap integral.
	// @details Rapidly determines whether the atoms don't overlap, in which case 0 is returned.
	core::Real atomic_overlap_integral( core::Real const radius1, core::Real const radius2, numeric::xyzVector< core::Real > const &pos1, numeric::xyzVector< core::Real > const &pos2 ) const;
	*/

	/// @brief Only for debugging!  This function returns a pose with a grid of copper atoms representing the voxels.
	/// @details If only_buried_voxels is true, then the voxel grid only has copper atoms at the "true" positions; otherwise, the full grid
	/// is dumped to the pose.
	core::pose::PoseOP visualize_voxel_grid( bool const only_buried_voxels ) const;

	/// @brief Given the indices of a voxel (keeping in mind that indices are ZERO-based), get the 3D coordinates.
	/// @details Note: This returns the 3D coordinates of the CENTRE of the voxel.
	numeric::xyzVector< core::Real > get_voxel_coordinates( core::Size const index_x, core::Size const index_y, core::Size const index_z ) const;

	/// @brief Given a point in 3D space, determine which voxel it lies in and return the indices of that voxel.
	/// @details Returns false if the point lies outside of the voxel grid.  Returns true otherwise and sets
	/// output_indices to the indices of the voxel grid cell.  If return value is false, the out-of-range indices are set to the
	/// values of the edge of the voxel grid (0 or voxel_data_dimensions_[n]-1, where n is 1, 2, or 3 for x, y, and z respectively),
	/// and all other indices are set appropriately.
	/// @note The xyzVector is deliberately copied in instead of being passed by reference because I'd otherwise need to allocate
	/// a new xyzVector and copy values that I want to modify as part of the calculation.
	bool get_indices_of_voxel_from_coordinates( numeric::xyzVector<core::Real> point, utility::fixedsizearray1< core::Size, 3 > &output_indices ) const;

	/// @brief Get the total buried volume (including volume that cannot be reached by any rotamer), in cubic Angstroms.
	core::Real total_buried_volume() const;

	/// @brief Get the buried volume EXCLUDING volume that cannot be reached by any rotamer, in cubic Angstroms.
	core::Real reachable_buried_volume() const;

private: //Private functions

	/// @brief Given a pose, find the lower-left and upper-right corners of the bounding volume, with padding.
	void compute_bounding_box( core::pose::Pose const &pose, core::Real const &padding, numeric::xyzVector< core::Real> &lower_left_coords, numeric::xyzVector< core::Real> &upper_right_coords ) const;

	/// @brief Given lower-left, upper-right, and voxel size data, resize the voxel grid appropriately.
	/// @details This also updates upper_right_coords to be the centre of the maximum voxel (pushing it a little beyond the old bounding box).
	void resize_voxel_grid( numeric::xyzVector< core::Real> const &lower_left_coords, numeric::xyzVector< core::Real> &upper_right_coords, core::Real const &voxel_size, utility::fixedsizearray1< core::Size, 3 > &voxel_data_dimensions, numeric::MathNTensor< bool, 3 > &voxel_data ) const;

	/// @brief Given a pose and a voxel grid that is already set up for the pose, compute whether each point is buried or not.
	void compute_burial( core::pose::Pose const &pose, utility::fixedsizearray1< core::Size, 3 > const &voxel_data_dimensions, numeric::MathNTensor< bool, 3 > &voxel_data, core::Size & buried_volume) const;

	/// @brief Given a pose and a point in space, determine whether the point is buried by the side-chain
	/// cones method.
	bool is_buried( numeric::xyzVector<core::Real> const &coords, utility::vector1< numeric::xyzVector<core::Real> > const &conevects, utility::vector1< numeric::xyzVector<core::Real> > const &conebases, utility::vector1< bool > const &skip_list ) const;

	/// @brief Given a residue index, get (a) the same residue index back in the asymmetric case, (b) the same residue index back in the
	/// symmetric case if this residue is independent, and (c) the residue index on which this residue depends in the symmetric case if
	/// this residue is dependent.
	core::Size get_current_res_symmetric( core::Size const index_in, core::conformation::symmetry::SymmetryInfoCOP symminfo ) const;

private:

	/// @brief The length, width, or height of a voxel, in Angstroms.
	/// @details Defaults to 0.5 A.  Smaller is not necessarily better.
	core::Real voxel_size_;

	/// @brief The voxel size divided by 2.  Cached to avoid repeated computation.
	core::Real half_voxel_size_;

	/// @brief The volume of a voxel, in cubic Angstroms.  Cached to avoid repeated computation.
	core::Real voxel_volume_;

	/// @brief The padding around the pose for the voxel grid.
	core::Real padding_;

	/// @brief The x, y, and z coordinates of the lower-left corner of the volume considered.
	numeric::xyzVector < core::Real > lower_left_coords_;

	/// @brief The x, y, and z coordinates of the upper-right corner of the volume considered.
	numeric::xyzVector < core::Real > upper_right_coords_;

	/// @brief The actual volumetric data.
	numeric::MathNTensor< bool, 3 > voxel_data_;

	/// @brief The dimensions, in cells, of the voxel data.
	utility::fixedsizearray1< core::Size, 3 > voxel_data_dimensions_;

	/// @brief The total buried volume (including volume that cannot be reached by any rotamer).
	core::Size total_buried_volume_;

	/// @brief The buried volume EXCLUDING volume that cannot be reached by any rotamer.
	core::Size reachable_buried_volume_;

	/// @brief Has this object been set up from a pose (i.e. has the voxel grid been initialized and the buried volume computed)?
	/// @details This step must precede rotamer setup.
	bool pose_setup_complete_;

	/// @brief Has the prune_voxels_for_fixed_residues() function been called?
	/// @details This must precede rotamer setup.
	bool non_packable_volume_pruned_;

	/// @brief Has this object been used to set up rotmers?
	/// @details This step must follow setup from pose, and can only be performed once.
	bool rotamer_setup_complete_;

	/// @brief The cutoff value for the dot product of a cone vector and a cone base-test point vector below which we declare the test point not to be
	/// within the cone.  Effectively, this is the cone width.  Lower values make broader cones.  Default 0.1.  Can range from 1.0 (infinitely thin
	/// cone) to -1.0 (full spherical volume), with 0.0 represeting all points on one side of the plane perpendicular to the cone vector.
	core::Real cone_dotproduct_cutoff_;

	/// @brief The square of the cutoff value for the distance from the cone base at which we are considered no longer to be within the cone.  Defaults
	/// to 64.0 (for a cutoff distance of 8.0).
	core::Real cone_distance_cutoff_sq_;

	/// @brief The minimum number of cones that must contain a particular voxel for that voxel to be considered "buried".
	/// @details Defaults to 6.
	core::Size containing_cones_cutoff_;

};


} //core
} //guidance_scoreterms
} //pack
} //voids_penalty_energy

#endif //INCLUDED_core_pack_guidance_scoreterms_voids_penalty_energy_VoidsPenaltyVoxelGrid_hh

