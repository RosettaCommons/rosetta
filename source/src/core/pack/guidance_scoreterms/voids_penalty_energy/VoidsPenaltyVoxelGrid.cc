// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyVoxelGrid.cc
/// @brief A 3D boolean array used for identifying core voxels in the VoidsPenaltyEnergy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyVoxelGrid.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>

#include <numeric/constants.hh>

#include <basic/Tracer.hh>


//Default values:
#define DEFAULT_VOXEL_SIZE 0.5
#define DEFAULT_PADDING 1.0
#define DEFAULT_CONE_DOTPRODUCT_CUTOFF 0.1
#define DEFAULT_CONE_DIST_CUTOFF_SQUARED 64.0
#define DEFAULT_CONTAINING_CONES_CUTOFF 6

static basic::Tracer TR( "core.pack.voids_penalty_energy.VoidsPenaltyVoxelGrid" );

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace voids_penalty_energy {

/// @brief Default constructor.
VoidsPenaltyVoxelGrid::VoidsPenaltyVoxelGrid():
	utility::pointer::ReferenceCount(),
	voxel_size_(DEFAULT_VOXEL_SIZE),
	half_voxel_size_(DEFAULT_VOXEL_SIZE/2.0),
	voxel_volume_( pow(DEFAULT_VOXEL_SIZE, 3) ),
	padding_(DEFAULT_PADDING),
	lower_left_coords_( 0.0, 0.0, 0.0 ),
	upper_right_coords_( 0.0, 0.0, 0.0 ),
	voxel_data_(),
	voxel_data_dimensions_(1),
	total_buried_volume_(0),
	reachable_buried_volume_(0),
	pose_setup_complete_(false),
	non_packable_volume_pruned_(false),
	rotamer_setup_complete_(false),
	cone_dotproduct_cutoff_(DEFAULT_CONE_DOTPRODUCT_CUTOFF),
	cone_distance_cutoff_sq_(DEFAULT_CONE_DIST_CUTOFF_SQUARED),
	containing_cones_cutoff_(DEFAULT_CONTAINING_CONES_CUTOFF)
{}

/// @brief Default destructor.
VoidsPenaltyVoxelGrid::~VoidsPenaltyVoxelGrid(){}

/// @brief Copy constructor.
VoidsPenaltyVoxelGrid::VoidsPenaltyVoxelGrid( VoidsPenaltyVoxelGrid const &src ) :
	voxel_size_(src.voxel_size_),
	half_voxel_size_(src.half_voxel_size_),
	voxel_volume_(src.voxel_volume_),
	padding_(src.padding_),
	lower_left_coords_(src.lower_left_coords_),
	upper_right_coords_(src.upper_right_coords_),
	voxel_data_(src.voxel_data_),
	voxel_data_dimensions_(src.voxel_data_dimensions_),
	total_buried_volume_(src.total_buried_volume_),
	reachable_buried_volume_(src.reachable_buried_volume_),
	pose_setup_complete_(src.pose_setup_complete_),
	non_packable_volume_pruned_(src.non_packable_volume_pruned_),
	rotamer_setup_complete_(src.rotamer_setup_complete_),
	cone_dotproduct_cutoff_(src.cone_dotproduct_cutoff_),
	cone_distance_cutoff_sq_(src.cone_distance_cutoff_sq_),
	containing_cones_cutoff_(src.containing_cones_cutoff_)
{}


/// @brief Clone operator: copy this object and return an owning pointer to the copy.
VoidsPenaltyVoxelGridOP
VoidsPenaltyVoxelGrid::clone() const {
	return VoidsPenaltyVoxelGridOP( new VoidsPenaltyVoxelGrid( *this ) );
}

/// @brief Sets the voxel size and padding.
/// @details Calls reset(), necessitating subsequent initialization from pose.  Computes half_voxel_size_ and voxel_volume_.
void
VoidsPenaltyVoxelGrid::set_voxel_size_and_padding(
	core::Real const &size_in,
	core::Real const &padding_in
) {
	static std::string const errmsg( "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::set_voxel_size_and_padding():  " );
	runtime_assert_string_msg( size_in > 0.0, errmsg + "The voxel size must be greater than zero." );
	runtime_assert_string_msg( padding_in > 0.0, errmsg + "The volume padding must be greater than zero." );
	voxel_size_ = size_in;
	half_voxel_size_ = voxel_size_ / 2.0;
	voxel_volume_ = pow(voxel_size_, 3);
	padding_ = padding_in;
	reset();
}

/// @brief Set the cone dotproduct cutoff and the cone distance cutoff.
/// @detals The cone dotproduct cutoff is the cutoff value for the dot product of a cone vector and a cone base-test point vector below which
/// we declare the test point not to be within the cone.  Effectively, this is the cone width.  Lower values make broader cones.  Default 0.1.
/// Can range from 1.0 (infinitely thin cone) to -1.0 (full spherical volume), with 0.0 represeting all points on one side of the plane perpendicular
/// to the cone vector.  The cone distance cutoff is the distance from the cone base that a test point must lie within to be in the cone.  This is
/// automatically converted to (and stored as) the square of the distance by this function for efficient comparison operations.
/// @note Calls reset(), necessitating subsequent initialization from pose.
void
VoidsPenaltyVoxelGrid::set_cone_parameters(
	core::Real const &cone_dotproduct_cutoff_in,
	core::Real const cone_distance_cutoff_in,
	core::Size const containing_cones_cutoff_in
) {
	cone_dotproduct_cutoff_ = cone_dotproduct_cutoff_in;
	cone_distance_cutoff_sq_ = pow( cone_distance_cutoff_in, 2);
	containing_cones_cutoff_ = containing_cones_cutoff_in;
	reset();
}


/// @brief Reset this object.
/// @details Clears stored volumetric data.  Does not clear settings, such as voxel size or cone dotproduct cutoff.
void
VoidsPenaltyVoxelGrid::reset() {
	lower_left_coords_ = numeric::xyzVector < core::Real >(0.0, 0.0, 0.0);
	upper_right_coords_ = numeric::xyzVector < core::Real >(0.0, 0.0, 0.0);
	voxel_data_dimensions_ = utility::fixedsizearray1< core::Size, 3 >( 1 );
	voxel_data_ = numeric::MathNTensor< bool, 3 >( voxel_data_dimensions_, false );
	total_buried_volume_ = 0;
	reachable_buried_volume_ = 0;
	pose_setup_complete_ = false;
	non_packable_volume_pruned_=false;
	rotamer_setup_complete_ = false;
}

/// @brief Set up voxel grid from a pose and determine which voxels are buried.
void
VoidsPenaltyVoxelGrid::set_up_voxel_grid_and_compute_burial(
	core::pose::Pose const &pose
) {
	runtime_assert_string_msg( !pose_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::set_up_voxel_grid_and_compute_burial(): Setup from a pose has already occurred.  Please reset this object before attempting this again." );
	runtime_assert_string_msg( !rotamer_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::set_up_voxel_grid_and_compute_burial(): Rotamer setup has already occurred.  Please reset this object before attempting setup from a pose again." );
	pose_setup_complete_ = true; //Mark setup as complete.
	compute_bounding_box( pose, padding_, lower_left_coords_, upper_right_coords_ );
	resize_voxel_grid( lower_left_coords_, upper_right_coords_, voxel_size_, voxel_data_dimensions_, voxel_data_ );
	compute_burial( pose, voxel_data_dimensions_, voxel_data_, total_buried_volume_);
}

/// @brief Remove voxels from voxel grid for positions that cannot repack.
/// @details set_up_voxel_grid_and_compute_burial() must be called first.
/// @note If symminfo == nullptr, then we assume this is the asymmetric case.  Otherwise, we
/// use the symmetry information to aid the pruning.
void
VoidsPenaltyVoxelGrid::prune_voxels_for_fixed_residues(
	core::pose::Pose const &pose,
	core::pack::rotamer_set::RotamerSets const &rotsets,
	core::conformation::symmetry::SymmetryInfoCOP symminfo,
	core::conformation::symmetry::SymmetricConformationCOP symmconf
) {
	runtime_assert_string_msg( pose_setup_complete_ , "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::prune_voxels_for_fixed_residues(): The voxel grid has not yet been set up from a pose." );
	runtime_assert_string_msg( !non_packable_volume_pruned_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::prune_voxels_for_fixed_residues():  Voxels for non-packable positions were already pruned.");

	non_packable_volume_pruned_ = true;

	TR << "Pruning voxels for non-packable residues.  Total volume before pruning: " << total_buried_volume_ << " cubic Angstroms." << std::endl;

	// Get mirror symmetry info, if available.
	core::conformation::symmetry::MirrorSymmetricConformationCOP mirrorsymmconf( nullptr );
	if ( symminfo != nullptr ) {
		debug_assert( symmconf != nullptr ); //Needs to be true.
		mirrorsymmconf = utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::MirrorSymmetricConformation const >( symmconf );
	}

	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		core::Size const current_res( get_current_res_symmetric( i, symminfo ) );
		if ( !rotsets.has_rotamer_set_for_residue(current_res) && (symminfo == nullptr || symminfo->bb_is_independent(current_res) ) && !pose.residue_type(current_res).is_virtual_residue() && !pose.residue_type(current_res).is_inverted_virtual_residue() ) {
			if ( symminfo == nullptr ) { //The asymmetric case is simplest:
				compute_volume_of_this_buried_rotamer( pose.conformation().residue_cop(current_res), voxel_data_, total_buried_volume_ ); //The volume calculation is discarded; only the effects on voxel_data_ and total_buried_volume_ matter here.
			} else { //The symmetric case is more complex.s
				compute_volume_of_this_buried_rotamer( pose.conformation().residue_cop(current_res), voxel_data_, total_buried_volume_ ); //The volume calculation is discarded; only the effects on voxel_data_ and total_buried_volume_ matter here.
				for ( core::Size k(1), kmax( symminfo->num_bb_clones()); k<=kmax; ++k ) {
					core::Size const symmpos( symminfo->bb_clones(current_res)[k] );
					compute_volume_of_this_buried_rotamer( pose.conformation().residue_cop(symmpos), voxel_data_, total_buried_volume_ ); //The volume calculation is discarded; only the effects on voxel_data_ and total_buried_volume_ matter here.  Here, we do it for the symmetry copies.
				}
			}
		}
	}

	TR << "Total volume after pruning: " << total_buried_volume_ << " cubic Angstroms." << std::endl;
}

/// @brief Given a rotamer set and a place to store rotamer volume information, computes volumes of buried parts of rotamers.
/// @details Requires that set_up_voxel_grd_and_compute_burial was called first. Input is a RotamersSet object.
/// Output is a vector of maps of ResidueOP to floats.  The float values are the total volume of each rotamer
/// that lies within the buried volume.
/// @note If symminfo == nullptr, then we assume this is the asymmetric case (in which case symmconf can also be nullptr).  Otherwise, we
/// use the symmetry information to simplify the calculation, returning rotamer volumes multiplied by the number of symmetry copies ONLY
/// for the independent positions.
void
VoidsPenaltyVoxelGrid::compute_volumes_of_buried_rotamers(
	core::pose::Pose const &pose,
	core::pack::rotamer_set::RotamerSets const &rotsets,
	utility::vector1< std::map< core::conformation::ResidueCOP, core::Real > > &rotamer_volume_maps,
	core::conformation::symmetry::SymmetryInfoCOP symminfo,
	core::conformation::symmetry::SymmetricConformationCOP symmconf
) {
	runtime_assert_string_msg( pose_setup_complete_ && non_packable_volume_pruned_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::compute_volumes_of_buried_rotamers(): The voxel grid has not yet been set up from a pose, or non-packable volume has not yet been pruned." );
	debug_assert( rotsets.total_residue() > 0 );
	rotamer_volume_maps.clear();
	rotamer_volume_maps.reserve( rotsets.total_residue() );

	numeric::MathNTensor< bool, 3 > untouched_voxels(voxel_data_);
	core::Size untouched_voxel_count( total_buried_volume_ );

	// Get mirror symmetry info, if available.
	core::conformation::symmetry::MirrorSymmetricConformationCOP mirrorsymmconf( nullptr );
	if ( symminfo != nullptr ) {
		debug_assert( symmconf != nullptr ); //Needs to be true.
		mirrorsymmconf = utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::MirrorSymmetricConformation const >( symmconf );
	}

	for ( core::Size i(1), imax(rotsets.total_residue()); i<=imax; ++i ) {
		std::map< core::conformation::ResidueCOP, core::Real > curmap;
		if ( rotsets.has_rotamer_set_for_residue(i) && (symminfo == nullptr || symminfo->bb_is_independent(i)) ) { //There is no rotamer set for those positions that aren't packable.
			core::pack::rotamer_set::RotamerSetCOP curset( rotsets.rotamer_set_for_residue(i) );
			for ( core::Size j(1), jmax(curset->num_rotamers() + 1); j<=jmax; ++j ) {
				core::conformation::ResidueCOP curres( j == jmax ? pose.conformation().residue_cop(i) : curset->rotamer(j) );
				if ( j == jmax ) {
					if ( curmap.count(curres) ) continue; //The current residue might already be processed.
				} else {
					debug_assert( !curmap.count(curres)  );
				}
				if ( symminfo != nullptr ) {
					core::Real accumulator( compute_volume_of_this_buried_rotamer( curres, untouched_voxels, untouched_voxel_count ) );
					for ( core::Size k(1), kmax(symminfo->num_bb_clones()); k<=kmax; ++k ) {
						core::Size const symmpos( symminfo->bb_clones(curres->seqpos())[k] );
						core::conformation::ResidueOP rotamer_clone( mirrorsymmconf == nullptr ? curres->clone() : curres->clone_flipping_chirality( *( symmconf->residue_type_set_for_conf( curres->type().mode() ) ) ) );
						for ( core::Size ia(1), iamax(rotamer_clone->natoms()); ia<=iamax; ++ia ) {
							rotamer_clone->set_xyz(ia, symmconf->apply_transformation_norecompute( curres->xyz(ia), curres->seqpos(), symmpos ) );
						}
						accumulator += compute_volume_of_this_buried_rotamer( rotamer_clone, untouched_voxels, untouched_voxel_count );
					}

					curmap[curres] = accumulator;
				} else {
					curmap[curres] = compute_volume_of_this_buried_rotamer( curres, untouched_voxels, untouched_voxel_count );
				}
			}
		}
		rotamer_volume_maps.push_back(curmap);
	}

	debug_assert( untouched_voxel_count <= total_buried_volume_ ); //Should be true -- can't overshoot past zero.
	reachable_buried_volume_ = total_buried_volume_ - untouched_voxel_count; //Update the amount of buried volume that can be touched by at least one rotamer.

	rotamer_setup_complete_ = true;
}

/// @brief Given a pose, compute the total volume of its buried residues.
/// @details Requires that set_up_voxel_grid_and_compute_burial() was called first.
core::Real
VoidsPenaltyVoxelGrid::compute_total_volume_of_current_residues(
	core::pose::Pose const &pose
) const {
	runtime_assert_string_msg( pose_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::compute_total_volume_of_current_residues(): The voxel grid has not yet been set up from a pose." );
	core::Size untouched_voxel_count( total_buried_volume_ );

	// Get the symmetry information, iff this is a symmetric pose.
	/*core::conformation::symmetry::SymmetryInfoCOP synninfo( nullptr );
	core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
	if( symmconf != nullptr ) {
	symminfo = symmconf->Symmetry_Info();
	}*/

	numeric::MathNTensor< bool, 3 > untouched_voxels(voxel_data_);
	core::Real accumulator(0.0);
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) { //Loop through current residues in the pose
		accumulator += compute_volume_of_this_buried_rotamer( pose.conformation().residue_cop(i), untouched_voxels, untouched_voxel_count );
	}

	return accumulator;
}

/// @brief Given a rotamer, compute the volume of the buried atoms.
/// @details This is a bit approximate.  The algorithm is as follows:
/// - Deterine the bounding box of the residue.
/// - Loop through every voxel in the bounding box.
/// - Check whether any atom in the residue overlaps the voxel.  If it does,
/// add 1 voxel volume to the total volume.
/// - Return the total volume.
core::Real
VoidsPenaltyVoxelGrid::compute_volume_of_this_buried_rotamer(
	core::conformation::ResidueCOP rot,
	numeric::MathNTensor< bool, 3 > &untouched_voxels,
	core::Size &untouched_voxel_count
) const {
	//utility::vector1< core::Size > buried_indices;
	//buried_indices.reserve(rot->natoms());

	core::Size rotamer_volume(0); //The number of voxels that this rotamer touches.

	utility::fixedsizearray1< core::Size, 3 > bounding_box_start, bounding_box_end; //Coordinates of corners of bounding volume for this residue.
	get_rotamer_bounding_box( rot, bounding_box_start, bounding_box_end );

	//Delete the following -- for debugging only.
	//TR << rot->name3() << rot->seqpos() << " bounding box: [" << bounding_box_start[1] << ", " << bounding_box_start[2] << ", " << bounding_box_start[1] << "] to [" << bounding_box_end[1] << ", " << bounding_box_end[2] << ", " << bounding_box_end[3] << "]." << std::endl;

	//Loop through voxels that overlap this atom:
	utility::fixedsizearray1< core::Size, 3 > curindex;
	for ( core::Size i(bounding_box_start[1]); i <= bounding_box_end[1]; ++i ) {
		for ( core::Size j(bounding_box_start[2]); j <= bounding_box_end[2]; ++j ) {
			for ( core::Size k(bounding_box_start[3]) ; k <= bounding_box_end[3]; ++k ) {
				curindex[1] = i; curindex[2] = j; curindex[3] = k;
				if ( !voxel_data_( curindex ) ) continue; //Only count buried.
				for ( core::Size ia(1), iamax( rot->natoms() ); ia<=iamax; ++ia ) {
					if ( atom_overlaps_voxel( curindex, rot->xyz( ia ), rot->atom_type( ia ).lj_radius() ) ) {
						if ( untouched_voxels(curindex) ) {
							--untouched_voxel_count;
							untouched_voxels(curindex) = false;
						}
						++rotamer_volume;
						break; //Break out of loop over atoms, but NOT loop over voxels.
					}
				} //loop over atoms
			} //loop over voxels -- z
		} //loop over voxels -- y
	} //loop over voxels -- x

	//Delete the following -- for debugging only.
	//TR << rot->name3() << rot->seqpos() << " volume: " << rotamer_volume << " voxels." << std::endl;

	return static_cast<core::Real>(rotamer_volume)*voxel_volume_;
}

/// @brief Given a rotamer, determine the bounding box (in voxel grid indices) of that rotamer.
/// @details Requires that set_up_voxel_grid_and_compute_burial() was called first.
/// @note If a rotamer exceeds the bounding box of the overall voxel grid, the edge indices are returned.  Outputs are bounding_box_start and bounding_box_end.
void
VoidsPenaltyVoxelGrid::get_rotamer_bounding_box(
	core::conformation::ResidueCOP rot,
	utility::fixedsizearray1< core::Size, 3 > &bounding_box_start,
	utility::fixedsizearray1< core::Size, 3 > &bounding_box_end
) const {
	runtime_assert_string_msg( pose_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::get_rotamer_bounding_box(): The voxel grid has not yet been initialized from a pose." );
	runtime_assert_string_msg( rot->natoms() > 0, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::get_rotamer_bounding_box(): The number of atoms in the rotamer is zero!" );

	utility::fixedsizearray1< core::Size, 3 > bb1, bb2; //Temporary storage for voxel indicies of bounding box of atom.
	numeric::xyzVector< core::Real > xyz1, xyz2; //Temporary storage of 3D coordinates of bounding box corners of atom.

	for ( core::Size ia(1), iamax(rot->natoms()); ia<=iamax; ++ia ) {
		core::Size const atom_rad( rot->atom_type(ia).lj_radius() );
		xyz1 = rot->xyz(ia);
		xyz2 = xyz1;
		xyz1.x( xyz1.x() - atom_rad ); xyz1.y( xyz1.y() - atom_rad ); xyz1.z( xyz1.z() - atom_rad );
		xyz2.x( xyz2.x() + atom_rad ); xyz2.y( xyz2.y() + atom_rad ); xyz2.z( xyz2.z() + atom_rad );
		get_indices_of_voxel_from_coordinates( xyz1, bb1 );
		get_indices_of_voxel_from_coordinates( xyz2, bb2 );
		if ( ia==1 || bb1[1] < bounding_box_start[1] ) bounding_box_start[1] = bb1[1];
		if ( ia==1 || bb1[2] < bounding_box_start[2] ) bounding_box_start[2] = bb1[2];
		if ( ia==1 || bb1[3] < bounding_box_start[3] ) bounding_box_start[3] = bb1[3];
		if ( ia==1 || bb2[1] > bounding_box_end[1] ) bounding_box_end[1] = bb2[1];
		if ( ia==1 || bb2[2] > bounding_box_end[2] ) bounding_box_end[2] = bb2[2];
		if ( ia==1 || bb2[3] > bounding_box_end[3] ) bounding_box_end[3] = bb2[3];
	}
}

/// @brief Given a voxel index, the coordinates of an atom centre, and an atom radius, determine
/// whether the atom overlaps the voxel.  Returns true if it does and false otherwise.
/// @note This works in three steps.  First, we determine whether the bonuding box of the atom
/// and the voxel overlap, and return false if they don't.  Next, we determine whether the atom
/// centre is in the voxel, and return true if it is.  Finally, we determine whether the voxel
/// centre to atom centre distance is less than the atomic radius, and return true if and only if
/// it is.
bool
VoidsPenaltyVoxelGrid::atom_overlaps_voxel(
	utility::fixedsizearray1< core::Size, 3 > const &voxel_index,
	numeric::xyzVector< core::Real > const &atom_xyz,
	core::Real const atomic_radius
) const {
	runtime_assert_string_msg( pose_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::atom_overlaps_voxel(): The voxel grid has not yet been initialized from a pose." );

	numeric::xyzVector<core::Real> const voxel_centre( get_voxel_coordinates( voxel_index[1], voxel_index[2], voxel_index[3] ) );
	numeric::xyzVector<core::Real> voxel_bounds_min;
	numeric::xyzVector<core::Real> voxel_bounds_max;
	core::Real const radval( atomic_radius + half_voxel_size_ );

	voxel_bounds_min.x( voxel_centre.x() - radval );
	voxel_bounds_min.y( voxel_centre.y() - radval );
	voxel_bounds_min.z( voxel_centre.z() - radval );
	voxel_bounds_max.x( voxel_centre.x() + radval );
	voxel_bounds_max.y( voxel_centre.y() + radval );
	voxel_bounds_max.z( voxel_centre.z() + radval );

	//Delete the following -- for debugging only.
	//TR << voxel_index[1] << " " << voxel_index[2] << " " << voxel_index[3] << " | " << voxel_bounds_min.x() << " " << voxel_bounds_min.y() << " " << voxel_bounds_min.z() << " | " << voxel_bounds_max.x() << " " << voxel_bounds_max.y() << " " << voxel_bounds_max.z() << " | " << atom_xyz.x() << " " << atom_xyz.y() << " " << atom_xyz.z() << std::endl;

	// Do the bounding box of the atom and the voxel overlap?
	if ( voxel_bounds_min.x() < atom_xyz.x() && atom_xyz.x() < voxel_bounds_max.x() &&
			voxel_bounds_min.y() < atom_xyz.y() && atom_xyz.y() < voxel_bounds_max.y() &&
			voxel_bounds_min.z() < atom_xyz.z() && atom_xyz.z() < voxel_bounds_max.z()
			) {
		//TR << "Atom_xyz=" << atom_xyz.x() << "," << atom_xyz.y() << "," << atom_xyz.z() <<  " with radius " << atomic_radius << " is in bounding box " << voxel_bounds_min.x() << "," << voxel_bounds_min.y() << "," << voxel_bounds_min.z() << " -- " << voxel_bounds_max.x() << "," << voxel_bounds_max.y() << "," << voxel_bounds_max.z() << " with indices " << voxel_index[1] << "," << voxel_index[2] << "," << voxel_index[3] << std::endl;

		//Contract the voxel again.
		voxel_bounds_min.x( voxel_bounds_min.x() + atomic_radius );
		voxel_bounds_min.y( voxel_bounds_min.y() + atomic_radius );
		voxel_bounds_min.z( voxel_bounds_min.z() + atomic_radius );
		voxel_bounds_max.x( voxel_bounds_max.x() - atomic_radius );
		voxel_bounds_max.y( voxel_bounds_max.y() - atomic_radius );
		voxel_bounds_max.z( voxel_bounds_max.z() - atomic_radius );

		// Is the atom centre inside the voxel?  If so, return true.
		if ( voxel_bounds_min.x() < atom_xyz.x() && atom_xyz.x() < voxel_bounds_max.x() &&
				voxel_bounds_min.y() < atom_xyz.y() && atom_xyz.y() < voxel_bounds_max.y() &&
				voxel_bounds_min.z() < atom_xyz.z() && atom_xyz.z() < voxel_bounds_max.z()
				) {
			return true;
		}

		// If we reach here, then the atom bounding box and voxel bounding box overlap, but the atom is not inside the voxel.
		// So we return true iff the atom centre - voxel centre distance is less than the atom radius.
		if ( voxel_centre.distance_squared( atom_xyz ) < pow(atomic_radius, 2) ) {
			return true;
		}

	}

	// If we reach here, then the bounding box of the atom and the voxel don't overlap.
	return false;
}


/*
// @brief Given two atomic radii and their positions, calculate and return the overlap integral.
// @details Rapidly determines whether the atoms don't overlap, in which case 0 is returned.
core::Real
VoidsPenaltyVoxelGrid::atomic_overlap_integral(
core::Real const radius1,
core::Real const radius2,
numeric::xyzVector< core::Real > const &pos1,
numeric::xyzVector< core::Real > const &pos2
) const {
core::Real const dist_sq( pos1.distance_squared(pos2) );
if ( dist_sq > pow(radius1+radius2, 2) ) return 0.0;
core::Real const dist( std::sqrt(dist_sq) );
// The following expression comes from Wolfram Mathworld -- I didn't bother to derive it myself. --VKM
// See http://mathworld.wolfram.com/Sphere-SphereIntersection.html.
return numeric::constants::d::pi*pow( radius1 + radius2 - dist, 2)*( dist_sq + 2 * dist * radius2 - 3 * pow( radius2, 2 ) + 2 * dist * radius1 + 6 * radius1 * radius2 - 3 * pow( radius1, 2 ) ) / (12.0 * dist);
}
*/


/// @brief Only for debugging!  This function returns a pose with a grid of copper atoms representing the voxels.
/// @details If only_buried_voxels is true, then the voxel grid only has copper atoms at the "true" positions; otherwise, the full grid
/// is dumped to the pose.
core::pose::PoseOP
VoidsPenaltyVoxelGrid::visualize_voxel_grid(
	bool const only_buried_voxels
) const {
	runtime_assert_string_msg( pose_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::visualize_voxel_grid(): The voxel grid has not yet been initialized from a pose." );

	core::pose::PoseOP pose( new core::pose::Pose );

	core::chemical::ResidueTypeCOP copper_type( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t )->get_representative_type_name3("CU") );

	for ( core::Size i(0); i<voxel_data_dimensions_[1]; ++i ) {
		for ( core::Size j(0); j<voxel_data_dimensions_[2]; ++j ) {
			for ( core::Size k(0); k<voxel_data_dimensions_[3]; ++k ) {
				if ( !only_buried_voxels || voxel_data_(i,j,k) ) {
					numeric::xyzVector< core::Real > const curpos( get_voxel_coordinates(i, j, k) );
					core::conformation::ResidueOP rsd( new core::conformation::Residue( copper_type, true ));
					rsd->set_xyz( 1, curpos );
					pose->append_residue_by_jump( *rsd, 1 );
				}
			}
		}
	}
	return pose;
}

/// @brief Given the indices of a voxel (keeping in mind that indices are ZERO-based), get the 3D coordinates.
/// @details Note: This returns the 3D coordinates of the CENTRE of the voxel.
numeric::xyzVector< core::Real >
VoidsPenaltyVoxelGrid::get_voxel_coordinates(
	core::Size const index_x,
	core::Size const index_y,
	core::Size const index_z
) const {
	runtime_assert_string_msg( pose_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::get_voxel_coordinates(): The voxel grid has not yet been initialized from a pose." );
	debug_assert( index_x < voxel_data_dimensions_[1] );
	debug_assert( index_y < voxel_data_dimensions_[2] );
	debug_assert( index_z < voxel_data_dimensions_[3] );
	return numeric::xyzVector< core::Real > (
		static_cast<core::Real>(index_x) * voxel_size_ + lower_left_coords_.x(),
		static_cast<core::Real>(index_y) * voxel_size_ + lower_left_coords_.y(),
		static_cast<core::Real>(index_z) * voxel_size_ + lower_left_coords_.z()
	);
}

/// @brief Given a point in 3D space, determine which voxel it lies in and return the indices of that voxel.
/// @details Returns false if the point lies outside of the voxel grid.  Returns true otherwise and sets
/// output_indices to the indices of the voxel grid cell.  If return value is false, the out-of-range indices are set to the
/// values of the edge of the voxel grid (0 or voxel_data_dimensions_[n]-1, where n is 1, 2, or 3 for x, y, and z respectively),
/// and all other indices are set appropriately.
/// @note The xyzVector is deliberately copied in instead of being passed by reference because I'd otherwise need to allocate
/// a new xyzVector and copy values that I want to modify as part of the calculation.
bool
VoidsPenaltyVoxelGrid::get_indices_of_voxel_from_coordinates(
	numeric::xyzVector<core::Real> point,
	utility::fixedsizearray1< core::Size, 3 > &output_indices
) const {
	bool in_box(true); //return value
	runtime_assert_string_msg( pose_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::get_indices_of_voxel_from_coordinates(): The voxel grid has not yet been initialized from a pose." );
	if (
			point.x() > ( upper_right_coords_.x() + half_voxel_size_ ) ||
			point.y() > ( upper_right_coords_.y() + half_voxel_size_ ) ||
			point.z() > ( upper_right_coords_.z() + half_voxel_size_ ) ||
			point.x() < ( lower_left_coords_.x() - half_voxel_size_ ) ||
			point.y() < ( lower_left_coords_.y() - half_voxel_size_ ) ||
			point.z() < ( lower_left_coords_.z() - half_voxel_size_ )
			) {
		in_box = false;
	}

	point -= lower_left_coords_; //Subtract off the lower-left corner of the box.  The position point is now the relative position within the box, starting from 0 at the edge.
	point.x( point.x() + half_voxel_size_ ); //Need to subtract off half a voxel more.  The lower-left coordintes are the coordinate of the CENTRE of the lower-left voxel.
	point.y( point.y() + half_voxel_size_ );
	point.z( point.z() + half_voxel_size_ );
	point /= voxel_size_; //Divide by voxel size.  The position point is now the fractional voxel position within the box.

	output_indices[1] = static_cast< core::Size >( std::min( static_cast< signed long >(voxel_data_dimensions_[1]) - static_cast< signed long >(1), std::max( static_cast< signed long >(0), static_cast< signed long >( std::floor( point.x() ) ) ) ) );
	output_indices[2] = static_cast< core::Size >( std::min( static_cast< signed long >(voxel_data_dimensions_[2]) - static_cast< signed long >(1), std::max( static_cast< signed long >(0), static_cast< signed long >( std::floor( point.y() ) ) ) ) );
	output_indices[3] = static_cast< core::Size >( std::min( static_cast< signed long >(voxel_data_dimensions_[3]) - static_cast< signed long >(1), std::max( static_cast< signed long >(0), static_cast< signed long >( std::floor( point.z() ) ) ) ) );

	return in_box;
}

/// @brief Get the total buried volume (including volume that cannot be reached by any rotamer), in cubic Angstroms.
core::Real
VoidsPenaltyVoxelGrid::total_buried_volume() const {
	runtime_assert_string_msg( pose_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::total_buried_volume(): The voxel grid has not yet been initialized from a pose." );
	return static_cast<core::Real>( total_buried_volume_) * voxel_volume_;
}

/// @brief Get the buried volume EXCLUDING volume that cannot be reached by any rotamer, in cubic Angstroms.
core::Real
VoidsPenaltyVoxelGrid::reachable_buried_volume() const {
	runtime_assert_string_msg( rotamer_setup_complete_, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::reachable_buried_volume(): The reachable volume calculation requires that rotamers are first analyzed using VoidsPenaltyVoxelGrid::compute_volumes_of_buried_rotamers()." );
	return static_cast<core::Real>( reachable_buried_volume_) * voxel_volume_;
}

/// @brief Given a pose, find the lower-left and upper-right corners of the bounding volume, with padding.
void
VoidsPenaltyVoxelGrid::compute_bounding_box(
	core::pose::Pose const &pose,
	core::Real const &padding,
	numeric::xyzVector< core::Real> &lower_left_coords,
	numeric::xyzVector< core::Real> &upper_right_coords
) const {
	TR << "Computing pose bounding box with padding of " << padding << " Angstroms..." << std::endl;

	runtime_assert_string_msg( pose.total_residue() > 0, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::compute_bounding_box():  The pose has no residues!" );
	lower_left_coords = pose.residue(1).xyz(1);
	upper_right_coords = lower_left_coords;
	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) { //Loop through all residues
		for ( core::Size ia(1), iamax(pose.residue_type(ir).natoms()); ia<=iamax; ++ia ) {
			core::id::AtomID const atomid(ia, ir);
			numeric::xyzVector< core::Real > const & curxyz( pose.xyz( atomid ) );
			if ( upper_right_coords.x() < curxyz.x() ) upper_right_coords.x( curxyz.x() );
			if ( upper_right_coords.y() < curxyz.y() ) upper_right_coords.y( curxyz.y() );
			if ( upper_right_coords.z() < curxyz.z() ) upper_right_coords.z( curxyz.z() );
			if ( lower_left_coords.x() > curxyz.x() ) lower_left_coords.x( curxyz.x() );
			if ( lower_left_coords.y() > curxyz.y() ) lower_left_coords.y( curxyz.y() );
			if ( lower_left_coords.z() > curxyz.z() ) lower_left_coords.z( curxyz.z() );
		}
	}
	lower_left_coords.x( lower_left_coords.x() - padding );
	lower_left_coords.y( lower_left_coords.y() - padding );
	lower_left_coords.z( lower_left_coords.z() - padding );
	upper_right_coords.x( upper_right_coords.x() + padding );
	upper_right_coords.y( upper_right_coords.y() + padding );
	upper_right_coords.z( upper_right_coords.z() + padding );

	TR << "Done.  Bounding box is from [" << lower_left_coords.x() << ", " << lower_left_coords.y() << ", " << lower_left_coords.z() << "] to [" << upper_right_coords.x() << ", " << upper_right_coords.y() << ", " << upper_right_coords.z() << "]." << std::endl;
}

/// @brief Given lower-left, upper-right, and voxel size data, resize the voxel grid appropriately.
void
VoidsPenaltyVoxelGrid::resize_voxel_grid(
	numeric::xyzVector< core::Real> const &lower_left_coords,
	numeric::xyzVector< core::Real> &upper_right_coords,
	core::Real const &voxel_size,
	utility::fixedsizearray1< core::Size, 3 > &voxel_data_dimensions,
	numeric::MathNTensor< bool, 3 > &voxel_data
) const {
	TR << "Setting up voxel grid..." << std::endl;

	voxel_data_dimensions[1] = static_cast< core::Size >( std::ceil( 1+(upper_right_coords.x() - lower_left_coords.x()) / voxel_size ) );
	voxel_data_dimensions[2] = static_cast< core::Size >( std::ceil( 1+(upper_right_coords.y() - lower_left_coords.y()) / voxel_size ) );
	voxel_data_dimensions[3] = static_cast< core::Size >( std::ceil( 1+(upper_right_coords.z() - lower_left_coords.z()) / voxel_size ) );

	//Update positions of upper-right coordinates.
	upper_right_coords.x( static_cast< core::Real >(voxel_data_dimensions[1]-1)*voxel_size_ + lower_left_coords.x() );
	upper_right_coords.y( static_cast< core::Real >(voxel_data_dimensions[2]-1)*voxel_size_ + lower_left_coords.y() );
	upper_right_coords.z( static_cast< core::Real >(voxel_data_dimensions[3]-1)*voxel_size_ + lower_left_coords.z() );

	voxel_data = numeric::MathNTensor< bool, 3 >( voxel_data_dimensions, false ); //Weirdly, in debug mode, this seems to be the very slow step.

	if ( TR.visible() ) {
		TR << "Updated upper-right coordinates.  Bounding box is now [" << lower_left_coords.x() << ", " << lower_left_coords.y() << ", " << lower_left_coords.z() << "] to [" << upper_right_coords.x() << ", " << upper_right_coords.y() << ", " << upper_right_coords.z() << "]." << std::endl;
		TR << "Set up a " << voxel_data_dimensions[1] << "x" << voxel_data_dimensions[2] << "x" << voxel_data_dimensions[3] << " voxel grid.  (" << voxel_data_dimensions[1]*voxel_data_dimensions[2]*voxel_data_dimensions[3] << " voxels total).";
	}
}

/// @brief Given a pose and a voxel grid that is already set up for the pose, compute whether each point is buried or not.
void
VoidsPenaltyVoxelGrid::compute_burial(
	core::pose::Pose const &pose,
	utility::fixedsizearray1< core::Size, 3 > const &voxel_data_dimensions,
	numeric::MathNTensor< bool, 3 > &voxel_data,
	core::Size & buried_volume
) const {
	TR << "Computing voxel burial by method of sidechain neighbours..." << std::endl;
	runtime_assert_string_msg( pose.total_residue() > 0, "Error in core::pack::voids_penalty_energy::VoidsPenaltyVoxelGrid::compute_burial():  The pose has no residues!"  );

	//Iterate over all residues and set up conevects ad conebases:
	utility::vector1< numeric::xyzVector< core::Real > > conevects( pose.total_residue() );
	utility::vector1< numeric::xyzVector< core::Real > > conebases( pose.total_residue() );
	utility::vector1< bool > skip_list( pose.total_residue(), false );
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		if ( !pose.residue_type(i).is_polymer() ) {
			skip_list[i] = true;
			continue; //Skip ligands
		}
		core::id::AtomID at2( pose.residue_type(i).aa() == core::chemical::aa_gly ? pose.residue_type(i).atom_index("2HA") : pose.residue_type(i).first_sidechain_atom(), i );
		core::id::AtomID at1( pose.residue_type(i).icoor(at2.atomno()).stub_atom1().atomno(), i );
		conevects[i] = numeric::xyzVector< core::Real >( (pose.xyz(at2) - pose.xyz(at1)).normalize() );
		conebases[i] = pose.xyz(at2);
	}

	//Iterate over all voxels:
	buried_volume = 0;
	for ( core::Size i(0); i<voxel_data_dimensions[1]; ++i ) {
		for ( core::Size j(0); j<voxel_data_dimensions[2]; ++j ) {
			for ( core::Size k(0); k<voxel_data_dimensions[3]; ++k ) {
				if ( is_buried( get_voxel_coordinates(i, j, k), conevects, conebases, skip_list ) ) {
					voxel_data(i,j,k) = true;
					++buried_volume;
				}
			}
		}
	}
}

/// @brief Given a pose and a point in space, determine whether the point is buried by the side-chain
/// cones method.
bool
VoidsPenaltyVoxelGrid::is_buried(
	numeric::xyzVector<core::Real> const &coords,
	utility::vector1< numeric::xyzVector<core::Real> > const &conevects,
	utility::vector1< numeric::xyzVector<core::Real> > const &conebases,
	utility::vector1< bool > const &skip_list
) const {
	debug_assert(conevects.size() == conebases.size()); //Should be true
	debug_assert(conevects.size() == skip_list.size()); //Should also be true.
	core::Size counter(0);
	for ( core::Size i(1), imax(conevects.size()); i<=imax; ++i ) {
		if ( skip_list[i] ) continue; //This residue is a ligand or is otherwise one to skip.
		if ( (coords - conebases[i]).length_squared() > cone_distance_cutoff_sq_ ) continue; //Skip most of the cone checks.
		if ( (coords - conebases[i]).normalize().dot_product( conevects[i] ) < cone_dotproduct_cutoff_ ) continue; //Skip cones pointing away.
		++counter;
		if ( counter >= containing_cones_cutoff_ ) return true;
	}

	return false; //TODO
}

/// @brief Given a residue index, get (a) the same residue index back in the asymmetric case, (b) the same residue index back in the
/// symmetric case if this residue is independent, and (c) the residue index on which this residue depends in the symmetric case if
/// this residue is dependent.
core::Size
VoidsPenaltyVoxelGrid::get_current_res_symmetric(
	core::Size const index_in,
	core::conformation::symmetry::SymmetryInfoCOP symminfo
) const {
	if ( symminfo == nullptr || symminfo->bb_is_independent(index_in) ) return index_in;
	return symminfo->bb_follows(index_in);
}


} //core
} //guidance_scoreterms
} //pack
} //voids_penalty_energy

