// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/protocols/generalized_kinematic_closure/util.hh
/// @brief  Headers for GeneralizedKIC utilities.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_generalized_kinematic_closure_util_hh
#define INCLUDED_protocols_generalized_kinematic_closure_util_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.fwd.hh>
#include <protocols/generalized_kinematic_closure/perturber/GeneralizedKICperturber.hh>
#include <protocols/generalized_kinematic_closure/filter/GeneralizedKICfilter.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Project Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace generalized_kinematic_closure {


	/// @brief Function to determine whether a value is in a list.
	bool is_in_list (core::Size const val, utility::vector1 < core::Size > const &list);


	/// @brief function to determine whether a residue index from the original pose is in the residue_map list.
	bool original_pose_residue_is_in_residue_map (core::Size const residue_index, utility::vector1 < std::pair<core::Size, core::Size> > const &residue_map);

	/// @brief Set the loop pose conformation based on a set of results from kinematic closure.
	/// @details
	/// @param[in,out] pose -- A pose consisting of the loop to be closed only.
	/// @param[in] atomlist -- A list of AtomIDs and corresponding xyz coordinates (though we don't use the latter) of the chain of atoms closed by KIC.  Note that the residue indices refer to the loop pose, not the original pose.
	/// @param[in] torsions -- The torsion angles values to set.
	/// @param[in] bondangles -- The bond angle values to set.
	/// @param[in] bondlengths -- The bond length values to set.
	void set_loop_pose (
		core::pose::Pose &pose,
		utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //I want this to be const access
		utility::vector1 < core::Real > const &t_ang,
		utility::vector1 < core::Real > const &b_ang,
		utility::vector1 < core::Real > const &b_len
	);
	
	/// @brief Set the loop pose conformation based on a set of results from kinematic closure.
	/// @details  This version ONLY sets mainchain torsions, and does no rebuilding of mainchain O or H atoms.
	/// @param[in,out] pose -- A pose consisting of the loop to be closed only.
	/// @param[in] atomlist -- A list of AtomIDs and corresponding xyz coordinates (though we don't use the latter) of the chain of atoms closed by KIC.  Note that the residue indices refer to the loop pose, not the original pose.
	/// @param[in] torsions -- The torsion angles values to set.
	void set_loop_pose (
		core::pose::Pose &pose,
		utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //I want this to be const access
		utility::vector1 < core::Real > const &t_ang
	);

	/// @brief Copy the atom positions of the residues in the loop pose to the original pose.
	///
	void copy_loop_pose_to_original (
		core::pose::Pose &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map,
		utility::vector1 < std::pair < core::Size, core::Size > > const &tail_residue_map
	);

	/// @brief Sets phi for an L-alpha amino acid, even if the nitrogen has a nonstandard connection.
	///
	void general_set_phi (
		core::pose::Pose &pose,
		core::Size const residue_index,
		core::Real const &phi_value
	);

	/// @brief Sets psi for an L-alpha amino acid, even if the carbonyl carbon has a nonstandard connection.
	///
	void general_set_psi (
		core::pose::Pose &pose,
		core::Size const residue_index,
		core::Real const &psi_value
	);

	/// @brief Given a residue_map vector of pairs, where each pair is < residue_index_in_perturbedloop_pose, residue_index_in_original_pose >,
	/// and a residue index in the perturbed loop, return the corresponding residue index in the original pose.
	core::Size get_original_pose_rsd ( core::Size const perturbedloop_rsd, utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map );

} //namespace generalized_kinematic_closure
} //namespace protocols

#endif //INCLUDED_protocols_generalized_kinematic_closure_util_hh
