// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/pose_manipulation/pose_manipulation.hh
/// @brief some general functions to manipulate poses. feel free to add your own
/// @brief if you add your own, please mention your name:)
/// @author Florian Richter, floric@u.washington.edu
/// @author Steven Lewis, smlewi@gmail.com (insert_pose_into_pose) domain insertion code


#ifndef INCLUDED_protocols_toolbox_pose_manipulation_hh
#define INCLUDED_protocols_toolbox_pose_manipulation_hh

// Unit headers
// Package headers

// Project headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>
//#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>


//#include <core/id/AtomID.fwd.hh>

// Utility Headers
//#include <utility/pointer/ReferenceCount.hh>
#ifdef WIN32
#include <string>
#endif

//Utility Headers

// C++ Headers

namespace protocols {
namespace toolbox {
namespace pose_manipulation {

/// @author Florian Richter( floric@u.washington.edu) , june 08
/// @brief puts in ala residues at the positions specified in the 'positions' input array
void
construct_poly_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
);

/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @brief puts in D-ala residues at the positions specified in the 'positions' input array
void
construct_poly_d_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
);

/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @brief puts in beta-3-ala residues at the positions specified in the 'positions' input array.
void
construct_poly_beta_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
);

/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @brief puts in D-beta-3-ala residues at the positions specified in the 'positions' input array
void
construct_poly_d_beta_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
);

/// @author Possu Huang ( possu@uw.edu)
/// @brief allows construction of a polymer of any residue type
void
construct_poly_uniq_restype_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	core::chemical::ResidueType const & restype,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
);

/// @author Nobuyasu Koga ( nobuyasu@uw.edu ), Oct 09
/// @brief puts in XXX residues at the positions specified in the 'positions' input array
void
construct_poly_XXX_pose(
	std::string const & aa,
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
);

/// @author Nobuyasu Koga ( nobuyasu@uw.edu ), Oct 09; Tom Linsky (tlinsky@uw.edu), Nov 2014
/// @brief puts in XXX residues at the positions from the given residue typeset specified in the 'positions' input array
void
construct_poly_XXX_pose(
	std::string const & aa,
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	core::chemical::ResidueTypeSetCOP restype_set,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
);

/// @author Florian Richter( floric@u.washington.edu) , aug 08
/// @brief deletes all nonprotein residues from a pose
void
remove_non_protein_residues(
	core::pose::Pose & pose
);


/// @author Florian Richter( floric@u.washington.edu), nov 11
/// @brief adds chainbreak residue types depending on fold tree jumps
void
add_chainbreaks_according_to_jumps( core::pose::Pose & pose );

/// @author Florian Richter( floric@u.washington.edu), nov 11
/// @brief removes chainbreak residue types depending on fold tree jumps
void
remove_chainbreaks_according_to_jumps( core::pose::Pose & pose );

/// @author Florian Richter( floric@u.washington.edu) , sep 08
/// @brief superimposes one pose onto the other at the positions specified and
/// @brief with the offset specified
core::Real
superimpose_pose_on_subset_CA(
	core::pose::Pose & pose,
	core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & positions,
	int const offset = 0
);

/// @author Brian Coventry (bcov@uw.edu)
/// @brief Repack a single residue
/// @detail Repacks a single residue using the command-line flags
///         A list of amino acid name1s may be provided if one wishes to design.
///         Leaving this "" or default implies no design
void
repack_this_residue(
	core::Size seq_pos,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP scorefxn,
	bool include_current = true,
	std::string name1s_if_design = "" );

/// @author Brian Coventry (bcov@uw.edu)
/// @brief Repack a subset of residues
/// @detail Repacks a multiple residues using the command-line flags
///         A list of amino acid name1s may be provided if one wishes to design.
///         Leaving this "" or default implies no design
void
repack_these_residues(
	core::select::residue_selector::ResidueSubset const & subset,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP scorefxn,
	bool include_current = true,
	std::string name1s_if_design = "" );

// 4 flavors of this function because vectors aren't easy to deal with in PyRosetta currently

/// @brief Performs a rigid body rotation then translation of a section
/// of a pose given by a ResidueSubset. By default this routine will
/// find the CA center of mass of the subset and use that as the center
/// of rotation.
void
rigid_body_move(
	numeric::xyzVector<core::Real> const & rotation_unit_vector,
	core::Real angle_deg,
	numeric::xyzVector<core::Real> const & translation_vector,
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	numeric::xyzVector<core::Real> center_of_rotation
	= numeric::xyzVector<core::Real>(std::numeric_limits<double>::quiet_NaN(), 0, 0));

/// @brief Performs a rigid body rotation then translation of a section
/// of a pose given by a ResidueSubset. By default this routine will
/// find the CA center of mass of the subset and use that as the center
/// of rotation.
void
rigid_body_move(
	numeric::xyzVector<core::Real> const & rotation_unit_vector,
	core::Real angle_deg,
	numeric::xyzVector<core::Real> const & translation_vector,
	core::Real translation_scalar,
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	numeric::xyzVector<core::Real> center_of_rotation
	= numeric::xyzVector<core::Real>(std::numeric_limits<double>::quiet_NaN(), 0, 0));

/// @brief Performs a rigid body rotation then translation of a section
/// of a pose given by a ResidueSubset. By default this routine will
/// find the CA center of mass of the subset and use that as the center
/// of rotation.
void
rigid_body_move(
	numeric::xyzMatrix<core::Real> rotation,
	numeric::xyzVector<core::Real> const & translation_vector,
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	numeric::xyzVector<core::Real> center_of_rotation
	= numeric::xyzVector<core::Real>(std::numeric_limits<double>::quiet_NaN(), 0, 0));

/// @brief Performs a rigid body rotation then translation of a section
/// of a pose given by a ResidueSubset. By default this routine will
/// find the CA center of mass of the subset and use that as the center
/// of rotation.
void
rigid_body_move(
	numeric::xyzMatrix<core::Real> rotation,
	numeric::xyzVector<core::Real> const & translation_vector,
	core::Real translation_scalar,
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	numeric::xyzVector<core::Real> center_of_rotation
	= numeric::xyzVector<core::Real>(std::numeric_limits<double>::quiet_NaN(), 0, 0));


} //pose_manipulation
} //toolbox
} //protocols


#endif // INCLUDED_protocols_toolbox_pose_manipulation_HH
