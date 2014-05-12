// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.hh
/// @brief
/// @detailed
///
///  @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSamplerUtil_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSamplerUtil_HH

#include <protocols/stepwise/sampling/rna/rigid_body/FloatingBaseClasses.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/stepwise/sampling/rna/util.hh>
#include <protocols/rotamer_sampler/rigid_body/EulerAngles.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh> //June 02, 2011
#include <core/pose/Pose.hh> //June 02, 2011

#define STANDARD_CENTROID_BIN_SIZE 1.0
#define STANDARD_EULER_ANGLE_BIN_SIZE 20
#define STANDARD_EULER_Z_BIN_SIZE 0.05

typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace rigid_body {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BaseBin
get_euler_stub_bin( numeric::xyzVector< core::Real > const & centroid, rotamer_sampler::rigid_body::EulerAngles const &  euler_angles );

core::kinematics::Stub
get_sugar_stub( core::conformation::Residue const & rsd, bool const is_prepend, bool const verbose = true );

int
DOF_bin_value( std::map< BaseBin, int, compare_base_bin > ::const_iterator const & base_bin_it, std::string const & DOF );

core::Real
DOF_bin_size( std::string const & DOF );

void
analyze_base_bin_map( std::map< BaseBin, int, compare_base_bin > const & base_bin_map, std::string const foldername );

void
analyze_base_bin_map( std::map< BaseBin, int, compare_base_bin > const & base_bin_map, std::string const & DOF_one, std::string const & DOF_two, std::string const foldername );

void
translate_then_rotate_pose( core::pose::Pose & pose, numeric::xyzVector< core::Real > const & vector, numeric::xyzMatrix< core::Real > const matrix, bool const verbose = false );

core::Real
get_max_centroid_to_atom_distance( utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, std::string const atom_name );

void
initialize_xyz_parameters( Distance & max_distance,
													 Distance & max_distance_squared,
													 int & centroid_bin_min,
													 int & centroid_bin_max,
													 utility::vector1< core::conformation::ResidueOP > const & moving_rsd_at_origin_list,
													 Size const gap_size_to_anchor );

utility::vector1<core::conformation::ResidueOP>
setup_residue_at_origin_list(
	core::pose::Pose const & pose,
	core::Size const & moving_res,
	bool const extra_chi,
	bool const use_phenix_geo
);

utility::vector1 < core::pose::PoseOP >
setup_pose_at_origin_list(
	pose::Pose const & pose,
	Size const & moving_res,
	bool const extra_chi,
	bool const use_phenix_geo
);

utility::vector1 < core::pose::PoseOP >
setup_pose_with_moving_residue_alternative_list(
	pose::Pose const & pose,
	Size const & moving_res,
	bool const extra_chi,
	bool const use_phenix_geo		);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
} //rigid_body
} //rna
} //sampling
} //stepwise
} //protocols

#endif


