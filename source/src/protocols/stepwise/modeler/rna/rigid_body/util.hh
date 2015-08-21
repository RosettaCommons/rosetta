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
/// @details
///
///  @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSamplerUtil_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSamplerUtil_HH

#include <protocols/stepwise/modeler/rna/rigid_body/FloatingBaseClasses.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/stepwise/sampler/rigid_body/EulerAngles.fwd.hh>

#define STANDARD_CENTROID_BIN_SIZE 1.0
#define STANDARD_EULER_ANGLE_BIN_SIZE 20
#define STANDARD_EULER_Z_BIN_SIZE 0.05

typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace rigid_body {

// following is rather RNA specific -- hope to deprecate/generalize soon
core::Real
get_max_centroid_to_atom_distance( utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, std::string const & atom_name );

void
initialize_xyz_parameters( core::Distance & max_distance,
	core::Distance & max_distance_squared,
	int & centroid_bin_min,
	int & centroid_bin_max,
	utility::vector1< core::conformation::ResidueOP > const & moving_rsd_at_origin_list,
	core::Size const gap_size_to_anchor );

utility::vector1 < core::pose::PoseOP >
setup_pose_with_moving_residue_alternative_list(
	core::pose::Pose const & pose,
	core::Size const & moving_res,
	bool const extra_chi,
	bool const use_phenix_geo  );

int
DOF_bin_value( std::map< BaseBin, int, compare_base_bin > ::const_iterator const & base_bin_it, std::string const & DOF );

core::Real
DOF_bin_size( std::string const & DOF );

// Base bin map was parin's way to assess rigid body sampling, but no longer really in use. Deprecate?
void
analyze_base_bin_map( std::map< BaseBin, int, compare_base_bin > const & base_bin_map, std::string const & foldername );

void
analyze_base_bin_map( std::map< BaseBin, int, compare_base_bin > const & base_bin_map, std::string const & DOF_one, std::string const & DOF_two, std::string const & foldername );


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
} //rigid_body
} //rna
} //modeler
} //stepwise
} //protocols

#endif


