// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_FloatingBaseSamplerUtil.hh
/// @brief
/// @detailed
///
///  @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_FloatingBaseSamplerUtil_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_FloatingBaseSamplerUtil_HH

#include <protocols/swa/rna/FloatingBaseClasses.hh>
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
#include <set>
#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh> //June 02, 2011
#include <core/pose/Pose.hh> //June 02, 2011

#define centroid_bin_size 1.0
#define euler_angle_bin_size 20
#define euler_z_bin_size 0.05

typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace swa {
namespace rna {

//typedef utility::pointer::owning_ptr< SugarModeling > SugarModelingOP;
//typedef utility::pointer::owning_ptr< SugarModeling const > SugarModelingCOP;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This is the version currently used for floating base sampling....slightly different from the old version of Base_centroid_screening which appear to StepWiseRNA_Sampling
//Should probably integrate this with Rhiju's class Jan 28, 2010. ***ALERT***RHIJU pointed out that the screening condition is slight different in his new class.

bool
is_base_stack( core::kinematics::Stub const & moving_res_base,
						  utility::vector1 < core::kinematics::Stub > const & other_residues_base_list,
				  	  core::Real const base_axis_CUTOFF,
	            core::Real const base_planarity_CUTOFF );

bool
is_base_pair( core::kinematics::Stub const & moving_res_base,
						 utility::vector1 < core::kinematics::Stub > const & other_residues_base_list,
				  	 core::Real const base_axis_CUTOFF,
	           core::Real const base_planarity_CUTOFF );

bool
is_strong_base_stack( core::kinematics::Stub const & moving_res_base, utility::vector1 < core::kinematics::Stub > const & other_residues_base_list );

bool
floating_base_centroid_screening( core::kinematics::Stub const & moving_res_base, utility::vector1 < core::kinematics::Stub > const & other_residues_base_list, core::Size const num_nucleotides, StepWiseRNA_CountStruct & count_data, bool const allow_base_pair_only_screen );


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Base_bin
get_euler_stub_bin( numeric::xyzVector< core::Real > const & centroid, Euler_angles const &  euler_angles );

core::kinematics::Stub
get_sugar_stub( core::conformation::Residue const & rsd, bool const is_prepend, bool const verbose = true );


int
DOF_bin_value( std::map< Base_bin, int, compare_base_bin > ::const_iterator const & base_bin_it, std::string const & DOF );

core::Real
DOF_bin_size( std::string const & DOF );

void
analyze_base_bin_map( std::map< Base_bin, int, compare_base_bin > const & base_bin_map, std::string const foldername );

void
analyze_base_bin_map( std::map< Base_bin, int, compare_base_bin > const & base_bin_map, std::string const & DOF_one, std::string const & DOF_two, std::string const foldername );


// Undefined, commenting out to fix PyRosetta build  void analyze_base_bin_map_old(std::map<Base_bin, int, compare_base_bin> const & base_bin_map, bool const is_dinucleotide);

void
translate_then_rotate_pose( core::pose::Pose & pose, numeric::xyzVector< core::Real > const & vector, numeric::xyzMatrix< core::Real > const matrix, bool const verbose = false );

void
set_to_origin( core::pose::Pose & pose, core::Size const seq_num, bool verbose = false );

Euler_angles
get_euler_angles( numeric::xyzMatrix< core::Real > const & coordinate_matrix );

void
convert_euler_to_coordinate_matrix( Euler_angles const & E, numeric::xyzMatrix< core::Real > & coordinate_matrix );

void
get_specific_atom_coordinate( std::string const & atom_name,
														 numeric::xyzVector< core::Real > & atom_pos,
										         core::conformation::Residue const & rsd_at_origin,
										         core::kinematics::Stub const & moving_res_base_stub );

core::Real
get_max_centroid_to_atom_distance( utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, std::string const atom_name );

//////////////////////////////////////////Sugar sugar, close_break closures function////////////////////////////

utility::vector1<core::conformation::ResidueOP>
setup_residue_at_origin_list(
	core::pose::Pose const & pose,
	core::Size const & moving_res,
	bool const extra_chi,
	bool const use_phenix_geo
);

bool
check_floating_base_chain_closable( core::Size const & reference_res,
																 core::pose::Pose const & pose,
																 utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
																 core::kinematics::Stub const & moving_res_base_stub,
																 bool const is_prepend,
																 core::Size const gap_size );

bool
check_floating_base_chain_closable( core::Size const & reference_res,
																		utility::vector1< core::pose::PoseOP >,
																		utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
																		core::kinematics::Stub const & moving_res_base_stub,
																		bool const is_prepend,
																		core::Size const gap_size );


void
set_base_coordinate_frame( core::pose::Pose & pose,
													core::Size const & seq_num,
													core::conformation::Residue const & rsd_at_origin,
													core::kinematics::Stub const & moving_res_base_stub );





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
}
}

#endif


