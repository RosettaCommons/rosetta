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

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/chemical/AA.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <set>
#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh> //June 02, 2011
#include <core/pose/Pose.hh> //June 02, 2011
#include <core/scoring/rna/RNA_Util.hh> //June 02, 2011

//#define centroid_bin_size 0.5
//#define euler_angle_bin_size 5
#define centroid_bin_size 1.0
#define euler_angle_bin_size 20
#define euler_z_bin_size 0.05


typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace swa {
namespace rna {

struct Anchor_ribose_stub{
	numeric::xyzVector<core::Real> origin;
  Matrix coordinate_matrix;
	Matrix invert_coordinate_matrix;
};

//Should make sure that alpha and gamma lies in the [-Pi:Pi] range.
struct Euler_angles{
	core::Real z; //z=cos(beta)
	core::Real alpha; //phi
	core::Real beta; //theta
	core::Real gamma; //psi
};


struct Base_bin{

	int centroid_x;
	int centroid_y;
	int centroid_z;
	int euler_alpha;
	int euler_z;
	int euler_gamma;

};


struct
compare_base_bin{

	//The expression comp(a,b), where comp is an object of this comparison class and a and b are key values, shall return true if a is to be placed at an earlier position than b in a strict weak ordering operation

  bool
	operator() (Base_bin const & first , Base_bin const & second) const {

		if(first.centroid_x != second.centroid_x) return (first.centroid_x<second.centroid_x); //x
		if(first.centroid_y != second.centroid_y) return (first.centroid_y<second.centroid_y); //y
		if(first.centroid_z != second.centroid_z) return (first.centroid_z<second.centroid_z) ; //z
		if(first.euler_alpha != second.euler_alpha) return (first.euler_alpha < second.euler_alpha);
		if(first.euler_gamma != second.euler_gamma) return (first.euler_gamma < second.euler_gamma);
		if(first.euler_z != second.euler_z) return (first.euler_z < second.euler_z);

		return false; //Equality case.
	}

};

struct
compare_int_pair{

	//The expression comp(a,b), where comp is an object of this comparison class and a and b are key values, shall return true if a is to be placed at an earlier position than b in a strict weak ordering operation


  bool
	operator() (std::pair<int, int> const & pair_one , std::pair<int, int> const & pair_two) const {

		if(pair_one.first != pair_two.first) return (pair_one.first<pair_two.first);
		if(pair_one.second != pair_two.second) return (pair_one.second<pair_two.second);

		return false; //Equality case.
	}

};

struct
compare_test{
  bool
	operator() (core::Real const & first, core::Real const & second) const {
		return first< second;
	}
};



//typedef utility::pointer::owning_ptr< FloatingBaseChainClosureJobParameter > FloatingBaseChainClosureJobParameterOP;
//typedef utility::pointer::owning_ptr< FloatingBaseChainClosureJobParameter const > FloatingBaseChainClosureJobParameterCOP;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This is the version currently used for floating base sampling....slightly different from the old version of Base_centroid_screening which appear to StepWiseRNA_Sampling
//Should probably integrate this with Rhiju's class Jan 28, 2010. ***ALERT***RHIJU pointed out that the screening condition is slight different in his new class.

bool
Is_base_stack(core::kinematics::Stub const & moving_res_base,
						  utility::vector1 < core::kinematics::Stub > const & other_residues_base_list,
				  	  core::Real const base_axis_CUTOFF,
	            core::Real const base_planarity_CUTOFF);

bool
Is_base_pair(core::kinematics::Stub const & moving_res_base,
						 utility::vector1 < core::kinematics::Stub > const & other_residues_base_list,
				  	 core::Real const base_axis_CUTOFF,
	           core::Real const base_planarity_CUTOFF);

bool
Is_strong_base_stack(core::kinematics::Stub const & moving_res_base, utility::vector1 < core::kinematics::Stub > const & other_residues_base_list);

// Undefined, commenting out to fix PyRosetta build  bool Is_medium_stack_base_and_medium_stack_base(core::kinematics::Stub const & moving_res_base, utility::vector1 < core::kinematics::Stub > const & other_residues_base_list);

bool
Base_centroid_screening(core::kinematics::Stub const & moving_res_base, utility::vector1 < core::kinematics::Stub > const & other_residues_base_list, core::Size const num_nucleotides, SillyCountStruct & count_data, bool const allow_base_pair_only_screen);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Base_bin
Get_euler_stub_bin(numeric::xyzVector<core::Real> const & centroid, Euler_angles const &  euler_angles);

core::kinematics::Stub
Get_ribose_stub(core::conformation::Residue const & rsd, bool const Is_prepend, bool const verbose=true);


int
DOF_bin_value(std::map<Base_bin , int , compare_base_bin>::const_iterator const & base_bin_it, std::string const & DOF);

core::Real
DOF_bin_size(std::string const & DOF);

void
Analyze_base_bin_map(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, std::string const foldername);

void
Analyze_base_bin_map(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, std::string const & DOF_one, std::string const & DOF_two, std::string const foldername);


// Undefined, commenting out to fix PyRosetta build  void Analyze_base_bin_map_old(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, bool const Is_dinucleotide);

void
translate_then_rotate_pose(core::pose::Pose & pose, numeric::xyzVector<core::Real> const & vector, numeric::xyzMatrix< core::Real > const matrix, bool const verbose=false);

void
set_to_origin(core::pose::Pose & pose, core::Size const seq_num, bool verbose=false);

Euler_angles
Get_euler_angles( numeric::xyzMatrix< core::Real > const & coordinate_matrix);

void
convert_euler_to_coordinate_matrix(Euler_angles const & E, numeric::xyzMatrix< core::Real > & coordinate_matrix);

void
get_specific_atom_coordinate(std::string const & atom_name,
														 numeric::xyzVector<core::Real> & atom_pos,
										         core::conformation::Residue const & rsd_at_origin,
										         core::kinematics::Stub const & moving_res_base_stub);

core::Real
get_max_centroid_to_atom_distance(utility::vector1 <core::conformation::ResidueOP> const & rsd_at_origin_list, std::string const atom_name);

//////////////////////////////////////////Ribose sugar , close_break closures function////////////////////////////

utility::vector1 <core::conformation::ResidueOP>
setup_residue_at_origin_list(core::pose::Pose const & pose, core::Size const & moving_res, bool const extra_anti_chi_rotamer, bool const extra_syn_chi_rotamer,std::string const pose_name);

bool
check_floating_base_chain_closable(core::Size const & reference_res,
																 core::pose::Pose const & pose,
																 utility::vector1 <core::conformation::ResidueOP> const & rsd_at_origin_list,
																 core::kinematics::Stub const & moving_res_base_stub,
																 bool const Is_prepend,
																 core::Size const gap_size);

bool
check_floating_base_chain_closable(core::Size const & reference_res,
																 utility::vector1< pose_data_struct2 >,
																 utility::vector1 <core::conformation::ResidueOP> const & rsd_at_origin_list,
																 core::kinematics::Stub const & moving_res_base_stub,
																 bool const Is_prepend,
																 core::Size const gap_size);


void
set_base_coordinate_frame(core::pose::Pose & pose,
													core::Size const & seq_num,
													core::conformation::Residue const & rsd_at_origin,
													core::kinematics::Stub const & moving_res_base_stub);





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
}
}

#endif


