// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @details
///
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_rna_RNA_BaseCentroidChecker_HH
#define INCLUDED_protocols_stepwise_rna_RNA_BaseCentroidChecker_HH

#include <core/types.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/rna/RNA_CentroidInfo.fwd.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <ObjexxFCL/FArray1D.hh> /* for some reason, can't get away with just .fwd.hh */
#include <ObjexxFCL/FArray2D.hh> /* for some reason, can't get away with just .fwd.hh */


#include <string>
#include <map>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {

class RNA_BaseCentroidChecker: public utility::pointer::ReferenceCount {
public:

	// Constructor
	RNA_BaseCentroidChecker( core::pose::Pose const & pose, working_parameters::StepWiseWorkingParametersCOP & working_parameters,
		bool const tether_jump = false );

	virtual ~RNA_BaseCentroidChecker();

	bool
	update_base_stub_list_and_check_centroid_interaction( core::pose::Pose const & pose, StepWiseRNA_CountStruct & count_data );

	bool
	update_base_stub_list_and_check_that_terminal_res_are_unstacked( core::pose::Pose const & pose, bool const reinitialize = false );

	bool
	check_that_terminal_res_are_unstacked( bool const verbose = false );

	bool
	check_centroid_interaction( core::kinematics::Stub const &  moving_res_base_stub, StepWiseRNA_CountStruct & count_data );

	void
	set_allow_base_pair_only_screen( bool const setting ){ allow_base_pair_only_screen_ = setting; }

	bool const &
	allow_base_pair_only_screen() const{ return allow_base_pair_only_screen_; }

	void
	set_floating_base( bool const setting ){ floating_base_ = setting; }

	bool
	found_centroid_interaction() const { return found_centroid_interaction_; }

private:

	void
	Initialize_is_virtual_base( core::pose::Pose const & pose, bool const verbose = false );

	void
	Initialize_base_stub_list( core::pose::Pose const & pose, bool const verbose = false );

	void
	Initialize_terminal_res( core::pose::Pose const & pose );

	bool
	check_base_stack( core::kinematics::Stub const & moving_residue_base_stub,
		core::kinematics::Stub const & other_base_stub,
		core::Real const base_axis_CUTOFF,
		core::Real const base_planarity_CUTOFF,
		bool const verbose = false ) const;

	bool
	check_base_stack( Size const & pos1, Size const & pos2, bool const verbose  = false  );

	bool
	check_base_stack( core::kinematics::Stub const & moving_res_base,
		core::Real const base_axis_CUTOFF,
		core::Real const base_planarity_CUTOFF ) const;

	bool
	check_base_pair( core::kinematics::Stub const & moving_residue_base_stub,
		core::kinematics::Stub const & other_base_stub,
		core::Real const base_axis_CUTOFF,
		core::Real const base_planarity_CUTOFF ) const;

	bool
	check_base_pair( core::kinematics::Stub const & moving_residue_base_stub,
		core::Real const base_axis_CUTOFF,
		core::Real const base_planarity_CUTOFF ) const;

	bool
	check_base_stack( core::kinematics::Stub const & moving_residue_base_stub,
		core::kinematics::Stub const & other_base_stub,
		bool const verbose = false ) const;

	bool
	check_centroid_interaction_floating_base( core::kinematics::Stub const &  moving_res_base_stub,
		StepWiseRNA_CountStruct & count_data ) const;

	bool
	check_centroid_interaction( StepWiseRNA_CountStruct & count_data );

	void
	update_base_stub_list( core::pose::Pose const & pose );

	bool
	is_strong_base_stack( core::kinematics::Stub const & moving_res_base ) const;

	bool
	is_medium_base_stack_and_medium_base_pair( core::kinematics::Stub const & moving_res_base ) const;

private:

	working_parameters::StepWiseWorkingParametersCOP working_parameters_;
	core::scoring::rna::RNA_CentroidInfoOP rna_centroid_info_;

	core::Real const base_stack_dist_cutoff_;
	core::Real const base_stack_z_offset_max_;
	core::Real const base_stack_z_offset_min_;
	core::Real const base_stack_axis_cutoff_;
	core::Real const base_stack_planarity_cutoff_;
	core::Real const base_pair_dist_min_;
	core::Real const base_pair_dist_max_;
	core::Real const base_pair_z_offset_cutoff_;
	core::Real const base_pair_axis_cutoff_;
	core::Real const base_pair_planarity_cutoff_;
	core::Real const base_pair_rho_min_;
	core::Real const base_pair_rho_max_;
	bool allow_base_pair_only_screen_;
	bool floating_base_;
	bool found_centroid_interaction_;
	bool tether_jump_;

	utility::vector1 < core::Size > moving_residues_;
	utility::vector1 < core::Size > fixed_residues_;
	utility::vector1 < core::kinematics::Stub > base_stub_list_;

	utility::vector1< core::Size > terminal_res_;
	ObjexxFCL::FArray1D < bool > is_terminal_res_;
	ObjexxFCL::FArray1D < bool > is_fixed_res_;
	ObjexxFCL::FArray1D < bool > is_moving_res_;
	ObjexxFCL::FArray1D < bool > is_virtual_base_; //Parin Mar 6
	ObjexxFCL::FArray2D < bool > stacked_on_terminal_res_in_original_pose_;

};

} //checker
} //rna
} //modeler
} //stepwise
} //protocols

#endif
