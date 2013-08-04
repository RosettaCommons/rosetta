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
/// @detailed
///
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_BaseCentroidScreener_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_BaseCentroidScreener_HH

#include <core/types.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/rna/RNA_CentroidInfo.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <ObjexxFCL/FArray1D.hh> /* for some reason, can't get away with just .fwd.hh */
#include <ObjexxFCL/FArray2D.hh> /* for some reason, can't get away with just .fwd.hh */


#include <string>
#include <map>

namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_BaseCentroidScreener: public utility::pointer::ReferenceCount {
	public:

	// Constructor
		StepWiseRNA_BaseCentroidScreener( core::pose::Pose const & pose, StepWiseRNA_JobParametersCOP & job_parameters );

		virtual ~StepWiseRNA_BaseCentroidScreener();

		bool
		Update_base_stub_list_and_Check_centroid_interaction( core::pose::Pose const & pose, SillyCountStruct & count_data );

		// Undefined, commenting out to fix PyRosetta build
		// bool non_adjacent_and_stack_base(core::pose::Pose const & pose,  Size const & pos1, Size const & pos2, bool const verbose = false  );

		bool
		Update_base_stub_list_and_Check_that_terminal_res_are_unstacked( core::pose::Pose const & pose, bool const reinitialize = false );

		bool
		Check_that_terminal_res_are_unstacked( bool const verbose = false );

		// Undefined, commenting out to fix PyRosetta build  utility::vector1< core::Size > const & moving_residues() const;

	private:

		void
		Initialize_is_virtual_base( core::pose::Pose const & pose, bool const verbose = false );

		void
		Initialize_base_stub_list( core::pose::Pose const & pose, bool const verbose = false );

		void
		Initialize_terminal_res( core::pose::Pose const & pose );

		bool
		check_stack_base( core::kinematics::Stub const & rebuild_residue_base_stub, core::kinematics::Stub const & base_stub, bool const verbose = false ) const;

		bool
		check_stack_base( Size const & pos1, Size const & pos2, bool const verbose  = false  );

		bool
		check_base_pairing( core::kinematics::Stub const & rebuild_residue_base_stub, core::kinematics::Stub const & base_stub   ) const;

		bool
		Check_centroid_interaction( SillyCountStruct & count_data ) const;

		void
		Update_base_stub_list( core::pose::Pose const & pose );


	private:

		StepWiseRNA_JobParametersCOP job_parameters_;
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

}
} //swa
} // protocols

#endif

