// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_PoseSetup.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_PoseSetup_hh
#define INCLUDED_protocols_swa_rna_StepWiseRNA_PoseSetup_hh

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
//#include <protocols/rna/RNA_StructureParameters.hh> //how about using a fwd file?
#include <core/types.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <string>
#include <map>

namespace protocols {
namespace swa {
namespace rna {

//	typedef std::map< std::string, core::pose::PoseOP > PoseList;

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseRNA_PoseSetup: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseRNA_PoseSetup( Size const & moving_res,
		std::string const & desired_sequence,
		utility::vector1< std::string > const & input_tags,
		utility::vector1< std::string > const & silent_files_in,
		utility::vector1< core::Size > const & input_res,
		utility::vector1< core::Size > const & input_res2,
		utility::vector1< core::Size > const & cutpoint_open,
		Size const & cutpoint_closed
	);

	//destructor -- necessary?
	~StepWiseRNA_PoseSetup();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

	virtual std::string get_name() const;

	StepWiseRNA_JobParametersOP & job_parameters();

	void
	set_fixed_res(utility::vector1 < core::Size > const & fixed_res );

	void
	set_terminal_res(utility::vector1 < core::Size > const & terminal_res );

	void
	set_bulge_res(utility::vector1 < core::Size > const & bulge_res );

private:

	void
	Import_pose( Size const & i, core::pose::Pose & import_pose) const;

	void
	check_moving_res_in_chain( Size const & start_chain, Size const & end_chain,
														 Size const & num_chains, Size & which_chain_has_moving_res  );

	void
	figure_out_working_sequence_and_mapping();

	void
	figure_out_jump_partners();

	void
	figure_out_cuts();


	void
	make_pose( core::pose::Pose & pose );

	void
	read_input_pose_and_copy_dofs( core::pose::Pose & pose );


	void
	figure_out_Prepend_Internal( core::pose::Pose const & pose );

	void
	figure_out_partition_definition( core::pose::Pose const & pose );

	void
	figure_out_gap_size_and_five_prime_chain_break_res();

	void
	reroot_fold_tree( core::pose::Pose & pose );

	void
	apply_cutpoint_variants( core::pose::Pose & pose ) const;

	void
	check_close_chain_break( core::pose::Pose const & pose ) const;

	void
	apply_bulge_variants( core::pose::Pose & pose ) const;

	void
	slice_native();

	utility::vector1< Size > const
	apply_full_to_sub_mapping( utility::vector1< Size > & res_vector) const;

	void
	apply_virtual_phosphate_variants( core::pose::Pose & pose ) const;


private:

	Size const moving_res_;
	std::string const desired_sequence_;
	core::chemical::ResidueTypeSetCAP rsd_set_;

	utility::vector1< std::string > const input_tags_;
	utility::vector1< std::string > const silent_files_in_;
	utility::vector1< utility::vector1< Size > > input_res_vectors_;
	utility::vector1< Size > const cutpoint_open_;
	utility::vector1< Size > fixed_res_;
	utility::vector1< Size > terminal_res_;
	utility::vector1< Size > bulge_res_;
	Size const cutpoint_closed_;
	ObjexxFCL::FArray1D< bool > is_cutpoint_;

	utility::vector1< std::pair< core::Size, core::Size > > jump_partners_;
	utility::vector1< core::Size > cuts_;

	StepWiseRNA_JobParametersOP job_parameters_;

	bool const virtualize_5prime_phosphates_;
};

}
} //swa
} // protocols

#endif
