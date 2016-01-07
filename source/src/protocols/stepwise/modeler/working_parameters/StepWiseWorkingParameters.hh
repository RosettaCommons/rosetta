// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_working_parameters_StepWiseWorkingParameters_HH
#define INCLUDED_protocols_stepwise_modeler_working_parameters_StepWiseWorkingParameters_HH

#include <protocols/stepwise/modeler/working_parameters/StepWiseBasicWorkingParameters.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.fwd.hh>

#define  BRIDGE_RES 123
#define  MOVING_RES 999

///////////////////////////////////////////////////////////////////////////////////////
//
// WorkingParameters was the original object for storing parameters for StepWise Assembly
//  Jobs. Now re-created 'on the fly' for each move in StepWise MonteCarlo.
///
// It (or a large fraction of it) may be largely deprecated with the
//  creation of the core::pose:FullModelInfo object.
//
// Currently in the midst of unifying RNA and Protein WorkingParameters objects, and
//  potentially deprecating some or all of it.
//
///////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace modeler {
namespace working_parameters {

class StepWiseWorkingParameters: public StepWiseBasicWorkingParameters {

public:

	//constructor
	StepWiseWorkingParameters();

	//destructor
	~StepWiseWorkingParameters();

public:

	/////////////////////////
	// RNA stuff
	/////////////////////////
	bool const & output_extra_RMSDs() const { return output_extra_RMSDs_; }
	void set_output_extra_RMSDs( bool const setting ){ output_extra_RMSDs_ = setting; }

	bool const & is_simple_full_length_job_params() const { return is_simple_full_length_job_params_; }
	void set_is_simple_full_length_job_params( bool const setting ){ is_simple_full_length_job_params_ = setting; }

	Size gap_size_to_anchor() const; // derives from internal parameters.

	virtual Size actually_moving_res() const; // need to unify with protein.
	Size working_reference_res() const;

	Size const & five_prime_chain_break_res() const { return five_prime_chain_break_res_; }
	void set_five_prime_chain_break_res( Size const & setting ) { five_prime_chain_break_res_ = setting; }

	bool const & add_virt_res_as_root() const { return add_virt_res_as_root_; }
	void set_add_virt_res_as_root( bool const setting ){ add_virt_res_as_root_ = setting; }

	bool const & floating_base() const{ return floating_base_;}
	void set_floating_base( bool const setting ){ floating_base_ = setting; }

	Size const & floating_base_anchor_res() const{ return floating_base_anchor_res_;}
	void set_floating_base_anchor_res( Size const setting ){ floating_base_anchor_res_ = setting; }
	Size working_floating_base_anchor_res() const;

	bool const & rebuild_bulge_mode() const{ return rebuild_bulge_mode_;}
	void set_rebuild_bulge_mode( bool const setting ){ rebuild_bulge_mode_ = setting; }

	bool const & sample_both_sugar_base_rotamer() const{ return sample_both_sugar_base_rotamer_;}
	void set_sample_both_sugar_base_rotamer( bool const setting ){ sample_both_sugar_base_rotamer_ = setting; }

	void set_is_prepend_map( std::map< core::Size, bool > const & setting ) { is_prepend_map_ = setting; }
	std::map< core::Size, bool > const & is_prepend_map() const { return is_prepend_map_; }

	utility::vector1< core::Size > const & native_alignment() const { return native_alignment_; }
	void set_native_alignment( utility::vector1< core::Size > const & setting ) { native_alignment_ = setting; }

	utility::vector1< core::Size > const & working_native_alignment() const { return working_native_alignment_; }
	void set_working_native_alignment( utility::vector1< core::Size > const & setting ) { working_native_alignment_ = setting; }

	utility::vector1< core::Size > const & working_best_alignment() const { return working_best_alignment_; }
	void set_working_best_alignment( utility::vector1< core::Size > const & setting ) { working_best_alignment_ = setting; }

	utility::vector1< core::Size > const & global_sample_res_list() const { return global_sample_res_list_; }
	void set_global_sample_res_list( utility::vector1< core::Size > const & setting ); // updates working list, too
	utility::vector1< core::Size > const & working_global_sample_res_list() const { return working_global_sample_res_list_; }

	utility::vector1< core::Size > const & force_syn_chi_res_list() const { return force_syn_chi_res_list_; }
	void set_force_syn_chi_res_list( utility::vector1< core::Size > const & setting );  // updates working list, too
	utility::vector1< core::Size > const & working_force_syn_chi_res_list() const { return working_force_syn_chi_res_list_; }

	utility::vector1< core::Size > const & force_anti_chi_res_list() const { return force_anti_chi_res_list_; }
	void set_force_anti_chi_res_list( utility::vector1< core::Size > const & setting );  // updates working list, too
	utility::vector1< core::Size > const & working_force_anti_chi_res_list() const { return working_force_anti_chi_res_list_; }

	utility::vector1< core::Size > const &  terminal_res() const { return terminal_res_; }
	void set_terminal_res( utility::vector1< core::Size > const & terminal_res );
	utility::vector1< core::Size > const &  working_terminal_res() const { return working_terminal_res_; }

	utility::vector1< core::Size > const &  block_stack_above_res() const { return block_stack_above_res_; }
	void set_block_stack_above_res( utility::vector1< core::Size > const & block_stack_above_res );
	utility::vector1< core::Size > const &  working_block_stack_above_res() const { return working_block_stack_above_res_; }

	utility::vector1< core::Size > const &  block_stack_below_res() const { return block_stack_below_res_; }
	void set_block_stack_below_res( utility::vector1< core::Size > const & block_stack_below_res );
	utility::vector1< core::Size > const &  working_block_stack_below_res() const { return working_block_stack_below_res_; }

	utility::vector1< core::Size > const & force_north_sugar_list() const { return force_north_sugar_list_; }
	void set_force_north_sugar_list( utility::vector1< core::Size > const & setting );  // updates working list, too
	utility::vector1< core::Size > const & working_force_north_sugar_list() const { return working_force_north_sugar_list_; }

	utility::vector1< core::Size > const & force_south_sugar_list() const { return force_south_sugar_list_; }
	void set_force_south_sugar_list( utility::vector1< core::Size > const & setting );  // updates working list, too
	utility::vector1< core::Size > const & working_force_south_sugar_list() const { return working_force_south_sugar_list_; }

	utility::vector1< core::Size > const & protonated_H1_adenosine_list() const { return protonated_H1_adenosine_list_; }
	void set_protonated_H1_adenosine_list( utility::vector1< core::Size > const & setting );  // updates working list, too
	utility::vector1< core::Size > const & working_protonated_H1_adenosine_list() const { return working_protonated_H1_adenosine_list_; }

	core::kinematics::FoldTree const & fold_tree() const;
	void set_fold_tree( core::kinematics::FoldTree const & setting );

	utility::vector1< utility::vector1< core::Size > > const & input_res_vectors() const { return input_res_vectors_; }
	void set_input_res_vectors( utility::vector1< utility::vector1< Size > > const & setting ) { input_res_vectors_ = setting; }

	/////////////////////////
	// protein stuff
	/////////////////////////
	utility::vector1< core::Size > const & working_bridge_res() const { return working_bridge_res_; }
	void set_working_bridge_res( utility::vector1< Size > const & setting ){ working_bridge_res_ = setting; }

	utility::vector1< core::Size > const &  working_superimpose_res() const { return working_superimpose_res_; }
	void set_working_superimpose_res( utility::vector1< core::Size > const & working_superimpose_res ) { working_superimpose_res_ = working_superimpose_res; }

	ObjexxFCL::FArray1D< core::Size > const & is_moving_res() const { return is_moving_res_; }
	void set_is_moving_res( ObjexxFCL::FArray1D< core::Size > const & setting ) { is_moving_res_ = setting; }

	utility::vector1< bool > const is_pre_proline() const;

private:

	/////////////////////////
	// RNA stuff
	/////////////////////////
	bool output_extra_RMSDs_; //Used in StepWiseRNA_output_Data.cc
	bool is_simple_full_length_job_params_;

	Size five_prime_chain_break_res_;

	std::map< core::Size, bool > is_prepend_map_;

	bool add_virt_res_as_root_;
	bool floating_base_;
	Size floating_base_anchor_res_;
	bool rebuild_bulge_mode_;
	bool sample_both_sugar_base_rotamer_;

	utility::vector1< core::Size > native_alignment_;
	utility::vector1< core::Size > working_native_alignment_;
	utility::vector1< core::Size > working_best_alignment_;

	utility::vector1< core::Size > global_sample_res_list_;
	utility::vector1< core::Size > working_global_sample_res_list_;

	utility::vector1< core::Size > terminal_res_;
	utility::vector1< core::Size > working_terminal_res_;

	utility::vector1< core::Size > block_stack_above_res_;
	utility::vector1< core::Size > working_block_stack_above_res_;

	utility::vector1< core::Size > block_stack_below_res_;
	utility::vector1< core::Size > working_block_stack_below_res_;

	utility::vector1< core::Size > force_syn_chi_res_list_;
	utility::vector1< core::Size > working_force_syn_chi_res_list_;

	utility::vector1< core::Size > force_anti_chi_res_list_;
	utility::vector1< core::Size > working_force_anti_chi_res_list_;

	utility::vector1< core::Size > force_north_sugar_list_;
	utility::vector1< core::Size > working_force_north_sugar_list_;

	utility::vector1< core::Size > force_south_sugar_list_;
	utility::vector1< core::Size > working_force_south_sugar_list_;

	utility::vector1< core::Size > protonated_H1_adenosine_list_;
	utility::vector1< core::Size > working_protonated_H1_adenosine_list_;

	core::kinematics::FoldTree fold_tree_;
	utility::vector1< utility::vector1< core::Size > > input_res_vectors_;

	/////////////////////////
	// protein stuff
	/////////////////////////
	ObjexxFCL::FArray1D< core::Size > is_moving_res_;
	utility::vector1< core::Size > working_bridge_res_;
	utility::vector1< core::Size > working_superimpose_res_;

};

} //working_parameters
} //modeler
} //stepwise
} //protocols

#endif
