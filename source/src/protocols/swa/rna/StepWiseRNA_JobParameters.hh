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


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_JobParameters_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_JobParameters_HH

#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>


#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

#include <string>
#include <map>

namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_JobParameters: public utility::pointer::ReferenceCount {
	public:

		//constructor!
    StepWiseRNA_JobParameters();

	  //destructor -- necessary?
    virtual ~StepWiseRNA_JobParameters();

		bool const & output_extra_RMSDs() const;
		bool const & Is_simple_full_length_job_params() const;

		std::string const & full_sequence() const;
		std::string const & working_sequence() const;

		Size const & moving_res() const;
		Size const & working_moving_res() const;
		utility::vector1< core::Size > const & working_moving_res_list() const;
		utility::vector1< core::Size > working_res_list() const;
		Size working_reference_res() const;

		Size const & working_moving_suite() const;
		utility::vector1< core::Size > const & working_moving_suite_list() const;

		Size actually_moving_res() const;

		utility::vector1< core::Size >  const & is_working_res() const;
		std::map< core::Size, core::Size > & full_to_sub();
		std::map< core::Size, core::Size > & sub_to_full();
		std::map< core::Size, core::Size > const & const_full_to_sub() const;
		std::map< core::Size, core::Size > const & const_sub_to_full() const;
		core::kinematics::FoldTree const & fold_tree() const;
		std::map< core::Size, bool > const & Is_prepend_map() const;


		utility::vector1< std::pair < core::Size, core::Size > > const & chain_boundaries() const;
		//Size const & which_chain_has_moving_res() const;
		Size const & gap_size() const;
		Size const & five_prime_chain_break_res() const;

		bool const & Is_prepend() const;
		bool const & Is_internal() const;
		bool const & add_virt_res_as_root() const;

		ObjexxFCL::FArray1D < bool > const & partition_definition() const;

		utility::vector1< core::Size > const &  working_fixed_res() const;
		utility::vector1< core::Size > const &  rmsd_res_list() const;
		utility::vector1< core::Size > const &  working_terminal_res() const;
		utility::vector1< core::Size > const &  working_moving_partition_pos() const;

		utility::vector1< utility::vector1< core::Size > > const & input_res_vectors() const;
		utility::vector1< core::Size > const & cutpoint_closed_list() const;
		utility::vector1< core::Size > const & working_best_alignment() const;


		utility::vector1< core::Size > const & native_alignment() const;
		utility::vector1< core::Size > const & working_native_alignment() const;


		utility::vector1< core::Size > const & global_sample_res_list() const;
		utility::vector1< core::Size > const & working_global_sample_res_list() const;

		utility::vector1< core::Size > const & force_syn_chi_res_list() const;
		utility::vector1< core::Size > const & working_force_syn_chi_res_list() const;

		utility::vector1< core::Size > const & force_north_ribose_list() const;
		utility::vector1< core::Size > const & working_force_north_ribose_list() const;

		utility::vector1< core::Size > const & force_south_ribose_list() const;
		utility::vector1< core::Size > const & working_force_south_ribose_list() const;

		utility::vector1< core::Size > const & protonated_H1_adenosine_list() const;
		utility::vector1< core::Size > const & working_protonated_H1_adenosine_list() const;

		void set_output_extra_RMSDs( bool const & setting );
		void set_Is_simple_full_length_job_params( bool const & setting );

		void set_full_sequence( std::string const & setting );
		//void set_working_sequence( std::string const & setting );
		void set_moving_res( Size const & setting );
		void set_working_moving_res_list( utility::vector1< core::Size > const & setting );

		//void set_working_moving_res( Size const & setting );
		//void set_working_moving_suite_list( utility::vector1< core::Size > const & setting );
		//void set_working_moving_suite( Size const & setting );

		void set_is_working_res( utility::vector1< core::Size > const & setting );
		void set_full_to_sub( std::map< core::Size, core::Size > const & setting );

		void set_fold_tree( core::kinematics::FoldTree const & setting );
		void set_Is_prepend_map( std::map< core::Size, bool > const & setting );


		void set_chain_boundaries( utility::vector1< std::pair < core::Size, core::Size > > const & setting );
		//void set_which_chain_has_moving_res( Size const & setting );
		void set_gap_size( Size const & setting );
		void set_five_prime_chain_break_res( Size const & setting );

		void set_Is_prepend( bool const & setting );
		void set_Is_internal( bool const & setting );
		void set_partition_definition( ObjexxFCL::FArray1D < bool > const & setting );

		void set_working_native_pose( core::pose::PoseOP & pose );
		void set_working_native_pose( core::pose::PoseCOP pose );
		void set_working_fixed_res(	utility::vector1< core::Size > const & working_fixed_res );
		void set_rmsd_res_list(	utility::vector1< core::Size > const & rmsd_res_list );
		void set_working_terminal_res(	utility::vector1< core::Size > const & working_terminal_res );
		void set_working_moving_partition_pos(	utility::vector1< core::Size > const & working_moving_partition_pos );
		void set_input_res_vectors(	utility::vector1< utility::vector1< Size > > const & setting );
		void set_cutpoint_closed_list( utility::vector1< core::Size >  const & setting );
		void set_working_best_alignment( utility::vector1< core::Size > const & setting );

		void set_native_alignment( utility::vector1< core::Size > const & setting );
		void set_working_native_alignment( utility::vector1< core::Size > const & setting );

		void set_global_sample_res_list( utility::vector1< core::Size > const & setting );
		void set_force_syn_chi_res_list( utility::vector1< core::Size > const & setting );
		void set_force_north_ribose_list( utility::vector1< core::Size > const & setting );
		void set_force_south_ribose_list( utility::vector1< core::Size > const & setting );
		void set_protonated_H1_adenosine_list( utility::vector1< core::Size > const & setting );
		void set_add_virt_res_as_root( bool const setting ){ add_virt_res_as_root_ = setting; }


		core::pose::PoseCOP	working_native_pose() const;



	private:

		void update_working_moving_suite();
		void update_working_sequence();

		std::map< core::Size, core::Size > create_sub_to_full_map( std::map< core::Size, core::Size > const & full_to_sub ) const;

	private:

		bool output_extra_RMSDs_; //Used in StepWiseRNA_Output_Data.cc
		bool Is_simple_full_length_job_params_;

		std::string full_sequence_;
		std::string working_sequence_;

		Size moving_res_;
		utility::vector1< core::Size > working_moving_res_list_;
		utility::vector1< core::Size > working_moving_suite_list_;
		Size working_moving_res_;
		Size working_moving_suite_;

		utility::vector1< core::Size >  is_working_res_;
		std::map< core::Size, core::Size > full_to_sub_;
		std::map< core::Size, core::Size > sub_to_full_;
		std::map< core::Size, bool > Is_prepend_map_;

		utility::vector1< std::pair < core::Size, core::Size > > chain_boundaries_;
		//Size which_chain_has_moving_res_;
		Size gap_size_;
		Size five_prime_chain_break_res_;

		bool Is_prepend_;
		bool Is_internal_;
		bool add_virt_res_as_root_;

		ObjexxFCL::FArray1D < bool > partition_definition_;

		core::pose::PoseCOP working_native_pose_;

		utility::vector1< core::Size > working_fixed_res_;
		utility::vector1< core::Size > rmsd_res_list_;
		utility::vector1< core::Size > working_terminal_res_;
		utility::vector1< core::Size > working_moving_partition_pos_;
		core::kinematics::FoldTree fold_tree_;
		utility::vector1< utility::vector1< core::Size > > input_res_vectors_;
		utility::vector1< core::Size > cutpoint_closed_list_;
		utility::vector1< core::Size > working_best_alignment_;

		utility::vector1< core::Size > native_alignment_;
		utility::vector1< core::Size > working_native_alignment_;

		utility::vector1< core::Size > global_sample_res_list_;
		utility::vector1< core::Size > working_global_sample_res_list_;

		utility::vector1< core::Size > force_syn_chi_res_list_;
		utility::vector1< core::Size > working_force_syn_chi_res_list_;

		utility::vector1< core::Size > force_north_ribose_list_;
		utility::vector1< core::Size > working_force_north_ribose_list_;

		utility::vector1< core::Size > force_south_ribose_list_;
		utility::vector1< core::Size > working_force_south_ribose_list_;


		utility::vector1< core::Size > protonated_H1_adenosine_list_;
		utility::vector1< core::Size > working_protonated_H1_adenosine_list_;


  };

}
} //swa
} // protocols

#endif

