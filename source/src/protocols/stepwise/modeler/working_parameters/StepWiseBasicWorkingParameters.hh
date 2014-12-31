// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/working_parameters/StepWiseBasicWorkingParameters.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_StepWiseBasicWorkingParameters_HH
#define INCLUDED_protocols_stepwise_modeler_StepWiseBasicWorkingParameters_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseBasicWorkingParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <map>

#define GAP_SIZE_DUMMY 999 //encodes for an 'infinite' gap between residues, e.g., if on different chains.

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

	class StepWiseBasicWorkingParameters: public utility::pointer::ReferenceCount {

	public:

		//constructor
		StepWiseBasicWorkingParameters();

		//destructor
		~StepWiseBasicWorkingParameters();

	public:

		std::string const & sequence() const { return full_sequence(); }
		void set_sequence( std::string const & setting ){ set_full_sequence( setting ); }

		std::string const & full_sequence() const { return full_sequence_; }
		void set_full_sequence( std::string const & setting ){ full_sequence_ = setting; update_working_sequence(); }
		std::string const & working_sequence() const{ return working_sequence_; }

		Size const & moving_res() const { return moving_res_; }
		void set_moving_res( Size const & setting ) { moving_res_ = setting; }
		Size const & working_moving_res() const { return working_moving_res_; }
		Size const & working_moving_suite() const { return working_moving_suite_; }

		Size actually_moving_res() const; // figured out on the fly!

		utility::vector1< core::Size > working_res_list() const; // figured out on the fly!
		//		void set_working_res_list( utility::vector1< core::Size > const & setting ); // is this in use?

		utility::vector1< core::Size > const & working_moving_res_list() const { return working_moving_res_list_; }
		void set_working_moving_res_list( utility::vector1< core::Size > const & setting ); // updates working_moving_suite!

		utility::vector1< core::Size > const & working_moving_suite_list() const { return working_moving_suite_list_; }
		//		void set_working_moving_suite_list( utility::vector1< core::Size > const & setting ) { working_moving_suite_list_ = setting; }

		utility::vector1< core::Size >  const & is_working_res() const { return is_working_res_; }
		void set_is_working_res( utility::vector1< core::Size > const & setting ){ is_working_res_ = setting; update_working_sequence(); }

		std::map< core::Size, core::Size > & full_to_sub() { return full_to_sub_; }
		void set_full_to_sub( std::map< core::Size, core::Size > const & setting ); // update sub_to_full as well.

		std::map< core::Size, core::Size > & sub_to_full() { return sub_to_full_; }
		std::map< core::Size, core::Size > const & const_full_to_sub() const { return full_to_sub_; }
		std::map< core::Size, core::Size > const & const_sub_to_full() const { return sub_to_full_; }
		core::Size full_to_sub( Size const i ) const {
			if ( full_to_sub_.find( i ) != full_to_sub_.end() ) return full_to_sub_.find( i )->second;
			return 0;
		}
		utility::vector1< Size > apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector);
		Size apply_full_to_sub_mapping( Size const res ) const; // unify with full_to_sub;
		std::map< core::Size, core::Size >
		create_sub_to_full_map( std::map< core::Size, core::Size > const & full_to_sub ) const;

		utility::vector1< core::Size > const & cutpoint_closed_list() const { return cutpoint_closed_list_; }
		void set_cutpoint_closed_list( utility::vector1< core::Size >  const & setting ) { cutpoint_closed_list_ = setting; }

		utility::vector1< core::Size > const & cutpoint_open_list() const { return cutpoint_open_list_; }
		void set_cutpoint_open_list( utility::vector1< core::Size >  const & setting ) { cutpoint_open_list_ = setting; }

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries() const { return chain_boundaries_; }
		void set_chain_boundaries( utility::vector1< std::pair< core::Size, core::Size > > const & setting ) { chain_boundaries_ = setting; }

		// Undefined, commenting out to fix PyRosetta build  Size const & first_chain_break_res() const;

		ObjexxFCL::FArray1D< bool > const & partition_definition() const { return partition_definition_; }
		void set_partition_definition( ObjexxFCL::FArray1D < bool > const & setting ) { partition_definition_ = setting; }
		void set_partition_definition( utility::vector1< Size > const & partition_definition_vector ); // convert vector1 to FArray1D

		bool const & is_prepend() const { return is_prepend_; }
		void set_is_prepend( bool const & setting ) { is_prepend_ = setting; update_working_moving_suite(); }

		bool const & is_internal() const { return is_internal_; }
		void set_is_internal( bool const & setting ) { is_internal_ = setting; }

		Size const & gap_size() const { return gap_size_; }
		void set_gap_size( Size const & setting ) { gap_size_ = setting; }

		utility::vector1< core::Size > const &  fixed_res() const { return fixed_res_; }
		void set_fixed_res(	utility::vector1< core::Size > const & fixed_res ){ fixed_res_ = fixed_res; }

		utility::vector1< core::Size > const &  working_fixed_res() const { return working_fixed_res_; }
		void set_working_fixed_res(	utility::vector1< core::Size > const & working_fixed_res ) { 	working_fixed_res_ = working_fixed_res; }

		utility::vector1< core::Size > const &  calc_rms_res() const { return calc_rms_res_; }
		void set_calc_rms_res(	utility::vector1< core::Size > const & calc_rms_res ){ calc_rms_res_ = calc_rms_res; std::sort( calc_rms_res_.begin(), calc_rms_res_.end() ); }

		utility::vector1< core::Size > const &  working_calc_rms_res() const { return working_calc_rms_res_; }
		void set_working_calc_rms_res(	utility::vector1< core::Size > const & working_calc_rms_res ) { working_calc_rms_res_ = working_calc_rms_res; }

		utility::vector1< core::Size > const &  working_moving_partition_res() const { return working_moving_partition_res_; }
		void set_working_moving_partition_res(	utility::vector1< core::Size > const & working_moving_partition_res ) {		working_moving_partition_res_ = working_moving_partition_res; }

		core::pose::PoseCOP	working_native_pose() const;
		void set_working_native_pose( core::pose::PoseOP & pose );
		void set_working_native_pose( core::pose::PoseCOP pose );

	private:

		void update_working_moving_suite();
		void update_working_sequence();

	protected:

		std::string full_sequence_;
		std::string working_sequence_;

    Size moving_res_;
		Size working_moving_res_;
		Size working_moving_suite_;

		utility::vector1< core::Size > working_moving_res_list_;
		utility::vector1< core::Size > working_moving_suite_list_;

		utility::vector1< core::Size >  is_working_res_;
		std::map< core::Size, core::Size > full_to_sub_;
		std::map< core::Size, core::Size > sub_to_full_;

		utility::vector1< core::Size > cutpoint_closed_list_;
		utility::vector1< core::Size > cutpoint_open_list_;
		utility::vector1< std::pair< core::Size, core::Size > > chain_boundaries_;

		ObjexxFCL::FArray1D < bool > partition_definition_;

		bool is_prepend_;
		bool is_internal_;
		Size gap_size_;

		core::pose::PoseCOP working_native_pose_;

		utility::vector1< core::Size > fixed_res_;
		utility::vector1< core::Size > working_fixed_res_;
		utility::vector1< core::Size > calc_rms_res_;
		utility::vector1< core::Size > working_calc_rms_res_;
		utility::vector1< core::Size > terminal_res_;
		utility::vector1< core::Size > working_terminal_res_;
		utility::vector1< core::Size > working_moving_partition_res_;

	};

} //working_parameters
} //modeler
} //stepwise
} //protocols

#endif
