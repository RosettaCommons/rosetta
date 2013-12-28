// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @detailed
///
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_StepWiseJobParameters_HH
#define INCLUDED_protocols_stepwise_StepWiseJobParameters_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

#include <string>
#include <map>

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

	class StepWiseJobParameters: public utility::pointer::ReferenceCount {
	public:

		//constructor!
    StepWiseJobParameters();

	  //destructor -- necessary?
    virtual ~StepWiseJobParameters();

		std::string const & sequence() const;
		std::string const & working_sequence() const;
		utility::vector1< core::Size > const & working_res_list() const;

		utility::vector1< core::Size > const & working_moving_res_list() const;
		utility::vector1< core::Size > const & working_moving_suite_list() const;

		Size actually_moving_res() const;

		ObjexxFCL::FArray1D< core::Size > const & is_working_res() const;
		ObjexxFCL::FArray1D< core::Size > const & is_moving_res() const;
		std::map< core::Size, core::Size > & full_to_sub();
		std::map< core::Size, core::Size > & sub_to_full();

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries() const;
		Size const & which_chain_has_moving_res() const;
		Size const & gap_size() const;
		Size const & first_chain_break_res() const;

		bool const & is_prepend() const;
		bool const & is_internal() const;

		ObjexxFCL::FArray1D< bool > const & partition_definition() const;

		utility::vector1< core::Size > const &  working_fixed_res() const;
		utility::vector1< core::Size > const &  working_terminal_res() const;
		utility::vector1< core::Size > const &  working_superimpose_res() const;
		utility::vector1< core::Size > const &  working_calc_rms_res() const;
		utility::vector1< core::Size > const &  working_bridge_res() const;
		utility::vector1< core::Size > const &  moving_pos() const;

		void set_sequence( std::string const & setting );
		void set_working_sequence( std::string const & setting );
		void set_working_res_list( utility::vector1< core::Size > const & setting );

		void set_working_moving_res_list( utility::vector1< core::Size > const & setting );
		void set_working_moving_suite_list( utility::vector1< core::Size > const & setting );

		void set_is_working_res( ObjexxFCL::FArray1D< core::Size > const & setting );
		void set_is_moving_res( ObjexxFCL::FArray1D< core::Size > const & setting );
		void set_full_to_sub( std::map< core::Size, core::Size > const & setting );

		std::map< core::Size, core::Size >
		create_sub_to_full_map(std::map< core::Size, core::Size > const & full_to_sub) const;

		void set_chain_boundaries( utility::vector1< std::pair< core::Size, core::Size > > const & setting );
		void set_which_chain_has_moving_res( Size const & setting );
		void set_gap_size( Size const & setting );
		void set_first_chain_break_res( Size const & setting );

		void set_is_prepend( bool const & setting );
		void set_is_internal( bool const & setting );
		void set_partition_definition( ObjexxFCL::FArray1D< bool > const & setting );

		void set_working_native_pose( core::pose::PoseOP & pose );
		void set_working_fixed_res(	utility::vector1< core::Size > const & working_fixed_res);
		void set_working_terminal_res(	utility::vector1< core::Size > const & working_terminal_res);
		void set_working_superimpose_res(	utility::vector1< core::Size > const & working_superimpose_res);
		void set_working_calc_rms_res(	utility::vector1< core::Size > const & working_calc_rms_res);
		void set_working_bridge_res(	utility::vector1< core::Size > const & working_bridge_res);
		void set_moving_pos(	utility::vector1< core::Size > const & moving_pos);

		utility::vector1< bool >const is_pre_proline() const;

		core::pose::PoseCOP	working_native_pose() const;

		utility::vector1< Size >
		apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector);

		Size
		apply_full_to_sub_mapping( Size const res ) const;

	private:

		std::string sequence_;
		std::string working_sequence_;

		utility::vector1< core::Size > working_res_list_;

		utility::vector1< core::Size > working_moving_res_list_;
		utility::vector1< core::Size > working_moving_suite_list_;

		ObjexxFCL::FArray1D< core::Size > is_working_res_;
		ObjexxFCL::FArray1D< core::Size > is_moving_res_;
		std::map< core::Size, core::Size > full_to_sub_;
		std::map< core::Size, core::Size > sub_to_full_;

		utility::vector1< std::pair< core::Size, core::Size > > chain_boundaries_;
		Size which_chain_has_moving_res_;
		Size gap_size_;
		Size first_chain_break_res_;

		bool is_prepend_;
		bool is_internal_;

		ObjexxFCL::FArray1D< bool > partition_definition_;

		core::pose::PoseOP working_native_pose_;

		utility::vector1< core::Size > working_fixed_res_;
		utility::vector1< core::Size > working_terminal_res_;
		utility::vector1< core::Size > working_superimpose_res_;
		utility::vector1< core::Size > working_calc_rms_res_;
		utility::vector1< core::Size > working_bridge_res_;
		utility::vector1< core::Size > moving_pos_;

  };

} //protein
} //enumerate
} //stepwise
} //protocols

#endif

