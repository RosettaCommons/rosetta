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


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_JobParameters_hh
#define INCLUDED_protocols_swa_rna_StepWiseRNA_JobParameters_hh

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
    ~StepWiseRNA_JobParameters();

		std::string const & sequence() const;
		std::string const & working_sequence() const;

		Size const & moving_res() const;
		Size const & working_moving_res() const;
		Size const & working_moving_suite() const;

		Size actually_moving_res() const;

		ObjexxFCL::FArray1D< core::Size > const & is_working_res() const;
		std::map< core::Size, core::Size > & full_to_sub();

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries() const;
		Size const & which_chain_has_moving_res() const;
		Size const & gap_size() const;
		Size const & five_prime_chain_break_res() const;

		bool const & Is_prepend() const;
		bool const & Is_internal() const;

		ObjexxFCL::FArray1D< bool > const & partition_definition() const;

		utility::vector1< core::Size > const &  working_fixed_res() const;
		utility::vector1< core::Size > const &  working_terminal_res() const;
		utility::vector1< core::Size > const &  moving_pos() const;

		void set_sequence( std::string const & setting );
		void set_working_sequence( std::string const & setting );
		void set_moving_res( Size const & setting );
		void set_working_moving_res( Size const & setting );
		void set_working_moving_suite( Size const & setting );

		void set_is_working_res( ObjexxFCL::FArray1D< core::Size > const & setting );
		void set_full_to_sub( std::map< core::Size, core::Size > & setting );

		void set_chain_boundaries( utility::vector1< std::pair< core::Size, core::Size > > const & setting );
		void set_which_chain_has_moving_res( Size const & setting );
		void set_gap_size( Size const & setting );
		void set_five_prime_chain_break_res( Size const & setting );

		void set_Is_prepend( bool const & setting );
		void set_Is_internal( bool const & setting );
		void set_partition_definition( ObjexxFCL::FArray1D< bool > const & setting );

		void set_working_native_pose( core::pose::PoseOP & pose );
		void set_working_fixed_res(	utility::vector1< core::Size > const & working_fixed_res);
		void set_working_terminal_res(	utility::vector1< core::Size > const & working_terminal_res);
		void set_moving_pos(	utility::vector1< core::Size > const & moving_pos);


		core::pose::PoseCOP	working_native_pose() const;

	private:

		std::string sequence_;
		std::string working_sequence_;

		Size moving_res_;
		Size working_moving_res_;
		Size working_moving_suite_;

		ObjexxFCL::FArray1D< core::Size > is_working_res_;
		std::map< core::Size, core::Size > full_to_sub_;

		utility::vector1< std::pair< core::Size, core::Size > > chain_boundaries_;
		Size which_chain_has_moving_res_;
		Size gap_size_;
		Size five_prime_chain_break_res_;

		bool Is_prepend_;
		bool Is_internal_;

		ObjexxFCL::FArray1D< bool > partition_definition_;

		core::pose::PoseOP working_native_pose_;

		utility::vector1< core::Size > working_fixed_res_;
		utility::vector1< core::Size > working_terminal_res_;
		utility::vector1< core::Size > moving_pos_;

  };

}
} //swa
} // protocols

#endif

