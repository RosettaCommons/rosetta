// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/BaseCentroidScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_BaseCentroidScreener_HH
#define INCLUDED_protocols_stepwise_screener_BaseCentroidScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/BaseCentroidScreener.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>

namespace protocols {
namespace stepwise {
namespace screener {

	class BaseCentroidScreener: public StepWiseScreener {

	public:

		//constructor
		BaseCentroidScreener( sampling::rna::checker::RNA_BaseCentroidCheckerOP base_centroid_checker,
													core::pose::PoseOP screening_pose,
													bool const force_centroid_interaction = true );

		//constructor
		BaseCentroidScreener( sampling::rna::checker::RNA_BaseCentroidCheckerOP base_centroid_checker,
													core::kinematics::Stub const & moving_res_base_stub );

		//destructor
		~BaseCentroidScreener();

	public:

		std::string
		name() const { return "BaseCentroidScreener"; }

		StepWiseScreenerType
		type() const { return BASE_CENTROID; }

		virtual
		bool
		check_screen();

		void
		fast_forward( rotamer_sampler::RotamerBaseOP sampler );

	private:

		sampling::rna::checker::RNA_BaseCentroidCheckerOP base_centroid_checker_;
		core::pose::PoseOP screening_pose_;
		bool const force_centroid_interaction_;

		bool const using_stub_;
		core::kinematics::Stub const & moving_res_base_stub_;

		sampling::rna::StepWiseRNA_CountStruct count_data_;

	};

} //screener
} //stepwise
} //protocols

#endif
