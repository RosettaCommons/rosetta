// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/SugarModeling.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_SugarModeling_HH
#define INCLUDED_protocols_swa_rna_SugarModeling_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/rna/SugarModeling.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace swa {
namespace rna {

	class SugarModeling: public utility::pointer::ReferenceCount {

	public:

		SugarModeling( core::Size const input_moving_res, core::Size const input_reference_res );

		SugarModeling();

		~SugarModeling();

		void
		check_compatibility( core::Size const nres ) const;

		void
		set_base_and_pucker_state( core::pose::Pose const & pose, StepWiseRNA_JobParametersCOP const & JP );

	public:
		bool sample_sugar;
		core::Size moving_res;
		core::Size reference_res;
		bool is_prepend;
		core::Size bulge_res;
		core::Size bulge_suite;
		core::Size five_prime_chain_break;
		core::Size moving_res_pucker_state;
		core::Size bulge_res_pucker_state;
		core::Size moving_res_base_state;
		core::Size bulge_res_base_state;
		utility::vector1< core::pose::PoseOP > pose_list; //pose_data_list of possible sugar conformations.

	};


} //rna
} //swa
} //protocols

#endif
