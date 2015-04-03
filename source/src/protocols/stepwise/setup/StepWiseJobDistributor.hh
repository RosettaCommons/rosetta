// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/setup/StepWiseJobDistributor.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_setup_StepWiseJobDistributor_HH
#define INCLUDED_protocols_stepwise_setup_StepWiseJobDistributor_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/setup/StepWiseJobDistributor.fwd.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace setup {

	class StepWiseJobDistributor: public protocols::moves::Mover {

	public:

		StepWiseJobDistributor( stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
														std::string const silent_file,
														core::Size const nstruct ):
		stepwise_monte_carlo_( stepwise_monte_carlo ),
		silent_file_( silent_file ),
		nstruct_( nstruct ),
		superimpose_over_all_( false )
		{
		}


		virtual std::string get_name() const {
			return "StepWiseJobDistributor";
		}

		virtual
		void
		apply( core::pose::Pose & pose ) = 0;

		virtual
		void
		initialize( core::pose::Pose const & pose ) = 0;

		virtual
		bool
		has_another_job() = 0;

		void set_superimpose_over_all( bool const & setting ){ superimpose_over_all_ = setting; }
		bool superimpose_over_all() const { return superimpose_over_all_; }

	protected:

		stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo_;
		std::string const silent_file_;
		core::Size const nstruct_;
		bool superimpose_over_all_;
		core::pose::PoseCOP start_pose_;

	};

} //setup
} //stepwise
} //protocols

#endif
