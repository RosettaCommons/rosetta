// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/setup/StepWiseCSA_JobDistributor.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_setup_StepWiseCSA_JobDistributor_HH
#define INCLUDED_protocols_stepwise_setup_StepWiseCSA_JobDistributor_HH

#include <protocols/stepwise/setup/StepWiseJobDistributor.hh>
#include <protocols/stepwise/setup/StepWiseCSA_JobDistributor.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>

namespace protocols {
namespace stepwise {
namespace setup {

class StepWiseCSA_JobDistributor: public StepWiseJobDistributor {

public:

	StepWiseCSA_JobDistributor( stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
		std::string const & silent_file,
		core::Size const nstruct,
		core::Size const csa_bank_size,
		core::Real const csa_rmsd,
		bool const output_round_silent_files );

	//destructor
	~StepWiseCSA_JobDistributor();

public:

	virtual
	void
	apply( core::pose::Pose & pose );

	virtual
	void
	initialize( core::pose::Pose const & pose );

	virtual
	bool
	has_another_job();

private:

	Size get_updates( core::io::silent::SilentStructCOP s ) const;
	void set_updates( core::io::silent::SilentStructOP s, Size const updates ) const;

	void
	update_bank( core::pose::Pose & pose );

	void
	read_in_silent_file();

	void
	write_out_silent_file( std::string const & silent_file = "" );

	void
	write_out_round_silent_file();

	void
	put_lock_on_silent_file();

	void
	free_lock_on_silent_file();

	bool
	check_for_closeness( core::pose::Pose & pose_test,
		core::pose::Pose const & full_model_pose ) const;

private:

	std::string const lock_file;
	core::Size const csa_bank_size_;
	core::Size const total_updates_;
	core::Real const csa_rmsd_;
	bool const output_round_silent_files_;
	core::Size total_updates_so_far_;
	core::io::silent::SilentFileDataOP sfd_;


};

} //setup
} //stepwise
} //protocols

#endif
