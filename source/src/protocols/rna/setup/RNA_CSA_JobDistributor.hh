// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/setup/RNA_CSA_JobDistributor.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rna_setup_RNA_CSA_JobDistributor_HH
#define INCLUDED_protocols_rna_setup_RNA_CSA_JobDistributor_HH

#include <protocols/rna/setup/RNA_JobDistributor.hh>
#include <protocols/rna/setup/RNA_CSA_JobDistributor.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>

namespace protocols {
namespace rna {
namespace setup {

class RNA_CSA_JobDistributor: public RNA_JobDistributor {

public:

	RNA_CSA_JobDistributor( stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
		std::string const & silent_file,
		core::Size const nstruct,
		core::Size const csa_bank_size,
		core::Real const csa_rmsd,
		bool const output_round_silent_files,
		bool const annealing );

	RNA_CSA_JobDistributor( protocols::rna::denovo::RNA_FragmentMonteCarloOP rna_fragment_monte_carlo,
		std::string const & silent_file,
		core::Size const nstruct,
		core::Size const csa_bank_size,
		core::Real const csa_rmsd,
		bool const output_round_silent_files,
		bool const annealing );

	//destructor
	~RNA_CSA_JobDistributor();

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

	core::Real
	average_pairwise_distance() const;

private:

	std::string const lock_file;
	core::Size const csa_bank_size_;
	core::Size const total_updates_;
	core::Real csa_rmsd_;
	bool const output_round_silent_files_;
	core::Size total_updates_so_far_;
	core::io::silent::SilentFileDataOP sfd_;
	bool const annealing_;


};

} //setup
} //rna
} //protocols

#endif
