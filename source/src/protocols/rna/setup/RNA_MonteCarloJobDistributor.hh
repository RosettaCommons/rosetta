// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/setup/RNA_MonteCarloJobDistributor.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rna_setup_RNA_MonteCarloJobDistributor_HH
#define INCLUDED_protocols_rna_setup_RNA_MonteCarloJobDistributor_HH

#include <protocols/rna/setup/RNA_JobDistributor.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.fwd.hh>
#include <protocols/rna/setup/RNA_MonteCarloJobDistributor.fwd.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace rna {
namespace setup {

class RNA_MonteCarloJobDistributor: public RNA_JobDistributor {

public:

	//constructor
	RNA_MonteCarloJobDistributor( stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
		std::string const & silent_file,
		core::Size const nstruct );

	//constructor
	RNA_MonteCarloJobDistributor( rna::denovo::RNA_FragmentMonteCarloOP rna_fragment_monte_carlo,
		std::string const & silent_file,
		core::Size const nstruct );

	//destructor
	~RNA_MonteCarloJobDistributor();

public:

	virtual std::string get_name() const {
		return "RNA_MonteCarloJobDistributor";
	}

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

	void
	move_forward_to_next_model();

	bool
	get_out_tag();

private:

	core::Size count_;
	std::string out_tag_;

	bool init_tag_is_done_;
	std::map< std::string, bool > tag_is_done_;


};

} //setup
} //rna
} //protocols

#endif
