// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/setup/RNA_JobDistributor.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rna_setup_RNA_JobDistributor_HH
#define INCLUDED_protocols_rna_setup_RNA_JobDistributor_HH

#include <protocols/moves/Mover.hh>
#include <protocols/rna/setup/RNA_JobDistributor.fwd.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.fwd.hh>
//#include <protocols/rna/denovo/RNA_DeNovoProtocol.fwd.hh>

// AMW TODO: write a forward header
#include <utility/pointer/owning_ptr.hh>
namespace protocols { namespace rna { namespace denovo {
class RNA_FragmentMonteCarlo;
typedef utility::pointer::shared_ptr< RNA_FragmentMonteCarlo > RNA_FragmentMonteCarloOP;
} } }
#include <core/types.hh>

namespace protocols {
namespace rna {
namespace setup {

class RNA_JobDistributor: public protocols::moves::Mover {

public:

	RNA_JobDistributor( stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
		std::string const & silent_file,
		core::Size const nstruct ):
		stepwise_monte_carlo_( stepwise_monte_carlo ),
		rna_fragment_monte_carlo_( nullptr ),
		silent_file_( silent_file ),
		nstruct_( nstruct ),
		superimpose_over_all_( false )
	{
	}

	RNA_JobDistributor( rna::denovo::RNA_FragmentMonteCarloOP rna_fragment_monte_carlo,
		std::string const & silent_file,
		core::Size const nstruct ):
		stepwise_monte_carlo_( nullptr ),
		rna_fragment_monte_carlo_( rna_fragment_monte_carlo ),
		silent_file_( silent_file ),
		nstruct_( nstruct ),
		superimpose_over_all_( false )
	{
	}

	virtual std::string get_name() const {
		return "RNA_JobDistributor";
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

	// The True way to achieve this effect would be through a sum type
	// (in... some other language) or via a superclass. Upon unification with
	// jd3, the natural superclass would be Mover, but at the moment this
	// makes some assumptions. So we use two pointers and a runtime_assert instead.
	stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo_;
	rna::denovo::RNA_FragmentMonteCarloOP rna_fragment_monte_carlo_;
	std::string const silent_file_;
	core::Size const nstruct_;
	bool superimpose_over_all_;
	core::pose::PoseCOP start_pose_;

};

} //setup
} //rna
} //protocols

#endif
