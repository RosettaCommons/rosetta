// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/FARNA_Optimizer.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_FARNA_Optimizer_HH
#define INCLUDED_protocols_farna_FARNA_Optimizer_HH

#include <protocols/moves/Mover.hh>
#include <protocols/farna/FARNA_Optimizer.fwd.hh>
#include <protocols/farna/RNA_FragmentMonteCarlo.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace farna {

class FARNA_Optimizer: public protocols::moves::Mover {

public:

	//constructor
	FARNA_Optimizer( utility::vector1< core::pose::PoseOP > const & pose_list,
		core::scoring::ScoreFunctionCOP scorefxn,
		core::Size cycles = 0 );

	//destructor
	~FARNA_Optimizer();

public:


	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const{ return "FARNA_Optimizer"; }

	void set_cycles( core::Size const & setting ){ cycles_ = setting; }
	core::Size cycles() const { return cycles_; }

	RNA_FragmentMonteCarloCOP rna_fragment_monte_carlo() const { return rna_fragment_monte_carlo_; }

private:

	utility::vector1< core::pose::PoseOP > pose_list_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	Size cycles_;

	RNA_FragmentMonteCarloOP rna_fragment_monte_carlo_;

};

} //farna
} //protocols

#endif
