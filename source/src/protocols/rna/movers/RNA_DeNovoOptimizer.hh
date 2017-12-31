// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/RNA_DeNovoOptimizer.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_DeNovoOptimizer_HH
#define INCLUDED_protocols_farna_RNA_DeNovoOptimizer_HH

#include <protocols/moves/Mover.hh>
#include <protocols/rna/movers/RNA_DeNovoOptimizer.fwd.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.fwd.hh>
#include <core/import_pose/options/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace rna {
namespace movers {

class RNA_DeNovoOptimizer: public protocols::moves::Mover {

public:

	//constructor
	RNA_DeNovoOptimizer( utility::vector1< core::pose::PoseOP > const & pose_list,
		core::scoring::ScoreFunctionCOP scorefxn,
		core::Size cycles = 0 );

	//destructor
	~RNA_DeNovoOptimizer() override;

public:


	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override{ return "RNA_DeNovoOptimizer"; }

	void set_cycles( core::Size const & setting ){ cycles_ = setting; }
	core::Size cycles() const { return cycles_; }

	protocols::rna::denovo::RNA_FragmentMonteCarloCOP rna_fragment_monte_carlo() const { return rna_fragment_monte_carlo_; }

private:

	utility::vector1< core::pose::PoseOP > pose_list_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	Size cycles_;

	protocols::rna::denovo::RNA_FragmentMonteCarloOP rna_fragment_monte_carlo_;
	core::import_pose::options::RNA_FragmentMonteCarloOptionsOP options_;

};

} //movers
} //rna
} //protocols

#endif
