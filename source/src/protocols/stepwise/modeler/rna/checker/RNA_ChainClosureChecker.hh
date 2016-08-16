// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/checker/RNA_ChainClosureChecker.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_checker_RNA_ChainClosureChecker_HH
#define INCLUDED_protocols_stepwise_rna_checker_RNA_ChainClosureChecker_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosureChecker.fwd.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {

class RNA_ChainClosureChecker: public protocols::moves::Mover {

public:

	//constructor
	RNA_ChainClosureChecker( pose::Pose const & pose, Size const five_prime_res );


	//destructor
	~RNA_ChainClosureChecker();

public:

	void
	apply( core::pose::Pose & pose ){ copy_CCD_torsions( pose ); }

	std::string get_name() const{ return "RNA_ChainClosureChecker"; }

	void
	set_reinitialize_CCD_torsions( bool const & setting ){ reinitialize_CCD_torsions_ = setting; };

	//void
	//add_harmonic_chain_break_constraint( Size const five_prime_res );

	void
	copy_CCD_torsions( pose::Pose & pose ) const;

	void
	copy_CCD_torsions_general( pose::Pose & pose, Size const five_prime_res, Size const three_prime_res ) const;

	bool
	check_loop_closed( pose::Pose const & pose );

	bool
	chain_break_screening_general( pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & chainbreak_scorefxn, Size const five_ );

	//bool
	//chain_break_screening( pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & chainbreak_scorefxn );

	bool
	check_screen();

	bool
	check_screen( pose::Pose & pose );

	pose::Pose & pose(){ return chain_break_screening_pose_; }

	Size const & five_prime_res() const{ return five_prime_res_;}

private:

	pose::Pose chain_break_screening_pose_;
	Size const five_prime_res_;
	bool reinitialize_CCD_torsions_;
	bool verbose_;

	core::scoring::ScoreFunctionOP chain_break_scorefxn_;

	StepWiseRNA_CountStruct count_data_;

};

} //checker
} //rna
} //modeler
} //stepwise
} //protocols

#endif
