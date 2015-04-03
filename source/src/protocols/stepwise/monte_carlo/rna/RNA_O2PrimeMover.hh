// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_O2PrimeMover.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_monte_carlo_RNA_O2PrimeMover_hh
#define INCLUDED_protocols_stepwise_monte_carlo_RNA_O2PrimeMover_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_O2PrimeMover.fwd.hh>


namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class RNA_O2PrimeMover: public protocols::moves::Mover {
public:

	RNA_O2PrimeMover( core::scoring::ScoreFunctionOP scorefxn,
									 bool const sample_all_o2prime,
									 core::Real const sample_range_small,
									 core::Real const sample_range_large );

	//destructor -- necessary? -- YES destructors are necessary.
	~RNA_O2PrimeMover();

	// Undefinded, commenting out to fix PyRosetta build  void apply( core::pose::Pose & pose, Size const res_to_delete, protocols::stepwise::monte_carlo::MovingResidueCase const moving_residue_case  );

	/// @brief Apply the minimizer to one pose
	using protocols::moves::Mover::apply;
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

  void
	apply( core::pose::Pose & pose, std::string & move_type );

private:

	void
	sample_near_o2prime_torsion( core::pose::Pose & pose, Size const moving_res, core::Real const sample_range);

	Size
	get_random_o2prime_residue( core::pose::Pose & pose );

	Size
	get_random_o2prime_residue_near_moving_residue( core::pose::Pose & pose, utility::vector1< Size > const moving_res_list );

private:

	core::scoring::ScoreFunctionOP scorefxn_;
	bool const sample_all_o2prime_;
	core::Real const sample_range_small_, sample_range_large_;

};

} //rna
} //monte_carlo
} //stepwise
} //protocols

#endif
