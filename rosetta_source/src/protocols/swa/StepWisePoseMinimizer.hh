// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_PoseMinimizer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_StepWisePoseMinimizer_hh
#define INCLUDED_protocols_swa_StepWisePoseMinimizer_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>

namespace protocols {
namespace swa {

typedef std::map< std::string, core::pose::PoseOP > PoseList;

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWisePoseMinimizer: public protocols::moves::Mover {
public:

	//constructor!
	StepWisePoseMinimizer( PoseList & pose_list, utility::vector1< Size > const & moving_residues );

	//destructor -- necessary?
	~StepWisePoseMinimizer();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void set_silent_file( std::string const & setting );

	void
	set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	//		void
	//		set_constraint_set( core::scoring::constraints::ConstraintSetOP const & cst_set );

private:

	void limit_minimize( core::kinematics::MoveMap & mm, Size const & nres );

	PoseList pose_list_;
	utility::vector1< core::Size > const & moving_residues_;
	std::string silent_file_;

	core::scoring::ScoreFunctionOP fa_scorefxn_;

	//		core::scoring::constraints::ConstraintSetOP cst_set_;

};

} //swa
} // protocols

#endif
