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
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinMinimizer_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinMinimizer_HH

#include <protocols/stepwise/modeler/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinMinimizer.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <string>
#include <map>

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace protein {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseProteinMinimizer: public protocols::moves::Mover {
public:

	StepWiseProteinMinimizer( utility::vector1< pose::PoseOP > const & pose_list,
		utility::vector1< Size > const & moving_residues );

	//destructor -- necessary?
	~StepWiseProteinMinimizer();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

	void set_min_tolerance( core::Real const & setting );
	void set_min_type( std::string const & setting );
	void set_fixed_res( utility::vector1< core::Size > const & fixed_res );
	void set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res );

	core::Size const & num_pose_minimize() const { return num_pose_minimize_; }
	void set_num_pose_minimize( core::Size const & setting ){ num_pose_minimize_ = setting; }

	void
	set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	void
	set_move_takeoff_torsions( bool const value ){ move_takeoff_torsions_ = value; }

	void
	set_rescore_only( bool const & setting ){ rescore_only_ = setting; }

	void
	set_move_jumps_between_chains( bool const & setting ){ move_jumps_between_chains_ = setting; }

	void
	set_cartesian( bool const setting );
	//  void
	//  set_constraint_set( core::scoring::constraints::ConstraintSetOP const & cst_set );

	void set_use_coordinate_constraints( bool const & setting ){ use_coordinate_constraints_ = setting; }
	bool use_coordinate_constraints() const{ return use_coordinate_constraints_; }

	utility::vector1< pose::PoseOP >
	pose_list() const { return pose_list_; }

private:

	void
	initialize_parameters();

	void
	let_neighboring_chis_minimize(
		core::kinematics::MoveMap & mm,
		core::pose::Pose & pose );

	bool
	pose_has_chainbreak( core::pose::Pose const & pose );

	utility::vector1< core::Size > const & moving_residues_;
	utility::vector1< core::Size > fixed_res_;
	utility::vector1< core::Size > calc_rms_res_;
	bool move_takeoff_torsions_;
	bool rescore_only_;
	bool move_jumps_between_chains_;
	bool cartesian_;

	core::scoring::ScoreFunctionOP fa_scorefxn_;

	utility::vector1< pose::PoseOP > pose_list_;

	std::string min_type_;
	core::Real min_tolerance_;
	bool use_coordinate_constraints_;
	Size num_pose_minimize_;
};

} //protein
} //modeler
} //legacy
} //stepwise
} //protocols

#endif
