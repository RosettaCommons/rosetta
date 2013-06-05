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


#ifndef INCLUDED_protocols_swa_protein_StepWiseProteinPoseMinimizer_HH
#define INCLUDED_protocols_swa_protein_StepWiseProteinPoseMinimizer_HH

#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>

namespace protocols {
namespace swa {
namespace protein {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseProteinPoseMinimizer: public protocols::moves::Mover {
  public:

    //constructor!
		StepWiseProteinPoseMinimizer( core::io::silent::SilentFileDataOP sfd, utility::vector1< Size > const & moving_residues );

    StepWiseProteinPoseMinimizer( PoseList & pose_list, utility::vector1< Size > const & moving_residues );

    //destructor -- necessary?
    ~StepWiseProteinPoseMinimizer();

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

    void set_silent_file( std::string const & setting );
    void set_min_tolerance( core::Real const & setting );
    void set_min_type( std::string const & setting );
    void set_fixed_res( utility::vector1< core::Size > const & fixed_res );
    void set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res );

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
		//		void
		//		set_constraint_set( core::scoring::constraints::ConstraintSetOP const & cst_set );

		core::io::silent::SilentFileDataOP & silent_file_data();

  private:

		void
		initialize_parameters();

		void
		initialize_protein_input_silent_file_data_from_pose_list( PoseList & pose_list );

		void
		let_neighboring_chis_minimize(
																	core::kinematics::MoveMap & mm,
																	core::pose::Pose & pose );

		bool
		pose_has_chainbreak( core::pose::Pose const & pose );

    PoseList pose_list_;
    utility::vector1< core::Size > const & moving_residues_;
    utility::vector1< core::Size > fixed_res_;
    utility::vector1< core::Size > calc_rms_res_;
		bool move_takeoff_torsions_;
		bool rescore_only_;
		bool move_jumps_between_chains_;
		bool cartesian_;
    std::string silent_file_;

		core::scoring::ScoreFunctionOP fa_scorefxn_;

		core::io::silent::SilentFileDataOP sfd_, input_silent_file_data_;

		std::string min_type_;
		core::Real min_tolerance_;

  };

} //swa
} // protocols

}
#endif
