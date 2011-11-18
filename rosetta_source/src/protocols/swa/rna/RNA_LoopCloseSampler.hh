// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_LoopCloseSampler
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_LoopCloseSampler_HH
#define INCLUDED_protocols_swa_rna_LoopCloseSampler_HH

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class RNA_LoopCloseSampler: public protocols::moves::Mover {
  public:

    //constructor!
    RNA_LoopCloseSampler( Size const moving_suite, Size const chainbreak_suite );

    //destructor -- necessary?
    ~RNA_LoopCloseSampler();

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;



    utility::vector1< utility::vector1< Real > >  all_torsion_info(){ return all_torsion_info_;}

		void set_bin_size( Size const setting ){ bin_size_ = setting; }

		void set_rep_cutoff( Real const rep_cutoff ){ rep_cutoff_ = rep_cutoff; }
		void set_torsion_range( Real const torsion_range ){ torsion_range_ = torsion_range; }
		void set_torsion_increment( Real const torsion_increment ){ torsion_increment_ = torsion_increment; }

		void set_center_around_native( bool const setting ){ center_around_native_ = setting; }
		void set_calculate_jacobian( bool const setting ){ calculate_jacobian_ = setting; }
		void set_save_torsion_info( bool const setting ){ save_torsion_info_ = setting; }
		void set_just_output_score( bool const setting ){ just_output_score_ = setting; }

		void set_rbs_new_pair( utility::vector1< Real > const & rbs_new_pair ){ rbs_new_pair_ = rbs_new_pair; }

		void set_silent_file( std::string const & silent_file);

		void set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

		core::io::silent::SilentFileDataOP silent_file_data();

  private:

		bool
		torsion_angles_within_cutoffs( pose::Pose const & pose,
																	 Size const moving_suite,
																	 Size const chainbreak_suite,
																	 int const bin_size_,
																	 int const bins2_ );

		void
		calculate_perturbed_solutions( pose::Pose const & pose,
																	 Size const moving_suite,
																	 Size const chainbreak_suite,
																	 utility::vector1< Real > const & rbs,
																	 Real const perturbation_size,
																	 utility::vector1< utility::vector1< utility::vector1< Real > > >  & perturbed_solution_torsions );

		utility::vector1< Real >
		find_closest_solution(
													utility::vector1< utility::vector1< Real > >  const & perturbed_torsions,
													utility::vector1< Real > const & solution_torsions );

		void
		get_jacobian( utility::vector1< Real > const & solution_torsions,
									utility::vector1< utility::vector1< utility::vector1< Real > > >  const & perturbed_solution_torsions,
									Real const perturbation_size,
									utility::vector1< utility::vector1< Real > > & J );


		void initialize_rep_scorefxn();

		Real
		initialize_fa_rep( pose::Pose const & pose,
											 utility::vector1< Size > const & moving_suites,
											 scoring::ScoreFunctionOP rep_scorefxn );

		bool
		check_clash( pose::Pose & pose,
								 Real const & fa_rep_score_baseline,
								 Real const & rep_cutoff_,
								 scoring::ScoreFunctionOP rep_scorefxn );

		bool
		check_close_to_native( pose::Pose const & pose1, pose::Pose const & pose2, Size const suite1, Size const suite2, Real const angle_range );

  private:

		Size const moving_suite_, chainbreak_suite_;

		core::io::silent::SilentFileDataOP sfd_;
		core::scoring::ScoreFunctionOP scorefxn_, rep_scorefxn_;
		Size bin_size_;
		Real const epsilon_range_;
		Real rep_cutoff_, torsion_range_, torsion_increment_;
		bool center_around_native_, calculate_jacobian_, save_torsion_info_, just_output_score_;
		utility::vector1< Real > rbs_new_pair_;

		utility::vector1< Real > torsion_info_;
    utility::vector1< utility::vector1< Real > > all_torsion_info_;
		std::string silent_file_;

  };

}
} //swa
} // protocols

#endif
