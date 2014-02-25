// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinResidueSampler_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinResidueSampler_HH

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>
#include <string>
#include <map>

//Auto Headers
#include <core/id/TorsionID.fwd.hh>
namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	//	typedef std::map< std::string, core::pose::PoseOP > PoseList;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseProteinResidueSampler: public protocols::moves::Mover {
  public:

    //constructor!
		StepWiseProteinResidueSampler(
													 utility::vector1< Size > const & moving_residues,
													 utility::vector1< core::id::TorsionID > const & which_torsions,
													 utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists );


    //constructor!
		StepWiseProteinResidueSampler(
													 utility::vector1< Size > const & moving_residues,
													 utility::vector1< core::id::TorsionID > const & which_torsions,
													 utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists,
													 Size const which_jump,
													 utility::vector1< core::kinematics::Jump > const & jumps );

    //destructor -- necessary?
    ~StepWiseProteinResidueSampler();

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;


    void
		set_silent_file( std::string const & setting );

		void
		set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

    void
		set_use_green_packer( bool const & setting );

    void
		set_use_packer_instead_of_rotamer_trials( bool const & setting );

    void
		set_do_prepack( bool const & setting );

    void set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res );

		core::io::silent::SilentFileDataOP & silent_file_data();

  private:

		void
		sample_residues( core::pose::Pose & pose );

		void
		quick_output( core::pose::Pose & pose,
									std::string const & tag	 );

		void
		initialize_moving_residues_including_junction( Size const & nres );

		void
		initialize_green_packer( core::Size const & nres );

		void
		initialize_for_regular_packer( core::pose::Pose & pose );

		void
		figure_out_neighbors( core::pose::Pose & pose,
													utility::vector1< bool > & residues_allowed_to_be_packed );

		void
		reinstate_side_chain_angles( core::pose::Pose & pose, core::pose::Pose const & src_pose );

		void
		apply_regular_packer( core::pose::Pose & pose );

	private:

		utility::vector1< Size > const moving_residues_;
		utility::vector1< core::id::TorsionID > const which_torsions_;
		utility::vector1< utility::vector1< core::Real > > main_chain_torsion_set_lists_;
		//		PoseList pose_list_;
		core::scoring::ScoreFunctionOP scorefxn_;
		protocols::simple_moves::GreenPackerOP green_packer_;
		bool use_green_packer_;
		bool use_packer_instead_of_rotamer_trials_;
		bool pack_at_neighbors_only_;
		bool do_prepack_;
		core::pack::task::PackerTaskOP pack_task_;
		core::pose::PoseOP pose_prepack_;

		std::string silent_file_;

		core::io::silent::SilentFileDataOP sfd_;

		utility::vector1< core::Size > calc_rms_res_;

		Size const which_jump_;
		utility::vector1< core::kinematics::Jump > jumps_;

  };

} //protein
} //sampling
} //stepwise
} //protocols


#endif
