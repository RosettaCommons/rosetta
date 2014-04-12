// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_Packer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinPacker_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinPacker_HH

#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <string>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseProteinPacker: public protocols::moves::Mover {
  public:

    //constructor!
		StepWiseProteinPacker( utility::vector1< Size > const & moving_residues );

    //destructor -- necessary?
    ~StepWiseProteinPacker();

    /// @brief Apply the minimizer to one pose
    virtual void apply( pose::Pose & pose );

		virtual std::string get_name() const;

		void
		initialize( pose::Pose & pose );

		void
		set_scorefxn( scoring::ScoreFunctionOP const & scorefxn );

    void
		set_use_green_packer( bool const & setting );

    void
		set_use_packer_instead_of_rotamer_trials( bool const & setting );

    void
		set_rescore_only( bool const & setting ){ rescore_only_ = setting; }

		void
		do_prepack( pose::Pose & pose );

		void
		set_moving_res_list( utility::vector1< Size > const & setting ){ moving_res_list_ = setting; }

		void set_allow_virtual_side_chains( bool const & setting ){ allow_virtual_side_chains_ = setting; }
		bool allow_virtual_side_chains() const{ return allow_virtual_side_chains_; }

  private:

		void
		initialize_moving_residues_including_junction( Size const & nres );

		void
		initialize_green_packer( Size const & nres );

		void
		initialize_for_regular_packer( pose::Pose const & pose );

		void
		figure_out_neighbors( pose::Pose & pose,
													utility::vector1< bool > & residues_allowed_to_be_packed );


		void
		reinstate_side_chain_angles( pose::Pose & pose, pose::Pose const & src_pose );

		void
		apply_regular_packer( pose::Pose & pose, bool const pack_at_neighbors_only );

	private:

		utility::vector1< Size > const moving_residues_;
		scoring::ScoreFunctionOP scorefxn_;
		protocols::simple_moves::GreenPackerOP green_packer_;
		bool use_green_packer_;
		bool use_packer_instead_of_rotamer_trials_;
		bool pack_at_neighbors_only_;
		bool allow_virtual_side_chains_;
		bool rescore_only_;
		pack::task::PackerTaskOP pack_task_;
		pose::PoseOP pose_init_;

		utility::vector1< Size > moving_res_list_; // is this ever distinct from moving_residues_ [?]

  };

} //protein
} //sampling
} //stepwise
} //protocols

#endif
