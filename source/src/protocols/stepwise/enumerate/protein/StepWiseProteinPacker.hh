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

#include <protocols/stepwise/enumerate/protein/sample_generators/StepWisePoseSampleGenerator.fwd.hh>
#include <protocols/stepwise/enumerate/protein/PoseFilter.fwd.hh>
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

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

	//	typedef std::map< std::string, core::pose::PoseOP > PoseList;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseProteinPacker: public protocols::moves::Mover {
  public:

    //constructor!
		StepWiseProteinPacker(
													 utility::vector1< Size > const & moving_residues,
													 protocols::stepwise::enumerate::protein::sample_generators::StepWisePoseSampleGeneratorOP sample_generator );

    //destructor -- necessary?
    ~StepWiseProteinPacker();

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
		set_rescore_only( bool const & setting ){ rescore_only_ = setting; }

    void set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res );

		core::io::silent::SilentFileDataOP & silent_file_data();

		void
		set_pose_filter( protocols::stepwise::enumerate::protein::PoseFilterOP pose_filter );

  private:

		void
		sample_residues( core::pose::Pose & pose );


		void
		print_tag( std::string const & tag, Size const k );

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
		//		PoseList pose_list_;
		core::scoring::ScoreFunctionOP scorefxn_;
		protocols::simple_moves::GreenPackerOP green_packer_;
		bool use_green_packer_;
		bool use_packer_instead_of_rotamer_trials_;
		bool pack_at_neighbors_only_;
		bool rescore_only_;
		core::pack::task::PackerTaskOP pack_task_;
		core::pose::PoseOP pose_init_;

		std::string silent_file_;

		core::io::silent::SilentFileDataOP sfd_;

		utility::vector1< core::Size > calc_rms_res_;

		protocols::stepwise::enumerate::protein::sample_generators::StepWisePoseSampleGeneratorOP sample_generator_;

		protocols::stepwise::enumerate::protein::PoseFilterOP pose_filter_;

  };

} //protein
} //enumerate
} //stepwise
} //protocols

#endif
