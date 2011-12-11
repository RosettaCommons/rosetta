// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/LoopRlxMover.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#ifndef INCLUDED_protocols_antibody2_LoopRlxMover_hh
#define INCLUDED_protocols_antibody2_LoopRlxMover_hh

#include <protocols/antibody2/LoopRlxMover.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>  // Needs to be the full header so the scorefxn can default to NULL
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

namespace protocols {
namespace antibody2 {

	/// @brief Closes only one CDR onto a framework
	class LoopRlxMover : public protocols::moves::Mover {
	public:
		// default constructor
		// LoopRlxMover();
		// constructor with arguments
		LoopRlxMover( core::Size query_start, core::Size query_end  );

		void set_default();
		void setup_objects( core::pose::Pose & pose );

		virtual void apply( core::pose::Pose & pose_in );
		virtual std::string get_name() const;

		void setup_packer_task( core::pose::Pose & pose_in );

		/// @brief enable benchmark mode
		inline void enable_benchmark_mode( bool setting ) {
			benchmark_ = setting;
		}

	private:
		// Limits of query loop
		core::Size loop_start_;
		core::Size loop_end_;
		core::Size inner_cycles_, outer_cycles_;
		core::Real temperature_, gamma_;
		/// @brief benchmark flag
		bool benchmark_;

		protocols::moves::MonteCarloOP mc_;
		protocols::simple_moves::MinMoverOP loop_min_mover_;
		protocols::moves::SequenceMoverOP wiggle_loop_;

		protocols::loops::LoopOP one_loop_;

		// score functions
		core::scoring::ScoreFunctionOP highres_scorefxn_;
		core::kinematics::MoveMapOP movemap_;

		//packer task
		core::pack::task::TaskFactoryOP tf_;


	}; // class LoopRlxMover

} // antibody2
} // protocols




#endif
