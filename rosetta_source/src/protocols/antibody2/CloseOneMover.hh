// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/CloseOneMover.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody2_CloseOneMover_hh
#define INCLUDED_protocols_antibody2_CloseOneMover_hh

#include <protocols/antibody2/CloseOneMover.fwd.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/moves/PyMolMover.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

namespace protocols {
namespace antibody2 {

	/// @brief Closes only one CDR onto a framework
	class CloseOneMover : public protocols::moves::Mover {
	public:
		// default constructor
		CloseOneMover();

		// constructor with arguments
		CloseOneMover( core::Size query_start, core::Size query_end  );

		// default destructor
		~CloseOneMover();

		void set_default();
		virtual void apply( core::pose::Pose & pose_in );
		virtual std::string get_name() const;

		void set_pymol( protocols::moves::PyMolMoverOP pymol ) { pymol_ = pymol; }

		/// @brief enable benchmark mode
		inline void enable_benchmark_mode( bool setting ) {
			benchmark_ = setting;
		}

	private:
		// Limits of query loop
		core::Size loop_start_, cdr_loop_start_;
		core::Size loop_end_, cdr_loop_end_;
		core::Real allowed_separation_;
		core::Size flanking_residues_;

		/// @brief benchmark flag
		bool benchmark_;

		core::kinematics::MoveMapOP movemap_;
		protocols::moves::PyMolMoverOP pymol_;

	}; // class CloseOneMover


} // antibody2
} // protocols


#endif
