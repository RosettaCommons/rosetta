// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/CloseOneCDRLoop.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_CloseOneCDRLoop_hh
#define INCLUDED_protocols_antibody_CloseOneCDRLoop_hh

#include <protocols/antibody/CloseOneCDRLoop.fwd.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>

#include <core/kinematics/MoveMap.fwd.hh>

namespace protocols {
namespace antibody {

/// @brief Closes only one CDR onto a framework
class CloseOneCDRLoop : public protocols::moves::Mover {
public:
	// default constructor
	CloseOneCDRLoop();

	// constructor with arguments
	CloseOneCDRLoop( core::Size query_start, core::Size query_end  );

	// default destructor
	~CloseOneCDRLoop();

	void set_default();
	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;


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

}; // class CloseOneCDRLoop


} // antibody
} // protocols


#endif
