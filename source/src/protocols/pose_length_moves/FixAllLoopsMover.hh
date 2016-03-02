// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/pose_length_moves/FixAllLoopsMover.hh
/// @brief

#ifndef INCLUDED_protocols_pose_length_moves_FixAllLoopsMover_hh
#define INCLUDED_protocols_pose_length_moves_FixAllLoopsMover_hh


#include <protocols/moves/Mover.hh>

#include <protocols/pose_length_moves/FixAllLoopsMover.fwd.hh>

#include <core/indexed_structure_store/ABEGOHashedFragmentStore.fwd.hh>
#include <core/pose/Pose.fwd.hh>

//#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

// C++ Headers
#include <string>
#include <map>
// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ctime>


namespace protocols {
namespace pose_length_moves {
using namespace core;
using namespace std;
using utility::vector1;



class FixAllLoopsMover : public protocols::moves::Mover {
public:
	FixAllLoopsMover();
	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	moves::MoverOP clone() const { return moves::MoverOP( new FixAllLoopsMover( *this ) ); }
	virtual void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	protocols::loops::Loops get_loops(core::pose::Pose const & pose);
	int resAdjustmentStartLow_;
	int resAdjustmentStartHigh_;
	int resAdjustmentStopLow_;
	int resAdjustmentStopHigh_;
	int resAdjustmentStartLow_sheet_;
	int resAdjustmentStartHigh_sheet_;
	int resAdjustmentStopLow_sheet_;
	int resAdjustmentStopHigh_sheet_;
	Size loopLengthRangeLow_;
	Size loopLengthRangeHigh_;
	Real rmsThreshold_;
	Real max_vdw_change_;
	bool idealExtension_;
	bool ideal_;
	bool reject_failed_loops_;
	Size firstResidue_;
	Size lastResidue_;
	//time_t start_time_;

	core::indexed_structure_store::ABEGOHashedFragmentStore * ABEGOHashedFragmentStore_;
};


} // pose_length_moves
} // protocols

#endif
