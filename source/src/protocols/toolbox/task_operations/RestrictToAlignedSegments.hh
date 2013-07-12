// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToAlignedSegmentsOperation.hh
/// @author Sarel Fleishman sarel@weizmann.ac.il

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToAlignedSegments_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictToAlignedSegments_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToAlignedSegments.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>
#include <set>

#include <utility/vector1.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

///@details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RestrictToAlignedSegmentsOperation : public RestrictOperationsBase
{
public:
	typedef protocols::toolbox::task_operations::RestrictOperationsBase parent;

	RestrictToAlignedSegmentsOperation();

	virtual ~RestrictToAlignedSegmentsOperation();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( TagPtr );

	utility::vector1< core::Size > start_res() const{ return start_res_; }
	void start_res( utility::vector1< core::Size > const s ){ start_res_ = s; }
	utility::vector1< core::Size > stop_res() const{ return stop_res_; }
	void stop_res( utility::vector1< core::Size > const s ){ stop_res_ = s; }
	core::Size chain()const {return chain_;}
	void chain( core::Size const c ){ chain_ = c; }

	core::Real repack_shell() const{ return repack_shell_;}
	void repack_shell( core::Real const r ){ repack_shell_ = r; }
private:
	utility::vector1< core::pose::PoseOP > source_pose_;
	utility::vector1< core::Size > start_res_; // start and end will be parsed at apply time to determine the relevant residue numbers
	utility::vector1< core::Size > stop_res_;
	core::Size chain_; //dflt 1; which chain to look at. 0 means all chains
	core::Real repack_shell_; //dflt 6A; allow repack in residues surrounding the aligned segment on chain_
};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif //INCLUDED_devel_splice_RestrictToAlignedSegments_hh
