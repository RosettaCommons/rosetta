// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictIdentitiesAtAlignedPositionsOperation.hh
/// @brief  TaskOperation class that restricts the identities at positions aligned with those of another pdb with the positions observed in that pdb
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictIdentitiesAtAlignedPositions_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictIdentitiesAtAlignedPositions_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictIdentitiesAtAlignedPositions.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

///@details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RestrictIdentitiesAtAlignedPositionsOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	RestrictIdentitiesAtAlignedPositionsOperation();

	virtual ~RestrictIdentitiesAtAlignedPositionsOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( TagPtr );

	utility::vector1< core::Size > res_ids() const{ return res_ids_; }
	void res_ids( utility::vector1< core::Size > const s ){ res_ids_ = s; }
	void source_pose( std::string const s );
	void chain( core::Size const c ){ chain_ = c; }
	core::Size chain() const{ return chain_; }
	void design_only_target_residues( bool const b ){ design_only_target_residues_ = b; }
	bool design_only_target_residues() const{ return design_only_target_residues_; }
private:
	core::pose::PoseOP source_pose_;
	utility::vector1< core::Size > res_ids_; // start and end will be parsed at apply time to determine the relevant residue numbers
	core::Size chain_; //dflt 1; chain on which to search for aligned residues
	bool design_only_target_residues_; //dflt false; if true, designs only the target residues to the identities seen in the source_pose and repacks a 6A shell around. If false, sets the target residues to design, and does not change the packer tasks for other residues
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictIdentitiesAtAlignedPositionsOperation_HH
