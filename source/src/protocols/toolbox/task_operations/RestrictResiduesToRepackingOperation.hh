// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictResiduesToRepackingOperation.hh
/// @brief  TaskOperation class that restricts a vector of Size defined residues to repacking
///   when parsed, it takes in a string and splits by ","
/// @author Eva-Maria Strauch (evas01@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictResiduesToRepackingOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictResiduesToRepackingOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictResiduesToRepackingOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RestrictResiduesToRepackingOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	RestrictResiduesToRepackingOperation();

	RestrictResiduesToRepackingOperation( utility::vector1 < core::Size > residues );

	utility::vector1< core::Size > get_residues() const;

	void set_residues( utility::vector1 < core::Size > residues_vec );

	virtual ~RestrictResiduesToRepackingOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( TagCOP, DataMap & );

	std::string reference_pdb_id() const{ return reference_pdb_id_; }
	void reference_pdb_id( std::string const s ){ reference_pdb_id_ = s; }

private:
	std::string unparsed_residues_;
	utility::vector1 < core::Size > residues_;
	std::string reference_pdb_id_;

};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictResiduesToRepackingOperation_HH
