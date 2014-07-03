// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/AlignedThreadOperation.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_AlignedThreadOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_AlignedThreadOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/AlignedThreadOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <core/types.hh>

// Utility Headers

// C++ Headers
#include <string>

namespace protocols {
namespace toolbox {
namespace task_operations {

///@details this class is a TaskOperation to prevent repacking of residues not near an interface.
class AlignedThreadOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	AlignedThreadOperation();
	AlignedThreadOperation( std::string const seq );

	virtual ~AlignedThreadOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( TagCOP, DataMap & );

	std::string alignment_file() const{ return alignment_file_; }
	void alignment_file( std::string const s ){ alignment_file_ = s; }

	std::string query_name() const{ return query_name_; }
	void query_name( std::string const s ){ query_name_ = s; }

	void start_res( core::Size const s ){ start_res_ = s ; }
	core::Size start_res() const{ return start_res_;}
private:
	std::string alignment_file_; //dflt ""
	std::string query_name_; //dflt "" ; the name of the sequence to thread on -s structure
	core::Size start_res_; //dflt 1; where does the alignment start
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_AlignedThreadOperation_HH
