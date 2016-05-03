// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/DatabaseThread.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Assaf Alon assafalon@gmail.com

#ifndef INCLUDED_protocols_toolbox_task_operations_DatabaseThread_hh
#define INCLUDED_protocols_toolbox_task_operations_DatabaseThread_hh

// Unit Headers
#include <protocols/toolbox/task_operations/DatabaseThread.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to thread sequences from a database.
class DatabaseThread : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;
	DatabaseThread();
	virtual ~DatabaseThread();
	virtual TaskOperationOP clone() const;
	virtual void apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	virtual void parse_tag( TagCOP, DataMap & );

	//function declerations:
	core::Size find_length( core::pose::Pose const & pose ) const;
	std::string pick_sequence_from_database( core::pose::Pose const & pose ) const;
	void mark_designable(std::string & sequence, core::pose::Pose const &) const;
	void mark_leave_as_is(std::string & sequence, core::pose::Pose const &) const;

	//setters and getters for functions in parse my tag:
	std::string template_file() const{ return template_file_; } //getter for the template file name
	void template_file( std::string const s ){ template_file_ = s; } // setter for the template file name
	core::Size start_res() const {return start_res_;};
	void start_res( core::Size const s ) { start_res_ = s; };
	core::Size end_res() const { return end_res_; }; //getter
	void end_res( core::Size const s ) { end_res_ = s; }; //setter
	bool allow_design_around() const{ return allow_design_around_;}
	void allow_design_around( bool const b ){ allow_design_around_ = b ; }
	std::string database_fname() const{ return database_fname_; };//getter
	void database_fname( std::string const d ){ database_fname_ = d; };//setter
	utility::vector1<core::Size> designable() const {return designable_; };//getter
	void designable(utility::vector1<core::Size> const vector) {designable_=vector; };//setter
	utility::vector1<core::Size> leave_as_is() const {return leave_as_is_; };//getter
	void leave_as_is(utility::vector1<core::Size> const vector) {leave_as_is_=vector; };//setter
	void target_sequence( std::string const seq ) {target_sequence_ = seq; };//setter
	std::string target_sequence() const {return target_sequence_; };//getter

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "DatabaseThread"; }

private:
	std::string target_sequence_;
	std::string template_file_;
	std::string database_fname_;
	core::pose::PoseOP template_pose_;
	core::Size start_res_;
	core::Size end_res_;
	bool allow_design_around_; //dflt true; if false restricts rest of the pose to repakcing
	utility::vector1<core::Size> design_;
	utility::vector1<core::Size> revert_to_template_;
	utility::vector1<std::string> full_database_;
	utility::vector1<core::Size> designable_,leave_as_is_;
};
} //namespace protocols
} //namespace toolbox
} //namespace task_operations
#endif // INCLUDED_protocols_toolbox_TaskOperations_DatabaseThread_HH
