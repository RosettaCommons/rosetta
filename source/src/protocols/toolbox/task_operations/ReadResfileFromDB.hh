// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/ReadResfileFromDB.hh
/// @brief  read a refile indexed by the input structure tag from a supplied
///         relational database
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_toolbox_task_operations_ReadResfileFromDB_hh
#define INCLUDED_protocols_toolbox_task_operations_ReadResfileFromDB_hh

// Unit Headers
#include <protocols/toolbox/task_operations/ReadResfileFromDB.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

#include <utility/vector1.hh>
#include <string>


namespace protocols {
namespace toolbox {
namespace task_operations {

class ReadResfileFromDB : public core::pack::task::operation::TaskOperation{
public:
	typedef core::pack::task::operation::TaskOperation parent;

public:
	ReadResfileFromDB();

	ReadResfileFromDB(
		utility::sql_database::sessionOP db_session,
		std::string const & database_table);

	ReadResfileFromDB(ReadResfileFromDB const & src);

	virtual ~ReadResfileFromDB();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;

	void db_session(utility::sql_database::sessionOP db_session);

	void database_table(std::string const & database_table );
	std::string const & database_table() const;

	virtual void parse_tag(utility::tag::TagCOP, DataMap &);

private:
	std::string database_table_;
	utility::sql_database::sessionOP db_session_;
};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif // include guard
