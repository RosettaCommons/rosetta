// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loops_definers/LoopsDatabaseDefiner.hh
/// @brief  A loops definer is creates a serialized loops list
/// @author Matthew O'Meara (mattjomear@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_LoopsDatabaseDefiner_HH
#define INCLUDED_protocols_loops_loops_definers_LoopsDatabaseDefiner_HH


// Unit Headers
#include <protocols/loops/loops_definers/LoopsDefiner.hh>
#include <protocols/loops/loops_definers/LoopsDatabaseDefiner.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>

// Platform Headers
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>


// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// C++ Headers
#include <string>


namespace protocols {
namespace loops {
namespace loops_definers {


///@ details given a database table with the following schema, defining a single loop per row
///
///CREATE TABLE loops (
///	tag TEXT,
///	start INTEGER,
///	stop INTEGER,
///	cut INTEGER,
///	skip_rate REAL,
///	extended BOOLEAN);
///
/// return all loops associated with the job distributor input tag
/// Note: you can specify a different table using the 'database_table' field
///
/// Note: if you would like to query the database for loops
/// differently, you can either pre-query the table and store it, or
/// extend, subclass, or create a different LoopsDefiner class.
class LoopsDatabaseDefiner : public LoopsDefiner {
public:

	LoopsDatabaseDefiner();

	virtual
	~LoopsDatabaseDefiner();

	LoopsDatabaseDefiner(
		LoopsDatabaseDefiner const & src);

	/// @brief Create another loops definer of the type matching the most-derived
	/// version of the class.
	virtual
	LoopsDefinerOP
	clone() const;


	/// @brief Used to parse an xml-like tag to load parameters and properties.
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap const & data,
		core::pose::Pose const &);

	virtual
	SerializedLoopList
	apply(
		core::pose::Pose const &);

private:
	utility::sql_database::sessionOP db_session_;
	std::string database_table_;

};

// do not add any derived classes to this file, unless they are
// generalized abstract base classes and do not actually 'do any work'

} //namespace
} //namespace
} //namespace

#endif // include guard


