// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/BatchFeatures.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_BatchFeatures_HH
#define INCLUDED_protocols_features_BatchFeatures_HH

// Unit Headers
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/BatchFeatures.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>
#include <utility/exit.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>

namespace protocols {
namespace features {

class BatchFeatures : public utility::pointer::ReferenceCount {
public:
	BatchFeatures();

	BatchFeatures( BatchFeatures const & src );

	virtual ~BatchFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session, core::Size batch_id) const;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	/// @brief return sql statements that setup helpful indices on the tables
	std::string
	indices() const;

	/// @brief return sql statments that add indexes check invariants etc.
	//std::string cleanup() const;

	core::Size
	report_features(
		core::Size batch_id,
		core::Size protocol_id,
		std::string name,
		std::string description,
		utility::sql_database::sessionOP db_session);

};

} // features namespace
} // protocols namespace

#endif
