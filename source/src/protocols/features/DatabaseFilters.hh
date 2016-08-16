// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/DatabaseFilters.hh
/// @brief  report atom-atom pair geometry and scores to features Statistics Scientific Benchmark
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_features_DatabaseFilters_hh
#define INCLUDED_protocols_features_DatabaseFilters_hh

// Unit Headers
#include <protocols/features/DatabaseFilters.fwd.hh>
#include <protocols/features/FeaturesReporter.fwd.hh>

//External

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

typedef std::pair<bool, utility::vector1<StructureID> > WriteDeletePair;

DatabaseFilterOP get_DB_filter_ptr();

class DatabaseFilter : public utility::pointer::ReferenceCount {
public:
	DatabaseFilter(){};
	virtual WriteDeletePair  operator()(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_session,
		core::Size const & protocol_id
	)=0;
	virtual ~DatabaseFilter(){};

private:
	DatabaseFilter( DatabaseFilter const & src );
};

class TopCountOfEachInput : public DatabaseFilter{
public:
	TopCountOfEachInput();
	TopCountOfEachInput(utility::vector1<std::string> arguments);
	WriteDeletePair  operator()(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_session,
		core::Size const & protocol_id
	);
	virtual ~TopCountOfEachInput(){};
	core::Size count_;
	std::string score_term_;
private:
	TopCountOfEachInput( TopCountOfEachInput const & src);
};

class TopCountOfAllInputs : public DatabaseFilter{
public:
	TopCountOfAllInputs();
	TopCountOfAllInputs(utility::vector1<std::string> arguments);
	WriteDeletePair  operator()(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_session,
		core::Size const & protocol_id
	);
	virtual ~TopCountOfAllInputs(){};
	core::Size count_;
	std::string score_term_;
private:
	TopCountOfAllInputs( TopCountOfAllInputs const & src);
};

class TopPercentOfEachInput : public DatabaseFilter{
public:
	TopPercentOfEachInput(utility::vector1<std::string> arguments);
	WriteDeletePair  operator()(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_session,
		core::Size const & protocol_id
	);
	virtual ~TopPercentOfEachInput(){};
private:
	TopPercentOfEachInput( TopPercentOfEachInput const & src);
	TopCountOfEachInput top_count_of_each_input_;
};

class TopPercentOfAllInputs: public DatabaseFilter{
public:
	TopPercentOfAllInputs(utility::vector1<std::string> arguments);
	WriteDeletePair  operator()(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_session,
		core::Size const & protocol_id
	);
	virtual ~TopPercentOfAllInputs(){};
private:
	TopPercentOfAllInputs( TopPercentOfAllInputs const & src);
	TopCountOfAllInputs top_count_of_all_inputs_;
};

} // namespace
} // namespace

#endif // include guard
