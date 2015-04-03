// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/FeaturesReporter.cc
/// @brief  Base class to report geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/database/sql_utils.hh>

// Boost Headers
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

// External Headers
#include <cppdb/frontend.h>
#include <boost/lexical_cast.hpp>

// C++ Headers
#include <string>
#include <sstream>
#include <iostream>
#include <utility/vector0.hh>


namespace protocols {
namespace features {

/// @details Auto-generated virtual destructor
FeaturesReporter::~FeaturesReporter() {}

using basic::Tracer;
using basic::datacache::CacheableString;
using boost::char_separator;
using boost::tokenizer;
using core::Size;
using core::pose::Pose;
using cppdb::statement;
using cppdb::cppdb_error;
using protocols::filters::Filters_map;
using protocols::jd2::JobDistributor;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using std::endl;
using std::string;
using std::stringstream;
using utility::tag::TagCOP;
using basic::database::safely_prepare_statement;
using basic::database::safely_write_to_database;
using utility::trim;
using utility::vector1;
using utility::sql_database::sessionOP;

static Tracer TR("protocols.features.FeaturesReporter");

void
FeaturesReporter::write_schema_to_db(
	sessionOP db_session
) const {

  string schema_str(schema());
	try{
		basic::database::write_schema_to_database(schema_str,db_session);
	} catch (cppdb::cppdb_error error){
		stringstream err_msg;
		err_msg
			<< "Failed to write database schema for '" << type_name() << "'" << endl
			<< "Error Message:" << endl << error.what() << endl
			<< "Schema:" << endl << schema_str << endl;
		utility_exit_with_message(err_msg.str());
	}
}

/// @details Extract all features from the pose to the database.  The
///input parent_id can be used to connect the features with the rest
///of the database.  Usually the parent_id will be struct_id which
///uniquely identifies each structure in the database.  The feature
///reporter can return a parent_id.
Size
FeaturesReporter::report_features(
	Pose const & pose,
	StructureID parent_id,
	sessionOP db_session
){
	vector1<bool> relevant_residues(pose.total_residue(), true);
	return report_features(pose, relevant_residues, parent_id, db_session);
}

/// @details Extract all features from the pose to the database.  The
///input parent_id can be used to connect the features with the rest
///of the database.  Usually the parent_id will be struct_id which
///uniquely identifies each structure in the database.  The feature
///reporter can return a parent_id.
Size
FeaturesReporter::report_features(
	Pose const & /*pose*/,
	vector1< bool > const & /*relevant_residues*/,
	StructureID /*parent id*/,
	sessionOP /*db_session*/
){
	Size dummy_parent_id(0);
	return dummy_parent_id;
}


void
FeaturesReporter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	runtime_assert(tag->getName() == "feature");
}


string
FeaturesReporter::find_tag(
	Pose const & pose
) const {
	//silent files and pdbs set the name of the pose differently
	string name = "";
	if (pose.pdb_info()){
		name = pose.pdb_info()->name();
	}
	if (name == ""){
		if (pose.data().has(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG)) {
			name = static_cast< basic::datacache::CacheableString const & >
				(pose.data().get(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG)).str();
		} else {
			name = JobDistributor::get_instance()->current_output_name();
		}
	}
	return name;
}

void
FeaturesReporter::delete_records_from_table(
	string const & table_name,
	StructureID struct_id,
	sessionOP db_session){
	stringstream sql;
	sql << "DELETE FROM " << table_name << " WHERE struct_id = ?;";
	statement stmt(safely_prepare_statement(sql.str(), db_session));
	stmt.bind(1,struct_id);
	safely_write_to_database(stmt);
}

void
FeaturesReporter::set_relevant_residues_mode(
	RelevantResiduesMode::T setting
) {
	relevant_residues_mode_ = setting;
}

RelevantResiduesMode::T
FeaturesReporter::get_relevant_residues_mode() const {
	return relevant_residues_mode_;
}

bool
FeaturesReporter::check_relevant_residues(
	vector1<bool> const & relevant_residues,
	Size res1
) const {
	switch(relevant_residues_mode_){
		case RelevantResiduesMode::Exclusive:
			return relevant_residues[res1];
		case RelevantResiduesMode::Inclusive:
			return relevant_residues[res1];
		default:
			utility_exit_with_message("Unrecognized relevant_residues_mode: " +
					boost::lexical_cast<std::string>(relevant_residues_mode_));
	}
}

bool
FeaturesReporter::check_relevant_residues(
	vector1<bool> const & relevant_residues,
	Size res1,
	Size res2
) const {
	switch(relevant_residues_mode_){
		case RelevantResiduesMode::Exclusive:
			return relevant_residues[res1] && relevant_residues[res2];
		case RelevantResiduesMode::Inclusive:
			return relevant_residues[res1] || relevant_residues[res2];
		default:
			utility_exit_with_message("Unrecognized relevant_residues_mode: " +
					boost::lexical_cast<std::string>(relevant_residues_mode_));
	}
}

bool
FeaturesReporter::check_relevant_residues_range(
	vector1<bool> const & relevant_residues,
	Size begin,
	Size end
) const {
	switch(relevant_residues_mode_){
		case RelevantResiduesMode::Exclusive: {
			for(Size ii=begin; ii != end; ++ii){
				if (!relevant_residues[ii]) {
					return false;
				}
			}
			return true;
		}
		case RelevantResiduesMode::Inclusive: {
			for(Size ii=begin; ii != end; ++ii){
				if (relevant_residues[ii]) {
					return true;
				}
			}
			return false;
		}
		default:
			utility_exit_with_message("Unrecognized relevant_residues_mode: " +
					boost::lexical_cast<std::string>(relevant_residues_mode_));
	}
}


bool
FeaturesReporter::check_relevant_residues(
	vector1<bool> const & relevant_residues,
	vector1< Size > const & residues
) const {
	switch(relevant_residues_mode_){
		case RelevantResiduesMode::Exclusive: {
			for( vector1< Size >::const_iterator ii = residues.begin(), ii_end = residues.end(); ii != ii_end; ++ii){
				if (!relevant_residues[*ii]) {
					return false;
				}
			}
			return true;
		}
		case RelevantResiduesMode::Inclusive: {
			for( vector1< Size >::const_iterator ii = residues.begin(), ii_end = residues.end(); ii != ii_end; ++ii){
				if (relevant_residues[*ii]) {
					return true;
				}
			}
			return false;
		}
		default:
			utility_exit_with_message("Unrecognized relevant_residues_mode: " +
					boost::lexical_cast<std::string>(relevant_residues_mode_));
	}
}

} // namespace
} // namespace

