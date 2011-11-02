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

// Boost Headers
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#define foreach BOOST_FOREACH

// External Headers
#include <cppdb/frontend.h>
// AUTO-REMOVED #include <cppdb/errors.h>

// C++ Headers
#include <string>
#include <iostream>

#include <utility/vector0.hh>


namespace protocols {
namespace features {

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
using protocols::moves::DataMap;
using protocols::moves::Movers_map;
using std::endl;
using std::string;
using utility::tag::TagPtr;

using utility::trim;
using utility::vector1;
using utility::sql_database::sessionOP;

static Tracer TR("protocols.features.FeaturesReporter");

void
FeaturesReporter::write_schema_to_db(
	sessionOP db_session
) const {

  string schema_str(schema());

	char_separator< char > sep(";");
	tokenizer< char_separator< char > > tokens( schema_str, sep );
	foreach( string const & stmt_str, tokens){
		string trimmed_stmt_str(trim(stmt_str, " \n\t"));
		if(trimmed_stmt_str.size()){
			try{
				statement stmt = (*db_session) << trimmed_stmt_str + ";";
				stmt.exec();
			} catch (cppdb_error e) {
				TR.Error
					<< "ERROR reading schema for FeaturesReporter: " << type_name() << "\n"
					<< trimmed_stmt_str << endl;
				TR.Error << e.what() << endl;
				utility_exit();
			}
		}
	}
}

///@details Extract all features from the pose to the database.  The
///input parent_id can be used to connect the features with the rest
///of the database.  Usually the parent_id will be struct_id which
///uniquely identifies each structure in the database.  The feature
///reporter can return a parent_id.
Size
FeaturesReporter::report_features(
	Pose const & pose,
	Size parent_id,
	sessionOP db_session
){
	vector1<bool> relevant_residues(true, pose.total_residue());
	return report_features(pose, relevant_residues, parent_id,db_session);
}

///@details Extract all features from the pose to the database.  The
///input parent_id can be used to connect the features with the rest
///of the database.  Usually the parent_id will be struct_id which
///uniquely identifies each structure in the database.  The feature
///reporter can return a parent_id.
Size
FeaturesReporter::report_features(
	Pose const & /*pose*/,
	vector1< bool > const & /*relevant_residues*/,
	Size /*parent_id*/,
	sessionOP /*db_session*/
){
	Size dummy_parent_id(0);
	return dummy_parent_id;
}


void
FeaturesReporter::parse_my_tag(
	TagPtr const tag,
	DataMap & /*data*/,
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
			name = JobDistributor::get_instance()->current_job()->input_tag();
		}
	}
	return name;
}

} // namespace
} // namespace

