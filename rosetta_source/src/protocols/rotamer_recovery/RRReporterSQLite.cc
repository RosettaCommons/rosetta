// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRReporterSQLite.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RRReporterSQLite.hh>

// Project Headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/types.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

//C++ Headers
#include <ostream>
#include <string>

namespace protocols {
namespace rotamer_recovery {

using std::endl;
using std::ostream;
using std::string;
using basic::datacache::CacheableString;
using basic::Tracer;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pose::Pose;
using protocols::jd2::JobDistributor;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::sessionOP;
using cppdb::statement;

static Tracer TR("protocols.rotamer_recovery.rRReporterSQLite");

RRReporterSQLite::RRReporterSQLite() :
	output_level_( OutputLevel::full ),
	struct_id1_( 0 ),
	struct_id2_( 0 ),
	comparer_name_(),
	comparer_params_(),
	residues_considered_( 0 ),
	rotamers_recovered_( 0 ),
	db_session_()
{
	db_session_ =
			DatabaseSessionManager::get_instance()->get_session("rotamer_recovery.db3");
}

RRReporterSQLite::RRReporterSQLite(
	string const & database_fname,
	OutputLevel::e output_level /* = OutputLevel::full */
) :
	output_level_( output_level ),
	struct_id1_( 0 ),
	struct_id2_( 0 ),
	comparer_name_(),
	comparer_params_(),
	residues_considered_( 0 ),
	rotamers_recovered_( 0 ),
	db_session_()
{
	db_session_ =
		DatabaseSessionManager::get_instance()->get_session(database_fname);

	statement stmt = (*db_session_) << schema(output_level);
	stmt.exec();

}

RRReporterSQLite::RRReporterSQLite(
	sessionOP db_session,
	OutputLevel::e const output_level /* = OutputLevel::full */
) :
	output_level_( output_level ),
	struct_id1_( 0 ),
	struct_id2_( 0 ),
	comparer_name_(),
	comparer_params_(),
	residues_considered_( 0 ),
	rotamers_recovered_( 0 ),
	db_session_( db_session )
{}



RRReporterSQLite::RRReporterSQLite( RRReporterSQLite const & src ) :
	RRReporter(),
	output_level_( src.output_level_ ),
	struct_id1_( src.struct_id1_ ),
	struct_id2_( src.struct_id2_ ),
	comparer_name_( src.comparer_name_ ),
	comparer_params_( src.comparer_params_ ),
	residues_considered_( src.residues_considered_ ),
	rotamers_recovered_( src.rotamers_recovered_ ),
	db_session_( src.db_session_ )
{}

RRReporterSQLite::~RRReporterSQLite() {}

string
RRReporterSQLite::schema(
	RRReporterSQLite::OutputLevel::e output_level /* = OutputLevel::full */
) {

	switch( output_level ){

	case OutputLevel::full:
		return
			"CREATE TABLE IF NOT EXISTS rotamer_recovery (\n"
			"	struct1_name TEXT,\n"
			"	aa TEXT,\n"
			"	residue_type TEXT,\n"
			"	chain1 TEXT,\n"
			"	res1 INTEGER,\n"
			"	struct2_name TEXT,\n"
			"	chain2 TEXT,\n"
			"	res2 INTEGER,\n"
			"	comparer_name TEXT,\n"
			"	comparer_params TEXT,\n"
			"	score REAL,\n"
			"	recovered BOOLEAN,"
			"	PRIMARY KEY (struct1_name, chain1, res1, struct2_name, chain2, res2));\n";
		break;

	case OutputLevel::features:
		return
			"CREATE TABLE IF NOT EXISTS rotamer_recovery (\n"
			"	struct_id INTEGER,\n"
			"	resNum INTEGER,\n"
			"	divergence REAL,\n"
			"	PRIMARY KEY(struct_id, resNum));\n";
		break;

	case OutputLevel::none:
		return "";

	default:
		utility_exit_with_message("Unrecognized Output Level.");
	}
	return "";
}

void
RRReporterSQLite::set_output_level(
	OutputLevel::e const output_level
){
	output_level_ = output_level;
}

RRReporterSQLite::OutputLevel::e
RRReporterSQLite::get_output_level(
) const {
	return output_level_;
}

void
RRReporterSQLite::set_struct_id1(
	Size const struct_id1
){
	struct_id1_ = struct_id1;
}

Size
RRReporterSQLite::get_struct_id1(
) const {
	return struct_id1_;
}

void
RRReporterSQLite::set_struct_id2(
	Size const struct_id2
){
	struct_id1_ = struct_id2;
}

Size
RRReporterSQLite::get_struct_id2(
) const {
	return struct_id2_;
}

void
RRReporterSQLite::set_comparer_info(
	string const & comparer_name,
	string const & comparer_params
)	{
	comparer_name_ = comparer_name;
	comparer_params_ = comparer_params;
}

void
RRReporterSQLite::reset_recovery(){
	residues_considered_=0;
	rotamers_recovered_=0;
}

void
RRReporterSQLite::report_rotamer_recovery(
	Pose const & pose1,
	Pose const & pose2,
	Residue const & res1,
	Residue const & res2,
	Real const score,
	bool recovered
){
	switch (output_level_) {
	case OutputLevel::full:
		report_rotamer_recovery_full(
			pose1, pose2, res1, res2, score, recovered );
		break;
	case OutputLevel::features:
		report_rotamer_recovery_features(struct_id1_, res1, score);
		break;
	case OutputLevel::none:
		break;
	default:
		utility_exit_with_message( "Unknown RRReporterSQLite output level." );
	}

	residues_considered_++;
	rotamers_recovered_ += recovered;
}

void
RRReporterSQLite::report_rotamer_recovery_full(
	Pose const & pose1,
	Pose const & pose2,
	Residue const & res1,
	Residue const & res2,
	Real const score,
	bool recovered
){
	// Argh why isn't there a good way go get the name of a pose?
	//silent files and pdbs set the name of the pose differently
	string struct1_name = "No_Name_Found";
	if (pose1.pdb_info() && ( pose1.pdb_info()->name() != "" ) ){
		struct1_name = pose1.pdb_info()->name();
	} else if ( pose1.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		struct1_name = static_cast< CacheableString const & >
			( pose1.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	} else {
		struct1_name = JobDistributor::get_instance()->current_job()->input_tag();
	}

	string struct2_name = "No_Name_Found";
	if (pose2.pdb_info() && ( pose1.pdb_info()->name() != "" ) ){
		struct2_name = pose2.pdb_info()->name();
	} else if ( pose2.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		struct2_name = static_cast< CacheableString const & >
			( pose2.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	} else {
		struct2_name = JobDistributor::get_instance()->current_job()->input_tag();
	}

	while(true)
	{
		try
		{
			statement stmt = (*db_session_)
				<< "INSERT INTO rotamer_recovery VALUES (?,?,?,?,?,?,?,?,?,?,?,?);"
				<< struct1_name
				<< res1.name1()
				<< res1.type().name()
				<< res1.chain()
				<< res1.seqpos()
				<< struct2_name
				<< res2.chain()
				<< res2.seqpos()
				<< comparer_name_
				<< comparer_params_
				<< score
				<< recovered;
			stmt.exec();
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}
}

void
RRReporterSQLite::report_rotamer_recovery_features(
	Size const struct_id1,
	Residue const & res1,
	Real const score
){


	while(true)
	{
		try
		{
			statement stmt = (*db_session_)
				<< "INSERT INTO rotamer_recovery VALUES (?,?,?);"
				<< struct_id1
				<< res1.seqpos()
				<< score;
			stmt.exec();
			break;
		}catch(cppdb::cppdb_error &)
		{
			usleep(10);
			continue;
		}
	}
}

void
RRReporterSQLite::show( ostream & out ) const {

	out
		<< "Recovered " << rotamers_recovered_
		<< " at " << residues_considered_
		<< " for a recovery rate of " << recovery_rate() << "." << endl;
}

void
RRReporterSQLite::show( ) const {
	show( TR );
}

Real
RRReporterSQLite::recovery_rate() const {
	return Real( rotamers_recovered_ / residues_considered_ );
}

} // namespace
} // namespace
