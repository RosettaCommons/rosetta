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

//Basic Headers
// #include <basic/database/sql_utils.hh> 09_11_2013, commented by Doonam due to double declaration
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>


// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

//C++ Headers
#include <boost/assign/list_of.hpp>
#include <ostream>
#include <string>
#include <vector>
#include <utility/vector1.hh>


namespace protocols {
namespace rotamer_recovery {

using std::endl;
using std::ostream;
using std::string;
using basic::datacache::CacheableString;
using basic::database::safely_prepare_statement;
using basic::database::safely_write_to_database;
using basic::database::write_schema_to_database;
using basic::database::get_db_session;
using basic::Tracer;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pose::Pose;
using protocols::jd2::JobDistributor;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using protocols::features::StructureID;

static Tracer TR("protocols.rotamer_recovery.RRReporterSQLite");

RRReporterSQLite::RRReporterSQLite() :
	output_level_( OutputLevel::full ),
//	struct_id1_(""),
//	struct_id2_(""),
	protocol_name_(),
	protocol_params_(),
	comparer_name_(),
	comparer_params_(),
	residues_considered_( 0 ),
	rotamers_recovered_( 0 ),
	database_name_("rotamer_recovery.db3"),
	database_pq_schema_(""),
	db_session_()
{}

RRReporterSQLite::RRReporterSQLite(
	string const & database_name,
	string const & database_pq_schema /* = "" */,
	OutputLevel::e output_level /* = OutputLevel::full */
) :
	output_level_( output_level ),
//	struct_id1_(""),
//	struct_id2_(""),
	protocol_name_(),
	protocol_params_(),
	comparer_name_(),
	comparer_params_(),
	residues_considered_( 0 ),
	rotamers_recovered_( 0 ),
	database_name_(database_name),
	database_pq_schema_(database_pq_schema),
	db_session_()
{}

RRReporterSQLite::RRReporterSQLite(
	sessionOP db_session,
	OutputLevel::e const output_level /* = OutputLevel::full */
) :
	output_level_( output_level ),
//	struct_id1_(""),
//	struct_id2_(""),
	protocol_name_(),
	protocol_params_(),
	comparer_name_(),
	comparer_params_(),
	residues_considered_( 0 ),
	rotamers_recovered_( 0 ),
 	database_name_(),
 	database_pq_schema_(),
	db_session_( db_session )
{}



RRReporterSQLite::RRReporterSQLite( RRReporterSQLite const & src ) :
	RRReporter(),
	output_level_( src.output_level_ ),
	struct_id1_( src.struct_id1_ ),
	struct_id2_( src.struct_id2_ ),
	protocol_name_( src.protocol_name_ ),
	protocol_params_( src.protocol_params_ ),
	comparer_name_( src.comparer_name_ ),
	comparer_params_( src.comparer_params_ ),
	residues_considered_( src.residues_considered_ ),
	rotamers_recovered_( src.rotamers_recovered_ ),
	database_name_( src.database_name_ ),
	database_pq_schema_( src.database_pq_schema_ ),
	db_session_( src.db_session_ )
{}

RRReporterSQLite::~RRReporterSQLite() {}

void
RRReporterSQLite::write_schema_to_db(
	sessionOP db_session,
	RRReporterSQLite::OutputLevel::e output_level /* = OutputLevel::ful */
) const {
	switch(output_level){
	case OutputLevel::full:
		write_nchi_table_schema(db_session);
		write_rotamer_recovery_full_table_schema(db_session);
		break;
	case OutputLevel::features:
		write_rotamer_recovery_features_table_schema(db_session);
		break;
	case OutputLevel::none:
		break;
	default:
		utility_exit_with_message("Unrecognized Output Level.");
	}
}

void
RRReporterSQLite::write_nchi_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;
	using namespace basic::database;
	using namespace boost::assign;

	Column name3("name3", new DbText());
	Column nchi("nchi", new DbInteger());

	Columns primary_key_columns;
	primary_key_columns.push_back(name3);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("nchi", primary_key);
	table.add_column(nchi);
	table.write(db_session);


	// insert values
	string table_name("nchi");
	std::vector<string> column_names;
	column_names.push_back("name3");
	column_names.push_back("nchi");
	insert_or_ignore(table_name, column_names, list_of("'ARG'")("4"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'LYS'")("4"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'MET'")("3"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'GLN'")("3"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'GLU'")("3"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'TYR'")("2"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'ILE'")("2"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'ASP'")("2"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'TRP'")("2"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'PHE'")("2"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'HIS'")("2"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'ASN'")("2"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'THR'")("1"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'SER'")("1"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'PRO'")("1"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'CYS'")("1"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'VAL'")("1"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'LEU'")("1"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'ALA'")("0"), db_session);
	insert_or_ignore(table_name, column_names, list_of("'GLY'")("0"), db_session);
}

void
RRReporterSQLite::write_rotamer_recovery_full_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct1_name("struct1_name", new DbText());
	Column name1("name1", new DbText());
	Column name3("name3", new DbText());
	Column residue_type("residue_type", new DbText());
	Column chain1("chain1", new DbText());
	Column res1("res1", new DbInteger());
	Column struct2_name("struct_name2", new DbText());
	Column chain2("chain2", new DbText());
	Column res2("res2", new DbInteger());
	Column protocol_name("protocol_name", new DbText());
	Column protocol_params("protocol_params", new DbText());
	Column comparer_name("comparer_name", new DbText());
	Column comparer_params("comparer_params", new DbText());
	Column score("score", new DbReal());
	Column recovered("recovered", new DbInteger());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct1_name);
	primary_key_columns.push_back(chain1);
	primary_key_columns.push_back(res1);
	primary_key_columns.push_back(struct2_name);
	primary_key_columns.push_back(chain2);
	primary_key_columns.push_back(res2);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("rotamer_recovery", primary_key);
	table.add_column(name1);
	table.add_column(name3);
	table.add_column(residue_type);
	table.add_column(protocol_name);
	table.add_column(protocol_params);
	table.add_column(comparer_name);
	table.add_column(comparer_params);
	table.add_column(score);
	table.add_column(recovered);

	table.write(db_session);
}


void
RRReporterSQLite::write_rotamer_recovery_features_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column resNum("resNum", new DbInteger());
	Column divergence("divergence", new DbReal());
	Column recovered("recovered", new DbInteger());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(resNum);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("rotamer_recovery", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(divergence);
	table.add_column(recovered);

	table.write(db_session);
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
	StructureID const struct_id1
){
	struct_id1_ = struct_id1;
}

StructureID
RRReporterSQLite::get_struct_id1(
) const {
	return struct_id1_;
}

void
RRReporterSQLite::set_struct_id2(
	StructureID const struct_id2
){
	struct_id1_ = struct_id2;
}

StructureID
RRReporterSQLite::get_struct_id2(
) const {
	return struct_id2_;
}

void
RRReporterSQLite::set_protocol_info(
	string const & protocol_name,
	string const & protocol_params
)	{
	protocol_name_ = protocol_name;
	protocol_params_ = protocol_params;
}

void
RRReporterSQLite::set_comparer_info(
	string const & comparer_name,
	string const & comparer_params
)	{
	comparer_name_ = comparer_name;
	comparer_params_ = comparer_params;
}

sessionOP
RRReporterSQLite::db_session(){
	if(!db_session_){
		db_session_ =
			basic::database::get_db_session(
				database_name_, database_pq_schema_);

		write_schema_to_db(db_session_, get_output_level());
	}
	return db_session_;
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
		report_rotamer_recovery_features(struct_id1_, res1, score, recovered);
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

	std::string statement_string = "INSERT INTO rotamer_recovery VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement stmt(safely_prepare_statement(statement_string,db_session()));

	stmt.bind(1,struct1_name);
	stmt.bind(2,string(1,res1.name1()));
	stmt.bind(3,res1.name3());
	stmt.bind(4,res1.type().name());
	stmt.bind(5,res1.chain());
	stmt.bind(6,res1.seqpos());
	stmt.bind(7,struct2_name);
	stmt.bind(8,res2.chain());
	stmt.bind(9,res2.seqpos());
	stmt.bind(10,protocol_name_);
	stmt.bind(11,protocol_params_);
	stmt.bind(12,comparer_name_);
	stmt.bind(13,comparer_params_);
	stmt.bind(14,score);
	stmt.bind(15,recovered);
	safely_write_to_database(stmt);

}

void
RRReporterSQLite::report_rotamer_recovery_features(
	StructureID const struct_id1,
	Residue const & res1,
	Real const score,
	bool const recovered
){


	std::string statement_string = "INSERT INTO rotamer_recovery VALUES (?,?,?,?);";
	statement stmt(safely_prepare_statement(statement_string, db_session_));

	stmt.bind(1,struct_id1);
	stmt.bind(2, res1.seqpos());
	stmt.bind(3, score);
	stmt.bind(4, recovered);
	safely_write_to_database(stmt);

}

void
RRReporterSQLite::show( ostream & out ) const {

	out
		<< "Recovered " << rotamers_recovered_
		<< " at " << residues_considered_ << " residues considered"
		<< " for a recovery rate of " << recovery_rate() << "." << endl;
}

void
RRReporterSQLite::show( ) const {
	show( TR );
}

Real
RRReporterSQLite::recovery_rate() const {
	return Real(rotamers_recovered_) / Real(residues_considered_);
}

} // namespace
} // namespace
