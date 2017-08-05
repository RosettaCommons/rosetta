// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.cc
///
/// @brief

/// @author Tim Jacobs
// MPI headers

#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <core/types.hh>

#include <protocols/features/util.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/BatchFeatures.hh>
#include <basic/mpi/MessageListenerFactory.hh>
#include <basic/mpi/MessageListener.hh>
#include <basic/mpi/DbMoverMessageListener.hh>
#include <basic/mpi/util.hh>
#include <basic/mpi/mpi_enums.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/DbDataType.hh>

#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/string_util.hh>
#include <utility/mpi_util.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/Binary_Util.hh>

#include <cppdb/frontend.h>
#include <cppdb/errors.h>


#include <core/conformation/Residue.hh>

#include <string>
#include <sstream>
#include <utility>
#include <iostream>


namespace protocols {
namespace features {

using std::endl;
using std::string;
using std::stringstream;
using std::pair;
using std::map;
using cppdb::cppdb_error;
using utility::sql_database::sessionOP;
using basic::Tracer;
using core::Size;
using basic::mpi::request_data_from_head_node;
using basic::mpi::send_data_to_head_node;
using basic::mpi::DATABASE_PROTOCOL_AND_BATCH_ID_TAG;
using cppdb::statement;
using cppdb::result;

/// Macros are not properly caught and passed along by my #inclusion
/// cleanup script


// Static data for the serial case
Size static_protocol_id_ = 0;
bool protocol_table_initialized_ = false;
map<string, Size> static_batch_id_map_;

static THREAD_LOCAL Tracer TR("protocols.features.util");
// End static data

/// @brief write the given protocol and batch ids to the database. The protocol and batches
///features reporters will check for an existing entry with the same key, and write if one
///does not exist. Not recommended for parallel use as it is subject to race conditions (due
///to the nature of 'insert or ignore' type database writing)
void
set_protocol_and_batch_id(
	core::Size protocol_id,
	core::Size batch_id,
	string const & batch_name,
	string const & batch_description,
	utility::vector1<FeaturesReporterOP> features_reporters,
	sessionOP db_session
){
	ProtocolFeaturesOP protocol_features( new ProtocolFeatures() );
	BatchFeaturesOP batch_features( new BatchFeatures() );

	db_session->begin_transaction();
	protocol_id = protocol_features->report_features(protocol_id, db_session);
	db_session->force_commit_transaction();
	write_features_reporters_table(features_reporters, db_session);

	db_session->begin_transaction();
	batch_id = batch_features->report_features(batch_id, protocol_id, batch_name, batch_description, db_session);
	db_session->force_commit_transaction();
	write_batch_reports_table(features_reporters, batch_id, db_session);
}

/// @brief Get the protocol and batch ids or create them if they don't
///yet exist. For MPI protocols, only allow the head node to create
///protocol or batch ids and have the other nodes ask the head node
///for the info.
pair<Size, Size>
get_protocol_and_batch_id(
	string const & batch_name,
	string const & batch_description,
	utility::vector1<FeaturesReporterOP> features_reporters,
	sessionOP db_session
) {
	using namespace basic::mpi;

	int protocol_id = 0;
	int batch_id = 0;
	ProtocolFeaturesOP protocol_features( new ProtocolFeatures() );
	BatchFeaturesOP batch_features( new BatchFeatures() );

	int rank = 0;

#ifdef USEMPI
	MPI_Comm_rank( MPI_COMM_WORLD, (int*)(&rank) );
#endif

	//Send an batch_name to the head node, along with a
	//message_listening tag so that the messageListenerFactory of the
	//job distributor knows who to give the data to.

	string listener_data="";
	TR << "Requesting data on '" << batch_name << "' "
		<< "from the DATABASE_PROTOCOL_AND_BATCH_ID_TAG message listener on head node." << std::endl;

	//Some implementations of mpi don't allow self messaging. So, if we
	//are the head node, don't try to message yourself, just access the listener directly
	if ( rank != 0 ) {
		listener_data = request_data_from_head_node(DATABASE_PROTOCOL_AND_BATCH_ID_TAG, batch_name);
	} else {
		MessageListenerOP listener(
			MessageListenerFactory::get_instance()->get_listener(
			DATABASE_PROTOCOL_AND_BATCH_ID_TAG));
		listener->request(batch_name,listener_data);
	}

	//deserialize data into protocol_id and batch_id
	pair<Size, Size> ids = deserialize_db_listener_data(listener_data);
	protocol_id = ids.first;
	batch_id = ids.second;
	TR << "Received protocol_id='" << protocol_id << "' and batch_id='" << batch_id << "'" << std::endl;

	//only send ids if we recieve a protocol or batch id equal to 0.
	bool send_new_ids_to_head_node=false;

	//no protocol id set yet - have protocol features reporter generate one
	if ( protocol_id==0 ) {
		try {
			db_session->begin_transaction();
			protocol_id = protocol_features->report_features(protocol_id, db_session);
			db_session->force_commit_transaction();
			write_features_reporters_table(features_reporters, db_session);

			send_new_ids_to_head_node=true;
			runtime_assert(protocol_id!=0);
		} catch (cppdb_error const & error){
			stringstream err_msg;
			err_msg
				<< "Failed to set the protocol id for batch "
				<< "'" << batch_name << "'" << endl
				<< "Error Message:" << endl << error.what() << endl;
			utility_exit_with_message(err_msg.str());
		}
	}

	//no batch id set yet - have the batch features reporter generate one
	if ( batch_id==0 ) {
		try {
			db_session->begin_transaction();
			batch_id = batch_features->report_features(batch_id, protocol_id, batch_name, batch_description, db_session);
			db_session->force_commit_transaction();
			write_features_reporters_table(features_reporters, db_session);
			write_batch_reports_table(features_reporters, batch_id, db_session);

			send_new_ids_to_head_node=true;
			runtime_assert(batch_id!=0);
		} catch (cppdb_error const & error){
			stringstream err_msg;
			err_msg
				<< "Failed to set the batch id for batch '" << batch_name << "' "
				<< "with protocol_id '" << protocol_id << "'" << endl
				<< "Error Message:" << endl << error.what() << endl;
			utility_exit_with_message(err_msg.str());
		}
	}

	//If we generate any new ids, then send them to the head node.
	if ( send_new_ids_to_head_node ) {
		TR
			<< "protocol_id '" << protocol_id << "' "
			<< "and batch_id '" << batch_id << "', for batch '" << batch_name
			<< "' have been written to the database. Now sending info back to head node." <<std::endl;

		if ( rank != 0 ) {
			send_data_to_head_node(serialize_ids(protocol_id, batch_name, batch_id));
		} else {
			MessageListenerOP listener(
				MessageListenerFactory::get_instance()->get_listener(
				DATABASE_PROTOCOL_AND_BATCH_ID_TAG));
			listener->receive(serialize_ids(protocol_id, batch_name, batch_id));
		}
	}
	TR << "Done with protocol and batch id generation." <<std::endl;

	return pair<Size, Size>(protocol_id, batch_id);
}

/// @detail The 'features_reporters' table lists the type_names of the
/// all defined features reporters. The 'batch_reports' table link the
/// features reporters with each batch defined in the 'batches' table.
void
write_features_reporters_table(
	utility::vector1<FeaturesReporterOP> features_reporters,
	utility::sql_database::sessionOP db_session
){
	using namespace basic::database;
	using namespace basic::database::schema_generator;

	Schema features_reporters_schema(
		"features_reporters",
		PrimaryKey( Column("report_name", DbDataTypeOP( new DbTextKey() ))));

	features_reporters_schema.write(db_session);

	std::vector< std::string > column_names;
	column_names.emplace_back("report_name");
	std::vector< std::string > values;

	for ( FeaturesReporterOP const & reporter : features_reporters ) {
		string const report_name(reporter->type_name());

		values.clear();
		values.push_back( report_name );
		insert_or_ignore(
			"features_reporters", column_names,
			values, db_session
		);
	}

}

void
write_batch_reports_table(
	utility::vector1<FeaturesReporterOP> features_reporters,
	core::Size batch_id,
	utility::sql_database::sessionOP db_session
){
	using namespace basic::database;
	using namespace basic::database::schema_generator;

	Schema batch_reports("batch_reports");
	Column report_name("report_name", DbDataTypeOP( new DbTextKey() ));
	Column batch_id_col("batch_id", DbDataTypeOP( new DbInteger() ));

	batch_reports.add_foreign_key(
		ForeignKey(batch_id_col, "batches", "batch_id", true /*defer*/));
	batch_reports.add_foreign_key(
		ForeignKey(report_name, "features_reporters", "report_name", true /*defer*/));

	utility::vector1<Column> batch_reports_unique;
	batch_reports_unique.push_back(batch_id_col);
	batch_reports_unique.push_back(report_name);
	batch_reports.add_constraint( ConstraintOP( new UniqueConstraint(batch_reports_unique) ) );

	batch_reports.write(db_session);

	//Only report features/batch_id pairs that aren't already in the database
	string select_string =
		"SELECT *\n"
		"FROM\n"
		"\tbatch_reports\n"
		"WHERE\n"
		"\tbatch_id = ? AND\n"
		" report_name = ?;";
	statement select_stmt(safely_prepare_statement(select_string, db_session));

	string insert_string = "INSERT INTO batch_reports (batch_id, report_name) VALUES (?,?);";
	statement insert_stmt(safely_prepare_statement(insert_string, db_session));

	for ( FeaturesReporterOP const & reporter : features_reporters ) {
		string const report_name(reporter->type_name());
		select_stmt.bind(1,batch_id);
		select_stmt.bind(2,report_name);

		result res(safely_read_from_database(select_stmt));
		if ( !res.next() ) {
			insert_stmt.bind(1, batch_id);
			insert_stmt.bind(2, report_name);
			safely_write_to_database(insert_stmt);
		}
	}
}

pair<Size, Size>
deserialize_db_listener_data(
	string data
){
	utility::vector1< std::string > tokens = utility::split(data);
	if ( tokens.size() != 2 ) {
		utility_exit_with_message("failed to deserialize the message from master node. Message was: " + data + " You will get this message if trying to run ReportToDB mover in MPI mode with only on processor.");
	}
	int protocol_id=utility::string2int(tokens[1]);
	int batch_id=utility::string2int(tokens[2]);
	return pair<Size, Size>(protocol_id, batch_id);
}

string
serialize_ids(
	int protocol_id,
	string batch_name,
	Size batch_id
){
	return
		utility::to_string(protocol_id) + " " +
		batch_name + " " +
		utility::to_string(batch_id);
}

/// @detail look up the batch id given a struct id. Note this should
///only be used once the structure's table has been created, eg in an
///average features reporter's report_features function.
Size
get_batch_id(
	StructureID struct_id,
	sessionOP db_session
) {

	std::string const stmt_str(
		"SELECT batch_id FROM structures WHERE struct_id = ?;");
	statement stmt(basic::database::safely_prepare_statement(stmt_str, db_session));
	stmt.bind(1, struct_id);
	result res(basic::database::safely_read_from_database(stmt));

	if ( !res.next() ) {
		stringstream err_msg;
		err_msg
			<< "No batch_id found for struct_id '"
			<< struct_id << "'";
		utility_exit_with_message(err_msg.str());
	}
	Size batch_id;
	res >> batch_id;
	return batch_id;
}

utility::vector1<StructureID>
struct_ids_from_tag(
	sessionOP db_session,
	string const & tag
){
	std::string statement_string = "SELECT struct_id FROM structures WHERE tag=?;";
	statement stmt = basic::database::safely_prepare_statement(statement_string,db_session);
	stmt.bind(1,tag);
	result res = stmt.query();

	utility::vector1<StructureID> struct_ids;

	while ( res.next() ) {
		StructureID id;
		res >> id;
		struct_ids.push_back(id);
	}

	return struct_ids;
}

std::string serialize_residue_xyz_coords(core::conformation::Residue const & residue)
{
	//6bitencode and decode work best with arrays
	auto* coord_data = new core::Real[residue.natoms()*3];
	for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
		core::Size array_index = (atom_index - 1)*3;
		numeric::xyzVector<core::Real> xyz_coords(residue.xyz(atom_index));

		coord_data[array_index] = xyz_coords.x();
		coord_data[array_index + 1] = xyz_coords.y();
		coord_data[array_index + 2] = xyz_coords.z();
	}

	std::string output_data;
	core::Size memory_size = residue.natoms()*3*sizeof(core::Real);
	utility::encode6bit((unsigned char*)coord_data,memory_size,output_data);
	delete [] coord_data; //YOLO
	return output_data;

}

utility::vector1< numeric::xyzVector<core::Real> > deserialize_xyz_coords(std::string const & data, core::Size natoms)
{
	//natoms really needs to be correct
	auto* coord_data = new core::Real[natoms*3];
	//core::Size memory_size = natoms*3*sizeof(core::Real);  // unused ~Labonte
	utility::decode6bit((unsigned char*)coord_data,data);

	utility::vector1< numeric::xyzVector<core::Real> > xyz_vector;
	for ( core::Size atom_index = 1; atom_index <= natoms; ++atom_index ) {
		core::Size array_index = (atom_index - 1)*3;
		xyz_vector.push_back(numeric::xyzVector<core::Real>(coord_data[array_index],coord_data[array_index + 1],coord_data[array_index + 2]));
	}
	delete [] coord_data;
	return xyz_vector;
}

std::string
get_question_mark_string(core::Size const n){
	std::string result;
	if ( n==1 ) {
		result = "(?)";
		return result;
	} else {
		result = "(?";
		for ( core::Size i =2; i<=n; ++i ) {
			result += ",?";
		}
		result += ")";
		return result;
	}
}

} //namespace protocols
} //namespace features
