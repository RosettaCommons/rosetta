// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/BatchFeatures.hh>
#include <protocols/jd2/util.hh>
#include <basic/message_listening/MessageListenerFactory.hh>
#include <basic/message_listening/MessageListener.hh>
#include <basic/message_listening/DbMoverMessageListener.hh>
#include <basic/message_listening/util.hh>

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

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <core/conformation/Residue.hh>

#include <string>
#include <sstream>
#include <utility>
#include <iostream>


namespace protocols{
namespace features{

using std::endl;
using std::string;
using std::stringstream;
using std::pair;
using std::map;
using cppdb::cppdb_error;
using utility::sql_database::sessionOP;
using basic::Tracer;
using core::Size;
using basic::message_listening::request_data_from_head_node;
using basic::message_listening::send_data_to_head_node;
using basic::message_listening::DATABASE_PROTOCOL_AND_BATCH_ID_TAG;
using cppdb::statement;
using cppdb::result;


// Static data for the serial case
Size static_protocol_id_ = 0;
bool protocol_table_initialized_ = false;
map<string, Size> static_batch_id_map_;

static Tracer TR("protocols.features.util");
// End static data


///@brief Get the protocol and batch ids or create them if they don't
///yet exist. For MPI protocols, only allow the head node to create
///protocol or batch ids and have the other nodes ask the head node
///for the info.
pair<Size, Size>
get_protocol_and_batch_id(
	string identifier,
	sessionOP db_session
) {
	using namespace basic::message_listening;

	int protocol_id = 0;
	int batch_id = 0;
	ProtocolFeaturesOP protocol_features = new ProtocolFeatures();
	BatchFeaturesOP batch_features = new BatchFeatures();

#ifdef USEMPI

	int rank = 0;
	MPI_Comm_rank( MPI_COMM_WORLD, (int*)(&rank) );

	//Send an identifier to the head node, along with a
	//message_listening tag so that the messageListenerFactory of the
	//job distributor knows who to give the data to.

	//Some implementations of mpi don't allow self messaging. So, if we
	//are the head node, don't try to message yourself, just access the listener directly

	string listener_data="";
	//Set the max_batch_id if necessary.
	if(rank == 0)
	{
		DbMoverMessageListenerOP listener(utility::pointer::dynamic_pointer_cast<DbMoverMessageListener,MessageListener>(MessageListenerFactory::get_instance()->get_listener(DATABASE_PROTOCOL_AND_BATCH_ID_TAG)));
		if(!listener->max_batch_id_set())
		{
			std::string select_max = "SELECT MAX(batch_id) FROM batches;";
			cppdb::statement stmt(basic::database::safely_prepare_statement(select_max,db_session));
			cppdb::result res(basic::database::safely_read_from_database(stmt));

			Size max_batch_id(0);
			while(res.next())
			{
				res >> max_batch_id;
			}
			listener->set_max_batch_id(max_batch_id);

		}
	}

	if(rank != 0) {
		listener_data = request_data_from_head_node(DATABASE_PROTOCOL_AND_BATCH_ID_TAG, identifier);
		TR
			<< "Requesting data on '" << identifier << "' "
			<< "from the DATABASE_PROTOCOL_AND_BATCH_ID_TAG message listener on head node." << std::endl;

	} else {
		MessageListenerOP listener(
			MessageListenerFactory::get_instance()->get_listener(
				DATABASE_PROTOCOL_AND_BATCH_ID_TAG));
		listener->request(identifier,listener_data);

		TR
			<< "Requesting data on '" << identifier << "' "
			<< "from the DATABASE_PROTOCOL_AND_BATCH_ID_TAG message listener." << std::endl;
	}

	//deserialize data into protocol_id and batch_id
	pair<Size, Size> ids = deserialize_db_listener_data(listener_data);
	protocol_id = ids.first;
	batch_id = ids.second;
	TR << "Received protocol_id='" << protocol_id << "' and batch_id='" << batch_id << "'" << std::endl;

	//no protocol id set yet - create protocol and first batch
	if(protocol_id==0){
		try {
			protocol_id = protocol_features->report_features(protocol_id, db_session);
		} catch (cppdb_error error){
			stringstream err_msg;
			err_msg
				<< "Failed to set the protocol id for batch "
				<< "'" << identifier << "'" << endl
				<< "Error Message:" << endl << error.what() << endl;
			utility_exit_with_message(err_msg.str());
		}

		TR
			<< "Initialize the protocol_id='" << protocol_id << "' "
			<< "and tell it to the head node." << std::endl;

		if(rank != 0) {
			TR << "Send the protocol_id '" << protocol_id << "' to the head node." << std::endl;
			utility::send_string_to_node(0/*HEAD*/, serialize_ids(protocol_id, identifier, batch_id));
		} else {
			MessageListenerOP listener(
				MessageListenerFactory::get_instance()->get_listener(
					DATABASE_PROTOCOL_AND_BATCH_ID_TAG));
			listener->receive(serialize_ids(protocol_id, identifier, batch_id));
		}

	}
	TR << "done with protocol" <<std::endl;
	// setup the batch_id

	try {
		batch_features->report_features(batch_id, protocol_id, identifier, "", db_session);
	} catch (cppdb_error error){
		stringstream err_msg;
		err_msg
			<< "Failed to set the batch id for batch '" << identifier << "' "
			<< "with protocol_id '" << protocol_id << "'" << endl
			<< "Error Message:" << endl << error.what() << endl;
		utility_exit_with_message(err_msg.str());
	}

#endif
#ifndef USEMPI
	//serial case
	if(!protocol_table_initialized_){
		TR << "Initializing protocol table" << endl;
		try{
			static_protocol_id_ = protocol_features->report_features(0, db_session);
		} catch (cppdb_error error){
			stringstream err_msg;
			err_msg
				<< "Failed to set the protocol id for batch "
				<< "'" << identifier << "'" << endl
				<< "Error Message:" << endl << error.what() << endl;
			utility_exit_with_message(err_msg.str());
		}
		protocol_table_initialized_ = true;
	}
	protocol_id = static_protocol_id_;

	if(!static_batch_id_map_.count(identifier)){
		TR << "Initializing batch table" << endl;


		std::string select_max = "SELECT MAX(batch_id) FROM batches;";
		cppdb::statement stmt(basic::database::safely_prepare_statement(select_max,db_session));
		cppdb::result res(basic::database::safely_read_from_database(stmt));

		Size max_batch_id(0);
		while(res.next())
		{
			res >> max_batch_id;
		}


		for(
			std::map< std::string, Size >::const_iterator
				i = static_batch_id_map_.begin(), ie = static_batch_id_map_.end();
			i != ie; ++i){
			max_batch_id = std::max(max_batch_id, i->second);
		}
		batch_id = max_batch_id + 1;
		static_batch_id_map_[identifier] = batch_id;

		try {
			batch_features->report_features(
				batch_id, protocol_id, identifier, "", db_session);
		} catch (cppdb_error error){
			stringstream err_msg;
			err_msg
				<< "Failed to set the batch id for batch '" << identifier << "' "
				<< "with protocol_id '" << protocol_id << "'" << endl
				<< "Error Message:" << endl << error.what() << endl;
			utility_exit_with_message(err_msg.str());
		}
	}
	else{
		batch_id=static_batch_id_map_[identifier];
	}
#endif

	return pair<Size, Size>(protocol_id, batch_id);
}

pair<Size, Size>
deserialize_db_listener_data(
	string data
){
	utility::vector1< std::string > tokens = utility::split(data);
	if(tokens.size() != 2){
		utility_exit_with_message("failed to deserialize the message from master node. Message was: " + data + " You will get this message if trying to run ReportToDB mover in MPI mode with only on processor.");
	}
	int protocol_id=utility::string2int(tokens[1]);
	int batch_id=utility::string2int(tokens[2]);
	return pair<Size, Size>(protocol_id, batch_id);
}

string
serialize_ids(
	int protocol_id,
	string identifier,
	Size batch_id
){
	return
		utility::to_string(protocol_id) + " " +
		identifier + " " +
		utility::to_string(batch_id);
}

///@detail look up the batch id given a struct id. Note this should
///only be used once the structure's table has been created, eg in an
///average features reporter's report_features function.
Size
get_batch_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session
) {

	std::string const stmt_str(
		"SELECT batch_id FROM structures WHERE struct_id = ?;");
	statement stmt(basic::database::safely_prepare_statement(stmt_str, db_session));
	stmt.bind(1, struct_id);
	result res(basic::database::safely_read_from_database(stmt));

	if(!res.next()){
		stringstream err_msg;
		err_msg
			<< "No batch_id found for struct_id '"
			<< to_string(struct_id) << "'";
		utility_exit_with_message(err_msg.str());
	}
	Size batch_id;
	res >> batch_id;
	return batch_id;
}

std::string serialize_residue_xyz_coords(core::conformation::Residue const & residue)
{
	//6bitencode and decode work best with arrays
	core::Real* coord_data = new core::Real[residue.natoms()*3];
	for(core::Size atom_index = 1; atom_index <= residue.natoms();++atom_index)
	{
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
	core::Real* coord_data = new core::Real[natoms*3];
	core::Size memory_size = natoms*3*sizeof(core::Real);
	utility::decode6bit((unsigned char*)coord_data,data);

	utility::vector1< numeric::xyzVector<core::Real> > xyz_vector;
	for(core::Size atom_index = 1; atom_index <= natoms;++atom_index)
	{
		core::Size array_index = (atom_index - 1)*3;
		xyz_vector.push_back(numeric::xyzVector<core::Real>(coord_data[array_index],coord_data[array_index + 1],coord_data[array_index + 2]));
	}
	delete [] coord_data;
	return xyz_vector;
}

} //namespace protocols
} //namespace features
