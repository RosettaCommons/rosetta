
#ifdef USEMPI
#include <mpi.h>
#endif

#include <basic/message_listening/MessageListener.fwd.hh>
#include <basic/message_listening/util.hh>
#include <utility/mpi_util.hh>
#include <utility/assert.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <string>

namespace basic {
namespace message_listening {

using std::string;
using std::endl;
using utility::send_string_to_node;
using utility::receive_string_from_node;


static basic::Tracer TR("protocols.jd2.JobDistributor");

///@brief used for message passing to the
///MPIWorkPoolJobDistributor. This function will ask the head node for
///data.  The type of data returned is based on the type of listener
///created based on the listener_tags of the MessageListenerFactory

string
request_data_from_head_node(
	listener_tags MPI_ONLY( listener_tag ) ,
	string const & MPI_ONLY( data )
){

#ifdef USEMPI

    //send a message to the head node that tells jd2 to create a message listener
    TR << "Requesting data from head node" << endl;
    TR.flush();
    MPI_Send( &listener_tag, 1, MPI_INT, 0/*head node*/, 40 /*REQUEST_MESSAGE_TAG*/, MPI_COMM_WORLD );

    //send a string to be processed by the listener
    TR << "Sending " << data << " to head node" << std::endl;
    TR.flush();
    send_string_to_node(0/*head node*/, data);

    //receive a response from the head node listener
    return utility::receive_string_from_node(0/*head node*/);
#endif
#ifndef USEMPI
    utility_exit_with_message(
		"ERROR: You have tried to request a message from the head node but you are not in mpi mode (compile with extras=mpi)");
#endif

	return "";  // required for compilation on Windows
}

void
send_data_to_head_node(
	basic::message_listening::listener_tags MPI_ONLY( listener_tag ),
	std::string const & MPI_ONLY( data )
){
#ifdef USEMPI
    //send a message to the head node that tells jd2 to create a message listener
  MPI_Send( &listener_tag, 1, MPI_INT, 0/*head node*/, 50 /*RECEIVE_MESSAGE_TAG*/, MPI_COMM_WORLD );

    //send a string to be processed by the listener
    utility::send_string_to_node(0/*head node*/, data);
#endif
#ifndef USEMPI
        utility_exit_with_message("ERROR: You have tried to send a message to the head node but you are not in mpi mode (compile with extras=mpi)");
#endif
}


} // namespace
} // namespace
