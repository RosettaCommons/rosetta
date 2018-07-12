// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/PyMOLMover.cc
/// @brief  Send infromation to PyMOL
/// @author Sergey Lyskov

#ifndef INCLUDED_protocols_moves_PyMOLMover_CC
#define INCLUDED_protocols_moves_PyMOLMover_CC

/*
/// Workaround for:
/// external/boost_1_55_0/boost/bind/mem_fn_template.hpp:156:30: error: no matching function for call to 'get_pointer'
///   BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1);
/// Get pointer of owning_ptr: needed by boost::mem_fn
#ifdef __clang__
#include <memory>
namespace boost {
template<typename T>
inline T* get_pointer(const std::shared_ptr<T>& p) { return p.get(); }
}
#endif
*/

// protocol headers
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/PyMOLMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/Edge.hh>

// numeric headers
#include <numeric/random/uniform.hh>
#include <numeric/types.hh>

// basic headers
#include <basic/Tracer.hh>

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/io/zipstream.ipp>

#include <utility/py/PyAssert.hh>
#include <utility/tag/Tag.hh>

#include <utility/exit.hh>

#include <sstream>

// c++ headers
#include <ctime>
#include <algorithm>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf

// boost headers
//#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/copy.hpp>

//#ifndef WIN_PYROSETTA  // CL compiler got horribly confused if our numeric header got included after <winsock2.h>
//#endif


/*#if  !defined(WINDOWS) && !defined(WIN32)
#include <sys/time.h>
#endif
*/

#ifndef _WIN32
#include "pthread.h"
// XSD XRW Includes
#endif


#ifdef __native_client__
#include <errno.h>
#endif

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace moves {

static basic::Tracer TR( "protocols.moves.PyMOLMover" );


/// We using independent numeric::random::rg() which is not connected to RNG system because we do not want PyMOL to interfere with other Rosetta systems.
/// I.e creating and using PyMOL mover should not change any trajectories of the running program even in production mode.
numeric::random::uniform_RG_OP
getRG()
{
	static numeric::random::uniform_RG_OP RG = nullptr;

	if ( RG == nullptr ) {
		//RG = new numeric::random::mt19937_RG;
		RG = numeric::random::uniform_RG_OP( new numeric::random::standard_RG );
		RG->setSeed( time(nullptr) );
	}

	return RG;
}


UDPSocketClient::UDPSocketClient(std::string const & address, unsigned int port, unsigned int max_packet_size) : max_packet_size_(max_packet_size), sentCount_(0)
{
	if ( port < 1  or  port > 65535 ) utility_exit_with_message("UDPSocketClient::UDPSocketClient: Invalid value for port:" + std::to_string(port) + " expected value must be in [1, 65535] range!");
	if ( max_packet_size_ < 1500  or  max_packet_size_ > 65535 ) utility_exit_with_message("UDPSocketClient::UDPSocketClient: Invalid value for max_packet_size:" + std::to_string(max_packet_size_) + " expected value must be in [1500, 65535] range!");

#ifndef  __native_client__
	/*
#ifdef __APPLE__
	max_packet_size_ = 8192-512-2;  // ← MacOS X kernel can send only small packets even locally.
#else
	max_packet_size_ = 64512; // 1024*63
#endif */

	// moved to PyMOLMover::_default_max_packet_size_, was:  max_packet_size_ = 8192-512-2;  // Choosing max pocket size seems to be less obvious when we have cross OS link, we have to pick less common denominator

	//#ifndef WIN_PYROSETTA
	// generating random uuid by hands
	for ( unsigned short & i : uuid_.shorts_ ) i = (unsigned short) getRG()->getRandom()*65536;  //RG.random_range(0, 65536);

	memset(&socket_addr_, '\0', sizeof(sockaddr_in));

	socket_addr_.sin_family = AF_INET;     // host byte order
	socket_addr_.sin_port = htons(port);  // short, network byte order
	socket_addr_.sin_addr.s_addr = inet_addr( address.c_str() );

	socket_h_ = socket(AF_INET, SOCK_DGRAM, 0);

	if ( socket_h_ == -1 ) utility_exit_with_message("UDPSocketClient::UDPSocketClient: Could not open socket, got error:" + std::to_string(errno) + "!");

	//#endif
#endif
}


UDPSocketClient::UDPSocketClient( UDPSocketClient const & other ) :
	max_packet_size_( other.max_packet_size_ ),
	uuid_( other.uuid_ ),
	sentCount_( other.sentCount_ ),
	socket_addr_( other.socket_addr_ ),
	socket_h_( other.socket_h_ )
{
#ifndef  __native_client__
	// reinit connection using coppied info
	socket_h_ = socket(AF_INET, SOCK_DGRAM, 0);
#endif
}


UDPSocketClient & UDPSocketClient::operator= (UDPSocketClient const&o)
{
#ifndef  __native_client__

	max_packet_size_ = o.max_packet_size_;
	uuid_            = o.uuid_;
	sentCount_       = o.sentCount_;
	socket_addr_     = o.socket_addr_;
	socket_h_        = o.socket_h_;

	socket_h_ = socket(AF_INET, SOCK_DGRAM, 0);

#endif

	return *this;
}


UDPSocketClient::~UDPSocketClient()
{
	//#ifndef WIN_PYROSETTA
#ifndef  __native_client__
#ifdef WIN32
	closesocket(socket_h_);
#else
	close(socket_h_);
#endif
	//#endif
#endif
}

void UDPSocketClient::sendMessage(std::string msg)
{
#ifndef  __native_client__
	int count = 1;
	if ( msg.size() > max_packet_size_ )  { count = msg.size()/max_packet_size_ + 1; }

	for ( int i=0; i<count; ++i ) {
		unsigned int last = (i+1)*max_packet_size_;
		if ( last>msg.size() ) last = msg.size();
		sendRAWMessage(sentCount_, i, count, &msg[i*max_packet_size_], &msg[last] );
		if ( count > 1 ) {
#ifdef _WIN32
			Sleep(10); // Sleep function takes milliseconds.
#else
			timespec ts;  ts.tv_sec=0;  ts.tv_nsec=1000000; //time to sleep in nanoseconds, we want to take a nap for ~0.001sec
			nanosleep(&ts, nullptr);
#endif
		}
	}

	++sentCount_;
#endif
}

void UDPSocketClient::sendRAWMessage(int globalPacketID, int packetI, int packetCount, char * msg_begin, char *msg_end)
{
#ifndef  __native_client__
	std::string buf(msg_end - msg_begin + sizeof(uuid_.bytes_) + sizeof(short)*3, 0);
	int i = 0;

	memcpy(&buf[i], uuid_.bytes_, sizeof(uuid_.bytes_));  i+= sizeof(uuid_.bytes_);
	memcpy(&buf[i], &globalPacketID, 2);  i+=2;
	memcpy(&buf[i], &packetI, 2);  i+=2;
	memcpy(&buf[i], &packetCount, 2);  i+=2;
	memcpy(&buf[i], msg_begin, msg_end-msg_begin); // i+=msg_end-msg_begin; //THIS VALUE IS NEVER USED. if more after msg added, uncomment!
	//#ifndef WIN_PYROSETTA
	sendto(socket_h_, &buf[0], buf.size(), 0 , (struct sockaddr *)&socket_addr_, sizeof(struct sockaddr_in));
	//#endif
#endif
}

void
UDPSocketClient::show(std::ostream & output) const
{
#ifndef  __native_client__
	output << "max packet size: " << max_packet_size_ << std::endl;
	output << "sent count: " << sentCount_ << std::endl;
	output << "socket handle: " << socket_h_ << std::endl;

	output << "uuid shorts: ";
	for ( unsigned short i : uuid_.shorts_ ) {
		output << i << " ";
	}
	output << std::endl;

	// output << "uuid byte: ";
	// for ( char byte : uuid_.bytes_ ) {
	//  output << static_cast<int>(byte) << " ";
	// }
	// output << std::endl;

	output << "socket address family: " << socket_addr_.sin_family << std::endl;
	output << "socket address address: " << socket_addr_.sin_addr.s_addr << std::endl;
	output << "socket address port: " << socket_addr_.sin_port << std::endl;
	output << "socket max_packet_size: " << max_packet_size_ << std::endl;
#endif
}

std::ostream &
operator<<(std::ostream & output, UDPSocketClient const & client)
{
	client.show(output);
	return output;
}

/* -------------------------------------------------------------------------------------------------
PyMOLMover Class
---------------------------------------------------------------------------------------------- */

// std::string  const PyMOLMover::_default_address_         = "127.0.0.1";
// unsigned int const PyMOLMover::_default_port_            = 65000;
// unsigned int const PyMOLMover::_default_max_packet_size_ = 8192-512-2;  // Choosing max pocket size seems to be less obvious when we have cross OS link, we have to pick less common denominator

std::string PyMOLMover::default_address()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::run::PyMOLMover;
	return option[address].user() ? option[address].value() : "127.0.0.1";
}

unsigned int PyMOLMover::default_port()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::run::PyMOLMover;
	return option[port].user() ? option[port].value() : 65000;
}

unsigned int PyMOLMover::default_max_packet_size(std::string const & address)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::run::PyMOLMover;
	if ( option[max_packet_size].user() ) return option[max_packet_size].value();
	else {
		// Choosing max pocket size seems to be less obvious when we have cross OS link, we have to pick less common denominator
		// the 8192-512-2 is a max value that seems to work on FreeBSD kernels
		if ( address == "127.0.0.1" ) return 8192-512-2;
		else return 1500; // MTU for enthernet equpment
	}
}


/// @brief ctor
PyMOLMover::PyMOLMover(std::string const & address, unsigned int port, unsigned int max_packet_size) :
	link_( address, port, max_packet_size ? max_packet_size : default_max_packet_size(address) ),
	update_energy_(false),
	energy_type_(core::scoring::total_score),
	update_membrane_(false),
	keep_history_(false),
	update_interval_(0),
	last_packet_sent_time_(0),
	pymol_name_()
{}

/// @brief cctor
PyMOLMover::PyMOLMover( PyMOLMover const & ) = default;

PyMOLMover::~PyMOLMover() = default;

// XRW TEMP std::string PyMOLMover::get_name() const
// XRW TEMP {
// XRW TEMP  return "PyMOL_Mover";
// XRW TEMP }

void PyMOLMover::set_PyMOL_model_name( std::string name ){
	pymol_name_ = name;
}

std::string PyMOLMover::get_PyMOL_model_name(Pose const & pose) const
{
	if ( pymol_name_.size() ) {
		return pymol_name_;
	} else {
		core::pose::PDBInfoCOP info = pose.pdb_info();
		if ( info && info->name().size() ) {
			std::string n = info->name();
			for ( char & i : n ) if ( i == '/' ) i = '_';
			return n;
		} else {
			return "pose";
		}
	}
}


bool PyMOLMover::is_it_time()
{
	// First let's check if enough time have passes since last time we send info...
	//double t = clock() / CLOCKS_PER_SEC;
	double t = time(nullptr);
	//TR << "t=" << t << " cl="<< clock() << std::endl;
	if ( t - last_packet_sent_time_ < update_interval_ ) return false;
	last_packet_sent_time_ = t;
	return true;
}


void PyMOLMover::apply( Pose const & pose)
{
	TR.Trace << "PyMOLMover::apply( Pose const & pose) ..." << std::endl;

	if ( !is_it_time() ) return;
	TR.Trace << "PyMOL_Mover::apply It is time!" << std::endl;

	std::string name = get_PyMOL_model_name(pose);
	TR.Trace << "PyMOL_Mover::apply name:" << name << std::endl;

	// Check if pose is a membrane pose
	if ( pose.conformation().is_membrane() ) {
		update_membrane_ = true;
	}

	// Creating message
	std::ostringstream os;
	pose.dump_pdb(os);

	// Compressing message
	std::ostringstream zmsg;
	zlib_stream::zip_ostream zipper(zmsg, true);
	zipper << os.str();
	zipper.zflush_finalize();

	//std::string message = "PDB     X" + name + zmsg.str();
	std::string message = "PDB.gzipXX" + name + zmsg.str();
	message[8] = keep_history_;
	message[9] = name.size();

	//TR << "Sending message: " << message << std::endl << "Size:" << message.size() << std::endl;
	//TR << "Sending message, Size:" << message.size() << std::endl;
	//TR << "debug in the apply function line 335 the send that the link_ sent would be "<< message << std::endl;
	link_.sendMessage(message);

	if ( update_membrane_ ) send_membrane_planes(pose);
	if ( update_energy_ ) send_energy(pose, energy_type_);
}

void PyMOLMover::apply( Pose & pose)
{
	Pose const & p(pose);
	apply(p);
}

void PyMOLMover::print(std::string const & message)
{
	if ( !is_it_time() ) return;

	std::string msg =  std::string("Text    ") + char(keep_history_) + char(0) /* Place holder for name size = 0 */ + message;

	link_.sendMessage(msg);
}

void PyMOLMover::send_energy(Pose const &pose, core::scoring::ScoreType score_type)
{
#ifndef  __native_client__
	if ( !is_it_time() ) return;

	if ( pose.energies().energies_updated() ) {

		utility::vector1<core::Real> e(pose.size());
		core::Real min=1e100, max=1e-100;
		for ( unsigned int i=1; i<=e.size(); ++i ) {
			if ( score_type == core::scoring::total_score ) e[i] = pose.energies().residue_total_energy(i);
			else e[i] = pose.energies().residue_total_energies(i)[score_type];

			if ( min > e[i] ) min = e[i];
			if ( max < e[i] ) max = e[i];
		}
		// We not using send_RAW_Energies for efficiency reasons...
		std::string msg(8*e.size(), ' ');
		core::pose::PDBInfoCOP info = pose.pdb_info();
		for ( unsigned int i=1; i<=e.size(); ++i ) {
			char chain = ' ';
			char icode = ' ';
			int  res = i;
			if ( info ) {
				chain = info->chain(i);
				icode = info->icode(i);
				res = info->number(i);
			}
			char buf[256];
			//TR << "Energy is for the residue " << i << " is " << e[i] << std::endl;
			e[i] = (e[i]-min)*255. / (max-min+1e-100);
			sprintf(buf, "%c%4d%c%02x", chain, res, icode, int(e[i]));
			for ( int k=0; k<8; ++k ) msg[(i-1)*8+k] = buf[k];
		}

		// Compressing message
		std::ostringstream zmsg;
		zlib_stream::zip_ostream zipper(zmsg, true);
		zipper << msg;
		zipper.zflush_finalize();

		std::string name = get_PyMOL_model_name(pose);
		std::string sname = core::scoring::name_from_score_type(score_type);

		std::string message =  std::string("Ene.gzip") + char(keep_history_) \
			+ char(name.size()) + name \
			+ char(sname.size()) + sname + zmsg.str();

		//TR << "Sending message: " << message << std::endl << "Size:" << message.size() << std::endl;
		//TR << "Sending message, Size:" << message.size() << std::endl;

		link_.sendMessage(message);
	}
#endif
}

/// Send specified energy to PyMOL.
void PyMOLMover::send_energy(Pose const &pose, std::string const & stype)
{
	send_energy(pose, core::scoring::ScoreTypeManager::score_type_from_name(stype) );
}


void PyMOLMover::label_energy(core::pose::Pose const &input_pose , std::string score_type ="total_score"){
	//Displays <sigs> number of characters for each energy with labels on CA.
	//&&did not find the equivalent of _get_energies
	core::scoring::Energies energy=input_pose.energies();
	//(this connection).send_energy();
	send_energy( input_pose, score_type);
}


void PyMOLMover::send_RAW_Energies(Pose const &pose, std::string energyType, utility::vector1<int> const & energies)
{
#ifndef  __native_client__
	if ( !is_it_time() ) return;

	std::string msg(8*energies.size(), ' ');
	core::pose::PDBInfoCOP info = pose.pdb_info();
	for ( unsigned int i=1; i<=energies.size(); ++i ) {
		char chain = ' ';
		char icode = ' ';
		int  res = i;
		if ( info ) {
			chain = info->chain(i);
			icode = info->icode(i);
			res = info->number(i);
		}
		char buf[256];
		sprintf(buf, "%c%4d%c%02x", chain, res, icode, energies[i]);
		for ( int k=0; k<8; ++k ) msg[(i-1)*8+k] = buf[k];
	}
	//TR << msg << std::endl;

	// Compressing message
	std::ostringstream zmsg;
	zlib_stream::zip_ostream zipper(zmsg, true);
	zipper << msg;
	zipper.zflush_finalize();

	std::string name = get_PyMOL_model_name(pose);
	std::string sname = energyType;

	std::string message = std::string("Ene.gzip") + char(keep_history_) \
		+ char(name.size()) + name \
		+ char(sname.size()) + sname + zmsg.str();

	//TR << "Sending message: " << message << std::endl << "Size:" << message.size() << std::endl;
	//TR << "Sending message, Size:" << message.size() << std::endl;

	link_.sendMessage(message);
#endif
}


/// @brief Send Membrane Planes to PyMOL
/// @details If pose is a membrane pose
/// pymol viewer will build CGO planes from points specified
void PyMOLMover::send_membrane_planes( Pose const & pose ) {

#ifndef __native_client__

	using namespace core::scoring::methods;

	if ( !is_it_time() ) return;

	// Check the membrane planes can be visualized
	if ( !pose.conformation().is_membrane() ) return;

	// Grab a list of relevant residues and go
	// Compute radius of gyration of the pose
	utility::vector1< bool > relevant_residues;
	relevant_residues.resize( pose.size() );
	for ( Size i = 1; i < relevant_residues.size(); ++i ) {
		relevant_residues[i] = true;
	}

	// Compute Radius of Gyration
	RG_Energy_Fast rg_method;
	core::Real rg =  2*rg_method.calculate_rg_score( pose, relevant_residues );

	// Get the normal vector and center position
	core::conformation::Conformation const & conf( pose.conformation() );
	core::Vector normal( conf.membrane_info()->membrane_normal( conf ) );
	core::Vector center( conf.membrane_info()->membrane_center( conf ) );

	// Get the plane thickness used by Rosetta
	core::Real thickness( conf.membrane_info()->membrane_thickness() );

	// Encode the Center vector
	std::string center_msg = "";
	center_msg += utility::to_string( center.x() );
	center_msg += ",";
	center_msg += utility::to_string( center.y() );
	center_msg += ",";
	center_msg += utility::to_string( center.z() );

	// Encode normal vector
	std::string normal_msg = ",";
	normal_msg += utility::to_string( normal.x() );
	normal_msg += ",";
	normal_msg += utility::to_string( normal.y() );
	normal_msg += ",";
	normal_msg += utility::to_string( normal.z() );

	// Encode the thickness
	std::string thickness_msg = ",";
	thickness_msg += utility::to_string( thickness );

	// Encode the radius of gyration
	std::string rg_msg = ",";
	rg_msg += utility::to_string( rg );

	// Construct full message
	std::string msg = center_msg + normal_msg + thickness_msg + rg_msg;

	// Compressing message
	std::ostringstream zmsg;
	zlib_stream::zip_ostream zipper( zmsg, true );
	zipper << msg;
	zipper.zflush_finalize();

	std::string name = get_PyMOL_model_name(pose);
	std::string sname = "membrane_planes";

	std::string message =  std::string("Mem.gzip") + char(keep_history_) \
		+ char(name.size()) + name \
		+ char(sname.size()) + sname + zmsg.str();

	link_.sendMessage(message);

#endif

}


void PyMOLMover::send_colors(Pose const &pose, std::map<int, int> const & colors, X11Colors default_color)
{
#ifndef  __native_client__
	utility::vector1<int> energies( pose.size(), default_color);  // energies = [ X11Colors[default_color][0] ] * pose.size()

	for ( auto const & color : colors ) {
		PyAssert( color.first >=1 && color.first <= static_cast<int>(pose.size()),
			"PyMOLMover::send_colors residue index is out of range!");
		PyAssert( color.second >= XC_first_color && color.second <= XC_last_color,
			"PyMOLMover::send_colors color index is out of range!");

		energies[ color.first ] = color.second;  // for r in colors: energies[r-1] = X11Colors[ colors[r] ][0]
	}
	send_RAW_Energies(pose, "X11Colors", energies);  //self._send_RAW_Energies(pose, 'X11Colors', energies, autoscale=False)
#endif
}

void PyMOLMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Keep history:          " << ( ( keep_history_ ) ? ("True") : ("False") ) << std::endl;
	output << "Update energy:         " << ( ( update_energy_ ) ? ("True") : ("False") ) << std::endl;
	output << "Last packet sent time: " << last_packet_sent_time_ << std::endl;
	output << "Update interval:       " << update_interval_ << std::endl;
	output << "Link:                  " << link_ << std::endl;
}


/// Xiyao's Code
/*
###########################################################################
# port "Luxury" methods in python code into c++
*/


void
PyMOLMover::send_any( std::string ptype, core::pose::Pose const & pose, utility::vector1< std::string > data, core::Size  size /* = 6 */ ){

	/// ptype is a tag for the type of data
	/// pose is the pose in PyMOL (the size what matters)
	/// data is a vector of strings, same size as pose
	std::ostringstream to_send;
	to_send << std::to_string( size );

	core::pose::PDBInfoCOP info( pose.pdb_info() );

	for ( core::Size i=1; i<=pose.size(); ++i ) {
		std::ostringstream pdb;
		pdb << std::right << std::setw(6); // set formatting

		if ( info != nullptr  && info->nres() != 0 ) {
			std::string pdb_resi_info = std::to_string(info->number(i)) + info->icode(i) + info->chain(i);
			pdb << pdb_resi_info;
		} else {
			pdb << std::to_string(i) + " A";

		}

		//proper formatting so 0 doesn't cast to 0.0 and is less than size
		std::ostringstream dat;
		dat << std::right << std::setw(size);

		if ( data[i] == "0" && size > 2 ) {
			dat << data[i];
		} else {
			dat << data[i];
		}
		to_send << pdb.str() << dat.str();
	}

	std::string name = get_PyMOL_model_name(pose);


	//compressing message
	std::ostringstream msgz;
	zlib_stream::zip_ostream zipper( msgz, true);
	//TR << "the ptype is \n" << ptype << std::endl;
	//TR << "the to send message  is \n" << to_send.str() << std::endl;
	zipper << to_send.str();

	zipper.zflush_finalize();
	std::string message = ptype + char(keep_history_) \
		+ char(name.size()) + name \
		+ msgz.str();
	// ptype << (char)keep_history() << (char)name.size() << name << to_send.rdbuf();
	link_.sendMessage( message );
}



//H-bonds.
//Sends list of hydrogen bonds and displays them in PyMOL.
//Makes use of PyMOL's "distance" function
void
PyMOLMover::send_hbonds( core::pose::Pose const & pose){

	// Check that the energies are updated.
	if ( !pose.energies().energies_updated() ) {
		//TR << "PyMOL_Mover::send_hbonds: Energy is not updated! please score the pose first!" << std::endl;
	}

	// Get the H-bonds.
	core::scoring::hbonds::HBondSet hbset;
	core::scoring::hbonds::fill_hbond_set(pose,false, hbset);

	std::stringstream to_send;
	to_send << std::right << std::setw(5) << std::to_string(hbset.nhbonds());

	//Get the energies.
	core::pose::PDBInfoCOP info( pose.pdb_info() );
	utility::vector1< numeric::Real> energy_list;
	for ( core::Size i = 1; i <= hbset.nhbonds(); ++i ) {
		energy_list.push_back(hbset.hbond(i).energy());
	}

	if ( !energy_list.empty() ) {
		numeric::Real max_e = *std::max_element(energy_list.begin(), energy_list.end());
		numeric::Real min_e = *std::min_element(energy_list.begin(), energy_list.end());

		for ( core::Size i = 1; i <= hbset.nhbonds(); ++i ) {
			core::scoring::hbonds::HBond hb = hbset.hbond(i);
			std::string  accatm = pose.residue(hb.acc_res()).atom_name(hb.acc_atm());
			std::string  donatm = pose.residue(hb.don_res()).atom_name(hb.don_hatm());

			// Each H-bond is 6 + 4 + 6 + 4 + 2 = 22 chars.
			if ( info != nullptr && info->nres() != 0 ) {
				// acc codes for acceptor atom of H bonds
				std::string acc_res = std::to_string(info->number(hb.acc_res()));
				char acc_icode = info->icode(hb.acc_res());
				char acc_chain = info->chain(hb.acc_res());

				// don codes for donor atom of H bonds
				std::string don_res = std::to_string(info->number(hb.don_res()));
				char don_icode = info->icode(hb.don_res());
				char don_chain = info->chain(hb.don_res());
				std::stringstream ss1, ss2;
				std::string acc_info = acc_res + acc_icode + acc_chain;
				std::string don_info = don_res + don_icode + don_chain;
				ss1 << std::right << std::setw(6) << acc_info;
				ss2 << std::right << std::setw(6) << don_info;

				// Compressing energy value. Format specifier %02x would print at least 2 digit and prepend 0 if there's less than 2, x means number is int.
				char buf[256];
				sprintf(buf, "%02x", int((energy_list[i] - min_e)*255./(max_e - min_e)));
				to_send  << ss1.str() << accatm << ss2.str() << donatm << buf;
			} else {

				std::stringstream ss1, ss2, out_energy;
				ss1 << std::right << std::setw(5)<<std::to_string(hb.acc_res());
				ss2 << std::right << std::setw(5)<<std::to_string(hb.don_res());
				out_energy << std::hex << std::setw(5)<< std::to_string(hb.energy());
				to_send << ss1.str() << accatm << ss2.str() << donatm << out_energy.str();
			}
		}
		//TR << to_send.str() << std::endl;

		// Compressing message
		std::ostringstream zmsg;
		zlib_stream::zip_ostream zipper(zmsg, true);
		zipper << to_send.str();
		zipper.zflush_finalize();
		std::string name = get_PyMOL_model_name( pose);
		std::string message = std::string("hbd.gzip") \
			+ char(keep_history_) \
			+ char(name.size()) \
			+ name \
			+ zmsg.str();
		link_.sendMessage( message);
	} else {
		//TR<< "No H-bonds could be determined for your pose!"<< std::endl;
	}
}



///@brief Returns a list of energies of type energy_type from the pose.
utility::vector1< numeric::Real > get_energies( core::pose::Pose const & pose, core::scoring::ScoreType energy_type ) {

	utility::vector1< numeric::Real > energies;

	// Check to make sure the energies are available.
	if ( !pose.energies().energies_updated() ) {
		//TR << "PyMOL_Mover::send_specific_energy:\n Energy is not updated; please score the pose first!" << std::endl;
	}

	for ( int i = 1; i <= int(pose.size()); ++i ) {
		energies.push_back(pose.energies().residue_total_energies(i)[energy_type]);
	}

	return energies;
}



/*
Secondary-structure assignment using DSSP.
Sends the DSSP assignment for pose to PyMOL and shows as a cartoon.
Useful for when you are making moves to a pose that change secondary
structure, and you wish for PyMOL to display the changes.
*/

void
PyMOLMover::send_ss(core::pose::Pose &pose, std::string ss = ""){

	// Get ss.
	utility::vector1< std::string > ssv;
	if ( ss.empty() ) {
		core::scoring::dssp::Dssp dssp = core::scoring::dssp::Dssp(pose);
		dssp.insert_ss_into_pose(pose);
		std::string dssp_str = dssp.get_dssp_secstruct();
		//or std::string const secstruct = get_secstruct( pose );?
		//ssv.push_back(pose.secstruct());
		for ( char& c : dssp_str ) {
			std::string s(1, c);
			ssv.push_back(s);
		}
	} else {
		for ( char& c : ss ) {
			std::string s(1, c);
			ssv.push_back(s);
		}
	}
	send_any(" ss.gzip", pose, ssv, 1);
}



//Polar identity per residue.
//Colors polar residues red and nonpolar residues blue.

void
PyMOLMover::send_polars(core::pose::Pose const &pose) {

	// Send 1 or 0, if polar or not.
	utility::vector1< std::string > data;

	// Get a string coding whether each residue is polar or not
	for ( int i = 1; i <= int(pose.size()); ++i ) {
		data.push_back(std::to_string(int(pose.residue_type(i).is_polar())));
	}

	send_any("pol.gzip", pose, data, 1);
}



//MoveMap DOF info per residue in pose.

void
PyMOLMover::send_movemap( core::pose::Pose const &pose, core::kinematics::MoveMap const &movemap) {

	//Colors movable regions green and non-movable regions red.
	utility::vector1< std::string > data;
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		data.push_back( "00");
	}
	for ( auto & i : data ) {
		std::cout << i << ' ';
	}
	for ( core::Size i=1 ; i <=pose.size(); ++ i ) {
		int num = 11 + int(movemap.get_bb(i)) * 10 + int(movemap.get_chi(i));
		data[i] = std::to_string(num);
	}
	send_any("mm1.gzip", pose, data, 2);
}



/*
Colors the pose by fold tree information.
Cutpoints (e.g., the C-termini of protein chains) are colored red.
Jump points are colored orange. (Unfortunately, no indication of which
jump point connects to which jump point is given at this time.)
Loops are colored an assortment of colors other than red or orange.
All other residues are colored gray.
See also:
PyMOL_Mover.view_foldtree_diagram()
*/

void
PyMOLMover::send_foldtree(core::pose::Pose const &pose, core::kinematics::FoldTree const &foldtree) {

	//If not sent, use pose's FoldTree.
	utility::vector1< std::string > data;
	for ( core::Size i=1; i<= pose.size(); i++ ) {
		data.push_back("00");
	}

	//# Remove jump data for identification later.
	core::Size njump = foldtree.num_jump();
	std::string starts = std::string( njump,'0');
	std::string stops = std::string( njump,'0');

	//in_loops = []
	utility::vector1<int> in_loops;

	//List of start and stop of jump edges, i.e., loops.
	//for x in range(0, njump):
	for ( core::Size x =0; x< njump; x++ ) {
		//TR << "this is in for loop and x is now " << x << std::endl;
		core::kinematics::Edge loop = foldtree.jump_edge(x + 1);
		Size s1 = loop.start();
		Size s2 = loop.stop();
		if ( s1 < s2 ) {
			starts[x] = s1;
			stops[x] = s2;
		} else {
			starts[x] = s2;
			stops[x] = s1;
		}
	}

	//Each residue is either a jump point, a cutpoint, in a loop, or else in a regular edge of the fold tree.
	//for i in xrange(1, pose.total_residue() + 1):
	for ( core::Size i =1; i<= pose.size(); i++ ) {

		// Keeps the identity relative to the jump, not entry.
		if ( foldtree.is_jump_point(i) ) { // # Jump point residue: color 1.
			data[i]="1";

			// After this point, we will either be entering a loop or leaving a loop.
			for ( int j =0; j< int(starts.length()); j++ ) {

				//j in range(0 , len(starts)):
				if ( i == core::Size(starts[j]) ) {
					in_loops.push_back(j + 1);
				} else if ( i == core::Size(stops[j]) ) {
					in_loops.erase(std::remove(in_loops.begin(), in_loops.end(), j+1));
				}
			}
		} else if ( foldtree.is_cutpoint(i) ) {

			//Cutpoint residue: color 0
			data[i] = "0";
		} else if ( !in_loops.empty() ) {

			//Residue in a loop: color varies.
			/* Up to 7 loops can easily be viewed.
			colors for more loops.) */
			data[i] = std::to_string(2 + *std::max_element(in_loops.begin(),in_loops.end()));

		} else {
			//All other residues: color 2
			data[i] = "2";
		}
	}

	send_any("ft1.gzip", pose, data, 1);
}


/// End Xiyao's Code

std::ostream &
operator<<(std::ostream & output, PyMOLMover const & mover)
{
	mover.show(output);
	return output;
}

PyMOLObserver::PyMOLObserver():
	CacheableObserver(),
	type_( no_observer ) // We have to set the observer type specifically.
{
}

PyMOLObserver::PyMOLObserver(PyMOLObserver const & rval) :
	CacheableObserver( rval ),
	type_( rval.type_ ),
	pymol_(rval.pymol_)
	// Do NOT copy the *_event_link_s
{
}

PyMOLObserver::~PyMOLObserver() {
	detach_from();
}

PyMOLObserver &
PyMOLObserver::operator= (PyMOLObserver const &rval) {
	if ( this != &rval ) {
		core::pose::datacache::CacheableObserver::operator=( rval );
		type_ = rval.type_;
		pymol_ = rval.pymol_;
		// Do NOT copy the *_event_link_s
	}
	return *this;
}

core::pose::datacache::CacheableObserverOP
PyMOLObserver::clone() {
	return core::pose::datacache::CacheableObserverOP( new PyMOLObserver( *this ) );
}

core::pose::datacache::CacheableObserverOP
PyMOLObserver::create() {
	return core::pose::datacache::CacheableObserverOP( new PyMOLObserver );
}

void
PyMOLObserver::set_type( ObserverType setting ) {
	type_ = setting;
	// We don't have a pose, so wait until we get attached
}

void
PyMOLObserver::add_type( ObserverType setting ) {
	type_ = type_ | setting;
	// We don't have a pose, so wait until we get attached
}

void PyMOLObserver::attach(core::pose::Pose &p)
{
	attach_to(p);
}

void PyMOLObserver::detach(core::pose::Pose & /*p*/)
{
	detach_from();
}

bool
PyMOLObserver::is_attached() const {
	return general_event_link_.valid() || energy_event_link_.valid() || conformation_event_link_.valid();
}

void
PyMOLObserver::attach_impl(core::pose::Pose & pose) {
	general_event_link_.invalidate();
	energy_event_link_.invalidate();
	conformation_event_link_.invalidate();

	if ( type_ & general_observer ) {
		general_event_link_ = pose.attach_general_obs( &PyMOLObserver::generalEvent, this );
	}
	if ( type_ & energy_observer ) {
		energy_event_link_ = pose.attach_energy_obs( &PyMOLObserver::energyEvent, this );
	}
	if ( type_ & conformation_observer ) {
		conformation_event_link_ = pose.attach_conformation_obs( &PyMOLObserver::conformationEvent, this );
	}
}

void
PyMOLObserver::detach_impl() {
	general_event_link_.invalidate();
	energy_event_link_.invalidate();
	conformation_event_link_.invalidate();
}

PyMOLObserverOP
get_pymol_observer(core::pose::Pose & pose) {
	using namespace core::pose::datacache;

	if ( !pose.observer_cache().has( CacheableObserverType::PYMOL_OBSERVER ) ) {
		PyMOLObserverOP obs( new PyMOLObserver );
		pose.observer_cache().set( CacheableObserverType::PYMOL_OBSERVER, obs, /*autoattach*/ false );
	}
	CacheableObserverOP obs = pose.observer_cache().get_ptr( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER );
	return utility::pointer::dynamic_pointer_cast< PyMOLObserver >( obs );
}

PyMOLObserverOP AddPyMOLObserver(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
	PyMOLObserverOP o( get_pymol_observer(p) );
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	o->add_type( PyMOLObserver::general_observer );
	o->attach(p);
	return o;
}

PyMOLObserverOP AddPyMOLObserver_to_energies(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
	PyMOLObserverOP o( get_pymol_observer(p) );
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	o->add_type( PyMOLObserver::energy_observer );
	o->attach(p);
	return o;
}

PyMOLObserverOP AddPyMOLObserver_to_conformation(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
	PyMOLObserverOP o( get_pymol_observer(p) );
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	o->add_type( PyMOLObserver::conformation_observer );
	o->attach(p);
	return o;
}

/// @brief PyMOLMoverCreator interface, name of the mover
// XRW TEMP std::string PyMOLMover::mover_name() {
// XRW TEMP  return "PyMOLMover";
// XRW TEMP }

/// @brief PyMOLMoverCreator interface, returns a unique key name to be used in xml file
// XRW TEMP std::string PyMOLMoverCreator::keyname() const {
// XRW TEMP  return PyMOLMover::mover_name();
// XRW TEMP }

/// @brief PyMOLMoverCreator interface, return a new instance
// XRW TEMP protocols::moves::MoverOP PyMOLMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new PyMOLMover() );
// XRW TEMP }

/// @brief allows for the setting of certain variabel from the rosetta scripts interface, only keep history
void
PyMOLMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	keep_history(tag->getOption<bool>( "keep_history", keep_history_ ) );

	std::string address = tag->getOption<std::string>("address", default_address());
	auto port = tag->getOption<unsigned int>("port", default_port());
	auto max_packet_size = tag->getOption<unsigned int>("max_packet_size", default_max_packet_size(address));

	//TR << "Settin addres to: " << address << std::endl;
	//TR << "Settin port to: " << port << std::endl;

	link_ = UDPSocketClient(address, port, max_packet_size);
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PyMOLMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PyMOLMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PyMOLMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::moves::PyMOLMover( *this ) );
}

std::string PyMOLMover::get_name() const {
	return mover_name();
}

std::string PyMOLMover::mover_name() {
	return "PyMOLMover";
}

void PyMOLMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction port_range = integer_range_restriction( "port_range", 1, 65535);
	xsd.add_top_level_element(port_range);

	XMLSchemaRestriction max_packet_size_range = integer_range_restriction( "max_packet_size_range", 1500, 65535);
	xsd.add_top_level_element(max_packet_size_range);

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"keep_history", xsct_rosetta_bool,
		"Each call to the mover stores the pose in a new state/frame of "
		"an object in PyMOL rather than overwriting it. Frames can then be "
		"played back like a movie to visualize the flow of a protocol.",
		"0")
		+ XMLSchemaAttribute::attribute_w_default("address",         xs_string, "IP address of machine running PyMOL with PyMOL-RosettaServer.py script running.", default_address())
		+ XMLSchemaAttribute::attribute_w_default("port",            "port_range", "Port number to which UDP/IP connection should be made", std::to_string( default_port() ) )
		+ XMLSchemaAttribute::attribute_w_default("max_packet_size", "max_packet_size_range", "Max size of packets to send over UDP/IP connection. Default value is dependable from address: for 127.0.0.1 it will be set to 7678 and any other addresses to 1500.", std::to_string( default_max_packet_size(default_address()) ) )
		;

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"PyMOLMover will send a pose to an instance of the PyMOL "
		"molecular visualization software running on the local host. "
		"Each call of the mover overwrites the object in PyMOL. "
		"It is not a full featured as the version built in to PyRosetta "
		"but is extremely useful for visualizing the flow of a protocol "
		"or generating a frames for a movie of a protocol.", attlist );
}

std::string PyMOLMoverCreator::keyname() const {
	return PyMOLMover::mover_name();
}

protocols::moves::MoverOP
PyMOLMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PyMOLMover );
}

void PyMOLMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PyMOLMover::provide_xml_schema( xsd );
}


} // moves
} // protocols

#endif // INCLUDED_protocols_moves_PyMOLMover_CC
