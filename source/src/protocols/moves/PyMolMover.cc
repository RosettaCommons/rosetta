// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/PyMolMover.cc
/// @brief  Send infromation to PyMol
/// @author Sergey Lyskov

#ifndef INCLUDED_protocols_moves_PyMolMover_CC
#define INCLUDED_protocols_moves_PyMolMover_CC

// protocol headers
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/PyMolMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/io/zipstream.ipp>

#include <utility/PyAssert.hh>
#include <utility/tag/Tag.hh>

// numeric headers
#include <numeric/random/uniform.hh>

// basic headers
#include <basic/Tracer.hh>

// c++ headers
#include <time.h>

//#ifndef WIN_PYROSETTA  // CL compiler got horribly confused if our numeric header got included after <winsock2.h>
//#endif


/*#if  !defined(WINDOWS) && !defined(WIN32)
	#include <sys/time.h>
#endif
*/

#ifndef _WIN32
	#include "pthread.h"
#endif


namespace protocols {
namespace moves {

static basic::Tracer TR("protocols.moves.PyMolMover");

//static numeric::random::RandomGenerator RG(9636236);

/// We using independent RG which is not connected to RNG system because we do not want PyMOL to interfere with other Rosetta systems.
/// I.e creating and using PyMOL mover should not change any trajectories of the running program even in production mode.
numeric::random::uniform_RG_OP getRG()
{
	static numeric::random::uniform_RG_OP RG = 0;

	if( RG == 0 ) {
		//RG = new numeric::random::mt19937_RG;
		RG = new numeric::random::standard_RG;
		RG->setSeed( time(NULL) );
	}

	return RG;
}


UDPSocketClient::UDPSocketClient() : sentCount_(0)
{
#ifndef  __native_client__
    /*
	#ifdef __APPLE__
		max_packet_size_ = 8192-512-2;  // ← MacOS X kernel can send only small packets even locally.
	#else
		max_packet_size_ = 64512; // 1024*63
	#endif */

	max_packet_size_ = 8192-512-2;  // Choosing max pocket size seems to be less obvious when we have cross OS link, we have to pick less common denominator

	//#ifndef WIN_PYROSETTA
	// generating random uuid by hands
	for(unsigned int i=0; i<sizeof(uuid_.shorts_)/sizeof(uuid_.shorts_[0]); i++) uuid_.shorts_[i] = (unsigned short) getRG()->getRandom()*65536;  //RG.random_range(0, 65536);

	memset(&socket_addr_, '\0', sizeof(sockaddr_in));

	socket_addr_.sin_family = AF_INET;     // host byte order
	socket_addr_.sin_port = htons(65000);  // short, network byte order
	socket_addr_.sin_addr.s_addr = inet_addr("127.0.0.1");

	socket_h_ = socket(AF_INET, SOCK_DGRAM, 0);
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
	if( msg.size() > max_packet_size_ )  { count = msg.size()/max_packet_size_ + 1; }

	for(int i=0; i<count; i++) {
		unsigned int last = (i+1)*max_packet_size_;
		if( last>msg.size() ) last = msg.size();
		sendRAWMessage(sentCount_, i, count, &msg[i*max_packet_size_], &msg[last] );
		if(count > 1) {
			#ifdef _WIN32
				Sleep(10); // Sleep function takes milliseconds.
			#else
				timespec ts;  ts.tv_sec=0;  ts.tv_nsec=1000000; //time to sleep in nanoseconds, we want to take a nap for ~0.001sec
				nanosleep(&ts, NULL);
			#endif
		}
	}

	sentCount_++;
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
		memcpy(&buf[i], msg_begin, msg_end-msg_begin);  i+=msg_end-msg_begin;
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
	output << "socket handel: " << socket_h_ << std::endl;

	output << "uuid short: ";
	for ( unsigned int i(0); i < sizeof(uuid_.shorts_) / sizeof(uuid_.shorts_[0]); ++i ) {
		output << uuid_.shorts_[i] << " ";
	}
	output << std::endl;

	output << "uuid byte: ";
	for ( unsigned char i(0); i < sizeof(uuid_.bytes_) / sizeof(uuid_.bytes_[0]); ++i ) {
		output << static_cast<int>(uuid_.bytes_[i]) << " ";
	}
	output << std::endl;

	output << "socket address family: " << socket_addr_.sin_family << std::endl;
	output << "socket address port: " << socket_addr_.sin_port << std::endl;
	output << "socket address address: " << socket_addr_.sin_addr.s_addr << std::endl;
#endif
}

std::ostream &
operator<<(std::ostream & output, UDPSocketClient const & client)
{
	client.show(output);
	return output;
}

/* -------------------------------------------------------------------------------------------------
         PyMolMover Class
   ---------------------------------------------------------------------------------------------- */

/// @brief ctor
PyMolMover::PyMolMover() :
	update_energy_(false),
	energy_type_(core::scoring::total_score),
	keep_history_(false),
	update_interval_(0),
	last_packet_sent_time_(0),
	pymol_name_()
{}

/// @breif cctor
PyMolMover::PyMolMover( PyMolMover const & other ) :
	protocols::moves::Mover( other ),
	link_( other.link_ ),
	update_energy_( other.update_energy_ ),
	energy_type_( other.energy_type_ ),
	keep_history_( other.keep_history_ ),
	update_interval_( other.update_interval_ ),
	last_packet_sent_time_( other.last_packet_sent_time_ ),
	pymol_name_( other.pymol_name_ )
{}

PyMolMover::~PyMolMover()
{}

std::string PyMolMover::get_name() const
{
	return "PyMOL_Mover";
}

std::string PyMolMover::get_PyMol_model_name(Pose const & pose) const
{
	if( pymol_name_.size() ) return pymol_name_;
    else {
		core::pose::PDBInfoCOP info = pose.pdb_info();
		if( info ) {
			std::string n = info->name();
			for(unsigned int i=0; i<n.size(); i++)
				if( n[i] == '/' ) n[i] = '_';
			return n;
		}
		else return "pose";
	}
}


bool PyMolMover::is_it_time()
{
	// First let's check if enough time have passes since last time we send info...
	//double t = clock() / CLOCKS_PER_SEC;
	double t = time(NULL);
	//TR << "t=" << t << " cl="<< clock() << std::endl;
	if( t - last_packet_sent_time_ < update_interval_ ) return false;
	last_packet_sent_time_ = t;
	return true;
}


void PyMolMover::apply( Pose const & pose)
{
	TR.Trace << "PyMolMover::apply( Pose const & pose)..." << std::endl;

	if( !is_it_time() ) return;
	TR.Trace << "PyMOL_Mover::apply It is time!" << std::endl;

	std::string name = get_PyMol_model_name(pose);

	TR.Trace << "PyMOL_Mover::apply name:" << name << std::endl;

	// Creating message...
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

	link_.sendMessage(message);

	if ( update_energy_ ) send_energy(pose, energy_type_);
}

void PyMolMover::apply( Pose & pose)
{
	Pose const & p(pose);
	apply(p);
}

void PyMolMover::print(std::string const & message)
{
	if( !is_it_time() ) return;

	std::string msg =  std::string("Text    ") + char(keep_history_) + char(0) /* Place holder for name size = 0 */ + message;

	link_.sendMessage(msg);
}

void PyMolMover::send_RAW_Energies(Pose const &pose, std::string energyType, utility::vector1<int> const & energies)
{
#ifndef  __native_client__
	if( !is_it_time() ) return;

	std::string msg(7*energies.size(), ' ');
	core::pose::PDBInfoCOP info = pose.pdb_info();
	for(unsigned int i=1; i<=energies.size(); i++) {
		char chain = ' ';
		int  res = i;
		if(info) {
			chain = info->chain(i);
			res = info->number(i);
		}
		char buf[256];
		sprintf(buf, "%c%4d%02x", chain, res, energies[i]);
		for(int k=0; k<7; k++) msg[(i-1)*7+k] = buf[k];
	}
	//TR << msg << std::endl;

	// Compressing message
	std::ostringstream zmsg;
	zlib_stream::zip_ostream zipper(zmsg, true);
	zipper << msg;
	zipper.zflush_finalize();

	std::string name = get_PyMol_model_name(pose);
	std::string sname = energyType;

	std::string message =  std::string("Ene.gzip") + char(keep_history_) \
						  + char(name.size()) + name \
						  + char(sname.size()) + sname + zmsg.str();

	//TR << "Sending message: " << message << std::endl << "Size:" << message.size() << std::endl;
	//TR << "Sending message, Size:" << message.size() << std::endl;

	link_.sendMessage(message);
#endif
}

void PyMolMover::send_energy(Pose const &pose, core::scoring::ScoreType score_type)
{
#ifndef  __native_client__
	if( !is_it_time() ) return;

	if( pose.energies().energies_updated() ) {

		utility::vector1<core::Real> e(pose.total_residue());
		core::Real min=1e100, max=1e-100;
		for(unsigned int i=1; i<=e.size(); i++) {
			if( score_type == core::scoring::total_score ) e[i] = pose.energies().residue_total_energy(i);
			else e[i] = pose.energies().residue_total_energies(i)[score_type];

			if( min > e[i] ) min = e[i];
			if( max < e[i] ) max = e[i];
		}
		// We not using send_RAW_Energies for efficiency reasons...
		std::string msg(7*e.size(), ' ');
        core::pose::PDBInfoCOP info = pose.pdb_info();
		for(unsigned int i=1; i<=e.size(); i++) {
			char chain = ' ';
			int  res = i;
            if(info) {
                chain = info->chain(i);
                res = info->number(i);
            }
			char buf[256];
			e[i] = (e[i]-min)*255. / (max-min+1e-100);
			sprintf(buf, "%c%4d%02x", chain, res, int(e[i]));
			for(int k=0; k<7; k++) msg[(i-1)*7+k] = buf[k];
		}
        //TR << msg << std::endl;

		// Compressing message
		std::ostringstream zmsg;
		zlib_stream::zip_ostream zipper(zmsg, true);
		zipper << msg;
		zipper.zflush_finalize();

		std::string name = get_PyMol_model_name(pose);
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
void PyMolMover::send_energy(Pose const &pose, std::string const & stype)
{
	send_energy(pose, core::scoring::ScoreTypeManager::score_type_from_name(stype) );
}


void PyMolMover::send_colors(Pose const &pose, std::map<int, int> const & colors, X11Colors default_color)
{
#ifndef  __native_client__
utility::vector1<int> energies( pose.total_residue(), default_color);  // energies = [ X11Colors[default_color][0] ] * pose.total_residue()

	for(std::map<int, int>:: const_iterator i = colors.begin(); i!=colors.end(); ++i) {
		PyAssert( (*i).first >=1 && (*i).first <= static_cast<int>(pose.total_residue()),
				"PyMolMover::send_colors residue index is out of range!");
		PyAssert( (*i).second >= XC_first_color && (*i).second <= XC_last_color,
				"PyMolMover::send_colors color index is out of range!");

		energies[ (*i).first ] = (*i).second;  // for r in colors: energies[r-1] = X11Colors[ colors[r] ][0]
	}
	send_RAW_Energies(pose, "X11Colors", energies);  //self._send_RAW_Energies(pose, 'X11Colors', energies, autoscale=False)
#endif
}

void PyMolMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Keep history:          " << ( ( keep_history_ ) ? ("True") : ("False") ) << std::endl;
	output << "Update energy:         " << ( ( update_energy_ ) ? ("True") : ("False") ) << std::endl;
	output << "Last packet sent time: " << last_packet_sent_time_ << std::endl;
	output << "Update interval:       " << update_interval_ << std::endl;
	output << "Link:                  " << link_ << std::endl;
}

std::ostream &
operator<<(std::ostream & output, PyMolMover const & mover)
{
	mover.show(output);
	return output;
}



void PyMolObserver::attach(core::pose::Pose &p)
{
		p.attach_general_obs(&PyMolObserver::generalEvent, this);
}

void PyMolObserver::detach(core::pose::Pose &p)
{
		p.detach_general_obs(&PyMolObserver::generalEvent, this);
}


PyMolObserverOP AddPyMolObserver(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
	//Add options
    PyMolObserverOP o = new PyMolObserver;
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	p.attach_general_obs(&PyMolObserver::generalEvent, o);
	return o;
}

PyMolObserverOP AddPyMolObserver_to_energies(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
    PyMolObserverOP o = new PyMolObserver;
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	p.attach_energy_obs(&PyMolObserver::energyEvent, o);
	return o;
}

PyMolObserverOP AddPyMolObserver_to_conformation(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
    PyMolObserverOP o = new PyMolObserver;
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	p.attach_conformation_obs(&PyMolObserver::conformationEvent, o);
	return o;
}

/// @brief PyMolMoverCreator interface, name of the mover
std::string PyMolMoverCreator::mover_name() {
  return "PyMolMover";
}

/// @brief PyMolMoverCreator interface, returns a unique key name to be used in xml file
std::string PyMolMoverCreator::keyname() const {
  return PyMolMoverCreator::mover_name();
}

/// @brief PyMolMoverCreator interface, return a new instance
protocols::moves::MoverOP PyMolMoverCreator::create_mover() const {
  return new PyMolMover();
}

/// @brief allows for the setting of certain variabel from the rosetta scripts interface, only keep history
void
PyMolMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	keep_history(tag->getOption<bool>( "keep_history", keep_history_ ) );
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PyMolMover::fresh_instance() const
{
	return new PyMolMover;
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PyMolMover::clone() const
{
	return new protocols::moves::PyMolMover( *this );
}

} // moves
} // protocols

#endif // INCLUDED_protocols_moves_PyMolMover_CC
