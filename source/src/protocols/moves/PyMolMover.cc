// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/PyMolMover.cc
/// @brief  Send infromation to PyMol
/// @author Sergey Lyskov

#ifndef INCLUDED_protocols_moves_PyMolMover_CC
#define INCLUDED_protocols_moves_PyMolMover_CC

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
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/PyMolMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/io/zipstream.ipp>

#include <utility/py/PyAssert.hh>
#include <utility/tag/Tag.hh>

// numeric headers
#include <numeric/random/uniform.hh>

// basic headers
#include <basic/Tracer.hh>

// c++ headers
#include <ctime>

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

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.PyMolMover" );


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


UDPSocketClient::UDPSocketClient(std::string const & address, int port) : sentCount_(0)
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
	for ( unsigned short & i : uuid_.shorts_ ) i = (unsigned short) getRG()->getRandom()*65536;  //RG.random_range(0, 65536);

	memset(&socket_addr_, '\0', sizeof(sockaddr_in));

	socket_addr_.sin_family = AF_INET;     // host byte order
	socket_addr_.sin_port = htons(port);  // short, network byte order
	socket_addr_.sin_addr.s_addr = inet_addr( address.c_str() );

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
	if ( msg.size() > max_packet_size_ )  { count = msg.size()/max_packet_size_ + 1; }

	for ( int i=0; i<count; i++ ) {
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
	output << "socket handel: " << socket_h_ << std::endl;

	output << "uuid short: ";
	for ( unsigned short i : uuid_.shorts_ ) {
		output << i << " ";
	}
	output << std::endl;

	output << "uuid byte: ";
	for ( char byte : uuid_.bytes_ ) {
		output << static_cast<int>(byte) << " ";
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
PyMolMover::PyMolMover(std::string const & address, int port) :
	link_(address, port),
	update_energy_(false),
	energy_type_(core::scoring::total_score),
	update_membrane_(false),
	keep_history_(false),
	update_interval_(0),
	last_packet_sent_time_(0),
	pymol_name_()
{}

/// @brief cctor
PyMolMover::PyMolMover( PyMolMover const & ) = default;

PyMolMover::~PyMolMover() = default;

std::string PyMolMover::get_name() const
{
	return "PyMOL_Mover";
}

void PyMolMover::set_PyMol_model_name( std::string name ){
	pymol_name_ = name;
}

std::string PyMolMover::get_PyMol_model_name(Pose const & pose) const
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


bool PyMolMover::is_it_time()
{
	// First let's check if enough time have passes since last time we send info...
	//double t = clock() / CLOCKS_PER_SEC;
	double t = time(nullptr);
	//TR << "t=" << t << " cl="<< clock() << std::endl;
	if ( t - last_packet_sent_time_ < update_interval_ ) return false;
	last_packet_sent_time_ = t;
	return true;
}


void PyMolMover::apply( Pose const & pose)
{
	TR.Trace << "PyMolMover::apply( Pose const & pose) ..." << std::endl;

	if ( !is_it_time() ) return;
	TR.Trace << "PyMOL_Mover::apply It is time!" << std::endl;

	std::string name = get_PyMol_model_name(pose);
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

	link_.sendMessage(message);

	if ( update_membrane_ ) send_membrane_planes(pose);
	if ( update_energy_ ) send_energy(pose, energy_type_);
}

void PyMolMover::apply( Pose & pose)
{
	Pose const & p(pose);
	apply(p);
}

void PyMolMover::print(std::string const & message)
{
	if ( !is_it_time() ) return;

	std::string msg =  std::string("Text    ") + char(keep_history_) + char(0) /* Place holder for name size = 0 */ + message;

	link_.sendMessage(msg);
}

void PyMolMover::send_RAW_Energies(Pose const &pose, std::string energyType, utility::vector1<int> const & energies)
{
#ifndef  __native_client__
	if ( !is_it_time() ) return;

	std::string msg(8*energies.size(), ' ');
	core::pose::PDBInfoCOP info = pose.pdb_info();
	for ( unsigned int i=1; i<=energies.size(); i++ ) {
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
		for ( int k=0; k<8; k++ ) msg[(i-1)*8+k] = buf[k];
	}
	//TR << msg << std::endl;

	// Compressing message
	std::ostringstream zmsg;
	zlib_stream::zip_ostream zipper(zmsg, true);
	zipper << msg;
	zipper.zflush_finalize();

	std::string name = get_PyMol_model_name(pose);
	std::string sname = energyType;

	std::string message = std::string("Ene.gzip") + char(keep_history_) \
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
	if ( !is_it_time() ) return;

	if ( pose.energies().energies_updated() ) {

		utility::vector1<core::Real> e(pose.total_residue());
		core::Real min=1e100, max=1e-100;
		for ( unsigned int i=1; i<=e.size(); i++ ) {
			if ( score_type == core::scoring::total_score ) e[i] = pose.energies().residue_total_energy(i);
			else e[i] = pose.energies().residue_total_energies(i)[score_type];

			if ( min > e[i] ) min = e[i];
			if ( max < e[i] ) max = e[i];
		}
		// We not using send_RAW_Energies for efficiency reasons...
		std::string msg(8*e.size(), ' ');
		core::pose::PDBInfoCOP info = pose.pdb_info();
		for ( unsigned int i=1; i<=e.size(); i++ ) {
			char chain = ' ';
			char icode = ' ';
			int  res = i;
			if ( info ) {
				chain = info->chain(i);
				icode = info->icode(i);
				res = info->number(i);
			}
			char buf[256];
			e[i] = (e[i]-min)*255. / (max-min+1e-100);
			sprintf(buf, "%c%4d%c%02x", chain, res, icode, int(e[i]));
			for ( int k=0; k<8; k++ ) msg[(i-1)*8+k] = buf[k];
		}

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

/// @brief Send Membrane Planes to PyMol
/// @details If pose is a membrane pose
/// pymol viewer will build CGO planes from points specified
void PyMolMover::send_membrane_planes( Pose const & pose ) {

#ifndef __native_client__

	using namespace core::scoring::methods;

	if ( !is_it_time() ) return;

	// Check the membrane planes can be visualized
	if ( !pose.conformation().is_membrane() ) return;

	// Grab a list of relevant residues and go
	// Compute radius of gyration of the pose
	utility::vector1< bool > relevant_residues;
	relevant_residues.resize( pose.total_residue() );
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

	std::string name = get_PyMol_model_name(pose);
	std::string sname = "membrane_planes";

	std::string message =  std::string("Mem.gzip") + char(keep_history_) \
		+ char(name.size()) + name \
		+ char(sname.size()) + sname + zmsg.str();

	link_.sendMessage(message);

#endif

}


void PyMolMover::send_colors(Pose const &pose, std::map<int, int> const & colors, X11Colors default_color)
{
#ifndef  __native_client__
	utility::vector1<int> energies( pose.total_residue(), default_color);  // energies = [ X11Colors[default_color][0] ] * pose.total_residue()

	 for ( auto const & color : colors ) {
		PyAssert( color.first >=1 && color.first <= static_cast<int>(pose.total_residue()),
			"PyMolMover::send_colors residue index is out of range!");
		PyAssert( color.second >= XC_first_color && color.second <= XC_last_color,
			"PyMolMover::send_colors color index is out of range!");

		energies[ color.first ] = color.second;  // for r in colors: energies[r-1] = X11Colors[ colors[r] ][0]
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

PyMolObserver::PyMolObserver():
	CacheableObserver(),
	type_( no_observer ) // We have to set the observer type specifically.
{
}

PyMolObserver::PyMolObserver(PyMolObserver const & rval) :
	CacheableObserver( rval ),
	type_( rval.type_ ),
	pymol_(rval.pymol_)
	// Do NOT copy the *_event_link_s
{
}

PyMolObserver::~PyMolObserver() {
	detach_from();
}

PyMolObserver &
PyMolObserver::operator= (PyMolObserver const &rval) {
	if ( this != &rval ) {
		core::pose::datacache::CacheableObserver::operator=( rval );
		type_ = rval.type_;
		pymol_ = rval.pymol_;
		// Do NOT copy the *_event_link_s
	}
	return *this;
}

core::pose::datacache::CacheableObserverOP
PyMolObserver::clone() {
	return core::pose::datacache::CacheableObserverOP( new PyMolObserver( *this ) );
}

core::pose::datacache::CacheableObserverOP
PyMolObserver::create() {
	return core::pose::datacache::CacheableObserverOP( new PyMolObserver );
}

void
PyMolObserver::set_type( ObserverType setting ) {
	type_ = setting;
	// We don't have a pose, so wait until we get attached
}

void
PyMolObserver::add_type( ObserverType setting ) {
	type_ = type_ | setting;
	// We don't have a pose, so wait until we get attached
}

void PyMolObserver::attach(core::pose::Pose &p)
{
	attach_to(p);
}

void PyMolObserver::detach(core::pose::Pose & /*p*/)
{
	detach_from();
}

bool
PyMolObserver::is_attached() const {
	return general_event_link_.valid() || energy_event_link_.valid() || conformation_event_link_.valid();
}

void
PyMolObserver::attach_impl(core::pose::Pose & pose) {
	general_event_link_.invalidate();
	energy_event_link_.invalidate();
	conformation_event_link_.invalidate();

	if ( type_ & general_observer ) {
		general_event_link_ = pose.attach_general_obs( &PyMolObserver::generalEvent, this );
	}
	if ( type_ & energy_observer ) {
		energy_event_link_ = pose.attach_energy_obs( &PyMolObserver::energyEvent, this );
	}
	if ( type_ & conformation_observer ) {
		conformation_event_link_ = pose.attach_conformation_obs( &PyMolObserver::conformationEvent, this );
	}
}

void
PyMolObserver::detach_impl() {
	general_event_link_.invalidate();
	energy_event_link_.invalidate();
	conformation_event_link_.invalidate();
}

PyMolObserverOP
get_pymol_observer(core::pose::Pose & pose) {
	using namespace core::pose::datacache;

	if ( !pose.observer_cache().has( CacheableObserverType::PYMOL_OBSERVER ) ) {
		PyMolObserverOP obs( new PyMolObserver );
		pose.observer_cache().set( CacheableObserverType::PYMOL_OBSERVER, obs, /*autoattach*/ false );
	}
	CacheableObserverOP obs = pose.observer_cache().get_ptr( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER );
	return utility::pointer::dynamic_pointer_cast< PyMolObserver >( obs );
}

PyMolObserverOP AddPyMolObserver(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
	PyMolObserverOP o( get_pymol_observer(p) );
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	o->add_type( PyMolObserver::general_observer );
	o->attach(p);
	return o;
}

PyMolObserverOP AddPyMolObserver_to_energies(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
	PyMolObserverOP o( get_pymol_observer(p) );
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	o->add_type( PyMolObserver::energy_observer );
	o->attach(p);
	return o;
}

PyMolObserverOP AddPyMolObserver_to_conformation(core::pose::Pose &p, bool keep_history, core::Real update_interval)
{
	PyMolObserverOP o( get_pymol_observer(p) );
	o->pymol().keep_history(keep_history);
	o->pymol().update_interval(update_interval);
	o->add_type( PyMolObserver::conformation_observer );
	o->attach(p);
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
	return protocols::moves::MoverOP( new PyMolMover() );
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

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PyMolMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PyMolMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PyMolMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::moves::PyMolMover( *this ) );
}

} // moves
} // protocols

#endif // INCLUDED_protocols_moves_PyMolMover_CC
