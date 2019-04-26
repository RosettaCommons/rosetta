// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/cloud.hh
/// @brief: RosettaCloud integration
///
/// @author Sergey Lyskov

#include <protocols/network/cloud.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/Tracer.hh>

#include <basic/options/keys/cloud.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/io/base64.hh>
#include <utility/io/zipstream.ipp>
#include <utility/thread/backwards_thread_local.hh>

#include <utility/json_utilities.hh>
//#include <utility/string_util.hh>


#define CPPHTTPLIB_ZLIB_SUPPORT
#include <httplib.h>

#ifdef MULTI_THREADED

#include <condition_variable>
#include <chrono>
#include <deque>
#include <thread>

#endif

static basic::Tracer TR( "protocols.network.cloud" );

namespace protocols {
namespace network {


using std::string;
using StringUP = std::unique_ptr<std::string>;


struct NetworkQueue
{
	~NetworkQueue();

	void add(std::string const &file_name, StringUP &&data, bool append=false);
	void add(std::string const &file_name, std::string const &data, bool append=false);

	void add(std::string const &file_name, core::pose::Pose const &);

	static NetworkQueue & instance();

private:
	NetworkQueue();

	static int get_execution_summary_id();


	// bag of data that represent individual network operation
	struct WorkUnit {
		//WorkUnit(std::string const &name, bool append) : name_(name), append_(append) {}

		std::string name;     // file name under which data should be posted
		bool append = false;  // only relevant if payload is string

		// payload: could be a Pose object or raw bytes
		core::pose::PoseUP pose;
		StringUP data;
	};

	using QueueType = std::deque<WorkUnit>;

	void add(WorkUnit&&);

	void run();

	void process_work_unit(WorkUnit & unit);

	struct ServerAddress {
		std::string host;
		std::string path;
		int port;
	};
	static ServerAddress ui_server_address();

	static std::string basic_auth_string();

private:

#ifdef MULTI_THREADED
	static int const HIGH_WATER_MARK = 128; //  maximum number of WorkUnit's to store in network queue

	std::mutex queue_mutex_;
	std::condition_variable queue_condition_variable_;
	QueueType queue_;

	bool terminate_ = false;

	//std::thread thread_;

	int const NUMBER_OF_WORKER_THREADS = 8;
	std::vector<std::thread> threads_;  // std::vector instead of std::array in case we need to specify number of worker though command-line options in the future

	//std::atomic<int> load_{0}; // load appoximation for reporting during shutdown

#endif // MULTI_THREADED
};


NetworkQueue::NetworkQueue()
{
#ifdef MULTI_THREADED
	threads_.reserve(NUMBER_OF_WORKER_THREADS);
	for ( int i = 0; i < NUMBER_OF_WORKER_THREADS; ++i ) {
		threads_.push_back( std::thread(&NetworkQueue::run, this) );
	}

	//std::atexit(atexit_handler);
#endif // MULTI_THREADED
}


NetworkQueue::~NetworkQueue()
{
#ifdef MULTI_THREADED
	{
		std::lock_guard<std::mutex> _(queue_mutex_);

		// here and below use of std::cout is intentional because usually this destructor will be called _after_ main is terminated and Tracer became unavailable
		if ( queue_.size() ) std::cout << "NetworkQueue::~NetworkQueue(): flushing " << queue_.size() << " network buffers..." << std::endl;
	}

	while ( true ) {
		{
			std::lock_guard<std::mutex> _(queue_mutex_);
			if ( queue_.empty() ) break;
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}

	{
		std::lock_guard<std::mutex> _(queue_mutex_);
		terminate_ = true;
	}
	queue_condition_variable_.notify_all();

	//if( thread_.joinable() ) thread_.join();
	for ( unsigned int i = 0; i < threads_.size(); ++i ) {
		if ( threads_[i].joinable() ) threads_[i].join();
	}

	std::cout << "NetworkQueue::~NetworkQueue(): flushing network buffers... Done." << std::endl;
#endif // MULTI_THREADED
}


NetworkQueue & NetworkQueue::instance()
{
	static THREAD_LOCAL std::unique_ptr<NetworkQueue> queue(new NetworkQueue);  // note: thread_local here is IMPORTANT because otherwise we will get into delete-race-condition during ~NetworkQueue call (attempt to use regex after main termination)
	return *queue;
}


NetworkQueue::ServerAddress NetworkQueue::ui_server_address()
{
	//return {"localhost", "/api.execution/", 64082};
	//return {"ui.graylab.jhu.edu", "/api.execution/", 80};

	using namespace basic::options;
	using namespace basic::options::OptionKeys::cloud;

	return { option[host].value(), "/api.execution", option[port].value() };
}


std::string NetworkQueue::basic_auth_string()
{
	auto login_and_password = basic::options::option[ basic::options::OptionKeys::cloud::auth ]();

	if ( login_and_password.empty() ) return "";
	else return "Basic " + utility::io::base64_encode(reinterpret_cast<unsigned char const *>( login_and_password.c_str() ), login_and_password.size());
}



void NetworkQueue::add(WorkUnit && wu)
{
#ifdef MULTI_THREADED

	static bool block = basic::options::option[ basic::options::OptionKeys::cloud::block ]();

	while ( true ) {
		bool hwm = false;

		{
			std::lock_guard<std::mutex> _(queue_mutex_);

			if ( queue_.size() < HIGH_WATER_MARK ) {
				queue_.push_back( std::move(wu) );
				break;
			} else hwm = true;
		}

		if ( hwm and  block ) {
			std::this_thread::sleep_for(std::chrono::milliseconds(64));
			//std::cerr << "Network queue is full, waiting..." << std::endl;
		} else {
			//std::cerr << "Dropping payload due queue at HIGH_WATER_MARK!" << std::endl;
			break;
		}
	}

	queue_condition_variable_.notify_one();
#else
	// If we're not multithreaded, we run the network job immediately.
	process_work_unit( wu );
#endif // MULTI_THREADED
}


void NetworkQueue::add(std::string const &file_name, StringUP &&data, bool append)
{
	WorkUnit wu;

	wu.name = file_name;
	wu.data = std::move(data);
	wu.append = append;

	add( std::move(wu) );
}


void NetworkQueue::add(std::string const &file_name, std::string const &data, bool append)
{
	WorkUnit wu;

	wu.name = file_name;
	wu.data = StringUP(new std::string(data));
	wu.append = append;

	add( std::move(wu) );
}


void NetworkQueue::add(string const &file_name, core::pose::Pose const &pose)
{
	// if( not thread_.joinable() ) thread_ = std::thread(&NetworkQueue::run, this);

	// core::pose::PoseUP p( new core::pose::Pose() );
	// p->detached_copy(pose);

	// {
	//  std::lock_guard<std::mutex> _(queue_mutex_);
	//  queue_.push_back( WorkUnit{ std::move(p) } );
	// }
	// queue_condition_variable_.notify_one();

	WorkUnit wu;

	wu.name = file_name;

	wu.pose = core::pose::PoseUP(new core::pose::Pose());
	wu.pose->detached_copy(pose);

	add( std::move(wu) );
}


int NetworkQueue::get_execution_summary_id()
{
	auto address = ui_server_address();

	httplib::Client cli(address.host.c_str(), address.port);

	nlohmann::json args;

	args["name"] = "no-name";

	auto auth = basic_auth_string();
	if ( auth.empty() ) {
		std::cout << TR.bgRed << TR.Black << "WARNING" << TR.Reset << TR.Red << " You trying to use RosettaCloud integration without providing user-name:password thorough `--cloud:auth` command line option! To get your credentials please visit https://ui.graylab.jhu.edu/settings page." << TR.Reset << std::endl;
#ifdef MULTI_THREADED
		std::this_thread::sleep_for(std::chrono::seconds(5));
#endif
		return 0;
	}

	httplib::Headers headers { {"Authorization", auth}, {"User-Agent", "Rosetta/3.7"} };

	string options = "?";

	bool clean = basic::options::option[basic::options::OptionKeys::cloud::clean].value();
	string key = basic::options::option[basic::options::OptionKeys::cloud::key].value();

	if ( clean ) options += "clean&";
	if ( !key.empty() ) options += "key=" + key + '&';

	auto res = cli.Post( (address.path + "/execution_summary" + options).c_str(), headers, args.dump(), "application/json"); // "text/plain" "application/octet-stream" .. "/_test"

	int result = 0;

	if ( res ) {
		// std::cout << "status: " << res->status << std::endl;
		// std::cout << res->body << std::endl;

		if ( res->status == 200 ) {

			auto jr = json::parse(res->body);

			// std::cout << "parsed result: " << jr.dump(2).c_str();

			utility::extract_value_if_present(jr, "execution_summary_id", result);
		} else {
			TR.Error <<  TR.bgRed << TR.Black << "Could not acquire ExecutionSummary ID, server response with code:" << res->status << " and message: " << res->body << TR.Reset << std::endl;
		}
	}

	if ( result ) TR << "Please visit " << TR.Black << TR.bgGreen << "http://" << address.host << ":" << address.port << "/execution/es/" << result << TR.Reset << " to view files posted by this run. You can also use 'magic' link http://" << address.host << ":" << address.port << "/execution/es/0 that will bind to _the latest_ run on refresh..."<< std::endl;
	else TR.Error <<  TR.bgRed << TR.Black << "Could not acquire ExecutionSummary ID!" << TR.Reset << std::endl;

	return result;
}


StringUP pose_to_string(core::pose::Pose const &pose)
{
	std::ostringstream os;
	pose.dump_pdb(os);

	return StringUP( new string( os.str() ) );
}


void NetworkQueue::run()
{
#ifdef MULTI_THREADED
	//std::cout << "NetworkQueue::run..." << std::endl;

	//QueueType queue;
	bool terminate = false;

	while ( !terminate ) {
		WorkUnit unit;

		{
			std::unique_lock<std::mutex> lock(queue_mutex_);
			queue_condition_variable_.wait(lock, [this]{ return terminate_ or !queue_.empty(); } );
			//queue.swap(queue_);
			if ( !queue_.empty() ) {
				unit = std::move( queue_.front() );
				queue_.pop_front();
			}

			terminate = terminate_;
			//load_ = queue.size();
			if ( terminate ) break;
		}

		process_work_unit(unit);

	}

	//std::cout << "NetworkQueue::run... Done!" << std::endl;
#else
	utility_exit_with_message("Attempted to use NetworkQueue::run() in a non-multithreaded environment!");
#endif
}

void
NetworkQueue::process_work_unit(WorkUnit & unit) {

	static auto execution_summary_id = get_execution_summary_id();

	if ( execution_summary_id  and  (unit.pose or unit.data) ) {
		//std::cout << "WS:" << std::this_thread::get_id() << std::endl;
		//WorkUnit & unit( queue.front() );

		auto addr = ui_server_address();

		httplib::Client cli(addr.host.c_str(), addr.port);

		httplib::Headers headers { {"Authorization", basic_auth_string()}, {"User-Agent", "Rosetta/3.7"} };

		StringUP payload;

		if ( unit.pose ) {
			payload = pose_to_string(*unit.pose);
		} else payload = std::move(unit.data);

		if ( !payload ) utility_exit_with_message("NetworkQueue::run(): Got EMPTY payload!");

		if ( payload->size() > 1024 ) { // Compressing payload when it make sense
			std::ostringstream zmsg;
			zlib_stream::zip_ostream zipper(zmsg, true);
			zipper << *payload;
			zipper.zflush_finalize();

			payload = StringUP( new string( zmsg.str() ) );

			headers.insert( {"Content-Encoding", "gzip"} );
		}

		string maybe_append = unit.append ? "?append" : "";

		auto res = cli.Post( (addr.path + "/file/" + std::to_string(execution_summary_id) + '/' + unit.name + maybe_append).c_str(),
			headers, *payload, "application/octet-stream"); // "text/plain"
	}
}

/// estimate if it resonable to expect for Cloud integration to work with current build settings and supplied commandline options
bool is_cloud_integration_enabled()
{
	// evaluating this only once to reduce overhead
	static bool r = basic::options::option[ basic::options::OptionKeys::cloud::auth ].user();
	return r;
}


void post_file(std::string const &file_name, std::string const &data, bool append)
{
	if ( not is_cloud_integration_enabled() ) return;

	NetworkQueue & queue = NetworkQueue::instance();
	queue.add(file_name, data, append);
}


void post_decoy(std::string const &file_name, core::pose::Pose const &pose)
{
	if ( not is_cloud_integration_enabled() ) return;

	NetworkQueue & queue = NetworkQueue::instance();

#ifdef MULTI_THREADED

	queue.add(file_name, pose);

#else
	static bool warn = true;

	if ( warn ) {
		TR << TR.bgRed << TR.Black << "WARNING" << TR.Reset << TR.Red << " You using RosettaCloud integration with non-multithreaded build: all Pose operation will be done in main thread which will lead to performance degradation... Please consider build with `extras=cxx11thread` instead! (if you using Ninja build: you can use source/cmake/build_cxx11thread build variant)..." << TR.Reset << std::endl;
		warn = false;
	}

	queue.add(file_name, pose_to_string(pose) );

#endif // MULTI_THREADED

}


void post_decoy(core::pose::Pose const &pose)
{
	if ( not is_cloud_integration_enabled() ) return;

	string name;

	core::pose::PDBInfoCOP info = pose.pdb_info();
	if ( info && info->name().size() ) {
		name = info->name();
		for ( char & i : name ) if ( i == '/' ) i = '_';
	} else name ="pose";

	post_decoy(name, pose);
}


// std::pair<std::string, std::string> user_password(std::string const & user_and_password)
// {
//  auto d = user_and_password.find(':');
//  if( d == std::string::npos) return std::make_pair(user_and_password, string() );
//  else return std::make_pair(user_and_password.substr(0, d), user_and_password.substr(d+1) );
// }


} // namespace network
} // namespace protocols
