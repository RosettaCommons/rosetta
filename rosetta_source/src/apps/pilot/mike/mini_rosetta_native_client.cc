/// @file minirosetta.cc
/// This example demonstrates loading, running and scripting a very simple NaCl
/// module.  To load the NaCl module, the browser first looks for the
/// CreateModule() factory method (at the end of this file).  It calls
/// CreateModule() once to load the module code from your .nexe.  After the
/// .nexe code is loaded, CreateModule() is not called again.
///
/// Once the .nexe code is loaded, the browser than calls the CreateInstance()
/// method on the object returned by CreateModule().  It calls CreateInstance()
/// each time it encounters an <embed> tag that references your NaCl module.
///
/// When the browser encounters JavaScript that references your NaCl module, it
/// calls the GetInstanceObject() method on the object returned from
/// CreateInstance().  In this example, the returned object is a subclass of
/// ScriptableObject, which handles the scripting support.

#include <ppapi/cpp/dev/scriptable_object_deprecated.h>
#include <ppapi/cpp/instance.h>
#include <ppapi/cpp/module.h>
#include <ppapi/cpp/var.h>

#include <cstdio>
#include <string>

// libRosetta headers
#include <basic/options/option.hh>
#include <core/init.hh>
// C++ headers
#include <iostream>
#include <string>
#include <basic/Tracer.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//Auto Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/relax_main.hh>
#include <protocols/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <ObjexxFCL/string.functions.hh>

using namespace ObjexxFCL;

#include <pthread.h>

#include <ppapi/cpp/url_loader.h>
#include <ppapi/cpp/url_request_info.h>
#include <ppapi/cpp/url_response_info.h>
#include <ppapi/c/pp_stdint.h>
#include <ppapi/c/pp_errors.h>
#include <ppapi/cpp/resource.h>
#include <ppapi/cpp/completion_callback.h>
#include <ppapi/cpp/dev/file_io_dev.h>

const std::string basic_url = "http://ralph.bakerlab.org/rah/temp/minirosetta_database";

// GetURLHandler is used to download data from |url|. When download is
// finished or when an error occurs, it calls JavaScript function
// reportResult(url, result, success) (defined in geturl.html) and
// self-destroys.
//
// EXAMPLE USAGE:
// GetURLHandler* handler* = GetURLHandler::Create(instance,url);
// if (!handler->Start())
//   delete handler;
//
//
//class GetURLHandler {
// public:
//  // Creates instance of GetURLHandler on the heap.
//  // GetURLHandler objects shall be created only on the heap (they
//  // self-destroy when all data is in).
//  static GetURLHandler* Create(pp::Instance* instance_);
//
//	static GetURLHandler* get_handler();
//
//  // Initiates page (URL) download.
//  // Returns false in case of internal error, and self-destroys.
// 	bool request_url( const std::string & url );
//
//
//  bool Start( );
//
//	std::string get_data() const { return url_response_body_; }
//
//	bool is_done() const { return is_done_; }
//
// public:
//  static const int kBufferSize = 4096;
//
//  GetURLHandler(pp::Instance* instance_ );
//  ~GetURLHandler();
//
//  // Callback fo the pp::URLLoader::Open().
//  // Called by pp::URLLoader when response headers are received or when an
//  // error occurs (in response to the call of pp::URLLoader::Open()).
//  // Look at <ppapi/c/ppb_url_loader.h> and
//  // <ppapi/cpp/url_loader.h> for more information about pp::URLLoader.
//  void OnOpen(int32_t result);
//
//  // Callback fo the pp::URLLoader::ReadResponseBody().
//  // |result| contains the number of bytes read or an error code.
//  // Appends data from this->buffer_ to this->url_response_body_.
//  void OnRead(int32_t result);
//
//  // Reads the response body (asynchronously) into this->buffer_.
//  // OnRead() will be called when bytes are received or when an error occurs.
//  void ReadBody();
//
//  // Calls JS reportResult() defined in geturl.html
//  void ReportResult(const std::string& fname,
//                    const std::string& text,
//                    bool success);
//  // Calls JS reportResult() and self-destroy (aka delete this;).
//  void ReportResultAndDie(const std::string& fname,
//                          const std::string& text,
//                          bool success);
//	bool is_done_;
//	PP_Instance instance_id_;
//  std::string url_;  // URL to be downloaded.
//  pp::URLRequestInfo url_request_;
//  pp::URLLoader url_loader_;  // URLLoader provides an API to download URLs.
//  char buffer_[kBufferSize];  // buffer for pp::URLLoader::ReadResponseBody().
//  std::string url_response_body_;  // Contains downloaded data.
//  pp::CompletionCallbackFactory<GetURLHandler> cc_factory_;
//
//  GetURLHandler(const GetURLHandler&);
//  void operator=(const GetURLHandler&);
//
//	static GetURLHandler* urlhandler_instance_;
//
//};
//
//GetURLHandler* GetURLHandler::urlhandler_instance_ = NULL;
//
//
//bool IsError(int32_t result) {
//  return ((PP_OK != result) && (PP_ERROR_WOULDBLOCK != result));
//}
//
//
//GetURLHandler* GetURLHandler::Create(pp::Instance* instance) {
//  if(!urlhandler_instance_){ urlhandler_instance_ = new GetURLHandler(instance); }
//	return urlhandler_instance_;
//}
//
//GetURLHandler* GetURLHandler::get_handler() {
//	return urlhandler_instance_;
//}
//
//GetURLHandler::GetURLHandler(pp::Instance* instance)
//    : is_done_(true),
//			instance_id_(instance->pp_instance()),
//			url_(""),
//      url_request_(instance),
//      url_loader_(instance),
//      cc_factory_(this) {
//}
//
//GetURLHandler::~GetURLHandler() {
//}
//
//bool GetURLHandler::request_url( const std::string & url){
//	if( url_ != "" ) return false; // re're not ready.
//	url_ = url;
//	return true;
//}
//
//bool GetURLHandler::Start(){
//  if( url_ == "" ) return false; // check if we're supposed to get something.
//
//	is_done_ = false;
//  url_request_.SetURL(url_);
//  url_request_.SetMethod("GET");
//  pp::CompletionCallback cc = cc_factory_.NewCallback(&GetURLHandler::OnOpen);
//  int32_t res = url_loader_.Open(url_request_, cc);
//  if (PP_ERROR_WOULDBLOCK != res)
//    cc.Run(res);
//  return !IsError(res);
//}
//
//void GetURLHandler::OnOpen(int32_t result) {
//  if (result < 0)
//    ReportResultAndDie(url_, "pp::URLLoader::Open() failed", false);
//  else
//    ReadBody();
//}
//
//void GetURLHandler::OnRead(int32_t result) {
//  if (result < 0) {
//    ReportResultAndDie(url_,
//                       "pp::URLLoader::ReadResponseBody() result<0",
//                       false);
//
//  } else if (result != 0) {
//    int32_t num_bytes = result < kBufferSize ? result : sizeof(buffer_);
//    url_response_body_.reserve(url_response_body_.size() + num_bytes);
//    url_response_body_.insert(url_response_body_.end(),
//                              buffer_,
//                              buffer_ + num_bytes);
//    ReadBody();
//  } else {  // result == 0, end of stream
//    ReportResultAndDie(url_, url_response_body_, true);
//  }
//}
//
//void GetURLHandler::ReadBody() {
//  // Reads the response body (asynchronous) into this->buffer_.
//  // OnRead() will be called when bytes are received or when an error occurs.
//  // Look at <ppapi/c/dev/ppb_url_loader> for more details.
//  pp::CompletionCallback cc = cc_factory_.NewCallback(&GetURLHandler::OnRead);
//  int32_t res = url_loader_.ReadResponseBody(buffer_,
//                                             sizeof(buffer_),
//                                             cc);
//  if (PP_ERROR_WOULDBLOCK != res)
//    cc.Run(res);
//}
//
//void GetURLHandler::ReportResultAndDie(const std::string& fname,
//                                       const std::string& text,
//                                       bool success) {
//  ReportResult(fname, text, success);
//  delete this;
//}
//
//void GetURLHandler::ReportResult(const std::string& fname,
//                                 const std::string& text,
//                                 bool success) {
//  is_done_ = true;
//	if (success)
//    printf("GetURLHandler::ReportResult(Ok).\n");
//  else
//    printf("GetURLHandler::ReportResult(Err). %s\n", text.c_str());
//  fflush(stdout);
//  pp::Module* module = pp::Module::Get();
//  if (NULL == module)
//    return;
//
//  pp::Instance* instance = module->InstanceForPPInstance(instance_id_);
//  if (NULL == instance)
//    return;
//
//  pp::Var window = instance->GetWindowObject();
//  // calls JavaScript function reportResult(url, result, success)
//  // defined in geturl.html.
//  pp::Var exception;
//  window.Call("reportResult", fname, text, success, &exception);
//}
//


////class LoadFileFromURLStorage {
////	public:
////		LoadFileFromURLStorage( pp::Instance *instance ):
////			loader_( instance ),
////			pp::URLRequestInfo request;
////			done_(false),
////			did_open_(false)
////		{
////
////		}
////
////
////		// note that this is a serial implementation
////		std::string read_file( std::string filename ){
////			done_     = false;
////			did_open_ = false;
////    	result_ = "";
////
////			std::string the_url =  basic_url + "/" + filename;
////			the_url = "geturl_success.html";
////			request.SetURL( the_url );
////      request.SetMethod("GET");
////      request.SetFollowRedirects(true);
////
////			pp::CompletionCallback cc = NewCallback();
////      int32_t rv = loader_.Open(request, cc);
////      result_ = result_ + string_of( rv );
////			if (rv != PP_ERROR_WOULDBLOCK){
////  			result_ = "NOTWOULDBLOCK:" + result_;
////  			cc.Run(rv);
////  		}
////
////			// wait until done
////  		int count = 10;
////			while( !done_ && count > 0 ){ sleep(1); count --; };
////
////
////			// lets try it the blocking way.
////
//////			int32_t rv = loader_.Open(request, NULL );
//////      result_ = result_ + "," + string_of( rv );
//////		 	DidCompleteIO(rv);
////
////			if( done_ ){ result_ = "RRDONE:" + result_; }
////			return result_;
////		}
////
////	private:
////    pp::CompletionCallback NewCallback() {
////       return factory_.NewCallback(&LoadFileFromURLStorage::DidCompleteIO);
////    }
////
////		 void DidCompleteIO(int32_t result) {
////       if (result > 0) {
////         // buf_ now contains 'result' number of bytes from the URL.
////         ProcessBytes(buf_, result);
////         ReadMore();
////     		result_ = result_ + "E";
////       } else if (result == PP_OK && !did_open_) {
////         // Headers are available, and we can start reading the body.
////         did_open_ = true;
////       	ProcessResponseInfo(loader_.GetResponseInfo());
////       	ReadMore();
////     		result_ = result_ + "J";
////       } else {
////       	// Done reading (possibly with an error given by 'result').
////       	done_ = true;
////     		result_ = result_ + "K";
////			 }
////     		result_ = result_ + "X";
////     }
////     void ReadMore() {
////       pp::CompletionCallback cc = NewCallback();
////       int32_t rv = fio_.Read(offset_, buf_, sizeof(buf_), cc);
////       if (rv != PP_ERROR_WOULDBLOCK) cc.Run(rv);
////     }
////     void ProcessResponseInfo(const pp::URLResponseInfo& response_info) {
////       // Read response headers, etc.
////			 // nothing here - just raw read whatever
////     		result_ = result_ + "PREPROCESS";
////     }
////     void ProcessBytes(const char* bytes, int32_t length) {
////       // Do work ...
////     		result_ = result_ + "<" +  std::string( bytes, length ) + ">";
////		 }
////
////		pp::URLLoader loader_;
////    pp::CompletionCallbackFactory<LoadFileFromURLStorage> factory_;
////    pp::FileIO_Dev fio_;
////		char buf_[4096];
////		int offset_;
////		bool done_;
////    bool did_open_;
////		std::string result_;
////};



/// These are the method names as JavaScript sees them.  Add any methods for
/// your class here.
namespace ncmini {
// A method consists of a const char* for the method ID and the method's
// declaration and implementation.
// TODO(sdk_user): 1. Add the declarations of your method IDs.

const char* const kRunRosettaMethodId = "runRosetta";
const char* const kGetStdoutMethodId = "getStdout";


//////
//////
////// class RunRosettaSingleton {
//////  private:
//////   RunRosettaSingleton():
////// 		is_running_(false),
////// 		command_line_(""),
////// 		pp_instance_(NULL)
////// 	{
////// 	}
//////
////// 	~RunRosettaSingleton(){
////// 	}
//////  public:
////// 	static void* run( void * void_the_instance) {
////// 		using namespace core;
////// 		using namespace protocols;
////// 		using namespace protocols::jobdist;
////// 		using namespace protocols::moves;
////// 		using namespace basic::options;
//////
////// 		//pp::Instance* the_instance = (pp::Instance*) void_the_instance;
//////
////// 		RunRosettaSingleton *rosetta = RunRosettaSingleton::get_instance();
////// 		rosetta->is_running_ = true;
//////
////// 		rosetta->ss_cout << "RunRosetta Native Client version of rosetta\n";
////// 		rosetta->ss_cout << "Arguments: '" + rosetta->command_line() + "'\n";
////// 		rosetta->ss_cout << std::endl;
////// 		rosetta->ss_cout.flush();
//////
//////
//////
//////
//////
////// //	std::string the_url = "geturl_success.html";
////// //  pp::Var window = the_instance->GetWindowObject();
////// //  pp::Var exception;
////// //  window.Call("getURL", the_url, &exception);
//////
////// 		// test the webbased file loading mechanism
////// 		// http://ralph.bakerlab.org/rah/temp/minirosetta_database/scoring/weights/score12_full.wts
//////
////// 		// why this fails i have no idea.
////// 		//std::string the_url = "geturl_success.html";
////// 		//GetURLHandler *handler = GetURLHandler::get_handler();
////// 		//handler->request_url( the_url );
//////
////// 		//if ( (handler != NULL) && (!handler->is_done() ) ){
////// 	//		if (!handler->Start(the_url)) {
////// 	//			rosetta->ss_cout << "ERROR: GetURLHandler::Start failed" << std::endl;
////// 	//		}
////// //			//while( !handler->is_done() ){ sleep(1); }
////// //			rosetta->ss_cout << "DATA:" << handler->get_data() << std::endl;
////// //		} else {
////// //			rosetta->ss_cout << "ERROR: GetURLHandler::Init failed" << std::endl;
//////
////// 		//rosetta->handler = GetURLHandler::Create(rosetta->pp_instance(), the_url );
//////   	//GetURLHandler*  test_handler = new GetURLHandler(rosetta->pp_instance(), the_url);
////// 		//delete test_handler;
//////
//////
////// //		if (rosetta->handler != NULL) {
////// 			// Starts asynchronous download. When download is finished or when an
////// 			// error occurs, |handler| calls JavaScript function
////// 			// reportResult(url, result, success) (defined in geturl.html) and
////// 			// self-destroys.
////// /////			if (!rosetta->handler->Start()) {
////// /////				rosetta->ss_cout << "ERROR: GetURLHandler::Start failed" << std::endl;
////// /////			}
////// /////
//////
////// 			//while( !rosetta->handler->is_done() ){ sleep(1); }
//////
//////
////// 			//rosetta->ss_cout << "DATA:" << rosetta->handler->get_data() << std::endl;
////// 	//	} else {
////// 	//		//rosetta->ss_cout << "ERROR: GetURLHandler::Create failed" << std::endl;
////// 	//	}
////// //////
////// //////		//LoadFileFromURLStorage loader( instance );
////// //////		//rosetta->ss_cout << "LOADFILE: " << loader.read_file( "scoring/weights/score12_full.wts"  ) << std::endl;
////// //////
////// 		// build c-style arguments.
////// 		std::vector< std::string > tokens ( utility::split( rosetta->command_line() ) );
////// 		utility::vector1<std::string> tokens_v1;
////// 		basic::final_channel = &rosetta->ss_cout;
//////
////// 		tokens_v1.push_back( "mini_rosetta_native_client.nexe" );
////// 		for( std::vector< std::string >::const_iterator it=tokens.begin(); it != tokens.end(); ++ it ){
////// 			tokens_v1.push_back( * it );
////// 			(*basic::final_channel) << (*it) << std::endl;
////// 		}
//////
//////
////// 		rosetta->ss_cout << "start" << std::endl;
////// 		try{
////// 			core::init(tokens_v1);
////// 			relax::ClassicRelax::register_options();
////// 			jd2::register_options();
////// 			option.add_relevant( OptionKeys::in::file::fullatom );
////// 			option.add_relevant( OptionKeys::relax::fast );
////// 			relax::Relax_main( false );
////// 		}
////// 		catch( std::string ex ){
////// 			rosetta->ss_cout << "Exception occured: " << ex << std::endl;
////// 		}
////// 		catch( ... ){
////// 			rosetta->ss_cout << "Unknown exception occured" << std::endl;
////// 		}
//////
//////
////// 		rosetta->ss_cout << "DONE" << std::endl;
//////
////// 		rosetta->is_running_ = false;
//////
//////
////// 		pthread_exit(0);
////// 		return (void*)rosetta;
////// 	}
//////
////// public:
//////
////// 	static RunRosettaSingleton *get_instance(){
////// 		if( rosetta_instance_ == NULL ) rosetta_instance_ = new RunRosettaSingleton();
////// 		return rosetta_instance_;
////// 	}
//////
////// 	void set_pp_instance( pp::Instance* pp_instance ){
////// 		pp_instance_ = pp_instance;
////// 	}
//////
////// 	void set_command_line( std::string command_line ){
////// 		command_line_ = command_line;
////// 	}
//////
////// 	pp::Instance* pp_instance() const { return pp_instance_; }
//////
////// 	const std::string& command_line() const { return command_line_; }
//////
////// 	bool is_running() const { return is_running_; }
//////
////// 	std::string get_stdout() { return ss_cout.str(); }
//////
////// 	pthread_t tid;
//////
////// 		std::stringstream ss_cout;
////// private: // data
////// 		static RunRosettaSingleton *rosetta_instance_;
////// 		bool is_running_;
////// 		std::string command_line_;
////// 		pp::Instance* pp_instance_;
//////
////// 		std::stringstream ss_cerr;
////// 		GetURLHandler* handler;
////// };
//////
////// RunRosettaSingleton *RunRosettaSingleton::rosetta_instance_ = NULL;





std::string global_command_line;

std::stringstream ss_cout;
pthread_t tid;

pp::Var MarshallGetStdout(  pp::Instance* instance, const std::vector<pp::Var>& args ) {
//	RunRosettaSingleton *rosetta = RunRosettaSingleton::get_instance();

	// initiate any URL requests pending;
//	GetURLHandler* handler = GetURLHandler::Create(instance);
//	if ( (handler != NULL) && (handler->is_done() ) ){
//		if (!handler->Start()) {
//			rosetta->ss_cout << "ERROR: GetURLHandler::Start failed or has no URL to retrieve" << std::endl;
//		}
//		//while( !handler->is_done() ){ sleep(1); }
//		rosetta->ss_cout << "DATA:" << handler->get_data() << std::endl;
//	} else {
//		rosetta->ss_cout << "ERROR: GetURLHandler::Init failed" << std::endl;
//	}

  return pp::Var( std::string( ss_cout.str() ) );
}


void * run( void * data ){

	using namespace core;
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace basic::options;

	std::string command_line = global_command_line;

	//pp::Instance* the_instance = (pp::Instance*) void_the_instance;
	ss_cout << "RunRosetta Native Client version of rosetta\n";
	ss_cout << "Arguments: '" + command_line + "'\n";

	// build c-style arguments.
	std::vector< std::string > tokens ( utility::split( command_line ) );
	utility::vector1<std::string> tokens_v1;
	basic::final_channel = &ss_cout;

	tokens_v1.push_back( "mini_rosetta_native_client.nexe" );
	for( std::vector< std::string >::const_iterator it=tokens.begin(); it != tokens.end(); ++ it ){
		tokens_v1.push_back( * it );
		(*basic::final_channel) << (*it) << std::endl;
	}

	ss_cout << "start" << std::endl;
	try{
		relax::ClassicRelax::register_options();
		jd2::register_options();
		option.add_relevant( OptionKeys::in::file::fullatom );
		option.add_relevant( OptionKeys::relax::fast );
		core::init(tokens_v1);
		relax::Relax_main( false );
	}
	catch( std::string ex ){
		ss_cout << "Exception occured: " << ex << std::endl;
	}
	catch( ... ){
		ss_cout << "Unknown exception occured" << std::endl;
	}

	ss_cout << "DONE" << std::endl;

	pthread_exit(0);
	return NULL;
}




// TODO(sdk_user): 2. Implement the methods that correspond to your method IDs.
/// This is the module's function that invokes FortyTwo and converts the return
/// value from an int32_t to a pp::Var for return.
pp::Var MarshallRunRosetta(  pp::Instance* instance, const std::vector<pp::Var>& args ) {

	//RunRosettaSingleton *rosetta = RunRosettaSingleton::get_instance();
	// interpret args and convert to rosetta argument string here
	//GetURLHandler::Create(instance);

	//std::string the_url = "geturl_success.html";
  //pp::Var window = instance->GetWindowObject();
  //pp::Var exception;
  //window.Call("getURL", the_url, &exception);

//	// make a RunRosetta Singleton and load it up with its paramters.
//	if( !rosetta->is_running() ){
//		rosetta->set_pp_instance( instance );
//		rosetta->set_command_line(  args[0].AsString() );
//
//		// now create a thread
//		//pthread_create(&rosetta->tid, NULL, RunRosettaSingleton::run, (void*) instance );
//
//		RunRosettaSingleton::run( (void*) instance  );
//	}else{
//  	return pp::Var( std::string("Already running..") );
//	}
	global_command_line = args[0].AsString();
	pthread_create(&tid, NULL, run, NULL );
  return pp::Var(  std::string("The Output") );
}

// Note to the user: This glue code reflects the current state of affairs.  It
// may change.  In particular, interface elements marked as deprecated will
// disappear sometime in the near future and replaced with more elegant
// interfaces.  As of the time of this writing, the new interfaces are not
// available so we have to provide this code as it is written below.

/// This class exposes the scripting interface for this NaCl module.  The
/// HasMethod method is called by the browser when executing a method call on
/// the object.  The name of the JavaScript function (e.g. "fortyTwo") is
/// passed in the |method| paramter as a string pp::Var.  If HasMethod()
/// returns |true|, then the browser will call the Call() method to actually
/// invoke the method.
class MinirosettaScriptableObject : public pp::deprecated::ScriptableObject {
 public:
   explicit MinirosettaScriptableObject(pp::Instance* instance)
	       : instance_(instance) {}

	/// Called by the browser to decide whether @a method is provided by this
  /// plugin's scriptable interface.
  /// @param[in] method The name of the method
  /// @param[out] exception A pointer to an exception.  May be used to notify
  ///     the browser if an exception occurs.
  /// @return true iff @a method is one of the exposed method names.
  virtual bool HasMethod(const pp::Var& method, pp::Var* exception);

  /// Invoke the function associated with @a method.  The argument list passed
  /// in via JavaScript is marshalled into a vector of pp::Vars.  None of the
  /// functions in this example take arguments, so this vector is always empty.
  /// @param[in] method The name of the method to be invoked.
  /// @param[in] args The arguments to be passed to the method.
  /// @param[out] exception A pointer to an exception.  May be used to notify
  ///     the browser if an exception occurs.
  /// @return true iff @a method was called successfully.
  virtual pp::Var Call(const pp::Var& method,
                       const std::vector<pp::Var>& args,
                       pp::Var* exception);

 private:
	pp::Instance* instance_;
};

bool MinirosettaScriptableObject::HasMethod(const pp::Var& method,
                                             pp::Var* exception) {
  if (!method.is_string()) {
    return false;
  }
  std::string method_name = method.AsString();
  // TODO(sdk_user): 3. Make this function return true iff method_name is equal
  // to any of your method IDs.
  bool has_method = false;

	if( method_name == kRunRosettaMethodId) has_method = true;
	if( method_name == kGetStdoutMethodId) has_method = true;
	return has_method;
}

pp::Var MinirosettaScriptableObject::Call(const pp::Var& method,
                                           const std::vector<pp::Var>& args,
                                           pp::Var* exception) {
  if (!method.is_string()) {
    return pp::Var();
  }
  std::string method_name = method.AsString();
  // TODO(sdk_user): 4. Make this function call whatever method has method_name
  // as its method ID.
	if( method_name == kRunRosettaMethodId){
		return MarshallRunRosetta(instance_, args);
	}else
	if( method_name == kGetStdoutMethodId){
		return MarshallGetStdout(instance_, args);
	}
  return pp::Var();
}

/// The Instance class.  One of these exists for each instance of your NaCl
/// module on the web page.  The browser will ask the Module object to create
/// a new Instance for each occurence of the <embed> tag that has these
/// attributes:
///     type="application/x-nacl"
///     nexes="ARM: minirosetta_arm.nexe
///            ..."
/// The Instance can return a ScriptableObject representing itself.  When the
/// browser encounters JavaScript that wants to access the Instance, it calls
/// the GetInstanceObject() method.  All the scripting work is done though
/// the returned ScriptableObject.
class MinirosettaInstance : public pp::Instance {
 public:
  /// The constructor creates the plugin-side instance.
  /// @param[in] instance the handle to the browser-side plugin instance.
  explicit MinirosettaInstance(PP_Instance instance) : pp::Instance(instance)
  {}
  virtual ~MinirosettaInstance() {}

  /// The browser calls this function to get a handle, in form of a pp::Var,
  /// to the plugin-side scriptable object.  The pp::Var takes over ownership
  /// of said scriptable, meaning the browser can call its destructor.  The
  /// MinirosettaScriptableObject is the plugin-side representation of that
  /// scriptable object.
  /// @return The browser's handle to the plugin side instance.
  virtual pp::Var GetInstanceObject() {
    MinirosettaScriptableObject* hw_object =
        new MinirosettaScriptableObject(this);
    return pp::Var(this, hw_object);
  }
};

/// The Module class.  The browser calls the CreateInstance() method to create
/// an instance of your NaCl module on the web page.  The browser creates a new
/// instance for each <embed> tag with type="application/x-nacl".
class MinirosettaModule : public pp::Module {
 public:
  MinirosettaModule() : pp::Module() {}
  virtual ~MinirosettaModule() {}

  /// Create and return a MinirosettaInstance object.
  /// @param[in] instance The browser-side instance.
  /// @return the plugin-side instance.
  virtual pp::Instance* CreateInstance(PP_Instance instance) {
		return new MinirosettaInstance(instance);
  }

};

}  // namespace

namespace pp {
/// Factory function called by the browser when the module is first loaded.
/// The browser keeps a singleton of this module.  It calls the
/// CreateInstance() method on the object you return to make instances.  There
/// is one instance per <embed> tag on the page.  This is the main binding
/// point for your NaCl module with the browser.
Module* CreateModule() {
  return new ncmini::MinirosettaModule();
}
}  // namespace pp
