#include <ppapi/cpp/instance.h>
#include <ppapi/cpp/module.h>

#include <ppapi/cpp/var.h>

#include <cstdio>
#include <string>

// libRosetta headers
#include <basic/options/option.hh>
#include <core/init/init.hh>
// C++ headers
#include <iostream>
#include <string>
#include <basic/Tracer.hh>
#include <protocols/rpc/rpc.hh>

//Auto Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/inline_file_provider.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/string.functions.hh>

#include "geturl_handler.h"
#include "static_database.hh"

using namespace ObjexxFCL;

#include <pthread.h>

/// These are the method names as JavaScript sees them.  Add any methods for
/// your class here.
namespace nacl_rosetta {
// A method consists of a const char* for the method ID and the method's
// declaration and implementation.
// TODO(sdk_user): 1. Add the declarations of your method IDs.
 
protocols::rpc::BasicCmdLineInit basic_init;

class RosettaInstance;

class RosettaURLFileHandler: public utility::Inline_File_Provider_Hook {
 public:
  RosettaURLFileHandler( RosettaInstance *rinst):
    rinst_(rinst) 
  {
  }
  // this will be called by Rosetta when it encounters a file it needs from the server. 
  // This function will request the file on the main thread (via the RosettaInstance pointer ) 
  // and then idle and periodically check the finished transfer variable until the return_file_callback
  // has happened. Once done, it will return the file contents to the Rosetta File provider and
  // return. In other words this function will block until the file transfer has completed (or until a 
  // timeout has hit ). 
  virtual bool request_file( const std::string filename, std::string &result_data );
  
  // this will be called by GetURLHandler when it is finished downloading the data and hands over the
  // file contents.
  virtual void return_file_callback( const std::string &result_data, bool error ); 
  
 private: // data
  RosettaInstance *rinst_;

  // this three variables are used to make the requesting thread (the rosetta thread) idle until 
  // the file results have been delivered or an error was returned.   
  bool finished_transfer_;
  pthread_mutex_t finished_transfer_mutex_; // mutex for file tansfer thread waiting
  pthread_cond_t finished_transfer_cond_;   // condition variable for file transfer thread waiting
  
  std::string result_data_; // buffer for file transfer 
  bool        error_; // error status (true is error, false is no error )
};

typedef utility::pointer::owning_ptr< RosettaURLFileHandler > RosettaURLFileHandlerOP;
typedef utility::pointer::owning_ptr< RosettaURLFileHandler const > RosettaURLFileHandlerCOP;

class RosettaInstance : public pp::Instance {
 public:
    explicit RosettaInstance(PP_Instance instance) : 
      pp::Instance(instance),
      tid_(0),callback_factory_(this), message_(""),rurlhandler_(NULL)
    {
    }
    virtual ~RosettaInstance() {}

    // Entry point for calls from JavaScript - interprets message and launches rosetta threads
    virtual void HandleMessage(const pp::Var& var_message);

    // This method will set up a call into the main Rosetta Library. This is typically called from
    // the standalone function StartRosettaThread 
    void RosettaThread();

    // Helper function to allow other threads to use PostMessage. Current Google Pepper API
    // does not support PostMessage to be called from anything but the main thread
    void PostMessage_from_main_thread( const std::string &data_to_send ){
      pp::Module::Get()->core()->CallOnMainThread( 0, callback_factory().NewCallback( &RosettaInstance::PostStringToBrowser,  data_to_send ));
    }
    
    // Helper function to allow other threads to request file transfers. Current Google Pepper API
    // does not seem to support this from arbitrary sub threads 
    void RequestFile_from_main_thread( const std::string &data_to_send ){
      pp::Module::Get()->core()->CallOnMainThread( 0, callback_factory().NewCallback( &RosettaInstance::request_file_transfer,  data_to_send ));
    }
  
  private: // functions

    // This method will set up file Handler and insert it into the RosettaFileProvider as a hook
    // to allow synconous callbacks from Rosetta Thread to request files.
    void setup_file_handler();

    // This method is called from the worker thread using CallOnMainThread.
    // It is not static, and allows PostMessage to be called.
    void* PostStringToBrowser(int32_t result, std::string data_to_send); 

    // This method is called from the worker thread using CallOnMainThread.
    // It is not static, and allows PostMessage to be called.
    void* request_file_transfer(int32_t result, std::string url); 

    // Return the callback factory.
    pp::CompletionCallbackFactory<RosettaInstance>& callback_factory(); 

  private: // data
    pthread_t tid_;
    pp::CompletionCallbackFactory<RosettaInstance> callback_factory_;
    std::string message_;

    RosettaURLFileHandlerOP rurlhandler_;
};
  


bool 
RosettaURLFileHandler::request_file( const std::string filename, std::string &result_data )
{
  // All Rosetta Datafiles are stored in this form on the server
  std::string url = "/data/rosetta_database/" + filename; // + ".txt" ; 
  rinst_->PostMessage_from_main_thread( "Trying to get: of '" + filename + "'" ); 
  
  const int MAX_ATTEMPTS = 3;

  for( int attempt = 0; attempt < MAX_ATTEMPTS; attempt ++ ){
    
    finished_transfer_  = false; 
    error_  = false; 
    result_data_.clear();
   
    //rinst_->PostMessage_from_main_thread( url ); 
    rinst_->RequestFile_from_main_thread( url );

    // now wait until the request file thread (it runs on the main thread, while this here is the rosetta thread)
    // has received the file 
    pthread_mutex_lock(&finished_transfer_mutex_);
    while( !finished_transfer_ ){
      pthread_cond_wait(&finished_transfer_cond_, &finished_transfer_mutex_ ); // we need a timeout here in case the other thread never unlocks us.
    }
    pthread_mutex_unlock(&finished_transfer_mutex_); 

    // if no error set the result variable to the data obtained. 
    if(!error_){
      result_data = result_data_; 
      
      // clean up memory
      result_data_.clear() ;
      // and report success
      std::stringstream debugmessage;
      debugmessage << result_data.size() << std::endl; 
      rinst_->PostMessage_from_main_thread( "Success!'" + filename + "' " + debugmessage.str() );
      return true;
    } else {
      std::stringstream debugmessage;
      debugmessage << result_data.size() << std::endl; 
      rinst_->PostMessage_from_main_thread( "Results of '" + filename + "' " + debugmessage.str() );
    }
  }
  
  // clean up memory even after failure
  result_data_.clear() ;
  
  // signal failed attempt at retrieving the resoruce
  return false;
}
  
void 
RosettaURLFileHandler::return_file_callback( const std::string &result_data, bool error ){
    
    // set the results of the file request and the error status
    result_data_ = result_data;
    error_ = error;
   
    if(error){
      rinst_->PostMessage_from_main_thread( result_data_ );
    }
    // now release the waiting thread
    finished_transfer_ = true; // << sempahore much better 
    pthread_mutex_lock(&finished_transfer_mutex_);
    pthread_cond_signal(&finished_transfer_cond_);
    pthread_mutex_unlock(&finished_transfer_mutex_); 

}



















void* StartRosettaThread(void *rosettaInst) {
  RosettaInstance *inst = static_cast<RosettaInstance*>(rosettaInst);
  inst->RosettaThread();
  return NULL;
}


void RosettaInstance::setup_file_handler(){
  if( rurlhandler_ == NULL){
    // set up the file handler
    rurlhandler_ = new RosettaURLFileHandler(this);  
    utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
    provider->add_file_provider_hook( rurlhandler_ ); 
  
    // also add statically compiled data files
    for(int i = 0; i < static_database_size; i++){
      provider->add_input_file( std::string( static_database[i][0]), std::string( static_database[i][1]) ); 
    }
      
    provider->add_black_listed_file( "" ); 
    provider->add_black_listed_file( ".gz" ); 
    provider->add_black_listed_file( "chemical/orbital_type_sets/fa_standard/.gz" ); 
    provider->add_black_listed_file( "chemical/residue_type_sets/fa_standard" );
    provider->add_black_listed_file( "chemical/residue_type_sets/fa_standard.gz" );
    provider->add_black_listed_file( "score12" );
    provider->add_black_listed_file( "standard" );

  }
}



// This method is called from the worker thread using CallOnMainThread.
// It is not static, and allows PostMessage to be called.
void* 
RosettaInstance::PostStringToBrowser(int32_t result, std::string data_to_send) {
  PostMessage(pp::Var(data_to_send));
  return 0;
}

// This method is called from the worker thread using CallOnMainThread.
// It is not static, and allows PostMessage to be called.
void* 
RosettaInstance::request_file_transfer(int32_t result, std::string url) {
  // create the URL handler which will perform the actual download. 
  GetURLHandler* handler = GetURLHandler::Create(this, url, &(*rurlhandler_) );
  if (handler != NULL) {
    handler->Start();
  }
  return 0;
}


// Return the callback factory.
pp::CompletionCallbackFactory<RosettaInstance>& 
RosettaInstance::callback_factory() {
  return callback_factory_;
}


void RosettaInstance::RosettaThread(){
  protocols::rpc::JSON_RPCOP rpc;     // rpc object - a Rosetta specific Object that manages RPC style calls into Rosetta
  utility::json_spirit::Object root;  // json object to be filled with all the information to be sent back. 
  
  try{
    // create the JSON_rpc object and initialize it with the json message from the caller
    rpc = new protocols::rpc::JSON_RPC( message_, true, &basic_init );

    // its own catch block so we can be sure to get the tracer output is also captured assuming the above constructor returns ok. 
    try{
      rpc->run();
    
      root.push_back( utility::json_spirit::Pair( "output", rpc->tracer() ) );
      utility::json_spirit::Object energies;  
      energies.push_back( utility::json_spirit::Pair( "score" ,  rpc->get_fa_score() ) ); 
      energies.push_back( utility::json_spirit::Pair( "irms", rpc->get_irms() ) ); 

      protocols::rpc::pose_energies_to_json( rpc->outputpose(), energies );
      root.push_back( utility::json_spirit::Pair( "energies",  energies ) );                              // rosetta energy values
      
      // Grab the output pose as a PDB file format and write it into the "pdbdata" field of the jsno return message 
      std::stringstream pdbdatastream;
      rpc->outputpose().dump_pdb( pdbdatastream );
      root.push_back( utility::json_spirit::Pair( "pdbdata", pdbdatastream.str() ) ); // the PDB data itself of course
      
      // Capture the CPU time spent
      root.push_back( utility::json_spirit::Pair( "cputime", (int) rpc->runtime() ) );
    }
    catch( utility::excn::EXCN_Msg_Exception &excn ){
      root.push_back( utility::json_spirit::Pair( "error", excn.msg() ) );
    }

  }
  catch( utility::excn::EXCN_Msg_Exception &excn ){
    // we're here because the constructuor of JSON_RPC failed - capture the message of the Exception
    root.push_back( utility::json_spirit::Pair( "error",  excn.msg() ) );
  }
  catch ( std::string &s ){ 
    root.push_back( utility::json_spirit::Pair( "error",  "Exception occured during Rosetta Execution:" + s  ) );
  }
  catch ( ... ){ 
    root.push_back( utility::json_spirit::Pair( "error",  "Unknown Exception occured during Rosetta Execution"  ) );
  }

  // write out json response and send back stringified
  std::stringstream sstr;
  write( root, sstr ); // stringify

  // return data. This will be unpacked and interpreted by the JavaScript caller.
  PostMessage_from_main_thread( sstr.str() );
  PostMessage_from_main_thread( "DONE!" ); 
}

void RosettaInstance::HandleMessage(const pp::Var& var_message) {
  if (!var_message.is_string()) {
    PostMessage( pp::Var( "ERROR; RosettaInstance only accepts strings as messages.\n" ) );
    return;
  }
  message_ = var_message.AsString();
  
  // ensure the file handler is ready to go.
  setup_file_handler();

  // Start main Rosetta thread 
  if( pthread_create( &tid_, NULL,  StartRosettaThread, this)) {
   // return Error message to JavaScript incase we failed to launch the thread
   PostMessage( pp::Var( "ERROR; pthread_create() failed.\n" ) );
   return;
  } 
  
}

class RosettaModule : public pp::Module {
 public:
  RosettaModule() : pp::Module() {
    
  }
  virtual ~RosettaModule() {}

  /// Create and return a RosettaInstance object.
  /// @param[in] instance a handle to a plug-in instance.
  /// @return a newly created RosettaInstance.
  /// @note The browser is responsible for calling @a delete when done.
  virtual pp::Instance* CreateInstance(PP_Instance instance) {
    return new RosettaInstance(instance);
  }
};
}  // namespace nacl_rosetta


namespace pp {
/// Factory function called by the browser when the module is first loaded.
/// The browser keeps a singleton of this module.  It calls the
/// CreateInstance() method on the object you return to make instances.  There
/// is one instance per <embed> tag on the page.  This is the main binding
/// point for your NaCl module with the browser.
/// @return new RosettaModule.
/// @note The browser is responsible for deleting returned @a Module.
  Module* CreateModule() {
    return new nacl_rosetta::RosettaModule();
  }
}  // namespace pp

