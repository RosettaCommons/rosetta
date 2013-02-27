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
#include <ObjexxFCL/string.functions.hh>

using namespace ObjexxFCL;

#include <pthread.h>

/// These are the method names as JavaScript sees them.  Add any methods for
/// your class here.
namespace nacl_rosetta {
// A method consists of a const char* for the method ID and the method's
// declaration and implementation.
// TODO(sdk_user): 1. Add the declarations of your method IDs.

const char* const kRunRosettaMethodId = "runRosetta";
const char* const kGetStdoutMethodId = "getStdout";

std::string global_command_line;
std::stringstream ss_cout;
pthread_t tid;

pp::Var MarshallGetStdout() {
  return pp::Var( std::string( ss_cout.str() ) );
}


void * run( void * data ){

	using namespace core;
	using namespace basic::options;

	std::string command_line = global_command_line;

	//pp::Instance* the_instance = (pp::Instance*) void_the_instance;
	ss_cout << "RunRosetta Native Client version of rosetta\n";
	ss_cout << "Arguments: '" + command_line + "'\n";

	// build c-style arguments.
	utility::vector1< std::string > tokens ( utility::split( command_line ) );
	utility::vector1<std::string> tokens_v1;
	basic::final_channel = &ss_cout;

	tokens_v1.push_back( "mini_rosetta_native_client.nexe" );
	for( std::vector< std::string >::const_iterator it=tokens.begin(); it != tokens.end(); ++ it ){
		tokens_v1.push_back( * it );
		(*basic::final_channel) << (*it) << std::endl;
	}

	ss_cout << "start" << std::endl;
	try{
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
pp::Var MarshallRunRosetta( const std::string& text ) { 
	//global_command_line = args[0].AsString();
	//pthread_create(&tid, NULL, run, NULL );
  return pp::Var(  std::string("The Rosetta Output") + text );
}

class RosettaInstance : public pp::Instance {
 public:
  explicit RosettaInstance(PP_Instance instance) : pp::Instance(instance) {
    printf("RosettaInstance.\n");
  }
  virtual ~RosettaInstance() {}

  virtual void HandleMessage(const pp::Var& var_message);
};

void RosettaInstance::HandleMessage(const pp::Var& var_message) {
  if (!var_message.is_string()) {
    return;
  }
  std::string message = var_message.AsString();
  pp::Var return_var;
	if( message == kRunRosettaMethodId){
		return_var = MarshallRunRosetta( "TestText" );
	}else
	if( message == kGetStdoutMethodId){
		return_var = MarshallGetStdout();
	}
  PostMessage(return_var);
}

class RosettaModule : public pp::Module {
 public:
  RosettaModule() : pp::Module() {
    printf("Got here.\n");
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

