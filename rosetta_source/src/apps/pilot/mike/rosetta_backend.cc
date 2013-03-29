// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/DockDesignParser.hh>
#include <protocols/frag_picker/VallChunk.hh>

#include <utility/pointer/owning_ptr.hh>

#include <core/svn_version.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <basic/options/option.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>

#include <basic/Tracer.hh>

#include <devel/init.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/match/Hit.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rpc/rpc.hh>
#include <protocols/loophash/LoopHashRelaxProtocol.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashRelaxProtocol.hh>

#include <protocols/jd2/DockDesignParser.hh>

// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/batch_relax.OptionKeys.gen.hh>
#include <basic/options/keys/rbe.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>


#include <utility/json_spirit/json_spirit_value.h>
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_tools.hh>

#include <utility/vector1.hh>
#include <utility/inline_file_provider.hh>
#include <numeric/random/random.hh>

#include <curl/curl.h>


static basic::Tracer TR("main");
static numeric::random::RandomGenerator RG(7293464);

class CurlGet {
 public:
   CurlGet() { }

  // This is the writer call back function used by curl
  static int writer(char *data, size_t size, size_t nmemb,
                    std::string *buffer)
  {
    // What we will return
    int result = 0;

    // Is there anything in the buffer?
    if (buffer != NULL)
    {
      // Append the data to the buffer
      buffer->append(data, size * nmemb);

      // How much did we write?
      result = size * nmemb;
    }

    return result;
  }

 private:
  char *getErrorBuffer() { return &errorBuffer[0]; }

  std::string * getbuffer() { return &buffer; }

 public:
  std::string get( const std::string &url ){
    std::cout << "Retrieving " << url << std::endl;

    // Our curl objects
    CURL *curl;
    CURLcode result;

    // Create our curl handle
    curl = curl_easy_init();

    if (curl)
    {
      // Now set up all of the curl options
      curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, getErrorBuffer() );
      curl_easy_setopt(curl, CURLOPT_URL, url.c_str() );
      curl_easy_setopt(curl, CURLOPT_HEADER, 0);
      curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);
      curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, this->writer);
      curl_easy_setopt(curl, CURLOPT_WRITEDATA, getbuffer() );

      // Attempt to retrieve the remote page
      result = curl_easy_perform(curl);

      // Always cleanup
      curl_easy_cleanup(curl);

      // Did we succeed?
      if (result == CURLE_OK)
      {
        return buffer;
      }
      else
      {
        std::stringstream ss;
        ss << "Error: [" << result << "] - " << errorBuffer;
        throw std::string( ss.str() );
      }
    }
  }

 private:

  // Write any errors in here
  char errorBuffer[CURL_ERROR_SIZE];

  // Write all expected data in here
  std::string buffer;
};


class CurlPost {
 public:
   CurlPost() { }

  // This is the writer call back function used by curl
  static int writer(char *data, size_t size, size_t nmemb,
                    std::string *buffer)
  {
    std::cout << "Post writer function callback" << std::endl;
    // What we will return
    int result = 0;
    //std::cout << "Buffer: " << *buffer << std::endl;
    // Is there anything in the buffer?
    if (buffer != NULL)
    {
      // Append the data to the buffer
      buffer->append(data, size * nmemb);

      // How much did we write?
      result = size * nmemb;
    }

    return result;
  }

 private:
  char *getErrorBuffer() { return &errorBuffer[0]; }

  std::string * getreadbuffer() { return &readbuffer; }

  std::string * getwritebuffer() { return &writebuffer; }

 public:
  std::string post( const std::string &url, const std::string &data, const std::string &fields ){
    std::cout << "Retrieving " << url << std::endl;

    // Our curl objects
    CURL *curl;
    CURLcode result;

    // Create our curl handle
    curl = curl_easy_init();

    if (curl)
    {
      // fill read buffer
      readbuffer = data;

      // Now set up all of the curl options
      curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, getErrorBuffer() );
      curl_easy_setopt(curl, CURLOPT_URL, url.c_str() );
      curl_easy_setopt(curl, CURLOPT_POSTFIELDS, fields.c_str() );
      curl_easy_setopt(curl, CURLOPT_HEADER, 0);
      curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);
      curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, this->writer);
      curl_easy_setopt(curl, CURLOPT_WRITEDATA, getwritebuffer() );

      // Attempt to retrieve the remote page
      result = curl_easy_perform(curl);

      // Always cleanup
      curl_easy_cleanup(curl);

      // Did we succeed?
      if (result == CURLE_OK)
      {
        return writebuffer;
      }
      else
      {
        std::stringstream ss;
        ss << "Error: [" << result << "] - " << errorBuffer;
        throw std::string( ss.str() );
      }
    }
  }

 private:

  // Write any errors in here
  char errorBuffer[CURL_ERROR_SIZE];

  // Write all expected input and out data in here
  std::string readbuffer;
  std::string writebuffer;
};




// Read to go protocols - these are slow to init:
protocols::loophash::LoopHashLibraryOP loop_hash_library;





class ServerInfo {
 public:
    ServerInfo( std::string server_url, std::string server_port, core::Real poll_frequency ):
      server_url_(server_url),
      server_port_(server_port),
      poll_frequency_(poll_frequency)
    {

    }

   const std::string url_gettask()   const { return full_url() + "/task/get";      };
   const std::string url_putresult() const { return full_url() + "/structure/put"; };
   const std::string full_url() const { return server_url_ + ":" + server_port_; }
   const core::Real poll_frequency() const { return poll_frequency_; }
 private:
   std::string server_url_;
   std::string server_port_;
   core::Real poll_frequency_;
};


class RosettaJob {
 public:
   RosettaJob( const ServerInfo &serverinfo, protocols::rpc::BasicInit *basic_init):
    serverinfo_(serverinfo),
    initialized_(false),
    basic_init_(basic_init)
   {}

   bool request_job_from_server( )
   {
      runtime_assert( !initialized_ )

      CurlPost cg;
      std::string data;
      std::cout << "URL: " << serverinfo_.url_gettask() << std::endl;
      try {
        data = cg.post( serverinfo_.url_gettask() , "", "lease_time=100"  );
      } catch ( std::string error ){
        std::cerr << "ERROR(" << serverinfo_.url_gettask() << std::endl;
        std::cerr << "):" << error << std::endl;
        return false;
      }

      if( data == std::string("[]") ){
        std::cout << "No work on server" << std::endl;
        return false;
      }

      // break down the input data
      //std::cout <<  data << std::endl;

      utility::json_spirit::mObject parsed_json;
      std::string payload;
      try{
        parsed_json = utility::json_spirit::read_mObject( data );
        payload = get_string( parsed_json, "payload" );
        taskname_ = get_string( parsed_json, "name" );
      }
      catch( utility::excn::EXCN_Msg_Exception &excn ){
        TR.Error << "EXCEPTION: " << excn.msg() << std::endl; // print the exception message to the Error stream.
      }

      //std::cout << payload << std::endl;
      utility::json_spirit::mObject parsed_payload;

      try{
        parsed_payload  = utility::json_spirit::read_mObject( payload );
      }
      catch( utility::excn::EXCN_Msg_Exception &excn ){
        std::cout << "Error extracting payload" << std::endl;
        TR.Error << "EXCEPTION: " << excn.msg() << std::endl; // print the exception message to the Error stream.
        return false;
      }

      hash_      = get_string_or_empty( parsed_payload, "hash");
      key_       = get_string_or_empty( parsed_payload, "key");
      user_id_   = get_string_or_empty( parsed_payload, "user_id");
      operation_ = get_string_or_empty( parsed_payload, "operation");

      std::cout << "NEW JOB: Hash: " << hash_ << " KEY: " << key_ << " USER: " << user_id_ << " TASKNAME: " << taskname_ << std::endl;

      std::string job_data_string;
      try{
        job_data_string = get_string( parsed_payload, "job_data" );
      }
      catch( utility::excn::EXCN_Msg_Exception &excn ){
        std::cout << "Error extracting job_data from RPC request" << std::endl;
        TR.Error << "EXCEPTION: " << excn.msg() << std::endl; // print the exception message to the Error stream.
        return false;
      }

      //std::cout << "Inputdata: " << job_data_string << std::endl;
      std::cout << "Inputdata: " << job_data_string.size() << " bytes received" << std::endl;

      try{
        rpc = new protocols::rpc::JSON_RPC( job_data_string, false, basic_init_ );
        // intrpret the actual RPC contents.
      }
      catch ( utility::excn::EXCN_Base& excn ) {
        return_results_to_server( true );
        return false;
      }

      // since a null operation would leave the pose unchanged:
      initialized_ = true;
      return true;
    }


    bool return_results_to_server( bool error = false ){

      utility::json_spirit::Object energies;
      utility::json_spirit::Object root;


      // set the output values
      root.push_back( utility::json_spirit::Pair( "parental_key",  key_ ) );   // the key and hash become parental_  key and hash
      root.push_back( utility::json_spirit::Pair( "parental_hash",  hash_ ) ); // so we can keep track of the geneaology of structures
      root.push_back( utility::json_spirit::Pair( "user_id",  user_id_ ) ); // so we can keep track of the geneaology of structures
      root.push_back( utility::json_spirit::Pair( "operation",  operation_ ) ); // so we can keep track of the geneaology of structures
      std::cout << "Finished TaskName: " << taskname_ << std::endl;
      root.push_back( utility::json_spirit::Pair( "taskname",  taskname_ ) ); // so we can keep track of the geneaology of structures
      std::stringstream rosetta_version;
      rosetta_version << "Mini-Rosetta version " << core::minirosetta_svn_version() << " from " << core::minirosetta_svn_url();
      root.push_back( utility::json_spirit::Pair( "workerinfo",  rosetta_version.str() ) );


      if( !error ){
        // do some basic measurements

        energies.push_back( utility::json_spirit::Pair( "score" ,  rpc->get_fa_score() ) );
        energies.push_back( utility::json_spirit::Pair( "irms", rpc->get_irms() ) );

        protocols::rpc::pose_energies_to_json( rpc->outputpose(), energies );

        // Now send back results to server.
        std::stringstream pdbdatastream;
        rpc->outputpose().dump_pdb( pdbdatastream );

        // DEBUG CODE
//        std::ofstream ofile("b4.pdb");
//        inputpose_.dump_pdb( ofile );
//        ofile.close();
//        std::ofstream ofile2("af.pdb");
//        outputpose_.dump_pdb( ofile2 );
//        ofile2.close();
        // ENDOF DEBUG CODE

        root.push_back( utility::json_spirit::Pair( "pdbdata", pdbdatastream.str() ) ); // the PDB data itself of course
        root.push_back( utility::json_spirit::Pair( "energies",  energies ) );                              // rosetta energy values
        root.push_back( utility::json_spirit::Pair( "cputime", (int) rpc->runtime() ) );

     }else{

        // set the output values
        root.push_back( utility::json_spirit::Pair( "pdbdata",  "" ) ); // the PDB data itself of course
        root.push_back( utility::json_spirit::Pair( "cputime",  0 ) );
      }

      root.push_back( utility::json_spirit::Pair( "error",  0 ) );

      std::cout << "STDERROR:" <<  rpc->tracer() << std::endl;
      root.push_back( utility::json_spirit::Pair( "stderr", rpc->tracer() ) ); // stderr output for debugging

      // add energy info etc other goodies here
      std::stringstream sstr;
      write( root, sstr );
      std::string output_json = "output=" + sstr.str();
      //std::cout << output_json << std::endl;
      CurlPost cg;
      try {
        std::string return_data = cg.post( serverinfo_.url_putresult() , "", output_json );
        std::cout << "POST " << serverinfo_.url_putresult() << " : " << return_data << std::endl;
      } catch ( std::string error ){
        std::cerr << "ERROR returning results:" << error << std::endl;
      }

    }


    bool run_and_return_to_server(){
      try{
        rpc->run();
      }
      catch ( utility::excn::EXCN_Base& excn ) {
        return_results_to_server( true );
        return false;
      }
      return_results_to_server();
    }

 private:

  // input stuff (coming from the server)
  core::pose::Pose inputpose_;
  std::string taskname_;
  std::string hash_;
  std::string key_;
  std::string user_id_;
  std::string operation_;

  bool initialized_;

  protocols::rpc::JSON_RPCOP rpc;

  ServerInfo serverinfo_;
  protocols::rpc::BasicInit *basic_init_;
};

class RosettaBackend {
 public:
    RosettaBackend( const ServerInfo &serverinfo, protocols::rpc::BasicInit *basic_init ) :
      serverinfo_(serverinfo),
      basic_init_(basic_init)
    {}

 private:

 public:

    void run(){
      do{
        RosettaJob newjob(serverinfo_, basic_init_ );

        core::Size wait_count = 0;

        while( !newjob.request_job_from_server() ){
          core::Real waittime;
          //waittime =  std::min( (double)10.0f, 0.5f* pow( (float)1.3, (float) wait_count )); // in seconds
          waittime =  serverinfo_.poll_frequency();
          std::cout << "No work. Waiting " << waittime << " seconds before retrying." << std::endl;
          sleep( waittime );
          wait_count ++;
        };

        newjob.run_and_return_to_server();

        std::cerr << "Returning results to server (mainloop)" << std::endl;


      } while (true);

    };


 private:
  ServerInfo serverinfo_;
  protocols::rpc::BasicInit *basic_init_;
};



int
main( int argc, char * argv [] )
{
    try {
 	using namespace core;
 	using namespace protocols;
 	using namespace protocols::loophash;
 	using namespace protocols::jd2;
 	using namespace basic::options;
 	using namespace basic::options::OptionKeys;
 	using io::silent::SilentStructFactory;
 	using io::silent::SilentStructOP;


 	// initialize core
 	devel::init(argc, argv);

  BasicCmdLineInit basic_init( argc, argv );

  evaluation::PoseEvaluatorsOP evaluators_( new protocols::evaluation::PoseEvaluators() );
  evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluators_);

  // initialize all the protocols
  // these are slow so only do them once at start up not when the job is requested for lower latency.
	if( option[lh::db_path].user() ){
    utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
    loop_hash_library = new LoopHashLibrary( loop_sizes );
    loop_hash_library->load_db();
  }

  // set up the application server information package
  ServerInfo server( option[rbe::server_url](), option[rbe::server_port](), option[rbe::poll_frequency ]() );
  std::cout << "server: " << server.url_gettask() << std::endl;
  RosettaBackend backend( server, &basic_init );


  backend.run();

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}

