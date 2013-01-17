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
#include <protocols/loophash/LoopHashRelaxProtocol.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

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

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
#include <utility/inline_file_provider.hh>
#include <numeric/random/random.hh>

#include <curl/curl.h>
#include "json/json.h"






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




void pose_energies_to_json( core::pose::Pose const & pose, Json::Value &root ) {
  using namespace core::pose::datacache;

  core::scoring::EnergyMap const emap = pose.energies().total_energies();
	core::scoring::EnergyMap const wts  = pose.energies().weights();

  core::scoring::EnergyMap::const_iterator emap_iter, wts_iter;
	for ( emap_iter = emap.begin(), wts_iter = wts.begin();
			emap_iter != emap.end() && wts_iter!= wts.end();
			++emap_iter && ++wts_iter
	) {

		// only grab scores that have non-zero weights.
		if ( *wts_iter != 0.0 ) {
			core::scoring::ScoreType sc_type
				= core::scoring::ScoreType( emap_iter - emap.begin() + 1 );
			std::string name = core::scoring::name_from_score_type( sc_type );

		  root[ name ] = Json::Value( (*emap_iter) * (*wts_iter) );

		} // if ( *wts_iter != 0.0 )
	} // for ( emap_iter ...)

}



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
   RosettaJob (): initialized_(false)
   {}

   bool request_job_from_server( const ServerInfo & server_info )
   {
      runtime_assert( !initialized_ )

      CurlPost cg;
      std::string data;
      std::cout << "URL: " << server_info.url_gettask() << std::endl;
      try {
        data = cg.post( server_info.url_gettask() , "", "lease_time=100"  );
      } catch ( std::string error ){
        std::cerr << "ERROR(" << server_info.url_gettask() << std::endl;
        std::cerr << "):" << error << std::endl;
        return false;
      }

      if( data == std::string("[]") ){
        std::cout << "No work on server" << std::endl;
        return false;
      }

      output_capture_clear();

      // break down the input data
      //std::cout <<  data << std::endl;

      Json::Value values;
      Json::Reader reader;
      reader.parse( data, values);

      std::string payload = values.get("payload","").asString();
      taskname_ = values.get("name","").asString();

      if( taskname_ == "" ){
        std::cout << "Error connecting to server (server didn't reply with JSON object or taskname was not set) " << std::endl;
        return false;
      }

      Json::Value payload_values;
      Json::Reader payload_reader;
      reader.parse( payload, payload_values);

      hash_      = payload_values.get("hash","").asString();
      key_       = payload_values.get("key","").asString();
      user_id_   = payload_values.get("user_id","").asString();
      operation_ = payload_values.get("operation","").asString();

      std::cout << "NEW JOB: Hash: " << hash_ << " KEY: " << key_ << " USER: " << user_id_ << " TASKNAME: " << taskname_ << std::endl;
      std::string job_data_string = payload_values.get("job_data","").asString();
      // std::cout << "Inputdata: " << job_data_string << std::endl;

      Json::Reader job_data_reader;
      job_data_reader.parse( job_data_string , job_data_);

      // extract the command name - this will be used to fork the code.
      command_ = job_data_.get("command","").asString();
      std::cout << "Rosetta command: " << command_ << std::endl;

      xmlscript_ = job_data_.get("xmlscript","").asString();
      std::cout << "XML script: " << xmlscript_ << std::endl;

      // extract the pdbdata (the starting structure)
      Json::Value pdbdata = job_data_.get("pdbdata","");
      std::string pdbdata_string = pdbdata.asString();

      output_capture_clear();
      output_capture_start();
      
      // create the pose
      try {
        core::import_pose::pose_from_pdbstring( inputpose_, pdbdata_string );
        
        // Load any files that were given as part of this job.
        load_new_set_of_virtual_files(  job_data_.get("user_files","") );

        // Load any flags that were given as part of this job.
        std::cout << "Initializing options: " << command_ << std::endl;

        if ( job_data_.isMember("user_flags") ){
          load_new_set_of_user_flags( job_data_.get("user_flags","[]") );
        }
       
        if ( job_data_.isMember("flags_file") ){
          std::string flags_file = job_data_.get("flags_file","").asString();
          std::cerr << "HERE:" << flags_file << std::endl; 
          load_user_flag_file( flags_file );
        }
      }
      catch ( utility::excn::EXCN_Base& excn ) {
				excn.show( TR.Error );
        // return error or empty pose
        return_error_to_server( server_info );
        return false;
      }
      catch ( std::string err ){ 
				TR.Error << err << std::endl;
        // return error or empty pose
        return_error_to_server( server_info );
        return false;
      }
      outputpose_ = inputpose_; 
      // since a null operation would leave the pose unchanged:
      error_ = 0;
      initialized_ = true;
      return true;
    }

    bool return_error_to_server( const ServerInfo & server_info ){
      std::cout << __FILE__ << " " << __LINE__ << std::endl; 
      // stop stdout capturing
      output_capture_stop();
      
      // do some basic measurements
      Json::Value energies;
      energies[ "score" ] = 0;
      energies[ "irms" ] = 0;

      // assemble the json structure needed.
      Json::Value root;

      // set the output values
      root["parental_key"] = Json::Value( key_ );   // the key and hash become parental_  key and hash
      root["parental_hash"] = Json::Value( hash_ ); // so we can keep track of the geneaology of structures
      root["user_id"] = Json::Value( user_id_ ); // so we can keep track of the geneaology of structures
      root["operation"] = Json::Value( operation_ ); // so we can keep track of the geneaology of structures
      root["taskname"] = Json::Value( taskname_ ); // so we can keep track of the geneaology of structures
      root["pdbdata"] = Json::Value( "" ); // the PDB data itself of course
      root["workerinfo"] = Json::Value(  "Rosetta Backend 0.1" );
      root["stderr"] = Json::Value( tracer_output_stream_.str() ); 
      std::cout << "STDERROR:" <<  root["stderr"] << std::endl;
      root["energies"] = energies;
      root["error"] =  Json::Value(1);
      root["cputime"] =  Json::Value(0);

      output_capture_clear();
      // add energy info etc other goodies here

      Json::FastWriter writer;
      std::string output_json = "output=" + writer.write( root );

      CurlPost cg;
      try {
        std::string return_data = cg.post( server_info.url_putresult() , "", output_json );
        std::cout << "POST " << server_info.url_putresult() << " : " << return_data << std::endl;
      } catch ( std::string error ){
        std::cerr << "ERROR returning results:" << error << std::endl;
      }
    }

    bool return_results_to_server( const ServerInfo & server_info ){

      std::cout << __FILE__ << " " << __LINE__ << std::endl; 
      // do some basic measurements
      Json::Value energies;
	    core::scoring::ScoreFunctionOP fascorefxn = core::scoring::getScoreFunction();
      energies[ "score" ] = (*fascorefxn)(outputpose_);
      core::scoring::calpha_superimpose_pose( outputpose_, inputpose_ );
      energies[ "irms" ] = core::scoring::rmsd_with_super(inputpose_, outputpose_  , core::scoring::is_protein_backbone_including_O );

      pose_energies_to_json( outputpose_, energies );

      // Now send back results to server.
      std::stringstream pdbdatastream;
      outputpose_.dump_pdb( pdbdatastream );
        
      // DEBUG CODE
      std::ofstream ofile("b4.pdb"); 
      inputpose_.dump_pdb( ofile );
      ofile.close();
      std::ofstream ofile2("af.pdb"); 
      outputpose_.dump_pdb( ofile2 );
      ofile2.close(); 
      // ENDOF DEBUG CODE

      // assemble the json structure needed.
      Json::Value root;

      // set the output values
      root["parental_key"] = Json::Value( key_ );   // the key and hash become parental_  key and hash
      root["parental_hash"] = Json::Value( hash_ ); // so we can keep track of the geneaology of structures
      root["user_id"] = Json::Value( user_id_ ); // so we can keep track of the geneaology of structures
      root["operation"] = Json::Value( operation_ ); // so we can keep track of the geneaology of structures
      std::cout << "Finished TaskName: " << taskname_ << std::endl;
      root["taskname"] = Json::Value( taskname_ ); // so we can keep track of the geneaology of structures
      root["pdbdata"] = Json::Value( pdbdatastream.str() ); // the PDB data itself of course
      root["error"] = Json::Value(error_);
      if( error_ ){ 
        root["pdbdata"] =  Json::Value("");  // if there was a fatal error dont return a pdb (to save space)
      }

      std::stringstream rosetta_version;
      rosetta_version << "Mini-Rosetta version " << core::minirosetta_svn_version() << " from " << core::minirosetta_svn_url();

      root["workerinfo"] = Json::Value( rosetta_version.str() ); 
      root["stderr"] = Json::Value( tracer_output_stream_.str() ); // stderr output for debugging
      root["energies"] = energies;                                 // rosetta energy values
      root["cputime"] =  Json::Value(Json::Int64(endtime_ - starttime_));

      // add energy info etc other goodies here

      Json::FastWriter writer;
      std::string output_json = "output=" + writer.write( root );
      //std::cout << writer.write( root ) << std::endl;
      // now do a POST on the server to hand back the data.

      CurlPost cg;
      try {
        std::string return_data = cg.post( server_info.url_putresult() , "", output_json );
        std::cout << "POST " << server_info.url_putresult() << " : " << return_data << std::endl;
      } catch ( std::string error ){
        std::cerr << "ERROR returning results:" << error << std::endl;
      }

    }


    


    void run(){
      using namespace basic::options;
      using namespace basic::options::OptionKeys;

      // start capturing Tracer outputs
      output_capture_start();
      starttime_ = time(NULL);

      std::cerr << "RANDOM CHECK: " << RG.uniform() << std::endl; 
      try{
        std::cout << "Executing: " << command_ << std::endl;
        TR << "Executing: " << command_ << std::endl;
        if( command_ == "score" ){
        }
        if( command_ == "relax" ){
          protocols::moves::MoverOP protocol = protocols::relax::generate_relax_from_cmd();
          protocol->apply( outputpose_ );
        }
        if( command_ == "loophash" ){
          TR << "LOOPHASH!!  " << rosetta_script_ << std::endl;
          protocols::loophash::LoopHashRelaxProtocolOP lh_protocol = new protocols::loophash::LoopHashRelaxProtocol( loop_hash_library );
          lh_protocol->manual_call( outputpose_ );
        }
        if( command_ == "enzdes" ){
        }
        if( command_ == "xmlscript" ){
          utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
          provider->add_input_file( "script.xml", xmlscript_ );
          protocols::jd2::DockDesignParser ddp;
          protocols::jd2::JobCOP job;
          protocols::moves::MoverOP protocol;
          ddp.generate_mover_from_pose( job, outputpose_ , protocol, true, "script.xml" );
          protocol->apply( outputpose_ );
        }
      } 
      catch( utility::excn::EXCN_Msg_Exception e ){
        TR << e.msg() << std::endl;
        error_ = 1;
      }
      catch(...) {
        TR << "UNKNOWN ERROR" << std::endl;
        error_ = 1;
      }
      
      // stop stdout capturing
      output_capture_stop();
      endtime_ = time(NULL);
    }

 private:

  void output_capture_start(){
      basic::set_new_final_channel( &tracer_output_stream_ );
  }

  void output_capture_stop(){
      basic::set_default_final_channel();
  }

  void output_capture_clear(){
    tracer_output_stream_.str(std::string());
  }

  void load_user_flag_file( const std::string &flags_file ){
    std::stringstream options_file;
    options_file << flags_file << std::endl; // ensure the last line has a CR/lF at the end. 
    basic::options::option.load_options_from_stream( options_file );
  }

  void load_new_set_of_user_flags(  const Json::Value &json_user_flags ){
    std::vector < std::string > user_flags;

    if( !json_user_flags.isObject() ){
      std::cerr << "User_flags is not an object - flags not loaded";
      throw 1;
    }
 
    for( Json::Value::const_iterator it = json_user_flags.begin(); 
         it != json_user_flags.end(); ++it ){
      std::cout << it.memberName() << "   " << (*it)["value"].asString() << std::endl;
      std::string new_flag(  std::string(it.memberName()) + " " + (*it)["value"].asString() );
      user_flags.push_back( new_flag );
      std::cout << "User flag conversion: " << new_flag << std::endl;
    }

    basic::options::option.load( user_flags, false);
  }

  void load_new_set_of_virtual_files( Json::Value newfiles, bool clear_previous = true ){
    utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
    if( clear_previous ) provider->clear_input_files();
    for( int i=0; i < newfiles.size(); ++i ){
      std::string filename = newfiles[i].get("filename","").asString();
      std::string contents = newfiles[i].get("content","").asString();
      TR << "Loading virtual file: " << filename << " " << contents.size() << " bytes" << std::endl;
      provider->add_input_file( filename, contents );
    }
  }

  // input stuff (coming from the server)
  core::pose::Pose inputpose_;
  std::string taskname_;
  std::string hash_;
  std::string key_;
  std::string user_id_;
  std::string operation_;
  Json::Value job_data_;

  std::string rosetta_script_;
  std::string command_;
  std::string xmlscript_;
  bool initialized_;

  // output stuff
  core::pose::Pose outputpose_;
  int              error_; 
  std::stringstream tracer_output_stream_;
  long             starttime_;
  long             endtime_;

  
};

class RosettaBackend {
 public:
    RosettaBackend( ServerInfo &serverinfo, int argc, char * argv [] ) :
      serverinfo_(serverinfo),
      argc_(argc),
      argv_(argv)
    {}

 private:

 public:

    void run(){
      do{
        RosettaJob newjob;
        core::Size wait_count = 0;

        using namespace basic::options;
        using namespace basic::options::OptionKeys;
        utility::options::OptionCollection &option_collection  = initialize();
        option_collection.load( argc_, argv_, false);

        while( !newjob.request_job_from_server( serverinfo_ ) ){
          core::Real waittime =  std::min( (double)10.0f, 0.5f* pow( (float)1.3, (float) wait_count )); // in seconds

          // HACK!
          waittime =  serverinfo_.poll_frequency();

          std::cout << "No work. Waiting " << waittime << " seconds before retrying." << std::endl;
          sleep( waittime );
          wait_count ++;
        };
        
        newjob.run();
        
        std::cerr << "Returning results to server (mainloop)" << std::endl;

        // now return the results to the server
        try {
          newjob.return_results_to_server( serverinfo_ );
        }
        catch ( utility::excn::EXCN_Base& excn ) {
          excn.show( TR.Error );
          // return error or empty pose
          newjob.return_error_to_server( serverinfo_ );
        }
        catch ( std::string err ){ 
          TR.Error << err << std::endl;
          // return error or empty pose
          newjob.return_error_to_server( serverinfo_ );
        }
      } while (true);

    };


 private:
  ServerInfo serverinfo_;
  int argc_;
  char ** argv_ ;

};

int
main( int argc, char * argv [] )
{
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
  RosettaBackend backend( server, argc, argv );


  backend.run();

 	return 0;
}

