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
#include <protocols/frag_picker/VallChunk.hh>

#include <utility/pointer/owning_ptr.hh>

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

#include <curl/curl.h>
#include "json/json.h"



static basic::Tracer TR("main");

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
        std::cout << "No work on serer" << std::endl;
        return false;
      }

      // break down the input data
      //std::cout <<  data << std::endl;
      
      Json::Value values;
      Json::Reader reader;
      reader.parse( data, values);
      
      std::string payload = values.get("payload","").asString();
      taskname_ = values.get("name","defaulttaskname").asString();
      

      Json::Value payload_values;
      Json::Reader payload_reader;
      reader.parse( payload, payload_values);
      
      hash_      = payload_values.get("hash","").asString();
      key_       = payload_values.get("key","").asString();
      operation_ = payload_values.get("operation","").asString();
     
      std::cout << "NEW JOB: Hash: " << hash_ << " KEY: " << key_ << " TASKNAME: " << taskname_ << std::endl;
      std::string inputdata_string = payload_values.get("pdbdata","").asString();
      //std::cout << "Inputdata: " << inputdata_string << std::endl;

      Json::Value inputdata_values;
      Json::Reader inputdata_reader;
      inputdata_reader.parse( inputdata_string , inputdata_values);
      
      Json::Reader rosetta_script_reader;
      Json::Value rosetta_script_value;
      //rosetta_script_ = inputdata_values.get("rosetta_script","").asString();
      rosetta_script_ = "";
      std::cout << "Rosetta script string: " << rosetta_script_ << std::endl;
      
      command_ = inputdata_values.get("command","").asString();
      std::cout << "Rosetta command: " << command_ << std::endl;
      
      Json::Value pdbdata = inputdata_values.get("pdbdata","");
      std::string pdbdata_string = pdbdata.asString();
      //std::cout << pdbdata_string << std::endl;
  
      // create the pose 
      try {
        core::import_pose::pose_from_pdbstring( inputpose_, pdbdata_string );
      } 
      catch ( utility::excn::EXCN_Base& excn ) {
				excn.show( std::cerr );
        // return error or empty pose
        outputpose_ = inputpose_;
        return_error_to_server( server_info );
        return false;
      }
      // since a null operation would leave the pose unchanged: 
      outputpose_ = inputpose_;
      
      initialized_ = true;

      return true;
    }
   
    bool return_error_to_server( const ServerInfo & server_info ){ 
      // do some basic measurements
      Json::Value energies;
      energies[ "score" ] = 0; 
      energies[ "irms" ] = 0; 

      // assemble the json structure needed.
      Json::Value root;
      
      // set the output values
      root["parental_key"] = Json::Value( key_ );   // the key and hash become parental_  key and hash
      root["parental_hash"] = Json::Value( hash_ ); // so we can keep track of the geneaology of structures
      root["operation"] = Json::Value( operation_ ); // so we can keep track of the geneaology of structures
      root["taskname"] = Json::Value( taskname_ ); // so we can keep track of the geneaology of structures
      root["pdbdata"] = Json::Value( "" ); // the PDB data itself of course
      root["workerinfo"] = Json::Value(  "Rosetta Backend 0.1" ); 
      root["energies"] = energies; 
      
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
      
      outputpose_.dump_pdb( "output.pdb" );

      // assemble the json structure needed.
      Json::Value root;
      
      // set the output values
      root["parental_key"] = Json::Value( key_ );   // the key and hash become parental_  key and hash
      root["parental_hash"] = Json::Value( hash_ ); // so we can keep track of the geneaology of structures
      root["operation"] = Json::Value( operation_ ); // so we can keep track of the geneaology of structures
      std::cout << "Finished TaskName: " << taskname_ << std::endl; 
      root["taskname"] = Json::Value( taskname_ ); // so we can keep track of the geneaology of structures
      root["pdbdata"] = Json::Value( pdbdatastream.str() ); // the PDB data itself of course
      root["workerinfo"] = Json::Value(  "Rosetta Backend 0.1" ); 
      root["energies"] = energies; 
      
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
  





    }

 private:

  core::pose::Pose inputpose_;
  core::pose::Pose outputpose_;
  std::string taskname_; 
  std::string hash_; 
  std::string key_; 
  std::string operation_; 
  std::string rosetta_script_;
  std::string command_;
  bool initialized_;


};

class RosettaBackend {
 public:
    RosettaBackend( ServerInfo &serverinfo) : serverinfo_(serverinfo) 
    {}

 private:

 public:

    void run(){
      do{ 
        RosettaJob newjob;
        core::Size wait_count = 0;
        
        while( !newjob.request_job_from_server( serverinfo_ ) ){
          core::Real waittime =  std::min( (double)10.0f, 0.5f* pow( (float)1.3, (float) wait_count )); // in seconds  
          
          // HACK!
          waittime = 1.5;
          
          std::cout << "No work. Waiting " << waittime << " seconds before retrying." << std::endl;
          sleep( waittime );
          wait_count ++;
        };

        // execute the run
        newjob.run();

        // now return the results to the server
        newjob.return_results_to_server( serverinfo_ );
      } while (true);

    };


 private:
  ServerInfo serverinfo_;
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
  RosettaBackend backend( server );


  backend.run();

 	return 0;
}

