// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rpc/rpc.cc
/// @brief
/// @author Mike Tyka

#include <basic/Tracer.hh>
#include <basic/options/option.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <protocols/rpc/rpc.hh>
#include <protocols/jd2/DockDesignParser.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/DockDesignParser.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loophash/LoopHashRelaxProtocol.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/assert.hh>
#include <utility/excn/Exceptions.hh>
#include <iostream>
#include <fstream>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/inline_file_provider.hh>
#include <boost/function.hpp>

#include <utility/json_spirit/json_spirit_value.h>
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_tools.hh>

namespace protocols {
namespace rpc {


  // Helper functions

void pose_energies_to_json( core::pose::Pose const & pose, utility::json_spirit::Object &json_energies ) {

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

      json_energies.push_back( utility::json_spirit::Pair( name, float(  (*emap_iter) * (*wts_iter) ) ) ); 
    } // if ( *wts_iter != 0.0 )
  } // for ( emap_iter ...)

}


  // RPC Class

  static basic::Tracer TR("rpc");
 
  using namespace utility::json_spirit;
        
  core::Real JSON_RPC::get_fa_score() {
    core::scoring::ScoreFunctionOP fascorefxn = core::scoring::getScoreFunction();
    return (*fascorefxn)(outputpose_); 
  }

  core::Real JSON_RPC::get_irms() const {
    //core::scoring::calpha_superimpose_pose( outputpose_, inputpose_ );
    return core::scoring::rmsd_with_super(inputpose_, outputpose_  , core::scoring::is_protein_backbone_including_O ) ;
  }


  void JSON_RPC::unpack( const std::string &msg ){
     
     output_capture_start();  // Make sure we're catching all the messages and errors 
    
     try{
       mObject parsed_json = read_mObject( msg );

       command_ = get_string(parsed_json, "command");
       std::cout << "Rosetta command: " << command_ << std::endl;
       xmlscript_ = "";
       if( has_value( parsed_json, "xmlscript" ) ){
         mObject parsed_xmlscript = get_mObject( parsed_json, "xmlscript" );
         xmlscript_ = get_string_or_empty( parsed_xmlscript, "contents" );
       }

       std::cout << "XML script: " << xmlscript_ << std::endl;

       std::string pdbdata_string = get_string(parsed_json, "pdbdata");
       // Load any flags that were given as part of this job.
       std::cout << "Initializing options: " << command_ << std::endl;
       
       if ( has_value(parsed_json, "user_flags") ){
         load_new_set_of_user_flags( get_mObject(parsed_json, "user_flags") );
       }
       
       if ( has_value(parsed_json, "flags_file") ){
         std::string flags_file = get_string(parsed_json, "flags_file");
         std::cerr << "Flags file: " << flags_file << std::endl; 
         load_user_flag_file( flags_file );
       }
       
       // Load any files that were given as part of this job.
       load_new_set_of_virtual_files(  get_Array(parsed_json, "user_files") );

       // finally load in the provided PDB
       core::import_pose::pose_from_pdbstring( inputpose_, pdbdata_string );
       
       outputpose_ = inputpose_; 
    }
    catch( ... ){
      output_capture_stop(); 
      throw;
    }
  }

  void JSON_RPC::run(){

      // start capturing Tracer outputs
      output_capture_start();
      starttime_ = time(NULL);

      try{
        // Leave some trace
        std::cout << "Executing: " << command_ << std::endl;
        TR << "Executing: " << command_ << std::endl;
       
        // scoring happens automatically at the end so this is a NOP
        if( command_ == "score" ){
        }
        
        
        // this covers a huge variety of possible operations and move movers are plugged into this system by now.
        if( command_ == "xmlscript" ){
          utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
          provider->add_input_file( "script.xml", xmlscript_ );
          protocols::jd2::DockDesignParser ddp;
          protocols::jd2::JobCOP job;
          protocols::moves::MoverOP protocol;
          ddp.generate_mover_from_pose( job, outputpose_ , protocol, true, "script.xml" );
          protocol->apply( outputpose_ );
        }
        
        // do a fast core::pose action - this mostly for providing access to simple core::.. functionality via Native Client.
        if( command_ == "poseop" ){
          

        }
        
        // here largely as examples - use XML scripts for this.
        if( command_ == "relax" ){
          //protocols::moves::MoverOP protocol = protocols::relax::generate_relax_from_cmd();
          //protocol->apply( outputpose_ );
        }
        if( command_ == "loophash" ){
          //protocols::loophash::LoopHashRelaxProtocolOP lh_protocol = new protocols::loophash::LoopHashRelaxProtocol( loop_hash_library );
          //lh_protocol->manual_call( outputpose_ );
        }
      } 
      catch( utility::excn::EXCN_Msg_Exception &excn ){
        std::cerr << "EXCEPTION: " << excn.msg() << std::endl; // print the exception message to the Error stream.
        TR.Error << "EXCEPTION: " << excn.msg() << std::endl; // print the exception message to the Error stream.
        endtime_ = time(NULL);
        output_capture_stop();  // Make sure we're catching the following message 
        throw; 
      }
      catch(...) {
        TR.Error << "UNKNOWN ERROR OCCURED DURING RUN" << std::endl;
        endtime_ = time(NULL);
        output_capture_stop();  // Make sure we're catching the following message 
        throw; 
      }
      
  }





  void JSON_RPC::output_capture_start(){
      basic::set_new_final_channel( &tracer_output_stream_ );
  }

  void JSON_RPC::output_capture_stop(){
      basic::set_default_final_channel();
  }

  void JSON_RPC::output_capture_clear(){
    tracer_output_stream_.str(std::string());
  }

  void JSON_RPC::load_user_flag_file( const std::string &flags_file ){
    std::stringstream options_file;
    options_file << flags_file << std::endl; // ensure the last line has a CR/lF at the end. 
    basic::options::option.load_options_from_stream( options_file );
  }

  void JSON_RPC::load_new_set_of_user_flags( const mObject &json_user_flags ){
    std::vector < std::string > user_flags;
    for( mObject::const_iterator it = json_user_flags.begin(); 
         it != json_user_flags.end(); ++it ){
      if( it->second.type() != obj_type ){
        throw utility::excn::EXCN_Msg_Exception("JSON error: expected an object for user_flag member:'" + it->first );
      }; 
      const mObject &flag = it->second.get_obj();      
      if( has_value( flag, "value" ) ){
        std::string flagname = it->first;
        std::string flagvalue = get_string( flag, "value" );
        std::cout << "User flag: " << flagname << " = " << flagvalue << std::endl;
        user_flags.push_back(  flagname + " " + flagvalue );
      }
    }

    basic::options::option.load( user_flags, false);
  }

  void JSON_RPC::load_new_set_of_virtual_files( const mArray &json_user_files , bool clear_previous ){
    utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
    if( clear_previous ) provider->clear_input_files();
    for( mArray::const_iterator it = json_user_files.begin(); 
         it != json_user_files.end(); ++it ){
      std::cout << it->type() << std::endl;
      if( it->type() != obj_type ){
        throw utility::excn::EXCN_Msg_Exception("JSON error: expected an object for user_file member:'");
      }; 
      const mObject &flag = it->get_obj();      
      if( !has_value( flag, "filename" ) ) {
        throw utility::excn::EXCN_Msg_Exception("JSON error: Syntax error in user_files field: 'filename' missing "); 
      }
      if( !has_value( flag, "contents" ) ) {
        TR.Error << "Warning: Filename " << get_string( flag, "filename" ) << " is missing content field - assuming empty file " << std::endl; 
      }
      
      std::string filename = get_string( flag, "filename" );
      std::string contents = get_string_or_empty( flag, "contents" );
      TR << "Loading virtual file: " << filename << " " << contents.size() << " bytes" << std::endl;
      provider->add_input_file( filename, contents );
    }
  }




}
}

