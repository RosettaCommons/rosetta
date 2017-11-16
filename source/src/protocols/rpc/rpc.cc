// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rpc/rpc.cc
/// @brief
/// @author Mike Tyka

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <protocols/rpc/rpc.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <core/init/init.hh>
#include <protocols/moves/Mover.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>

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
#include <utility/json_spirit/json_spirit_error_position.h>

#include <ObjexxFCL/string.functions.hh>

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

			json_energies.emplace_back( name, float(  (*emap_iter) * (*wts_iter) ) );
		} // if ( *wts_iter != 0.0 )
	} // for ( emap_iter ...)

}


// RPC Class

static basic::Tracer TR( "rpc" );
using namespace utility::json_spirit;


bool
BasicCmdLineInit::do_init(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::init::init(argc_, argv_);
	//utility::options::OptionCollection &option_collection  = initialize();
	//option_collection.load( argc_, argv_, false);
	return true;
}


JSON_RPC::JSON_RPC(const std::string &msg, bool capture_tracer, BasicInit *basic_init  ):
	capture_tracer_(capture_tracer),
	starttime_(0),
	endtime_(0),
	basic_init_( basic_init )
{
	msg_ = msg;
	unpack( msg_ );
}

JSON_RPC::JSON_RPC( JSON_RPC const & json_rpc) : ReferenceCount(json_rpc) {
	(*this) = json_rpc;
}

JSON_RPC const & JSON_RPC::operator= ( JSON_RPC const & json_rpc )
{
	msg_                    = json_rpc.msg_;
	pdbdata_string_         = json_rpc.pdbdata_string_;
	inputpose_              = json_rpc.inputpose_;
	outputpose_             = json_rpc.outputpose_;
	xmlscript_              = json_rpc.xmlscript_;
	command_                = json_rpc.command_;
	capture_tracer_         = json_rpc.capture_tracer_;
	starttime_              = json_rpc.starttime_;
	endtime_                = json_rpc.endtime_;
	parsed_json_            = json_rpc.parsed_json_ ;
	tracer_output_stream_   <<  json_rpc.tracer_output_stream_.str();
	return *this;
}

core::Real JSON_RPC::get_fa_score() {
	core::scoring::ScoreFunctionOP fascorefxn = core::scoring::get_score_function();
	return (*fascorefxn)(outputpose_);
}

core::Real JSON_RPC::get_irms() const {
	//core::scoring::calpha_superimpose_pose( outputpose_, inputpose_ );
	return core::scoring::rmsd_with_super(inputpose_, outputpose_  , core::scoring::is_protein_backbone_including_O ) ;
}


void JSON_RPC::unpack( const std::string &msg ){
	output_capture_start();  // Make sure we're catching all the messages and errors

	try{
		std::cout << "Reading the message" << std::endl;
		parsed_json_ = read_mObject( msg );

		if ( !has_value(parsed_json_, "command") ) {
			throw utility::excn::EXCN_Msg_Exception("RPC calls must provide command field");
		}
		command_ = get_string(parsed_json_, "command");
		std::cout << "Rosetta command: " << command_ << std::endl;
		xmlscript_ = "";
		if ( has_value( parsed_json_, "xmlscript" ) ) {
			mObject parsed_xmlscript = get_mObject( parsed_json_, "xmlscript" );
			xmlscript_ = get_string_or_empty( parsed_xmlscript, "content" );
			std::cout << "XML script: " << xmlscript_ << std::endl;
		}

		if ( !has_value(parsed_json_, "pdbdata") ) {
			throw utility::excn::EXCN_Msg_Exception("RPC calls must provide pdbdata field");
		}
		pdbdata_string_ = get_string(parsed_json_, "pdbdata");

		// do a re-init if user has set a basic_init functor
		if ( basic_init_ != nullptr ) basic_init_->do_init();

		// Load any flags that were given as part of this job.
		std::cout << "Initializing options: " << command_ << std::endl;

		if ( has_value(parsed_json_, "user_flags") ) {
			mObject parsed_user_flags = get_mObject(parsed_json_, "user_flags");
			load_new_set_of_user_flags( parsed_user_flags );
		}

		std::cout << "Loading flags file: " << std::endl;
		if ( has_value(parsed_json_, "flags_file") ) {
			std::string flags_file = get_string(parsed_json_, "flags_file");
			std::cerr << "Flags file: " << flags_file << std::endl;
			//load_user_flag_file( flags_file );
		}

		// Load any files that were given as part of this job.
		std::cout << "Loading user files: " << std::endl;
		if ( has_value(parsed_json_, "user_files") ) {
			load_new_set_of_virtual_files(  get_mArray(parsed_json_, "user_files") );
		}
	}
catch( utility::excn::EXCN_Msg_Exception &excn ){
	output_capture_stop();  // Make sure we're catching the following message
	throw;
}
catch( std::string &s ){
	// convert Exception into a rosetta exception
	output_capture_stop();
	throw utility::excn::EXCN_Msg_Exception( tracer() + s );
}
catch( utility::json_spirit::Error_position &ep ){
	// convert Exception into a rosetta exception
	output_capture_stop();
	throw utility::excn::EXCN_Msg_Exception( tracer() + ep.reason_ + " Line: " + ObjexxFCL::string_of( ep.line_ ) + " Col: " + ObjexxFCL::string_of( ep.column_ ) );
}
catch( ... ){
	output_capture_stop();
	throw utility::excn::EXCN_Msg_Exception( tracer() + "Unknown Exception happened during unpacking of JSON data" );
}
}

void JSON_RPC::run(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// start capturing Tracer outputs
	output_capture_start();
	starttime_ = time(nullptr);

	// finally load in the provided PDB
	std::cout << "Loading PDB file: " << std::endl;
	try{
		core::import_pose::pose_from_pdbstring( inputpose_, pdbdata_string_ );

		outputpose_ = inputpose_;

		// Leave some trace
		std::cout << "Executing: " << command_ << std::endl;
		TR << "Executing: " << command_ << std::endl;

		// scoring happens automatically at the end so this is a NOP
		if ( command_ == "score" ) {}


		// this covers a huge variety of possible operations and move movers are plugged into this system by now.
		if ( command_ == "xmlscript" ) {
			if ( xmlscript_ == "" ) {
				throw utility::excn::EXCN_Msg_Exception("RPC error: XML script is empty! " );

			}
			utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
			provider->add_input_file( "script.xml", xmlscript_ );
			protocols::rosetta_scripts::RosettaScriptsParser rsp;
			protocols::moves::MoverOP protocol;
			rsp.generate_mover_from_pose( outputpose_ , protocol, true, "script.xml" );
			protocol->apply( outputpose_ );
		}

		// do a fast core::pose action - this mostly for providing access to simple core::.. functionality via Native Client.
		if ( command_ == "poseop" ) {

			mArray poseops = get_mArray( parsed_json_, "poseop" );

			for ( mArray::const_iterator it = poseops.begin();
					it != poseops.end(); ++it ) {

				std::cout << it->type() << std::endl;
				if ( it->type() != obj_type ) {
					throw utility::excn::EXCN_Msg_Exception("JSON error: expected an object for poseop array member'");
				};

				const mObject &poseop_params = it->get_obj();

				std::string poseop = get_string(poseop_params, "op");

				if ( false ) ;
				else if ( poseop == "set_phi" )      outputpose_.set_phi(        get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_psi" )      outputpose_.set_psi(        get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_omega" )    outputpose_.set_omega(      get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_alpha" )    outputpose_.set_alpha(      get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_beta" )     outputpose_.set_beta(       get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_gamma" )    outputpose_.set_gamma(      get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_delta" )    outputpose_.set_delta(      get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_epsilon" )  outputpose_.set_epsilon(    get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_zeta" )     outputpose_.set_zeta(       get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));
				else if ( poseop == "set_chi" )      outputpose_.set_chi(        get_int( poseop_params, "chino" ), get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ) );
				else if ( poseop == "set_chi_nucl" ) outputpose_.set_chi(        get_int( poseop_params, "seqpos" ),  get_real( poseop_params, "value" ));


				//            bool empty() const;
				//            bool is_fullatom() const;
				//            bool is_centroid() const;
				//
				//            Size total_residue() const;
				//            Size n_residue() const;
				//            Size num_jump() const;
				//
				//            chemical::AA aa( Size const seqpos) const;
				//            char secstruct( Size const seqpos ) const;
				//            std::string secstruct() const;
				//            std::string sequence() const;
				//            std::string annotated_sequence( bool show_all_variants = false ) const;
				//            std::string chain_sequence( core::Size const chain_in ) const;
				//
				//            Real phi( Size const seqpos ) const;
				//            Real psi( Size const seqpos ) const;
				//            Real omega( Size const seqpos ) const;
				//            Real alpha( Size const pos ) const;
				//            Real beta( Size const seqpos ) const;
				//            Real gamma( Size const seqpos ) const;
				//            Real delta( Size const pos ) const;
				//            Real epsilon( Size const seqpos ) const;
				//            Real zeta( Size const seqpos ) const;
				//            Real chi( int const chino, Size const seqpos) const;
				//            Real chi( Size const seqpos ) const;

			}

		}

		// here largely as examples - use XML scripts for this.
		//        if( command_ == "relax" ){
		//          protocols::moves::MoverOP protocol = protocols::relax::generate_relax_from_cmd();
		//          protocol->apply( outputpose_ );
		//        }
		//        if( command_ == "loophash" ){
		//          protocols::loophash::LoopHashRelaxProtocolOP lh_protocol = new protocols::loophash::LoopHashRelaxProtocol( loop_hash_library );
		//          lh_protocol->manual_call( outputpose_ );
		//        }
	}
catch( utility::excn::EXCN_Msg_Exception &excn ){
	std::cerr << "EXCEPTION: " << excn.msg() << std::endl; // print the exception message to the Error stream.
	TR.Error << "EXCEPTION: " << excn.msg() << std::endl; // print the exception message to the Error stream.
	endtime_ = time(nullptr);
	output_capture_stop();  // Make sure we're catching the following message
	throw;
}
catch( std::string &s ){
	// convert Exception into a rosetta exception
	endtime_ = time(nullptr);
	output_capture_stop();
	throw utility::excn::EXCN_Msg_Exception( tracer() + s );
}
catch( utility::json_spirit::Error_position &ep ){
	// convert Exception into a rosetta exception
	endtime_ = time(nullptr);
	output_capture_stop();
	throw utility::excn::EXCN_Msg_Exception( tracer() + ep.reason_ + " Line: " + ObjexxFCL::string_of( ep.line_ ) + " Col: " + ObjexxFCL::string_of( ep.column_ ) );
}
catch( ... ){
	endtime_ = time(nullptr);
	output_capture_stop();
	throw utility::excn::EXCN_Msg_Exception( tracer() + "Unknown Exception happened during execution of json RPC call " );
}

	endtime_ = time(nullptr);
}


void JSON_RPC::output_capture_start(){
	basic::TracerImpl::set_new_final_stream( &tracer_output_stream_ );
}

void JSON_RPC::output_capture_stop(){
	basic::TracerImpl::set_default_final_stream();
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
	for ( auto const & json_user_flag : json_user_flags ) {
		if ( json_user_flag.second.type() != obj_type ) {
			throw utility::excn::EXCN_Msg_Exception("JSON error: expected an object for user_flag member:'" + json_user_flag.first );
		};
		const mObject &flag = json_user_flag.second.get_obj();
		if ( has_value( flag, "value" ) ) {
			std::string flagname = json_user_flag.first;
			std::string flagvalue = get_string( flag, "value" );
			std::cout << "User flag: " << flagname << " = " << flagvalue << std::endl;
			user_flags.push_back(  flagname + " " + flagvalue );
		}
	}

	basic::options::option.load( user_flags, false);
}

void JSON_RPC::load_new_set_of_virtual_files( const mArray &json_user_files , bool clear_previous ){
	std::cerr << __FILE__ << __LINE__ << std::endl;
	utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
	std::cerr << __FILE__ << __LINE__ << std::endl;
	if ( clear_previous ) provider->clear_input_files();
	std::cerr << __FILE__ << __LINE__ << std::endl;
	for ( auto const & json_user_file : json_user_files ) {
		std::cerr << __FILE__ << __LINE__ << std::endl;
		std::cout << json_user_file.type() << std::endl;
		std::cerr << __FILE__ << __LINE__ << std::endl;
		if ( json_user_file.type() != obj_type ) {
			throw utility::excn::EXCN_Msg_Exception("JSON error: expected an object for user_file member:'");
		};
		const mObject &flag = json_user_file.get_obj();
		if ( !has_value( flag, "filename" ) ) {
			throw utility::excn::EXCN_Msg_Exception("JSON error: Syntax error in user_files field: 'filename' missing ");
		}
		if ( !has_value( flag, "contents" ) ) {
			TR.Error << "Filename " << get_string( flag, "filename" ) << " is missing content field - assuming empty file " << std::endl;
		}

		std::string filename = get_string( flag, "filename" );
		std::string contents = get_string_or_empty( flag, "contents" );
		TR << "Loading virtual file: " << filename << " " << contents.size() << " bytes" << std::endl;
		provider->add_input_file( filename, contents );
	}
}


}
}
