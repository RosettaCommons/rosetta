// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rpc/JSON_RPC.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_rpc_JSON_RPC_hh
#define INCLUDED_protocols_rpc_JSON_RPC_hh

#include <protocols/rpc/rpc.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/json_spirit/json_spirit_value.h>
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/json_spirit/json_spirit_writer.h>

#include <string>
#include <vector>

#include <utility/vector1.hh>

namespace protocols {
namespace rpc {


void pose_energies_to_json( core::pose::Pose const & pose, utility::json_spirit::Object &json_energies );

// this is a virtual functor which is going to act basically as a callback or function pointer but cleaner.
class BasicInit {
public:
	BasicInit() {}

	virtual ~BasicInit(){};

	virtual bool do_init() { return false; };
};

// this is a functor that just clears the current options and sets
class BasicCmdLineInit: public protocols::rpc::BasicInit
{
public:
	BasicCmdLineInit()
	{
		argc_ = 0;
		empty_ = &std::string("")[0];
		argv_ = &empty_;
	}

	BasicCmdLineInit( int argc, char * argv [] ) :
		argc_(argc),
		empty_(&std::string("")[0]),
		argv_(argv)
	{}

	virtual bool do_init();

private:
	int argc_;
	char * empty_;
	char ** argv_ ;
};

// This is an interface to rosetta tailored for remote procedure calls through HTTP requests or via Native Client.
// The idea is that everything needed for a rosetta run is encoded in a JSON object (including an input PDB or PDBs,
// the command, any parameters, etc etc. ). This information is received via some channel and then used to execute
// the command. The results (another PDB or other data) is then sent back to the client.

class JSON_RPC : public utility::pointer::ReferenceCount {
public:

	JSON_RPC(std::string const &msg, bool capture_tracer = true,  BasicInit *basic_init = NULL );

	virtual ~JSON_RPC() {}

	JSON_RPC( JSON_RPC const & json_rpc);

	JSON_RPC const & operator = ( JSON_RPC const & json_rpc );

	virtual void run();

	long runtime() const { return endtime_ - starttime_; }

	std::string tracer() const { return tracer_output_stream_.str(); }

	const core::pose::Pose& outputpose() const { return outputpose_; }

	core::Real get_fa_score();

	core::Real get_irms() const;

private: // member functions
	void unpack( const std::string &msg );

	void output_capture_start();
	void output_capture_stop();
	void output_capture_clear();

	void load_user_flag_file( const std::string &flags_file );
	void load_new_set_of_user_flags(  const utility::json_spirit::mObject &json_user_flags );
	void load_new_set_of_virtual_files(  const utility::json_spirit::mArray &json_user_files, bool clear_previous = true );

private: // member variables
	std::string       msg_;
	std::string       pdbdata_string_;

	core::pose::Pose  inputpose_;
	core::pose::Pose  outputpose_;
	std::string       xmlscript_;
	std::string       command_;

	bool              capture_tracer_;
	std::stringstream tracer_output_stream_;
	long              starttime_;
	long              endtime_;

	utility::json_spirit::mObject parsed_json_;

	BasicInit         *basic_init_;
};

}
}

#endif

