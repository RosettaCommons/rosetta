// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @details
///
///
/// @author Oliver Lange
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


// Unit Headers
#include <protocols/simple_filters/ExternalEvaluator.hh>
#include <basic/options/option.hh>

// Package Headers

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <numeric/random/random.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <iostream>
// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/vector1.hh>


#ifdef  __native_client__
#define system(a) 1
#endif

// C++ headers

static THREAD_LOCAL basic::Tracer tr( "protocols.simple_filter.ExternalEvaluator" );


namespace protocols {
namespace simple_filters {

using namespace core;
using namespace std;


ExternalEvaluator::ExternalEvaluator( std::string tag, std::string command )
: evaluation::SingleValuePoseEvaluator<core::Real>( tag ),
	command_( command )
{
	// this probably shouldn't go on BOINC
#ifndef WIN32
#ifndef BOINC
#ifndef  __native_client__

	using namespace basic::options;
	if ( !option[ OptionKeys::out::path::scratch ].user() ) {
		tr.Warning << "******************************************************************************************************\n"
			<< "                       no scratch dir defined  \n" << "   use -out:path:scratch \n"
			<< "******************************************************************************************************\n" << endl;
		utility_exit_with_message(" define your scratch dir with -out:path:scratch " );
		scratch_dir_ = "/scratch/USERS/";
	} else {
		scratch_dir_ = option[ OptionKeys::out::path::scratch ]();
	}
	if ( !utility::file::file_exists( scratch_dir_.c_str() ) ) mkdir( scratch_dir_.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

	string const tmp_file_name ("_ExternalEvaluator_"+name( 1 ));

	string dir = "./";
	if ( option[ OptionKeys::out::file::silent ].user() ) {
		dir = option[ OptionKeys::out::file::silent ]();
	}

	//also use random number!

	string sub_work_dir;
	{// if processes write silent files to subdirectory ... get this path
		string dircmd( "echo `dirname "+dir+"` | sed s@/@_@g");
		FILE* get_dir = popen(dircmd.c_str(),"r");
		char buf[500];
		if ( fgets(buf, 500, get_dir) == NULL ) {
			tr.Error << "Read error!" << std::endl;
		}
		buf[strlen(buf)-1]='\0';
		sub_work_dir = string( buf )+"_"+tmp_file_name; //get rid of newline
		pclose( get_dir );
	}
	sub_work_dir = sub_work_dir+"_"+ ObjexxFCL::string_of( numeric::random::rg().random_range(0, 999999) );
	// set npes and rank based on whether we are using MPI or not
#ifdef USEMPI
	int rank_;
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &rank_ ) );
	sub_work_dir = "mpi_" + ObjexxFCL::string_of( rank_ ) + sub_work_dir;
#endif


	{ // get our working path
		string dircmd( "pwd | sed s@/@_@g");
		FILE* get_dir = popen(dircmd.c_str(),"r");
		char buf[500];
		if ( fgets(buf, 500, get_dir) == NULL ) {
			tr.Error << "Read error!" << std::endl;
		}
		buf[strlen(buf)-1]='\0';
		string tmp_id = string( buf );
		work_dir_ = tmp_id.substr( max(0, (int)tmp_id.size()-40) );
		pclose( get_dir );
	}


	work_dir_ = scratch_dir_+work_dir_;
	tr.Info << "create scratch space... : " << work_dir_ << std::endl;
	if ( !utility::file::file_exists( work_dir_.c_str() ) ) mkdir(work_dir_.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
	work_dir_ = work_dir_+"/"+sub_work_dir;
	tr.Info << "create scratch space... : " << work_dir_ << std::endl;
	if ( !utility::file::file_exists( work_dir_.c_str() ) ) mkdir(work_dir_.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
#endif
#endif
#endif
}


core::Real ExternalEvaluator::apply( core::pose::Pose& pose ) const {
	using namespace basic::options;

	if ( !applicable( pose ) ) return 99999;

	string const command_buf( "DIR="+work_dir_+"; cd $DIR; "+command_+"; cd -");

	string const pose_file_name( work_dir_+string("/__POSE.pdb") );
	string const result_file_name( work_dir_+string("/__RESULT") );
	std::ofstream pose_stream(pose_file_name.c_str() ); //make sure not the MPI FileBuf is used !
	if ( !pose_stream ) tr.Error << "can't write pose to file " << pose_file_name << std::endl;
	pose.dump_pdb( pose_stream, pose_file_name );
	pose_stream.close();
	tr.Info << "write pose : " << pose_file_name << endl;
	tr.Info << "execute command: "<< command_buf << endl;
	int ret(system(command_buf.c_str()));
	if ( ret ) {
		tr.Warning << "Applying the external evaluator failed!" << endl;
	}

	//execl("/bin/bash","bash","-c",command_buf.c_str(), (char *)0);

	core::Real result;
	{
		utility::io::izstream result_file( result_file_name );
		result_file >> result;
	}
	// AMW: cppcheck does not see that delete_result_cmd is used in the next line?
	std::string const delete_result_cmd( "rm -f "+work_dir_+"/__RESULT" );
	ret = system( delete_result_cmd.c_str() );
	if ( ret ) {
		tr.Warning << "Deleting the file '" << work_dir_  << "/__RESULT' failed!" << endl;
	}

	tr.Debug << "obtained result: " << result << endl;
	return result;
}

ExternalEvaluator::~ExternalEvaluator() {
	std::string command = "rm -Rf "+work_dir_;
	// AMW: cppcheck does not see that command is used in the next line?
	int const ret(system(command.c_str()));
	if ( ret ) {
		tr.Warning << "Deleting work directory '" << work_dir_ << "' failed!" << endl;
	}
	//clean up
}

// core::Real ExternalEvaluator::apply( core::pose::Pose& pose_in  ) const {
//  pose::Pose pose( pose_in );

//  runtime_assert( constraints_ );
//  pose.constraint_set( constraints_ );

//  ScoreFunction scfxn;
//  scfxn.set_weight( atom_pair_constraint, 1.0 );
//  return scfxn( pose );

// }


}
}
