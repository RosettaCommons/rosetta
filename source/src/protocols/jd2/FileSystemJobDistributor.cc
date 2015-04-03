// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/FileSystemJobDistributor.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Simplest class FileSystemJobDistributor
/// @author Andrew Leaver-Fay
/// @author Steven Lewis smlewi@gmail.com

// Unit headers
#include <protocols/jd2/FileSystemJobDistributor.hh>

//Package headers
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/exit.hh>

///C++ headers
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <utility/vector1.hh>


using basic::Warning;
static thread_local basic::Tracer TR( "protocols.jd2.FileSystemJobDistributor" );

namespace protocols {
namespace jd2 {

FileSystemJobDistributor::FileSystemJobDistributor() :
	JobDistributor(),
	extension_(".in_progress"),
	//	next_job_to_try_assigning_( 1 ),
	retry_count_( 0 )
{
	if ( basic::options::option[ basic::options::OptionKeys::out::path::pdb ].user() ){
		path_ = basic::options::option[ basic::options::OptionKeys::out::path::pdb ]().path();
	} else if( basic::options::option[ basic::options::OptionKeys::out::path::all ].user() ){
		path_ = basic::options::option[ basic::options::OptionKeys::out::path::all ]().path();
	}else{
		path_ = "";
	}

}

///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
FileSystemJobDistributor::~FileSystemJobDistributor()
{}

void FileSystemJobDistributor::restart() {
	//	next_job_to_try_assigning_ = 1;
	retry_count_ = 0;
	JobDistributor::restart();
}


core::Size get_min_nstruct_index_checkpoint_file(){
	using namespace std;
	using namespace basic::options;
	if ( option[ basic::options::OptionKeys::jd2::checkpoint_file ].user()) {
		std::string saved_nstruct_number_file = option[ basic::options::OptionKeys::jd2::checkpoint_file ]();
		ifstream f( saved_nstruct_number_file.c_str(), ios::in );
		if( !f.good() ) return 0;
		core::Size const begin = f.tellg();
		f.seekg( 0, ios::end );
		core::Size const end = f.tellg();
		if( end - begin == 0 ) return 0;
		f.seekg( 0, ios::beg );
		std::string line;
		getline( f, line );
		std::istringstream line_stream( line );
		core::Size nstruct_index;
		line_stream >> nstruct_index;
		return nstruct_index;
	}
	return 0;
}


core::Size
FileSystemJobDistributor::get_new_job_id()
{

	if( basic::options::option[ basic::options::OptionKeys::out::overwrite ].value()  &&
  		basic::options::option[ basic::options::OptionKeys::run::multiple_processes_writing_to_one_directory ].value() ){
		utility_exit_with_message("ambiguous, cannot have both -out::overwrite and -run::multiple_processes_writing_to_one_directory");
	}

	Jobs const & jobs( get_jobs() );
	JobOutputterOP outputter = job_outputter();
	core::Size next_job_to_try_assigning( current_job_id() + ( retry_count_ > 0 ? 0 : 1 ) );
	next_job_to_try_assigning = std::max( next_job_to_try_assigning, get_min_nstruct_index_checkpoint_file() );
	while ( next_job_to_try_assigning <= jobs.size() ) {
		if ( jobs[ next_job_to_try_assigning ]->bad() ) {
			++next_job_to_try_assigning;
		} else if ( outputter->job_has_completed( jobs[ next_job_to_try_assigning ] ) &&
				 !basic::options::option[ basic::options::OptionKeys::out::overwrite ].value() ) {
			++next_job_to_try_assigning;
		} else if ( basic::options::option[ basic::options::OptionKeys::run::multiple_processes_writing_to_one_directory ].value() ) {
            std::string const next_job_output_name( temporary_file_name( jobs[ next_job_to_try_assigning ] ) );
			if ( utility::file::file_exists( next_job_output_name) ) {
				++next_job_to_try_assigning;
			} else {
				/// create a temporary file in this directory to indicate that the job
				/// has started to other processes.  This is not in any way safe.
				/// nor does it guarantee that no two processes will attempt to perform the same job.
				/// If two processes do attempt the same job, then the one that completes last will
				/// have its data saved, and the other will be overwritten.
				utility::io::ozstream outfile( next_job_output_name ); //opens automatically
				outfile << "temp" << std::flush;
				outfile.close();
				break;
			}
		} else {
			break;
		}
	} //while

	if ( next_job_to_try_assigning <= jobs.size() ) {
		//	core::Size job_to_assign = next_job_to_try_assigning;
		//	++next_job_to_try_assigning;
		return next_job_to_try_assigning;
	}

	// indicate that no jobs remain
	return 0;
}

void
FileSystemJobDistributor::mark_current_job_id_for_repetition()
{
	using namespace basic::options;
// 	if( current_job_id() != next_job_to_try_assigning - 1 ){
// 		std::cerr << current_job_id() << "  " << next_job_to_try_assigning << std::endl;
// 		runtime_assert( current_job_id() == next_job_to_try_assigning - 1 );
// 	}
	++retry_count_;
	if( (int)retry_count_ <= (int)basic::options::option[ basic::options::OptionKeys::run::max_retry_job ].value() ) {
		// 		--next_job_to_try_assigning;
		if ( !option[ OptionKeys::run::write_failures ]() ){
			clear_current_job_output(); // DONT DO THIS WHEN WRITING FAILURES COS IT LOOSES THE PAIR DATA
		}
	}else{
		retry_count_=0; // reset
		TR << "Too many retries (max_retry_job = " << basic::options::option[ basic::options::OptionKeys::run::max_retry_job ].value() << ") " << std::endl;
	}
}

/// @details this function handles the FAIL_BAD_INPUT mover status by removing other jobs with the same input from consideration
void
FileSystemJobDistributor::remove_bad_inputs_from_job_list()
{
	//we should only fail on a job which was the first of its type - if it fails on nstruct=2 with BAD_INPUT then why did it not fail on nstruct=1?
	//runtime_assert( current_job()->nstruct_index() == 1 );
	if( current_job()->nstruct_index() != 1 ){
		Warning() << "A job reported bad input, but was not the first input of its type!  You should figure out why the first one passed if later ones failed!" << std::endl;
	}

	//	std::string const & current_input_tag(current_job()->inner_job()->input_tag());
	//std::string const & current_native_tag(current_job()->inner_job()->native_tag());

	TR << "job failed, reporting bad input; other jobs of same input will be canceled: "
		 << job_outputter()->output_name( current_job() ) << std::endl;
	mark_job_as_bad( current_job_id() );
	// this latter stuff requires that jobs come in sequence of input... some application might prefer to
	// reshuffle ...

// 	Jobs const & jobs( get_jobs() );
// 	core::Size next_job_to_try_assigning( current_job_id() );
// 	while(next_job_to_try_assigning <= jobs.size() && //MUST BE FIRST for c++ shortcut logical evaluation
// 				jobs[next_job_to_try_assigning]->inner_job()->input_tag() == current_input_tag) {
// 		//&&jobs[next_job_to_try_assigning]->inner_job()->native_tag() == current_native_tag){
// 		TR << "job canceled without trying due to previous bad input: "
// 			 << job_outputter()->output_name( jobs[next_job_to_try_assigning] ) << std::endl;
// 		mark_job_as_bad( next_job_to_try_assigning );
// 		++next_job_to_try_assigning;
// 	}

}

/// @details when multiple processes are writing to the same directory, then
/// after a job completes, the "in_progress" file for the job must be deleted.
void
FileSystemJobDistributor::current_job_finished()
{
	delete_in_progress_files();
	using namespace basic::options;
	if ( option[ basic::options::OptionKeys::jd2::checkpoint_file ].user() ) {
		std::string saved_nstruct_number_file = option[ basic::options::OptionKeys::jd2::checkpoint_file ]();
		std::ofstream f;
	    f.open( saved_nstruct_number_file.c_str(), std::ios::out );
	    if( !f.good() ) utility_exit_with_message( "Unable to open MC checkpointing file " + saved_nstruct_number_file );
		f<<(1+current_job()->nstruct_index());
	    f.close();
	}
}


void
FileSystemJobDistributor::go( protocols::moves::MoverOP mover )
{
	setup_system_signal_handler();
	go_main( mover );
	remove_system_signal_handler();
}


std::string const
FileSystemJobDistributor::temporary_file_name( JobCOP job ) const {
	return path_ + job_outputter()->output_name( job ) + extension_;
}


void
FileSystemJobDistributor::job_succeeded(core::pose::Pose & pose, core::Real run_time, std::string const & tag)
{
	JobDistributor::job_succeeded( pose, run_time, tag );
 	retry_count_ = 0;
	return;
}

void
FileSystemJobDistributor::job_failed( core::pose::Pose & pose, bool will_retry )
{
	using namespace basic::options;

	JobDistributor::job_failed( pose, will_retry );

	if ( option[ OptionKeys::run::write_failures ]() ){
		current_job()->set_status_prefix("C");
		job_succeeded( pose, 0, "" );
	}

	if( !will_retry ) retry_count_ = 0;
}

void FileSystemJobDistributor::handle_interrupt()
{
	delete_in_progress_files();
}

void FileSystemJobDistributor::delete_in_progress_files()
{
	if ( basic::options::option[ basic::options::OptionKeys::run::multiple_processes_writing_to_one_directory ].value() ) {
		Jobs const & jobs( get_jobs() );
		JobOutputterCOP outputter = job_outputter();
		std::string const output_name( temporary_file_name( jobs[ current_job_id() ] ) );
		utility::file::file_delete( output_name.c_str() );
	}
}


}//jd2
}//protocols
