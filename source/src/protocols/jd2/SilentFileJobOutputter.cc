// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/SilentFileJobOutputter.cc
/// @brief
/// @author Oliver Lange

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/pose/Pose.hh>
#include <boost/foreach.hpp>

///Utility headers
#include <utility/file/FileName.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

#include <string>
#include <algorithm>

static THREAD_LOCAL basic::Tracer tr( "protocols.jd2.SilentFileJobOutputter" );

namespace protocols {
namespace jd2 {

SilentFileJobOutputter::SilentFileJobOutputter() {
	set_defaults();
	read_done_jobs();
}

SilentFileJobOutputter::~SilentFileJobOutputter() {
	//DO NOT PUT THINGS HERE - it is not guarunteed to get called - use flush below instead.
}

void SilentFileJobOutputter::flush() {
	write_all_structs();
}

void SilentFileJobOutputter::write_all_structs() {
	using std::pair;
	using utility::vector1;
	using utility::file::FileName;
	using core::io::silent::SilentStructOP;

	typedef std::map< std::string, core::io::silent::SilentFileData > SFD_MAP;
	typedef vector1< std::pair< SilentStructOP, FileName > >::iterator iter;
	SFD_MAP sfds;
	// Only write structures if the user hasn't disabled it - otherwise it totally breaks
	// the user's expectation.
	if ( !bWriteNoStructures_ ) {
		tr.Debug << "writing " << saved_structs_.size() << " structs." << std::endl;
		for ( iter it = saved_structs_.begin(), end = saved_structs_.end();
				it != end; ++it
				) {
			//tr.Debug << "writing struct " << ss->decoy_tag() << std::endl;
			//tr.Debug << "writing struct " << (*it->first)->decoy_tag() << std::endl;
			//SilentStructOP ss = it->first;
			sfds[ it->second ].add_structure( (*it->first) );
		}
		for ( SFD_MAP::iterator it = sfds.begin(); it!=sfds.end(); ++it ) {
			it->second.write_all( it->first );
		}
	}
	// very important to clear after writing!
	saved_structs_.clear();

	tr.Debug << "currently have " << saved_structs_.size() << " structs."
		<< std::endl;
	tr.flush();
}

void SilentFileJobOutputter::set_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	forced_silent_struct_type_ = "any";
	silent_file_ = option[ out::file::silent ]();
	if ( silent_file_().compare("") == 0 ) {
		utility_exit_with_message("Please supply a file name with the -out:file:silent flag. If you specify a path the -out:path:all flag, the -out:file:silent file name will be relative to it.");
	}

	if ( ! silent_file_.absolute() ) {
		silent_file_.path( option[ out::path::all ]().path() + "/" + silent_file_.path() );
		//FileName takes care of platform-specific path seperator, i.e.,  "/" or "\" ...
	}
	tr.Debug << "SilentFileJobOutputter setup for file " << silent_file_ << std::endl;

#ifdef USEMPI
	int mpi_rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);/* get current process id */
	std::string name = silent_file_;

	//Unused variables, commented out to silence warnings
	// attach mpi rank to out files
	//size_t lastslash = name.find_last_of("/\\");
	//size_t lastdot   = name.find_last_of('.');
	silent_file_ = name;
#endif

	// dd_parser should use binary silent files as a default. but score-only
	// silent files are also an option
	// and protein-style too, if that's my choice? -JRP 2012
	if ( option[ basic::options::OptionKeys::jd2::dd_parser ]() ) {
		if ( option[ out::file::silent_struct_type ]() != "score" ) {
			option[ out::file::silent_struct_type ].def( "binary" );
		}
	}

	//default is 1
	n_to_buffer_ = option[ basic::options::OptionKeys::jd2::buffer_silent_output ]();
	random_flush_frequency_ = option[ basic::options::OptionKeys::jd2::buffer_flush_frequency ]();

	bWriteIntermediateFiles_ = (
		option[ run::intermediate_scorefiles ]() ||
		option[ run::intermediate_structures ]()
	);

	bWriteIntermediateStructures_ = option[ run::intermediate_structures ]();

	bWriteNoStructures_ = false;

	if ( option[ OptionKeys::out::file::scorefile ].user() ) {
		write_separate_scorefile_ = true;
	} else write_separate_scorefile_ = false;
}

void SilentFileJobOutputter::read_done_jobs() {
	core::io::silent::SilentFileData sfd;
	if ( utility::file::file_exists( silent_file_ ) ) {
		silent_file_tags_ = sfd.read_tags_fast( silent_file_ );
		BOOST_FOREACH ( std::string & tag, silent_file_tags_ ) {
			/// eliminate the FAILURE_ prefix so that jobs know to start from
			/// the 'next' nstruct on restart. This is important to avoid duplicate
			/// entries
			if ( tag.substr( 0, 8 ) == "FAILURE_" ) {
				tag = tag.substr( 8, 1000 );
			}
		} //foreach
	} //fi
}

void SilentFileJobOutputter::final_pose(
	JobOP job,
	core::pose::Pose const & pose,
	std::string const & tag
) {
	call_output_observers( pose, job );

	core::io::silent::SilentStructOP ss =
		dump_pose( silent_file_, job, pose,  bWriteNoStructures_, -1 /* copy_count */, tag);

	// only write a scorefile if specified by the user
	if ( write_separate_scorefile_ ) {

		// Adding Luki & Brian's JD2 Patch (@ralford 8/23/14)
		core::pose::Pose pose_copy( pose );
		ss->energies_into_pose( pose_copy );
		scorefile(job, pose_copy, "", (tag.empty() ? "" : std::string("_") + tag));

		//scorefile(job, pose, "", (tag.empty() ? "" : std::string("_") + tag));
	}
}

/// @brief this function is intended for saving mid-protocol poses; for example
/// the final centroid structure in a combined centroid/fullatom protocol.
/// --->these go to file silent_filename+tag
void SilentFileJobOutputter::other_pose(
	JobOP job,
	core::pose::Pose const & pose,
	std::string const & tag,
	int copy_count, /*default -1 */
	bool score_only /*default false*/
) {
	call_output_observers( pose, job );
	utility::file::FileName filename( silent_file_ );
	filename.base( silent_file_.base() +"_"+ tag );

	utility::file::FileName tagged_score_filename( scorefile_name() );
	tagged_score_filename.base( scorefile_name().base() +"_"+ tag );

	core::io::silent::SilentStructOP ss( NULL );
	if ( bWriteIntermediateFiles_ ) {
		if ( write_separate_scorefile_ && score_only ) {
			dump_pose( tagged_score_filename, job, pose, !bWriteIntermediateStructures_ || score_only , copy_count );
		} else {
			ss=dump_pose( filename, job, pose, !bWriteIntermediateStructures_ || score_only , copy_count );
		}
	}

	if ( write_separate_scorefile_ && ss ) {
		//write here to avoid 2x evaluated for same structure
		//core::io::silent::SilentFileData sfd;
		//sfd.write_silent_struct( *ss, tagged_score_filename, true );
		scorefile(job, pose, tag, "", basic::options::option[ basic::options::OptionKeys::run::other_pose_scorefile ].value());
	}


}

core::io::silent::SilentStructOP SilentFileJobOutputter::dump_pose(
	utility::file::FileName const & filename,
	JobCOP job,
	core::pose::Pose const & pose_in,
	bool bWriteScoreOnly,
	int copy_count,
	std::string const & suffix
) {
	PROF_START( basic::JD2_SILENT_OUTPUTTER );
	core::io::silent::SilentFileData sfd;

	using core::io::silent::SilentStructFactory;
	core::io::silent::SilentStructOP ss;
	if ( bWriteScoreOnly ) {
		ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("score");
	} else if ( forced_silent_struct_type_ != "any" ) {
		ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( forced_silent_struct_type_ );
	} else {
		ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( pose_in );
	}

	std::ostringstream tag;
	tag << output_name( job );
	if ( copy_count>=0 ) {
		core::Size const replica( current_replica() );
		if ( replica ) {
			tag << "_" << std::setfill('0') << std::setw(3) << replica;
		}
		tag << '_' << std::setfill('0') << std::setw(8) << copy_count;
	}

	if ( !suffix.empty() ) {
		tag << "_" << suffix;
	}

	ss->fill_struct( pose_in, tag.str() );
	add_job_data_to_ss( ss, job );

	core::pose::Pose pose( pose_in );
	evaluate( pose, tag.str(), *ss );


	add_silent_struct( ss, filename );
	tr.Debug << "adding struct " << ss->decoy_tag() << std::endl;
	tr.Debug << "have " << saved_structs_.size() << ", buffering " << n_to_buffer_ << std::endl;
	core::Real rand;
	if ( random_flush_frequency_ < 1.0 ) rand = numeric::random::rg().uniform();
	else rand=0.0;
	if ( saved_structs_.size() >= n_to_buffer_ && random_flush_frequency_ >= rand ) {
		write_all_structs();
	}

	tr.flush();
	PROF_STOP( basic::JD2_SILENT_OUTPUTTER );
	return ss;
}

void SilentFileJobOutputter::add_silent_struct(
	core::io::silent::SilentStructOP ss,
	utility::file::FileName const & fn
) {
	saved_structs_.push_back( std::make_pair( ss, fn ) );
}

/////////////////////////////////state of output functions/////////////////////////////////
bool SilentFileJobOutputter::job_has_completed( JobCOP job ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// did we complete the job later ?
	if ( job->completed() ) {
		return true;
	}

	// was the job completed beforehand ( in the silent file before the app even
	// started ) ?
	if ( option[ run::multiple_processes_writing_to_one_directory ].value() ) {
		read_done_jobs(); // refresh silent_file_tags_ for parallel processes
	}
	CompareTags predicate( output_name(job) );

	bool const already_written(
		find_if(silent_file_tags_.begin(), silent_file_tags_.end(), predicate) != silent_file_tags_.end()
	);

	using std::pair;
	using std::string;
	using utility::vector1;
	using utility::file::FileName;
	using core::io::silent::SilentStructOP;

	vector1< string > tags;
	typedef vector1< pair< SilentStructOP, FileName > >::const_iterator iter;
	for ( iter it = saved_structs_.begin(), end = saved_structs_.end(); it != end; ++it ) {
		tags.push_back( it->first->decoy_tag() );
	}

	bool const already_buffered(
		find_if( tags.begin(), tags.end(), predicate ) != tags.end()
	);

	return ( already_written || already_buffered );
}

/// @details
/// SilentFile tags should preserve the FULL NAME such that we don't end up with
/// duplicate tags. This will cause problems on BOINC if changed.
std::string SilentFileJobOutputter::output_name( JobCOP job ) {
	return affixed_numbered_name( job );
}


void
SilentFileJobOutputter::set_silent_file_name( utility::file::FileName name ){
	silent_file_ = name;
	read_done_jobs(); //safer to do this again
}

void
SilentFileJobOutputter::set_forced_silent_struct_type( std::string const& setting ) {
	forced_silent_struct_type_ = setting;
}

void
SilentFileJobOutputter::set_write_separate_scorefile( bool write_separate_scorefile ){
	write_separate_scorefile_ = write_separate_scorefile;
}

//CREATOR SECTION
std::string
SilentFileJobOutputterCreator::keyname() const
{
	return "SilentFileJobOutputter";
}

protocols::jd2::JobOutputterOP
SilentFileJobOutputterCreator::create_JobOutputter() const {
	return protocols::jd2::JobOutputterOP( new SilentFileJobOutputter );
}

} //jd2
} //protocols
