// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/jd2/ScoreOnlyJobOutputter.cc
/// @author Sam DeLuca

#include <protocols/jd2/ScoreOnlyJobOutputter.hh>
#include <protocols/jd2/ScoreOnlyJobOutputterCreator.hh>

#include <protocols/jd2/SilentFileJobOutputter.hh> // For CompareTags
#include <core/io/silent/SilentFileData.hh> // used in read_done_jobs()
#include <core/io/silent/SilentFileOptions.hh> // used in read_done_jobs()
#include <utility/file/file_sys_util.hh> // for file_exists()

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <protocols/jd2/Job.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace jd2 {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.ScoreOnlyJobOutputter" );

ScoreOnlyJobOutputter::ScoreOnlyJobOutputter(): FileJobOutputter()
{
	read_done_jobs();
}

//void ScoreOnlyJobOutputter::file(JobCOP,std::string const &)
//{}

void ScoreOnlyJobOutputter::final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag )
{
	call_output_observers( pose, job );
	scorefile(job, pose, tag);
}

void ScoreOnlyJobOutputter::other_pose(
	JobOP job,
	core::pose::Pose const & pose,
	std::string const & tag,
	int /*copy_count*/, /*default -1 */
	bool /*score_only*/ /*default false*/
) {
	call_output_observers( pose, job );
	if ( basic::options::option[ basic::options::OptionKeys::run::other_pose_to_scorefile ].value() ) {
		scorefile(job, pose, tag, "", basic::options::option[ basic::options::OptionKeys::run::other_pose_scorefile ].value());
	}
}

bool ScoreOnlyJobOutputter::job_has_completed(JobCOP job) {

	// Inelegant cut-paste job from SilentFileJobOutputter::job_has_completed()

	// Is the job already marked as done?
	if ( job->completed() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "Skipping job " << output_name(job) << " because it has been marked as already completed." << std::endl;
		}
		return true;
	}

	// Was the job completed before the app even started?
	if ( basic::options::option[ basic::options::OptionKeys::run::multiple_processes_writing_to_one_directory ].value() ) {
		read_done_jobs(); // refresh score_file_tags_ for parallel processes
	}
	CompareTags predicate( output_name(job) );

	bool const already_written(
		find_if(score_file_tags_.begin(), score_file_tags_.end(), predicate) != score_file_tags_.end()
	);

	if ( TR.Debug.visible() && already_written ) {
		TR.Debug << "Skipping job " << output_name(job) << " because it has been already written to disk." << std::endl;
	}

	return already_written;
}

void ScoreOnlyJobOutputter::read_done_jobs() {
	// Inelegant cut-paste job from SilentFileJobOutputter::job_has_completed()

	if ( utility::file::file_exists( scorefile_name() ) ) {
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd( opts );
		score_file_tags_ = sfd.read_tags_fast( scorefile_name() );
		for ( std::string & tag : score_file_tags_ ) {
			/// eliminate the FAILURE_ prefix so that jobs know to start from
			/// the 'next' nstruct on restart. This is important to avoid duplicate
			/// entries
			if ( tag.substr( 0, 8 ) == "FAILURE_" ) {
				tag = tag.substr( 8 ); //start at 8, go for as many possible charachers as possible.
			}
		} //foreach
	} else {
		score_file_tags_.clear();
	}
}


std::string ScoreOnlyJobOutputter::output_name(JobCOP job)
{
	return affixed_numbered_name(job);
}

//CREATOR SECTION
std::string
ScoreOnlyJobOutputterCreator::keyname() const
{
	return "ScoreOnlyJobOutputter";
}

protocols::jd2::JobOutputterOP
ScoreOnlyJobOutputterCreator::create_JobOutputter() const {
	return protocols::jd2::JobOutputterOP( new ScoreOnlyJobOutputter );
}

}
}
