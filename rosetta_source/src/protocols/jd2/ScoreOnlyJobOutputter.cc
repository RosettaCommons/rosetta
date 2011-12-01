// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/jd2/ScoreOnlyJobOutputter.cc
/// @author Sam DeLuca

#include <protocols/jd2/ScoreOnlyJobOutputter.hh>
#include <protocols/jd2/ScoreOnlyJobOutputterCreator.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <protocols/jd2/Job.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

ScoreOnlyJobOutputter::ScoreOnlyJobOutputter(): FileJobOutputter()
{}

//void ScoreOnlyJobOutputter::file(JobCOP,std::string const &)
//{}

void ScoreOnlyJobOutputter::final_pose(JobCOP job ,core::pose::Pose const& pose )
{
scorefile(job,pose);
}

void ScoreOnlyJobOutputter::other_pose(
  JobCOP job,
	core::pose::Pose const & pose,
	std::string const & tag,
	int copy_count, /*default -1 */
	bool score_only /*default false*/
) {
	if( basic::options::option[ basic::options::OptionKeys::run::other_pose_to_scorefile ].value() )
	{
		scorefile(job, pose, tag, basic::options::option[ basic::options::OptionKeys::run::other_pose_scorefile ].value());
	}
}

bool ScoreOnlyJobOutputter::job_has_completed(JobCOP)
{
	return false;
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
        return new ScoreOnlyJobOutputter;
}

}
}
