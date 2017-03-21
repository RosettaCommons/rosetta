// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/AtomTreeDiffJobInputter.cc
/// @brief
/// @author Gordon Lemmon (gordon.h.lemmon@vanderbilt.edu); Rocco Moretti (rmoretti@u.washington.edu)

///Unit headers
#include <protocols/jd2/AtomTreeDiffJobInputter.hh>
#include <protocols/jd2/AtomTreeDiffJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

///Project headers

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

///C++ headers
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.jd2.AtomTreeDiffJobInputter" );

namespace protocols {
namespace jd2 {

protocols::jd2::AtomTreeDiffJobInputter::AtomTreeDiffJobInputter()
// silent_files_( option[ in::file::silent ]() ) {
{
	tr.Debug << "Instantiate AtomTreeDiffJobInputter" << std::endl;
}

protocols::jd2::AtomTreeDiffJobInputter::~AtomTreeDiffJobInputter() = default;

/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
void protocols::jd2::AtomTreeDiffJobInputter::pose_from_job(
	core::pose::Pose & pose,
	JobOP job
) {
	tr.Debug << "AtomTreeDiffJobInputter::pose_from_job" << std::endl;

	if ( !job->inner_job()->get_pose() ) {
		//core::import_pose::pose_from_file( pose, job->inner_job()->input_tag() , core::import_pose::PDB_file);
		//  core::io::silent::SilentStructOP ss = atom_tree_diff_[ job->inner_job()->input_tag() ];
		tr.Debug << "filling pose from AtomTree (tag = " << job->input_tag() << ")" << std::endl;
		pose.clear();

		if ( atom_tree_diff_.has_ref_tag( job->input_tag() ) ) {
			atom_tree_diff_.read_pose(job->input_tag(), pose);
		} else if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			// A "native" pose for the diff reference point.
			// This used to be required but will be rarely used now.
			tr << "Using -in:file:native for reference pose with AtomTreeDiff input format.";
			core::pose::Pose native_pose;
			core::import_pose::pose_from_file( native_pose, basic::options::option[ basic::options::OptionKeys::in::file::native ]().name() , core::import_pose::PDB_file);
			atom_tree_diff_.read_pose(job->input_tag(), pose, native_pose);
		} else {
			utility_exit_with_message("Can't find the reference structure for job with input tag " + job->input_tag());
		}
		load_pose_into_job( pose, job );
	} else {
		pose.clear();

		pose = *(job->inner_job()->get_pose());
		tr.Debug << "filling pose from saved copy (tag = " << job->input_tag() << ")" << std::endl;
	}
}

/// @details this function determines what jobs exist -in:file:atom_tree_diff

void protocols::jd2::AtomTreeDiffJobInputter::fill_jobs( JobsContainer & jobs ){
	tr.Debug << "AtomTreeDiffJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	using std::string;
	using utility::file::FileName;
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< FileName > const atom_tree_diff_files( option[ in::file::atom_tree_diff ]() );

	runtime_assert( atom_tree_diff_files.size() > 0 );// getting here should be impossible

	if ( atom_tree_diff_files.size() != 1 ) {
		utility_exit_with_message( "Code currently deals with only one input atom_tree_diff file" );
	}
	FileName file_name= atom_tree_diff_files.at(1);
	tr.Debug << "reading " << file_name << std::endl;

	atom_tree_diff_.read_file( file_name );

	core::Size const nstruct( get_nstruct() );

	core::import_pose::atom_tree_diffs::TagScoreMap tags;

	if ( option[ in::file::tags ].user() ) {
		core::import_pose::atom_tree_diffs::TagScoreMap const & atd_tsm= atom_tree_diff_.get_tag_score_map();
		for ( string const & tag : option[ in::file::tags ]() ) {
			if ( ! atom_tree_diff_.has_tag(tag) ) {
				utility_exit_with_message("Input AtomTreeDiff file does not have tag "+tag);
			}
			auto tag_score_iter(atd_tsm.find(tag));
			debug_assert( tag_score_iter != atd_tsm.end() ); // Shouldn't happen - we've checked it was present
			tags[tag] = tag_score_iter->second;
		}
	} else {
		tags = atom_tree_diff_.get_tag_score_map();
	}

	if ( ! tags.size() ) {
		utility_exit_with_message("No valid input structures found for AtomTreeDiff file.");
	}

	core::import_pose::atom_tree_diffs::ScoresPairList const & all_scores( atom_tree_diff_.scores() );
	for ( core::import_pose::atom_tree_diffs::TagScorePair const & tag_score : tags ) {
		InnerJobOP ijob( new InnerJob( tag_score.first, nstruct ) );

		//second entry in tag_score is index to entry in ScoresPairList
		debug_assert( tag_score.first == all_scores[tag_score.second].first );
		core::import_pose::atom_tree_diffs::Scores scores= all_scores[tag_score.second].second;

		for ( core::Size j=1; j<=nstruct; ++j ) {
			JobOP job( new Job( ijob, j) );
			if ( basic::options::option[ basic::options::OptionKeys::in::file::keep_input_scores ] ) {
				for ( core::import_pose::atom_tree_diffs::ScorePair const & score : scores ) {
					job->add_string_real_pair(score.first, score.second);
				}
			}
			jobs.push_back( job );
		}
	}

} // fill_jobs

/// @brief Return the type of input source that the AtomTreeDiffJobInputter is currently
///  using.
/// @return Always <em>ATOM_TREE_FILE</em>.
JobInputterInputSource::Enum AtomTreeDiffJobInputter::input_source() const {
	return JobInputterInputSource::ATOM_TREE_FILE;
}

/// @brief Utility function to allow mutatable (!) access to the reference
/// poses of the undelying AtomTreeDiff object. Use with caution.
///
/// @details This exists to allow setup on stored reference poses for properties that
/// don't get saved/restored in PDB format, like covalent constraints for enzyme design.
utility::vector1< core::pose::PoseOP > const &
AtomTreeDiffJobInputter::all_ref_poses() const {
	return atom_tree_diff_.all_ref_poses();
}

//CREATOR SECTION
std::string
AtomTreeDiffJobInputterCreator::keyname() const
{
	return "AtomTreeDiffJobInputter";
}

protocols::jd2::JobInputterOP
AtomTreeDiffJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new AtomTreeDiffJobInputter );
}

} // jd2
} // protocols
