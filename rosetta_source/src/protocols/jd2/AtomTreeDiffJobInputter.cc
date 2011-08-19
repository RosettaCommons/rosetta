// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/AtomTreeDiffJobInputter.cc
/// @brief
/// @author James Thompson

///Unit headers
#include <protocols/jd2/AtomTreeDiffJobInputter.hh>
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

// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>


static basic::Tracer tr("protocols.jd2.AtomTreeDiffJobInputter");

namespace protocols {
namespace jd2 {

protocols::jd2::AtomTreeDiffJobInputter::AtomTreeDiffJobInputter()
	//	silent_files_( option[ in::file::silent ]() ) {
{
	tr.Debug << "Instantiate AtomTreeDiffJobInputter" << std::endl;
}

protocols::jd2::AtomTreeDiffJobInputter::~AtomTreeDiffJobInputter() {}

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
		//core::import_pose::pose_from_pdb( pose, job->inner_job()->input_tag() );
		//		core::io::silent::SilentStructOP ss = atom_tree_diff_[ job->inner_job()->input_tag() ];
		tr.Debug << "filling pose from AtomTree (tag = " << job->input_tag() << ")" << std::endl;
		pose.clear();

		if ( atom_tree_diff_.has_ref_tag( job->input_tag() ) ) {
			atom_tree_diff_.read_pose(job->input_tag(), pose);
		} else if( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			// A "native" pose for the diff reference point.
			// This used to be required but will be rarely used now.
			tr << "Using -in:file:native for reference pose with AtomTreeDiff input format.";
			core::pose::Pose native_pose;
			core::import_pose::pose_from_pdb( native_pose, basic::options::option[ basic::options::OptionKeys::in::file::native ]().name() );
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

void protocols::jd2::AtomTreeDiffJobInputter::fill_jobs( Jobs & jobs ){
	tr.Debug << "AtomTreeDiffJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	using std::string;
	using utility::file::FileName;
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< FileName > const atom_tree_diff_files( option[ in::file::atom_tree_diff ]() );

	runtime_assert( atom_tree_diff_files.size() > 0 );// getting here should be impossible

	if( atom_tree_diff_files.size() != 1 ){
		utility_exit_with_message( "Code currently deals with only one input atom_tree_diff file" );
	}
	FileName file_name= atom_tree_diff_files.at(1);
	tr.Debug << "reading " << file_name << std::endl;

		atom_tree_diff_.read_file( file_name );

	core::Size const nstruct(	get_nstruct()	);


	if ( option[ in::file::tags ].user() ) {
		utility::vector1< string > const tags=  option[ in::file::tags ]();
		utility::vector1< string >::const_iterator tag= tags.begin();
		utility::vector1< string >::const_iterator end= tags.end();
		for(; tag != end; ++tag){
			assert(atom_tree_diff_.has_tag(*tag));
			InnerJobOP ijob( new InnerJob( *tag, nstruct ) );
			for(core::Size i=1; i<=nstruct; ++i){
				jobs.push_back( JobOP( new Job( ijob, i) ) );
			}
		}
	}
	else{
		std::map< std::string, core::Size > const & tags= atom_tree_diff_.get_tag_score_map();
		std::map< std::string, core::Size > ::const_iterator tag= tags.begin();
		std::map< std::string, core::Size > ::const_iterator const end= tags.end();

		for(; tag != end; ++tag){
			InnerJobOP ijob( new InnerJob( tag->first, nstruct ) );
			core::import_pose::atom_tree_diffs::ScoresPairList::const_iterator pairs= atom_tree_diff_.scores().begin();
			core::import_pose::atom_tree_diffs::ScoresPairList::const_iterator const end= atom_tree_diff_.scores().end();
			for(; pairs != end; ++pairs){
				if(pairs->first == tag->first)
					break;
			}
			assert(pairs != end);

			core::import_pose::atom_tree_diffs::Scores scores= pairs->second;
			core::import_pose::atom_tree_diffs::Scores::const_iterator const & scores_begin= scores.begin();
			core::import_pose::atom_tree_diffs::Scores::const_iterator const & scores_end= scores.end();

			for(core::Size j=1; j<=nstruct; ++j){
				JobOP job = new Job( ijob, j);
				core::import_pose::atom_tree_diffs::Scores::const_iterator scores_itr= scores_begin;
				for(; scores_itr != scores_end; ++scores_itr){
					job->add_string_real_pair(scores_itr->first, scores_itr->second);
				}
				jobs.push_back( job );
			}
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

} // jd2
} // protocols
