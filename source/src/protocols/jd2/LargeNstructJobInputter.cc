// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/LargeNstructJobInputter.cc
/// @brief  A JobInputter for cases where the number of jobs is large enough to fill up memory.
/// @author Vikram K. Mulligan, Baker Laboratory (vmullig@uw.edu)

// Unit headers
#include <protocols/jd2/LargeNstructJobInputter.hh>
#include <protocols/jd2/LargeNstructJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <core/pose/symmetry/util.hh>

// Project headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/util.hh>
//#include <protocols/simple_moves/ExtendedPoseMover.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.jd2.LargeNstructJobInputter" );

namespace protocols {
namespace jd2 {

LargeNstructJobInputter::LargeNstructJobInputter() {
	tr.Debug << "Instantiate LargeNstructJobInputter" << std::endl;
}

/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
void LargeNstructJobInputter::pose_from_job( core::pose::Pose& pose, protocols::jd2::JobOP job) {
	//using protocols::simple_moves::ExtendedPoseMover;
	using std::string;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	tr.Debug << "LargeNstructJobInputter::pose_from_job" << std::endl;
	if ( !job->inner_job()->get_pose() ) {
		//fpd if pose is symmetric we need to change SymmetricConformation to Conformation
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			core::pose::symmetry::make_asymmetric_pose( pose );
		}

		// Creates an extended, idealized pose from the first sequence in the first
		// file in -in:file:fasta. Preserves current behavior for the TopologyBroker
		if ( option[OptionKeys::in::file::fasta].user() ) {
			string fasta = option[in::file::fasta]()[1];
			string sequence = core::sequence::read_fasta_file_str(fasta)[1];
			core::pose::make_pose_from_sequence(pose, sequence, core::chemical::CENTROID );
			//ExtendedPoseMover m(sequence);
			//m.apply(pose);
		}
	} else {
		pose = *(job->inner_job()->get_pose());
		tr.Debug << "filling pose from saved copy " << job->inner_job()->input_tag() << std::endl;
	}
}

/// @details this function determines what jobs exist.  Note that if the value of the -jd2:max_nstruct_in_memory flag
/// is less than the value specified with the -nstruct flag, this ONLY sets up the first max_nstruct_in_memory jobs.
void LargeNstructJobInputter::fill_jobs( protocols::jd2::JobsContainer & jobs ){
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::jd2;

	if ( tr.Debug.visible() ) tr.Debug << "LargeNstructJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	jobs.set_job_inputter( utility::pointer::dynamic_pointer_cast< protocols::jd2::JobInputter >(get_self_ptr()) );

	core::Size const nstruct( get_nstruct () );
	core::Size initial_nstruct(nstruct);
	core::Size const max_nstruct( basic::options::option[ max_nstruct_in_memory ].value() );

	if ( max_nstruct>0 && nstruct>max_nstruct ) {
		initial_nstruct=max_nstruct;
	}

	populate_next_n_jobs(jobs, 1, initial_nstruct, nstruct );

	if ( max_nstruct>0 && nstruct>max_nstruct ) {
		jobs.set_total_jobs( nstruct ); //Specify that there are more jobs than are currently loaded into memory.
		if ( tr.visible() ) {
			tr << "The user specified -nstruct " << nstruct << ", but the maximum allowed jobs in memory (-jd2:max_nstruct_in_memory flag) is " << max_nstruct << ".  ";
			tr << "Setting up only the first " << initial_nstruct << " jobs.  The job list will be updated once these jobs run." << std::endl;
		}
	}

	if ( tr.visible() ) tr.flush(); //Flush the tracer.
	if ( tr.Debug.visible() ) tr.Debug.flush(); //Flush the debug tracer.
} //fill_jobs

/// @brief This function is only called by certain JobInputters to update the jobs list after it has already been created.
/// @details An example case would be the LargeNstructJobInputter, which uses this function to load additional jobs after
/// the first N have started to come back.
void LargeNstructJobInputter::update_jobs_list( JobsContainerOP jobs ) {
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::jd2;

	core::Size jobs_in_memory(0); //Counter for jobs currently in memory.
	core::Size jobs_were_in_memory(0); //Counter for jobs in memory prior to deletions.
	core::Size const max_nstruct( basic::options::option[ max_nstruct_in_memory ].value() );
	core::Size const nstruct( get_nstruct () );
	core::Size const highest_job_was_in_memory( jobs->highest_job_index() );

	for ( core::Size i=1, imax=jobs->highest_job_index(); i<=imax; ++i ) {
		if ( jobs->has_job(i) ) {
			if ( tr.Debug.visible() ) tr.Debug << "Job " << i << " is in memory." << std::endl;
			++jobs_in_memory;
			++jobs_were_in_memory;
			if ( jobs->can_be_deleted(i) ) {
				jobs->erase(i); //Delete all jobs that can be deleted (since they've completed and unloaded their data) or which are bad.
				--jobs_in_memory;
				if ( tr.Debug.visible() ) tr.Debug << "Deleting job " << i << "." << std::endl;
			} else {
				if ( tr.Debug.visible() ) tr << "Keeping job " << i << "." << std::endl;
			}
		}
	}

	if ( tr.Debug.visible() ) tr.Debug << jobs_were_in_memory << " jobs were in memory; " << jobs_in_memory << " remain after pruning list of jobs that have completed and unloaded data or which are bad." << std::endl;

	if ( highest_job_was_in_memory < nstruct ) {
		//Figure out how many jobs to add:
		core::Size jobs_to_add( max_nstruct - jobs_in_memory );
		if ( (highest_job_was_in_memory + jobs_to_add) > nstruct ) jobs_to_add = nstruct - highest_job_was_in_memory;

		//Add the additional jobs:
		populate_next_n_jobs( (*jobs), highest_job_was_in_memory+1, jobs_to_add, nstruct );
	}

	if ( tr.Debug.visible() ) tr.Debug.flush(); //Flush the debug output.

	return;
} //LargeNstructJobInputter::update_jobs_list( JobsContainerOP jobs )

/// @brief Return the type of input source that the LargeNstructJobInputter is currently using
/// @return Always <em>POSE</em>.
protocols::jd2::JobInputterInputSource::Enum LargeNstructJobInputter::input_source() const {
	return protocols::jd2::JobInputterInputSource::POSE;
}

/// @brief Private function to add N jobs to the list of jobs.
///
void LargeNstructJobInputter::populate_next_n_jobs(
	protocols::jd2::JobsContainer & jobs,
	core::Size const first_job_index,
	core::Size const number_of_jobs_to_add,
	core::Size const total_jobs
) {
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::jd2;

	//note that we are not really using the second and third fields in this implementation
	protocols::jd2::InnerJobOP ijob( new protocols::jd2::InnerJob( basic::options::option[ generic_job_name ].value() , total_jobs ) );

	for ( core::Size index = first_job_index, index_max = number_of_jobs_to_add + first_job_index - 1 ; index <= index_max; ++index ) {
		jobs.push_back( protocols::jd2::JobOP( new protocols::jd2::Job( ijob, index ) ) );
		//if(tr.Debug.visible()) tr.Debug << "create job index " << index << std::endl;
	}
}

//CREATOR SECTION
std::string
LargeNstructJobInputterCreator::keyname() const
{
	return "LargeNstructJobInputter";
}

protocols::jd2::JobInputterOP
LargeNstructJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new LargeNstructJobInputter );
}

}// jd2
}// protocols
