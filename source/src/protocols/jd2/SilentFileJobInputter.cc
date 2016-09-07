// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/SilentFileJobInputter.cc
/// @brief
/// @author James Thompson

///Unit headers
#include <protocols/jd2/SilentFileJobInputter.hh>
#include <protocols/jd2/SilentFileJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

///C++ headers
#include <string>

// External headers
#include <boost/algorithm/string/predicate.hpp>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/symmetry/util.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.jd2.SilentFileJobInputter" );

namespace protocols {
namespace jd2 {

protocols::jd2::SilentFileJobInputter::SilentFileJobInputter()
// silent_files_( option[ in::file::silent ]() ) {
{
	tr.Debug << "Instantiate SilentFileJobInputter" << std::endl;
}

protocols::jd2::SilentFileJobInputter::~SilentFileJobInputter() = default;


/// @brief this function returns the SilentStruct that belongs to the given job
core::io::silent::SilentStruct const&
protocols::jd2::SilentFileJobInputter::struct_from_job( JobOP job ) {
	if ( !sfd_.has_tag( job->inner_job()->input_tag() ) ) {
		utility_exit_with_message(" job with input tag " + job->inner_job()->input_tag() +" can't find his input structure ");
	}
	return sfd_.get_structure( job->inner_job()->input_tag() );
}


/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
void protocols::jd2::SilentFileJobInputter::pose_from_job(
	core::pose::Pose & pose,
	JobOP job
) {
	tr.Debug << "SilentFileJobInputter::pose_from_job" << std::endl;

	if ( !job->inner_job()->get_pose() ) {
		//core::import_pose::pose_from_file( pose, job->inner_job()->input_tag() , core::import_pose::PDB_file);
		//  core::io::silent::SilentStructOP ss = sfd_[ job->inner_job()->input_tag() ];
		tr.Debug << "filling pose from SilentFile (tag = " << job->inner_job()->input_tag()
			<< ")" << std::endl;
		pose.clear();

		// kinda hacky fix for symmetry ... this should probably be added to Pose::clear()
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			core::pose::symmetry::make_asymmetric_pose( pose );
		}

		struct_from_job( job ).fill_pose( pose );
		tag_into_pose( pose, job->inner_job()->input_tag() );
		//   if ( !sfd_.has_tag( job->inner_job()->input_tag() ) ) {
		//    utility_exit_with_message(" job with input tag " + job->inner_job()->input_tag() +" can't find his input structure ");
		//   }
		//   sfd_.get_structure( job->inner_job()->input_tag() ).fill_pose( pose );
		// ss->fill_pose( pose );//, ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		//  load_pose_into_job( pose, job ); //this is a HUGE memory leak... not really useful
		//   making structure from silent-structs doesn't take that long...
	} else {
		pose.clear();

		// kinda hacky fix for symmetry ... this should probably be added to Pose::clear()
		// this will likely be overwritten in next assignment.... I move it further down OL 11/21/09
		//  if ( core::pose::symmetry::is_symmetric( pose ) ) {
		//    core::conformation::symmetry::make_asymmetric_pose( pose );
		//  }

		pose = *(job->inner_job()->get_pose());

		/// should this really be here? If there is a pose in the job-object it should have the right properties already, no?
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			core::pose::symmetry::make_asymmetric_pose( pose );
		}

		tr.Debug << "filling pose from saved copy (tag = " << job->input_tag()
			<< ")" << std::endl;
	}
}

/// @details this function determines what jobs exist from -in::file::silent,
/// -in::file::tags and -in::file::tagfile
void protocols::jd2::SilentFileJobInputter::fill_jobs( JobsContainer & jobs ){
	tr.Debug << "SilentFileJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	using std::string;
	using utility::file::FileName;
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< FileName > const silent_files( option[ in::file::silent ]() );

	for (const auto & silent_file : silent_files) {
		tr.Debug << "reading " << silent_file << std::endl;
		if ( option[ in::file::tags ].user() || option[ in::file::tagfile ].user() ) {
			utility::vector1< string > tags;
			if ( option[ in::file::tags ].user() ) {
				tags.append( option[ in::file::tags ]() );
			}
			if ( option[ in::file::tagfile ].user() ) {
				utility::io::izstream tag_file( option[ in::file::tagfile ]() );

				// Add the whitespace separated tag entries to the file
				std::copy( std::istream_iterator< std::string >( tag_file ), std::istream_iterator< std::string >(),
					std::back_inserter( tags ) );
			}
			sfd_.read_file( silent_file, tags );
		} else {
			sfd_.read_file( silent_file );
		}
	}

	core::Size const nstruct( get_nstruct() );

	using namespace core::io::silent;
	utility::vector1< InnerJobOP > inner_jobs;

	//save list of all inner_jobs first... this allows better sampling of jobs in case of unfinished runs:
	// input1_0001
	// input2_0001
	// ...
	// inputn_0001
	// input1_0002
	// input2_0002
	// ....
	tr.Debug << "reserve memory for InnerJob List " << sfd_.size() << std::endl;
	inner_jobs.reserve( sfd_.size() );
	tr.Debug << "fill list with " << sfd_.size() << " InnerJob Objects" << std::endl;
	for ( SilentFileData::iterator iter = sfd_.begin(), end = sfd_.end();
			iter != end; ++iter
			) {
		const std::string tag = iter->decoy_tag();

		// Optionally ignore failed simulations. Supporting protocols are not consistent
		// in their support of this option. Abrelax, for example, writes models from
		// failed simulations in centroid residue type set, despite the fact that
		// fullatom was requested. This can lead to issues during clustering, rescoring,
		// etc.
		bool failed_simulation = boost::starts_with(tag, "W_");
		if ( failed_simulation && option[OptionKeys::in::file::skip_failed_simulations]() ) {
			continue;
		}

		InnerJobOP ijob( new InnerJob( tag, nstruct ) );
		inner_jobs.push_back( ijob );
	}

	//tr.Debug << "reserve list for " << inner_jobs.size() * nstruct << " Job Objects" << std::endl;
	//jobs.reserve( nstruct * inner_jobs.size() );

	tr.Debug << "fill job list with... " << std::endl;
	for ( core::Size index = 1; index <= nstruct; ++index ) {
		for ( utility::vector1< InnerJobOP >::const_iterator ijob = inner_jobs.begin(), end = inner_jobs.end(); ijob != end; ++ijob ) {
			jobs.push_back( JobOP( new Job( *ijob, index ) ) );
			tr.Trace << "pushing " << (*ijob)->input_tag() << " nstruct index " << index << std::endl;
		} // loop over nstruct
	} // loop over inputs
} // fill_jobs

/// @brief Return the type of input source that the SilentFileJobInputter is currently
///  using.
/// @return Always <em>SILENT_FILE</em>.
JobInputterInputSource::Enum SilentFileJobInputter::input_source() const {
	return JobInputterInputSource::SILENT_FILE;
}

//CREATOR SECTION
std::string
SilentFileJobInputterCreator::keyname() const
{
	return "SilentFileJobInputter";
}

protocols::jd2::JobInputterOP
SilentFileJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new SilentFileJobInputter );
}

} // jd2
} // protocols
