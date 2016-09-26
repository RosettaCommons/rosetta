// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/SilentFilePoseInputter.cc
/// @brief
/// @author Andy Watkins (amw579@stanford.edu)

///Unit headers
#include <protocols/jd3/pose_inputters/SilentFilePoseInputter.hh>
#include <protocols/jd3/pose_inputters/SilentFilePoseInputterCreator.hh>

// Package headers
#include <protocols/jd3/PoseInputSource.hh>
#include <protocols/jd3/pose_inputters/PoseInputterFactory.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/OptionCollection.fwd.hh>
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

static THREAD_LOCAL basic::Tracer tr( "protocols.jd3.SilentFilePoseInputter" );

namespace protocols {
namespace jd3 {
namespace pose_inputters {

SilentFilePoseInputter::SilentFilePoseInputter()
{
	tr.Debug << "Instantiate SilentFilePoseInputter" << std::endl;
}

SilentFilePoseInputter::~SilentFilePoseInputter() {}

bool SilentFilePoseInputter::job_available_on_command_line() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return option[ in::file::silent ].user();
}

/*
/// @brief this function returns the SilentStruct that belongs to the given job
core::io::silent::SilentStruct const&
SilentFilePoseInputter::struct_from_job( JobOP job ) {
	if ( !sfd_.has_tag( job->inner_job()->input_tag() ) ) {
		utility_exit_with_message(" job with input tag " + job->inner_job()->input_tag() +" can't find his input structure ");
	}
	return sfd_.get_structure( job->inner_job()->input_tag() );
}
*/

/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
core::pose::PoseOP SilentFilePoseInputter::pose_from_input_source(
	PoseInputSource const & input_source,
	utility::options::OptionCollection const & options
) const {
	assert( input_source.string_string_map().find( "filename" ) != input_source.string_string_map().end() );
	assert( input_source.string_string_map().find( "tag" ) != input_source.string_string_map().end() );

	core::import_pose::ImportPoseOptions import_opts( options );
	return core::import_pose::pose_from_file(
		input_source.string_string_map().find( "filename" )->second,
		import_opts,
		import_opts.read_fold_tree(),
		core::import_pose::PDB_file );
}

/// @details this function determines what jobs exist from -in::file::silent,
/// -in::file::tags and -in::file::tagfile
void SilentFilePoseInputter::fill_jobs( JobsContainer & jobs ){
	tr.Debug << "SilentFilePoseInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	using std::string;
	using utility::file::FileName;
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< FileName > const silent_files( option[ in::file::silent ]() );

	for ( vector1< FileName >::const_iterator current_fn_ = silent_files.begin();
			current_fn_ != silent_files.end(); ++current_fn_
			) {
		tr.Debug << "reading " << *current_fn_ << std::endl;
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
			sfd_.read_file( *current_fn_, tags );
		} else {
			sfd_.read_file( *current_fn_ );
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

/// @brief returns the schema for the PDB element used in a job-definition file
/// including all options that govern how a PDB is loaded.
static void SilentFilePoseInputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {

}

static void SilentFilePoseInputter::list_options_read( utility::options::OptionKeyList & read_options ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// possibly relevant options:
	StructFileRepOptions::list_options_read( read_options );
	read_options
		+ in::file::tags
		+ in::file::user_tags
		+ in::file::silent
		+ basic::options::OptionKeys::in::file::lazy_silent
		+ basic::options::OptionKeys::in::file::force_silent_bitflip_on_read
		+ basic::options::OptionKeys::in::file::atom_tree_diff
		+ basic::options::OptionKeys::in::file::fullatom
		+ basic::options::OptionKeys::in::file::centroid
		+ basic::options::OptionKeys::in::file::silent_energy_cut
		+ basic::options::OptionKeys::in::file::silent_list
		+ basic::options::OptionKeys::in::file::silent_struct_type
		+ basic::options::OptionKeys::in::file::silent_read_through_errors
		+ basic::options::OptionKeys::in::file::silent_score_prefix
		+ basic::options::OptionKeys::in::file::silent_select_random
		+ basic::options::OptionKeys::in::file::silent_select_range_start
		+ basic::options::OptionKeys::in::file::silent_select_range_mul
		+ basic::options::OptionKeys::in::file::silent_select_range_len
		+ basic::options::OptionKeys::in::file::skip_failed_simulations
		+ basic::options::OptionKeys::in::file::silent_scores_wanted
		+ basic::options::OptionKeys::in::file::template_silent;
}

//CREATOR SECTION
std::string
SilentFilePoseInputterCreator::keyname() const
{
	return "SilentFilePoseInputter";
}

protocols::jd3::PoseInputterOP
SilentFilePoseInputterCreator::create_inputter() const {
	return protocols::jd3::PoseInputterOP( new SilentFilePoseInputter );
}

} // pose_inputters
} // jd3
} // protocols
