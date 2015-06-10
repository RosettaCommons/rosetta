// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/LazySilentFileJobInputter.cc
/// @brief
/// @author

#include <protocols/jd2/LazySilentFileJobInputter.hh>
#include <protocols/jd2/LazySilentFileJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

///C++ headers
#include <string>


#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers


static thread_local basic::Tracer tr( "protocols.jd2.LazySilentFileJobInputter" );

namespace protocols {
namespace jd2 {

  protocols::jd2::LazySilentFileJobInputter::LazySilentFileJobInputter()
  {
    tr.Debug << "instantiating LazySilentFileJobInputter" << std::endl;
  }

  protocols::jd2::LazySilentFileJobInputter::~LazySilentFileJobInputter() {}
  

  void
  protocols::jd2::LazySilentFileJobInputter::fill_jobs( JobsContainer & jobs ) {
    using namespace utility;
    using namespace core;
    using namespace basic::options;

    tr.Debug << "LazySilentFileJobInputter::fill_jobs" << std::endl;

    utility::vector1< file::FileName > const silent_files( option[ OptionKeys::in::file::silent ]() );
    utility::vector1<std::string> tags;

    // This code appears to suggest that it will work for multiple files - alas
    // it doesnt. It only saves the tag names not the file names so ALL the tags
    // will be attempted to be read from the first file.
    // I added this assert here to so that the user doesnt accidentally try
    // feedign through multiple files.
    runtime_assert( silent_files.size() <= 1 );

    for ( vector1< file::FileName >::const_iterator current_fn_ = silent_files.begin();
	  current_fn_ != silent_files.end(); ++current_fn_
	  ) {
      utility::vector1< std::string > filetags;
      //core::Size startindex = tags.size();
      sfd_.read_tags_fast( *current_fn_, filetags );

      for( core::Size jj = 1; jj <= filetags.size(); jj++ ) {
	      tags.push_back( filetags[jj] );
      }
    }

    core::Size const nstruct(	get_nstruct()	);

    using namespace core::io::silent;
    utility::vector1< InnerJobOP > inner_jobs;

    tr.Debug << "reserve memory for InnerJob List " << tags.size() << std::endl;
    inner_jobs.reserve( tags.size() );
    tr.Debug << "fill list with " << tags.size() << " InnerJob Objects" << std::endl;
    for ( core::Size iter = 1; iter <= tags.size(); iter++ ) {
      protocols::jd2::InnerJobOP ijob( new InnerJob( tags[iter], nstruct ) );
      inner_jobs.push_back( ijob );
    }
    
    //tr.Debug << "reserve list for " << inner_jobs.size() * nstruct << " Job Objects" << std::endl;
    //jobs.reserve( nstruct * inner_jobs.size() );
    
    tr.Debug << "fill job list with... " << std::endl;
    for ( core::Size index = 1; index <= nstruct; ++index ) {
      for ( utility::vector1< InnerJobOP >::const_iterator ijob = inner_jobs.begin(); ijob != inner_jobs.end(); ijob ++ ) {
	jobs.push_back( protocols::jd2::JobOP( new Job( *ijob, index ) ) );
	tr.Trace << "pushing " << (*ijob)->input_tag() << " nstruct index " << index	<< std::endl;
      } // loop over nstruct
    } // loop over inputs
  } // fill_jobs
  

  
  
  core::io::silent::SilentStruct const&
  protocols::jd2::LazySilentFileJobInputter::struct_from_job( JobOP job ) {
    using namespace basic::options;
    utility::vector1<std::string> tag_to_read;
    tag_to_read.push_back( job->inner_job()->input_tag());

    if ( !sfd_.read_file( (option[ OptionKeys::in::file::silent ]()[1]).name(), tag_to_read ) ) {
      utility_exit_with_message(" job with input tag " + job->inner_job()->input_tag() +" can't find his input structure ");
    }
    return sfd_.get_structure( job->inner_job()->input_tag() );
  }
  

  void
  protocols::jd2::LazySilentFileJobInputter::pose_from_job( core::pose::Pose & pose,
							    JobOP job
							    ) {
    tr.Debug << "LazySilentFileJobInputter::pose_from_job" << std::endl;

    if ( !job->inner_job()->get_pose() ) {
      tr.Debug << "filling pose from SilentFile (tag = " << job->inner_job()->input_tag()
	       << ")" << std::endl;
      pose.clear();
      struct_from_job( job ).fill_pose( pose );
      tag_into_pose( pose, job->inner_job()->input_tag() );
    } else {
      pose.clear();
      pose = *(job->inner_job()->get_pose());
      tr.Debug << "filling pose from saved copy (tag = " << job->input_tag()
	       << ")" << std::endl;
    }
  }

  JobInputterInputSource::Enum LazySilentFileJobInputter::input_source() const {
    return JobInputterInputSource::SILENT_FILE;
  }

//CREATOR SECTION
std::string
LazySilentFileJobInputterCreator::keyname() const
{
  return "LazySilentFileJobInputter";
}

protocols::jd2::JobInputterOP
LazySilentFileJobInputterCreator::create_JobInputter() const {
  return protocols::jd2::JobInputterOP( new LazySilentFileJobInputter );
}

} //jd2
} ///protocols
