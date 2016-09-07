// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/domain_assembly/DomainAssemblyJobInputter.cc
/// @brief  The DomainAssemblyJobInputter takes two poses, fuses them and creates jobs with the fused pose as starting pose
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

///Unit headers
#include <devel/domain_assembly/DomainAssemblyJobInputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

///Project headers
#include <devel/domain_assembly/domain_assembly.hh>
#include <devel/domain_assembly/domain_assembly_setup.hh>
#include <devel/domain_assembly/DomainAssemblyReader.hh>

#include <core/pose/Pose.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <utility/vector1.hh>

///Basic headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

///Option key headers
#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

///C++ headers
#include <string>


// option key includes

static THREAD_LOCAL basic::Tracer TR( "devel.domain_assembly.DomainAssemblyJobInputter" );

namespace devel {
namespace domain_assembly {

devel::domain_assembly::DomainAssemblyJobInputter::DomainAssemblyJobInputter(){
	TR << "Instantiate DomainAssemblyJobInputter" << std::endl;
}

devel::domain_assembly::DomainAssemblyJobInputter::~DomainAssemblyJobInputter() = default;

/// @details This function will first see if the pose already exists in the Job.  If not, it will read it into the pose reference, and hand a COP cloned from that pose to the Job. If the pose pre-exists it just copies the COP's pose into it.
void devel::domain_assembly::DomainAssemblyJobInputter::pose_from_job( core::pose::Pose & pose, protocols::jd2::JobOP job){
	TR << "DomainAssemblyJobInputter::pose_from_job" << std::endl;

	if ( !job->inner_job()->get_pose() ) {

		std::string option_filename( basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_setup_option_file ]() );
		utility::vector1< DomainInfo > domain_info;
		// this will utility_exit if things go wrong.
		TR << "reading DA setup file from " << option_filename << std::endl;
		parse_da_option_file( domain_info, option_filename );

		TR << "Processing domains" << std::endl;
		process_domains( domain_info );

		// go through domain_info and save domain and linker
		// definitions in the pose before they're lost.
		for ( utility::vector1< DomainInfo >::const_iterator di_it = domain_info.begin();
				di_it != domain_info.end();
				++di_it ) {
			TR << "domain info ("
				<< (int)(di_it - domain_info.begin()) << "): "
				<< di_it->get_domain_begin() << "-" << di_it->get_domain_end() << std::endl;
		}

		TR << "filling pose with connected domains" << std::endl;
		connect_domains( domain_info, pose );

		load_pose_into_job(pose, job); // save pose in inner job so we don't have to do this again.
	} else {
		TR << "filling pose from saved copy " << job->input_tag() << std::endl;
		pose = *(job->inner_job()->get_pose());
	}
}

/// @details this function determines what jobs exist
void devel::domain_assembly::DomainAssemblyJobInputter::fill_jobs( protocols::jd2::JobsContainer & jobs ){
	TR << "DomainAssemblyJobInputter::fill_jobs" << std::endl;

	jobs.clear();

	std::string option_filename( basic::options::option[  basic::options::OptionKeys::DomainAssembly::da_setup_option_file ]() );
	core::Size const nstruct( get_nstruct() );

	protocols::jd2::InnerJobOP ijob( new protocols::jd2::InnerJob( option_filename, nstruct ) );

	for ( core::Size index(1); index <= nstruct; ++index ) {
		jobs.push_back( protocols::jd2::JobOP( new protocols::jd2::Job( ijob, index ) ) );
	}

}//fill_jobs

/// @brief Return the type of input source that the DomainAssemblyJobInputter is currently
///  using.
/// @return PDB_FILE.
protocols::jd2::JobInputterInputSource::Enum DomainAssemblyJobInputter::input_source() const {
	// returning PDB_FILE here should be OK because the only thing
	// it does is remove stuff that makes the tag look like a file name.
	return protocols::jd2::JobInputterInputSource::PDB_FILE;
}


}//domain_assembly
}//devel
