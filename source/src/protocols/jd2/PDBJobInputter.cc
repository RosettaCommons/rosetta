// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/PDBJobInputter.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class PDBJobInputter
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/PDBJobInputter.hh>
#include <protocols/jd2/PDBJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

///Project headers

#include <core/pose/Pose.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/vector1.hh>

#include <core/util/metalloproteins_util.hh>

///C++ headers
#include <string>

#include <core/import_pose/import_pose.hh>


// option key includes

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.PDBJobInputter" );

namespace protocols {
namespace jd2 {

protocols::jd2::PDBJobInputter::PDBJobInputter(){
	TR << "Instantiate PDBJobInputter" << std::endl;
}

protocols::jd2::PDBJobInputter::~PDBJobInputter(){}

/// @details This function will first see if the pose already exists in the Job.  If not,
/// it will read it into the pose reference, and hand a COP cloned from that pose to the
/// Job. If the pose pre-exists it just copies the COP's pose into it.
void protocols::jd2::PDBJobInputter::pose_from_job( core::pose::Pose & pose, JobOP job){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "PDBJobInputter::pose_from_job" << std::endl;

	if ( !job->inner_job()->get_pose() ) {
		TR << "filling pose from PDB " << job->input_tag() << std::endl;
		core::import_pose::pose_from_pdb( pose, job->input_tag() );
		load_pose_into_job(pose, job);
	} else {
		TR << "filling pose from saved copy " << job->input_tag() << std::endl;
		pose = *(job->inner_job()->get_pose());

		if ( option[in::auto_setup_metals].user()
				&& (option[in::metals_distance_constraint_multiplier]() > 1.0e-10 ||
				option[in::metals_angle_constraint_multiplier]() > 1.0e-10) ) {

			//If the user has specified the auto_setup_metals option, we need to add back the metal constraints.
			TR << "setting up metal constraints for saved pose copy " << job->input_tag() << std::endl;
			core::util::auto_setup_all_metal_constraints( pose,
				option[in::metals_distance_constraint_multiplier](),
				option[in::metals_angle_constraint_multiplier]() );

		}
	}
}

/// @details this function determines what jobs exist from -s/-l
void protocols::jd2::PDBJobInputter::fill_jobs( JobsContainer & jobs ){
	TR << "PDBJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	utility::vector1< std::string > const inputs( basic::options::start_files() );
	core::Size const nstruct( get_nstruct() );

	for ( core::Size i(1); i <= inputs.size(); ++i ) {

		//protocols::jobdist::BasicJob = protocols::jd2::InnerJob
		//note that we are not really using the second and third fields in this implementation
		InnerJobOP ijob( new InnerJob( inputs[i], nstruct ) );

		for ( core::Size index(1); index <= nstruct; ++index ) {
			jobs.push_back( JobOP( new Job( ijob, index ) ) );
			TR.Debug << "pushing " << inputs[i] << " nstruct index " << index << std::endl;
		}//loop over nstruct
		TR << "pushed " << inputs[i] << " nstruct "
			<< ( ( int(nstruct)==1 ) ? "index " : "indices 1 - ")
			<< nstruct << std::endl;
	} //loop over inputs
} //fill_jobs

/// @brief Return the type of input source that the PDBJobInputter is currently
///  using.
/// @return Always <em>PDB_FILE</em>.
JobInputterInputSource::Enum PDBJobInputter::input_source() const {
	return JobInputterInputSource::PDB_FILE;
}

//CREATOR SECTION
std::string
PDBJobInputterCreator::keyname() const
{
	return "PDBJobInputter";
}

protocols::jd2::JobInputterOP
PDBJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new PDBJobInputter );
}

}//jd2
}//protocols
