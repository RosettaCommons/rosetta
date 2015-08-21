// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/GenericJobInputter.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class GenericJobInputter
/// @author Oliver Lange

// Unit headers
#include <protocols/comparative_modeling/GenericJobInputter.hh>
#include <protocols/comparative_modeling/GenericJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <core/pose/symmetry/util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/sequence/util.hh>
#include <protocols/simple_moves/ExtendedPoseMover.hh>

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


static thread_local basic::Tracer tr( "protocols.comparative_modeling.GenericJobInputter" );

namespace protocols {
namespace comparative_modeling {

GenericJobInputter::GenericJobInputter() {
	tr.Debug << "Instantiate GenericJobInputter" << std::endl;
}

/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
/// why does this code live in comparative_modeling?
void GenericJobInputter::pose_from_job( core::pose::Pose& pose, protocols::jd2::JobOP job) {
	using protocols::simple_moves::ExtendedPoseMover;
	using std::string;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	tr.Debug << "GenericJobInputter::pose_from_job" << std::endl;
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
			ExtendedPoseMover m(sequence);
			m.apply(pose);
		}
	} else {
		pose = *(job->inner_job()->get_pose());
		tr.Debug << "filling pose from saved copy " << job->inner_job()->input_tag() << std::endl;
	}
}

/// @details this function determines what jobs exist from -s/-l
void GenericJobInputter::fill_jobs( protocols::jd2::JobsContainer & jobs ){
	tr.Debug << "GenericJobInputter::fill_jobs" << std::endl;

	jobs.clear(); //should already be empty anyway

	core::Size const nstruct( get_nstruct () );

	//note that we are not really using the second and third fields in this implementation
	using basic::options::OptionKeys::jd2::generic_job_name; //This option defaults to 'S' for original behavior
	protocols::jd2::InnerJobOP ijob( new protocols::jd2::InnerJob( basic::options::option[ generic_job_name ].value() , nstruct ) );

	for ( core::Size index = 1; index <= nstruct; ++index ) {
		jobs.push_back( protocols::jd2::JobOP( new protocols::jd2::Job( ijob, index ) ) );
		tr.Trace << "create job index " << index << std::endl;
	}
}

/// @brief Return the type of input source that the GenericJobInputter is currently using
/// @return Always <em>POSE</em>.
protocols::jd2::JobInputterInputSource::Enum GenericJobInputter::input_source() const {
	return protocols::jd2::JobInputterInputSource::POSE;
}

//CREATOR SECTION
std::string
GenericJobInputterCreator::keyname() const
{
	return "GenericJobInputter";
}

protocols::jd2::JobInputterOP
GenericJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new GenericJobInputter );
}

}// comparative_modeling
}// protocols
