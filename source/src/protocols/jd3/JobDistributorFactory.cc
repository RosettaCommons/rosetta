// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobDistributorFactory.cc
/// @brief  JobDistributorFactory class; reads the command line to create a JobDistributor
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/JobDistributorFactory.hh>

// Package headers
#include <protocols/jd3/JobDistributor.hh>
//#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd3/job_distributors/MPIWorkPoolJobDistributor.hh>
// not yet #include <protocols/jd3/job_distributors/MultiThreadedJobDistributor.hh>
#include <protocols/jd3/job_distributors/VanillaJobDistributor.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>


//#include <protocols/jd3/job_distributors/FileSystemJobDistributor.hh>

namespace protocols {
namespace jd3 {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.JobDistributorFactory" );

JobDistributorOP
JobDistributorFactory::create_job_distributor()
{
	// TEMP!
#ifdef USEMPI
#ifdef SERIALIZATION
	// In future: some command line option?
	// if ( basic::options::option[ basic::options::OptionKeys::jd3::mpi_work_partition_job_distributor ].value() ) {
	// 	return JobDistributorOP( new job_distributors::MPIWorkPartitionJobDistributor );
	// } else {
	return JobDistributorOP( new job_distributors::MPIWorkPoolJobDistributor );
	//}
#else
	TR.Error << "JD3 cannot be used with MPI but without serialization; you must add \"extras=mpi,serlization\" to your ./scons.py build command" << std::endl;
	return nullptr;
#endif // SERIALIZATION

#ifdef MULTI_THREADED
	// not yet return JobDistributorOP( new job_distributors::MultiThreadedJobDistributor );
#endif

#else
	return JobDistributorOP( new job_distributors::VanillaJobDistributor );
#endif


}


} // namespace jd3
} // namespace protocols

