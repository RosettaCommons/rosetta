// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/internal_util.cc
/// @brief  Utilities for JD2, intended to be used within the JD2 system itself
/// @details For utilities which might be used outside JD2, see util.hh

// Unit Headers

#ifdef USEMPI
#include <mpi.h>
#endif

#include <protocols/jd2/util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/MPIMultiCommJobDistributor.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>
#include <basic/mpi/MessageListenerFactory.hh>

#include <protocols/evaluation/util.hh>

#include <core/pose/Pose.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>


#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>


#include <core/chemical/ChemicalManager.hh>
#include <protocols/moves/Mover.hh>

#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/mpi_util.hh>
#include <utility/assert.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

namespace protocols {
namespace jd2 {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.internal_util" );

void add_job_data_to_ss( core::io::silent::SilentStructOP ss, JobCOP job_op ) {
	using namespace core::pose;

	for ( auto iter = job_op->output_string_string_pairs_begin(),
			end = job_op->output_string_string_pairs_end();
			iter != end; ++iter
			) {
		ss->add_string_value(iter->first, iter->second );
	}

	for ( auto iter = job_op->output_string_real_pairs_begin(),
			end = job_op->output_string_real_pairs_end();
			iter != end; ++iter
			) {
		ss->add_energy( iter->first, iter->second, 1.0 );
	}
}

void register_options() {
	evaluation::register_options();
	using namespace basic::options::OptionKeys;
	OPT( in::file::silent );
	OPT( in::file::s );
	OPT( in::file::l );
	OPT( in::file::native );
	OPT( in::file::silent_read_through_errors );
	OPT( out::file::silent );
	OPT( out::file::scorefile );
	OPT( run::batches );
	OPT( run::archive );
}

JobOP get_current_job() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->job_inputter() ) {
		return jd->current_job();
	} else return nullptr;
}


#ifdef USEMPI
/// @brief returns communicator defined by the JobDistributor or MPI_COMM_WORLD
MPI_Comm const& current_mpi_comm() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd ) {
		MPIMultiCommJobDistributor* mpi_jd = dynamic_cast< MPIMultiCommJobDistributor* >( jd );
		if ( mpi_jd ) {
			return mpi_jd->current_mpi_comm();
		}
	}
	//return MPI_COMM_WORLD; //causes warning: returning reference to temporary
	//workaround to avoid warning ( MPI_COMM_WORLD is actually a macro )
	TR.Trace << "Requested jd2::current_mpi_comm() but apparently flag -run:n_replica was not set.";
	TR.Trace << "Returning MPI_COMM_WORLD" << std::endl;
	static MPI_Comm my_mpi_comm_world = MPI_COMM_NULL;
	MPI_Comm_dup( MPI_COMM_WORLD, &my_mpi_comm_world );
	return my_mpi_comm_world;
}
#endif

} // jd2
} // protocols
