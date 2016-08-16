// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobDistributor.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class
/// @author Oliver Lange

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
#include <core/io/silent/SilentStructFactory.hh>


#include <basic/Tracer.hh>

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


static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.JobDistributor" );

utility::vector1< utility::file::FileName >
input_pdb_files_from_command_line()
{
	using basic::options::option;
	using utility::vector1;
	using utility::file::FileName;
	using namespace basic::options::OptionKeys;

	// concatenate -s and -l flags together to get total list of PDB files
	vector1< FileName > pdb_file_names;
	if ( option[ in::file::s ].active() ) {
		pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)
	}

	vector1< FileName > list_file_names;
	if ( option[ in::file::l ].active() ) {
		list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)
	}
	if ( option[ in::file::list ].active() ) {
		vector1< FileName > better_list_file_names;
		better_list_file_names = option[in::file::list ]().vector(); // make a copy (-list)
		for ( vector1< FileName >::iterator i = better_list_file_names.begin(), i_end = better_list_file_names.end(); i != i_end; ++i ) {
			list_file_names.push_back(*i); // make a copy (-l)
		}
	}

	for ( vector1< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i ) {
		std::string filename( i->name() );
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while ( getline(data, line) ) {
			pdb_file_names.push_back( FileName(line) );
		}
		data.close();
	}
	return pdb_file_names;
}

//multithreaded case requires specia
///end parser interface, start Job Distributor interface/////////////
void output_intermediate_pose(
	core::pose::Pose const & pose,
	std::string const & stage_tag,
	int copy_count,
	bool score_only
) {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->job_outputter() && jd->current_job() ) {
		jd->job_outputter()->other_pose( jd->current_job(), pose, stage_tag, copy_count, score_only );
	} else {
		TR.Warning << "can't output intermediate pose if not running with  jobdistributor ( jd2 / 2008 )" << std::endl;
	}
}

std::string current_output_name() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->job_outputter() && jd->current_job() ) {
		return jd->job_outputter()->output_name( jd->current_job() );
	} else return "NoTag";
}

std::string current_batch() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd ) {
		return jd->get_current_batch();
	} else return "NoJD2";
}

bool jd2_used() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	return ( jd && jd->job_outputter() && jd->current_job()->inner_job() != JD2_BOGUS_JOB->inner_job());
}

std::string current_output_filename() {
	jd2::JobDistributor* jd  = jd2::JobDistributor::get_instance();
	if ( jd && jd->job_outputter() ) {
		JobOP job = jd->current_job();
		if ( job ) {
			return jd->job_outputter()->filename( job );
		}
	}
	return "JD2_OUTPUT_FILE_UNKNOWN"; //else
}

void
write_score_tracer( core::pose::Pose const& pose_in, std::string tracer_point ) {
	basic::Tracer tr_score("protocols.jd2.score", basic::t_info, true /*muted by default*/ );

	if ( !tr_score.visible() ) return;

	JobDistributor* jd = JobDistributor::get_instance();

	if ( !jd || !jd->job_outputter() ) {
		tr_score.Warning << "can't output intermediate pose if not running with  jobdistributor ( jd2 / 2008 )" << std::endl;
		return;
	}

	using core::io::silent::SilentStructFactory;
	core::io::silent::SilentStructOP ss;
	ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("score");
	JobCOP job( get_current_job() );
	std::string tag( jd->job_outputter()->output_name( job ) );
	ss->fill_struct( pose_in, tag );
	add_job_data_to_ss( ss, job );

	core::pose::Pose pose( pose_in );
	jd->job_outputter()->evaluate( pose, tag, *ss );

	ss->add_string_value("tracer_point", tracer_point );

	core::io::silent::SilentFileData sfd;
	if ( !basic::options::option[ basic::options::OptionKeys::out::file::silent_print_all_score_headers ]() ) {
		ss->print_header( tr_score );
	}
	sfd.write_silent_struct( *ss, tr_score, true /*write scores only*/ );
	tr_score.flush();
}

JobOP get_current_job() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->job_inputter() ) {
		return jd->current_job();
	} else return NULL;
}

core::pose::PoseCOP get_current_jobs_starting_pose() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	core::pose::PoseCOP pose( NULL );
	if ( jd && jd->job_outputter() && jd->job_inputter() && jd->current_job() ) {
		JobOP job = jd->current_job();
		core::pose::PoseOP aPose( new core::pose::Pose );
		jd->job_inputter()->pose_from_job( *aPose, job);
		pose = aPose;
	}
	return pose;
}

void add_job_data_to_ss( core::io::silent::SilentStructOP ss, JobCOP job_op ) {
	using namespace core::pose;

	typedef Job::StringStringPairs::const_iterator str_iter;
	for ( str_iter iter = job_op->output_string_string_pairs_begin(),
			end = job_op->output_string_string_pairs_end();
			iter != end; ++iter
			) {
		ss->add_string_value(iter->first, iter->second );
	}

	typedef Job::StringRealPairs::const_iterator real_iter;
	for ( real_iter iter = job_op->output_string_real_pairs_begin(),
			end = job_op->output_string_real_pairs_end();
			iter != end; ++iter
			) {
		ss->add_energy( iter->first, iter->second, 1.0 );
	}
}


void set_native_in_mover( protocols::moves::Mover &mover ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose );
		core::chemical::ResidueTypeSetCOP rsd_set;
		if ( option[ in::file::fullatom ]() ) {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		} else {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		}
		std::string native_pdb_file  = option[ in::file::native ]();
		core::import_pose::pose_from_file( *native_pose, *rsd_set, native_pdb_file , core::import_pose::PDB_file);
		mover.set_native_pose( native_pose );
	}
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

core::Size current_replica() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd ) {
		MPIMultiCommJobDistributor* mpi_jd = dynamic_cast< MPIMultiCommJobDistributor* >( jd );
		if ( mpi_jd ) {
			return mpi_jd->sub_rank()+1;
		}
	}
	return 0;
}

} // jd2
} // protocols
