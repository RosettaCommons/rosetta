// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobDistributor.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class
/// @author Oliver Lange

// Unit Headers
#include <protocols/jd2/util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/evaluation/util.hh>

#include <core/pose/Pose.hh>
#include <core/io/silent/SilentStruct.hh>

#include <basic/Tracer.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>


#include <core/chemical/ChemicalManager.hh>
#include <core/io/pdb/pose_io.hh>
#include <protocols/moves/Mover.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


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


static basic::Tracer TR("protocols.jd2.JobDistributor");

//multithreaded case requires specia
///end parser interface, start Job Distributor interface/////////////
void output_intermediate_pose(
	core::pose::Pose const & pose,
	std::string const & stage_tag
) {
  protocols::jd2::JobDistributor* jd
  	= protocols::jd2::JobDistributor::get_instance();
  if ( jd && jd->job_outputter() && jd->current_job() ) {
    jd->job_outputter()->other_pose( jd->current_job(), pose, stage_tag );
  } else {
    TR.Warning << "can't output intermediate pose if not running with  jobdistributor ( jd2 / 2008 )" << std::endl;
  }
}

std::string current_output_name() {
	protocols::jd2::JobDistributor* jd
  	= protocols::jd2::JobDistributor::get_instance();
	if ( jd && jd->job_outputter() && jd->current_job() ) {
		return jd->job_outputter()->output_name( jd->current_job() );
	} else return "NoTag";
}

bool jd2_used() {
	protocols::jd2::JobDistributor* jd
  	= protocols::jd2::JobDistributor::get_instance();
	return ( jd && jd->job_outputter() && jd->current_job() != JD2_BOGUS_JOB );
}


JobOP get_current_job() {
	protocols::jd2::JobDistributor* jd
  	= protocols::jd2::JobDistributor::get_instance();
	if ( jd && jd->job_inputter() ) {
		return jd->current_job();
	}
	else return NULL;
}

core::pose::PoseCOP get_current_jobs_starting_pose() {
	protocols::jd2::JobDistributor* jd
  	= protocols::jd2::JobDistributor::get_instance();
	core::pose::PoseCOP pose( NULL );
	if ( jd && jd->job_outputter() && jd->job_inputter() && jd->current_job() ) {
		JobOP job = jd->current_job();
		core::pose::PoseOP aPose = new core::pose::Pose;
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
		core::pose::PoseOP native_pose = new core::pose::Pose;
		core::chemical::ResidueTypeSetCAP rsd_set;
 		if ( option[ in::file::fullatom ]() ) {
 			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
 		} else {
 			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
 		}
		std::string native_pdb_file  = option[ in::file::native ]();
		core::import_pose::pose_from_pdb( *native_pose, *rsd_set, native_pdb_file );
		mover.set_native_pose( native_pose );
	}

}


} // jd2
} // protocols
