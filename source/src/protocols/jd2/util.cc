// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/util.cc
/// @brief  Utilities for JD2, which are somewhat safe to use outside of the JD2 system itself
/// @details For utilities internal to the jd2 system, see internal_util.hh

// Unit Headers

#include <protocols/jd2/util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/MPIMultiCommJobDistributor.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>
#include <protocols/jd2/internal_util.hh>
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
#include <utility/assert.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/options/keys/OptionKeyList.hh>

namespace protocols {
namespace jd2 {

static basic::Tracer TR( "protocols.jd2.util" );
static basic::Tracer tr_score("protocols.jd2.score", basic::t_info, true /*muted by default*/ );

bool jd2_used() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	return ( jd && jd->job_outputter() && jd->current_job() && jd->current_job()->inner_job() != JD2_BOGUS_JOB->inner_job());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string current_input_tag() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		return jd->current_job()->input_tag(); // Should defer to the inner job
	} else return "UnknownInput";

}

core::Size current_nstruct_index() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		return jd->current_job()->nstruct_index();
	} else return 0;

}

core::Size max_nstruct_index() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		return jd->current_job()->nstruct_max(); // Should defer to the inner job
	} else return basic::options::option[ basic::options::OptionKeys::out::nstruct ]();
}


std::string current_output_name() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->job_outputter() && jd->current_job() ) {
		return jd->job_outputter()->output_name( jd->current_job() );
	} else return "NoTag";
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

std::string current_batch() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd ) {
		return jd->get_current_batch();
	} else return "NoJD2";
}

core::Size current_batch_id() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd ) {
		return jd->current_batch_id();
	} else return 0;
}

core::Size current_replica() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd ) {
		auto* mpi_jd = dynamic_cast< MPIMultiCommJobDistributor* >( jd );
		if ( mpi_jd ) {
			return mpi_jd->sub_rank()+1;
		}
	}
	return 0;
}

core::pose::PoseCOP get_current_jobs_starting_pose() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	core::pose::PoseCOP pose( nullptr );
	if ( jd && jd->job_outputter() && jd->job_inputter() && jd->current_job() ) {
		JobOP job = jd->current_job();
		core::pose::PoseOP aPose( new core::pose::Pose );
		jd->job_inputter()->pose_from_job( *aPose, job);
		pose = aPose;
	}
	return pose;
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

/// @brief add output string
void add_string_to_current_job( std::string const & string_in ) {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		jd->current_job()->add_string( string_in );
	} else {
		// Do nothing ... for now
	}
}

/// @brief add output strings
void add_strings_to_current_job( std::list< std::string > const & strings) {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		jd->current_job()->add_strings( strings );
	} else {
		// Do nothing ... for now
	}
}

/// @brief add a string/string pair
void add_string_string_pair_to_current_job( std::string const & string1, std::string const & string2 ) {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		jd->current_job()->add_string_string_pair( string1, string2 );
	} else {
		// Do nothing ... for now
	}
}

/// @brief add a string/real pair
void add_string_real_pair_to_current_job( std::string const & string_in, core::Real const real_in ) {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		jd->current_job()->add_string_real_pair( string_in, real_in );
	} else {
		// Do nothing ... for now
	}
}

std::list< std::string > get_strings_from_current_job() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		return jd->current_job()->get_strings();
	} else {
		return {}; // return empty object
	}
}


std::map< std::string, std::string > get_string_string_pairs_from_current_job() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		return jd->current_job()->get_string_string_pairs();
	} else {
		return {}; // return empty object
	}
}

std::map< std::string, core::Real > get_string_real_pairs_from_current_job() {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		return jd->current_job()->get_string_real_pairs();
	} else {
		return {}; // return empty object
	}
}

void
add_current_job_data_to_ss( core::io::silent::SilentStructOP ss ) {
	JobDistributor* jd
		= JobDistributor::get_instance();
	if ( jd && jd->current_job() ) {
		protocols::jd2::add_job_data_to_ss( ss, jd->current_job() );
	} else {
		// Do nothing, currently.
	}
}

void
write_score_tracer( core::pose::Pose const& pose_in, std::string tracer_point ) {

	if ( !tr_score.visible() ) return;

	using core::io::silent::SilentStructFactory;
	core::io::silent::SilentStructOP ss;
	core::io::silent::SilentFileOptions opts;
	ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("score", opts);
	std::string tag( protocols::jd2::current_output_name() );
	ss->fill_struct( pose_in, tag );
	add_current_job_data_to_ss( ss );

	core::pose::Pose pose( pose_in );
	if ( protocols::jd2::jd2_used() ) {
		JobDistributor* jd = JobDistributor::get_instance();
		debug_assert( jd );
		debug_assert( jd->job_outputter() );
		jd->job_outputter()->evaluate( pose, tag, *ss );
	}

	ss->add_string_value("tracer_point", tracer_point );

	core::io::silent::SilentFileData sfd( opts );
	if ( !basic::options::option[ basic::options::OptionKeys::out::file::silent_print_all_score_headers ]() ) {
		ss->print_header( tr_score );
	}
	sfd.write_silent_struct( *ss, tr_score, true /*write scores only*/ );
	tr_score.flush();
}

//////////////////////////////////////////////////////////////////////////////

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
		for ( auto & better_list_file_name : better_list_file_names ) {
			list_file_names.push_back(better_list_file_name); // make a copy (-l)
		}
	}

	for ( auto & list_file_name : list_file_names ) {
		std::string filename( list_file_name.name() );
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

void set_native_in_mover( protocols::moves::Mover &mover ){
	set_native_in_mover( mover, basic::options::option );
}

void set_native_in_mover(
	protocols::moves::Mover & mover,
	utility::options::OptionCollection const & options
){
	using namespace basic::options::OptionKeys;

	if ( options[ in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose );
		core::chemical::ResidueTypeSetCOP rsd_set;
		if ( options[ in::file::fullatom ]() ) {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		} else {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		}
		std::string native_pdb_file  = options[ in::file::native ]();
		core::import_pose::pose_from_file( *native_pose, *rsd_set, native_pdb_file , core::import_pose::PDB_file);

		mover.set_native_pose( native_pose );
	}
}

void
options_for_set_native_in_mover(
	utility::options::OptionKeyList & opts
)
{
	using namespace basic::options::OptionKeys;
	opts
		+ in::file::native
		+ in::file::fullatom;
}


} // jd2
} // protocols
