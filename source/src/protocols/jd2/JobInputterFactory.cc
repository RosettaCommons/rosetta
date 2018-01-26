// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobInputterFactory.cc
/// @brief  JobInputterFactory
/// @author Steven Lewis smlewi@gmail.com, Brian Weitzner brian.weitzner@gmail.com

#include <protocols/jd2/JobInputterFactory.hh>
#include <protocols/jd2/JobInputterCreator.hh>
#include <protocols/jd2/JobInputter.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/make_rot_lib.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <basic/Tracer.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace jd2 {

static basic::Tracer TR( "protocols.jd2.JobInputterFactory" );

JobInputterFactory::JobInputterFactory() = default;

JobInputterFactory::~JobInputterFactory()= default;

/// @brief add a JobInputter prototype, using its default type name as the map key
void
JobInputterFactory::factory_register( JobInputterCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	std::string const job_inputter_type( creator->keyname() );
	if ( job_inputter_creator_map_.find( job_inputter_type ) != job_inputter_creator_map_.end() ) {
		utility_exit_with_message("JobInputterFactory::factory_register already has a mover creator with name \"" + job_inputter_type + "\".  Conflicting JobInputter names" );
	}
	job_inputter_creator_map_[ job_inputter_type ] = creator;
}

/// @details return new JobInputter by key lookup in map
JobInputterOP
JobInputterFactory::get_JobInputter_from_string( std::string const & job_inputter_type )
{
	//get pointer to Creator
	JobInputterMap::const_iterator iter( job_inputter_creator_map_.find( job_inputter_type ) );
	if ( iter != job_inputter_creator_map_.end() ) { //if Creator has an entry
		if ( ! iter->second ) { //if Creator inexplicably fails to exist, crash
			utility_exit_with_message( "Error: JobInputterCreatorOP for " + job_inputter_type + " is NULL, you should never have been able to get here!" );
		}

		//if creator exists, return a JobInputter from it (this is good)
		return iter->second->create_JobInputter();
	} else { //else, a non-existent JobInputter has been requested.  Print existing ones and exit.
		TR << "Available : ";
		for ( JobInputterMap::const_iterator mover_it = job_inputter_creator_map_.begin(); mover_it != job_inputter_creator_map_.end(); ++mover_it ) {
			TR << mover_it->first<<", ";
		}
		TR << std::endl;
		utility_exit_with_message( job_inputter_type + " is not known to the JobInputterFactory. Was it registered via a JobInputterRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return nullptr;
	}
}

/// @brief return new JobInputter from logic of option system plus compilation options.  All the logic for determining job input type lives here.
JobInputterOP
JobInputterFactory::get_new_JobInputter()
{
	using namespace basic::options;

	//initial copy of this code copied at XRW2 by SML+BDW from about SVN:46190 from JobDistributorFactory.cc

	if ( basic::options::option[ basic::options::OptionKeys::jd2::pose_input_stream ]() ) {
		return get_JobInputter_from_string( "PoseInputStreamJobInputter" );
	}

	//if ( basic::options::option[ basic::options::OptionKeys::jd2::resource_definition_files ].user() ) {
	// return get_JobInputter_from_string("JD2ResourceManagerJobInputter" );
	//}

	if ( basic::options::option[ basic::options::OptionKeys::in::file::screening_job_file].user() ) {
		return get_JobInputter_from_string("ScreeningJobInputter");
	}

	//PDB input block
	if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].user() || basic::options::option[ basic::options::OptionKeys::in::file::l ].user() || basic::options::option[ basic::options::OptionKeys::in::file::list ].user() || basic::options::option[ basic::options::OptionKeys::in::file::screening_list ].user() ) {
		if ( basic::options::option[ basic::options::OptionKeys::enzdes::parser_read_cloud_pdb ].user() ) {
			return get_JobInputter_from_string( "EnzdesJobInputter" );
		}
		if ( option[ OptionKeys::jd2::dd_parser ].user() && option[ OptionKeys::parser::patchdock ].user() ) {
			return get_JobInputter_from_string( "ParserJobInputter" );
		} else if ( option[ OptionKeys::jd2::grid_ensemble]() || option[ OptionKeys::jd2::seed_ensemble]() ||
				option[ OptionKeys::jd2::seed_ensemble_weights].user() || option[ basic::options::OptionKeys::jd2::seed_ensemble_weights_file].user() ) {
			return get_JobInputter_from_string( "EnsembleJobInputter");
		} else {
			return get_JobInputter_from_string( "PDBJobInputter" ); //SML override until we have other child classes
		}
		//silent file block
	} else if ( option[ OptionKeys::in::file::silent ].user() ) {
		if ( option[ OptionKeys::jd2::lazy_silent_file_reader ].user() ) {
			return get_JobInputter_from_string( "LazySilentFileJobInputter" );
		} else {
			return get_JobInputter_from_string( "SilentFileJobInputter" );
		}

	} else if ( option[OptionKeys::in::file::atom_tree_diff].user() ) {
		return get_JobInputter_from_string( "AtomTreeDiffJobInputter" );
	} else if ( option[ OptionKeys::in::file::template_pdb ].user() || option[ OptionKeys::in::file::template_silent ].user() ) {
		return get_JobInputter_from_string( "ThreadingJobInputter" );
	} else if ( option[OptionKeys::in::use_database].user() ) {
		return get_JobInputter_from_string( "DatabaseJobInputter" );
	} else if ( option[ OptionKeys::make_rot_lib::options_file ].user() ) {
		return get_JobInputter_from_string( "MakeRotLibJobInputter" );
	} else if (
			option[ OptionKeys::jd2::max_nstruct_in_memory ].value()!=0 &&
			option[ OptionKeys::jd2::max_nstruct_in_memory ].value() < option[ OptionKeys::out::nstruct ].value()
			) {
		return get_JobInputter_from_string( "LargeNstructJobInputter" ); //handles cases where the total job list won't fit into memory (e.g. big runs on a Blue Gene/Q machine)
	} else {
		return get_JobInputter_from_string( "GenericJobInputter" ); //handles -nstruct alone; works for abinitio with no structure input
	}
}

} //namespace jd2
} //namespace protocols
