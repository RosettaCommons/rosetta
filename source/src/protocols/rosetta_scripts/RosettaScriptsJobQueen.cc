// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/rosetta_scripts/RosettaScriptsJobQueen.cc
/// @brief  JD3 JobQueen for RosettaScripts class method implementation
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/rosetta_scripts/RosettaScriptsJobQueen.hh>

// protocol headers
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd3/JobDistributor.hh>
//#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/deallocation/ResourceDeallocationMessage.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>

#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.hh>
#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>


#include <basic/resource_manager/ResourceManager.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>



namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR("protocols.rosetta_scripts.RosettaScriptsJobQueen");

using namespace protocols::jd3;


RosettaScriptsJobQueen::RosettaScriptsJobQueen()
{

	//This is where you will list the options you will be using.  These options are command-line options
	// that can be specified at the <common> and <job> level.

	utility::options::OptionKeyList opts;

	/// Static Functions setup the OptionKeyList:
	//
	RosettaScriptsParser::list_options_read( opts );

	add_options( opts );
	parser_ = RosettaScriptsParserOP( new RosettaScriptsParser() );
	resource_manager_ = std::make_shared< basic::resource_manager::ResourceManager >();
	if ( basic::options::option[ basic::options::OptionKeys::jd3::resource_definition_files ].user() ) {

		for ( std::string const & fname : basic::options::option[ basic::options::OptionKeys::jd3::resource_definition_files ]  ) {
			std::string contents;
			try {
				contents = utility::file_contents( fname );
			} catch ( utility::excn::Exception const & e ) {
				std::ostringstream oss;
				oss << "Failed to open requested resource definition file named \"" << fname <<"\"" << std::endl;
				throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
			}
			TR << "Reading resource definitions from " << fname << std::endl;
			resource_manager_->read_resources_from_xml( fname, contents );
		}
		TR << "Read definitions for " << resource_manager_->n_resources_declared() << " resources" << std::endl;
	}
}

RosettaScriptsJobQueen::~RosettaScriptsJobQueen() {}

protocols::jd3::JobDigraphOP
RosettaScriptsJobQueen::initial_job_dag()
{
	// Ask the base class to read the job-definition file or the command line;
	// next, we'll construct the Tags for each of the preliminary job nodes
	// and write down which Resources depend on which preliminary job nodes
	JobDigraphOP sjq_digraph = parent::initial_job_dag();
	initialize_resources_for_preliminary_job_nodes();
	return sjq_digraph;
}


protocols::jd3::JobOP
RosettaScriptsJobQueen::complete_larval_job_maturation(
	protocols::jd3::LarvalJobCOP larval_job,
	utility::options::OptionCollectionCOP job_options,
	utility::vector1< protocols::jd3::JobResultCOP > const &
)
{

	TR << "Completing larval job maturation" << std::endl;

	protocols::jd3::standard::MoverAndPoseJobOP mature_job( new protocols::jd3::standard::MoverAndPoseJob );
	core::pose::PoseOP pose = pose_for_job( larval_job, *job_options );
	mature_job->pose( pose );

	bool modified_pose;

	parser_->set_recursion_limit( *job_options );
	protocols::jd3::JobOutputIndex index;
	index.primary_output_index = larval_job->nstruct_index();
	protocols::rosetta_scripts::ParsedProtocolOP mover_protocol = parser_->generate_mover_and_apply_to_pose(
		*pose, *job_options, modified_pose,
		larval_job->input_tag(),
		larval_job->job_tag_with_index_suffix( index ),
		resource_manager_ );

	mature_job->mover( mover_protocol );
	return mature_job;
}

/// @brief the RosettaScriptsJobQueen needs to tell remote nodes to deallocate Resources held by the ResourceManager that are no longer needed
std::list< deallocation::DeallocationMessageOP >
RosettaScriptsJobQueen::deallocation_messages()
{
	// TO DO: figure out whether any resources need to be deallocated
	std::list< deallocation::DeallocationMessageOP > messages = parent::deallocation_messages();
	messages.splice( messages.end(), deallocation_messages_to_send_ );
	return messages;
}

/// @brief A deallocation message first sent to the JobDistributor on this host originating from
/// a remote JobQueen
void
RosettaScriptsJobQueen::derived_process_deallocation_message( deallocation::DeallocationMessageOP message )
{
	using namespace deallocation;
	ResourceDeallocationMessageOP resource_message = utility::pointer::dynamic_pointer_cast< ResourceDeallocationMessage >( message );
	if ( resource_message ) {
		resource_manager_->deallocate_resource( resource_message->resource_name() );
		TR << "Deallocated resource \"" << resource_message->resource_name() << "\"; " << resource_manager_->n_resources_in_memory() << " resources currently reside in memory" << std::endl;
	}
}

void
RosettaScriptsJobQueen::note_preliminary_job_node_is_complete( core::Size pjn_index )
{
	for ( std::string const & resource_name : resources_for_pjn_[ pjn_index ] ) {
		debug_assert( pjns_requiring_resources_[ resource_name ].count( pjn_index ) );
		pjns_requiring_resources_[ resource_name ].erase( pjn_index );
		if ( pjns_requiring_resources_[ resource_name ].empty() ) {
			deallocation_messages_to_send_.push_back(
				std::make_shared< deallocation::ResourceDeallocationMessage >( resource_name ) );
			resource_manager_->deallocate_resource( resource_name );
		}
	}
	resources_for_pjn_[ pjn_index ].clear();
}

void
RosettaScriptsJobQueen::initialize_resources_for_preliminary_job_nodes()
{
	using utility::tag::TagCOP;
	using namespace protocols::jd3::standard;

	utility::vector1< PreliminaryLarvalJob > const & pjns( parent::preliminary_larval_jobs() );
	resources_for_pjn_.resize( pjns.size() );
	for ( Size ii = 1; ii <= pjns.size(); ++ii ) {
		utility::options::OptionCollectionCOP iiopts = pjns[ ii ].job_options;
		TagCOP rstag = parser_->create_tag_from_xml( (*iiopts)[ basic::options::OptionKeys::parser::protocol ], *iiopts );
		if ( rstag->hasTag( "RESOURCES" ) ) {
			TagCOP resources_tag = rstag->getTag( "RESOURCES" );
			for ( TagCOP const & resource_tag : resources_tag->getTags() ) {
				std::string resource_name( resource_tag->getOption< std::string >( "name" ) );
				pjns_requiring_resources_[ resource_name ].insert( ii );
				resources_for_pjn_[ ii ].push_back( resource_name );
			}
		}
	}
}

basic::resource_manager::ResourceManagerOP
RosettaScriptsJobQueen::resource_manager()
{
	return resource_manager_;
}


}
}
