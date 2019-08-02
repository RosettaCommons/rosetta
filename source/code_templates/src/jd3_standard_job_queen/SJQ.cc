// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

// Unit headers
#include <--path--/--class--.hh>

// protocol headers
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd3/JobDistributor.hh>
//#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/JobGenealogist.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/dag_node_managers/NodeManager.hh>
#include <protocols/jd3/dag_node_managers/SimpleNodeManager.hh>
#include <protocols/jd3/deallocation/ResourceDeallocationMessage.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/jobs/MoverJob.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/ConstDataMap.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>


#include <basic/resource_manager/ResourceManager.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


--namespace--

static basic::Tracer TR("--namespace_dot--.--class--");

using namespace protocols::jd3;
using namespace protocols::jd3::standard;

using namespace protocols::jd3::jobs;
using namespace protocols::jd3::job_results;
using namespace protocols::jd3::job_summaries;

--class--::--class--():
    StandardJobQueen()
{


	//This is where you will list the options you will be using.  These options are command-line options
	// that can be specified at the <common> and <job> level.

	utility::options::OptionKeyList opts;

	/// Static Functions setup the OptionKeyList:
	//
	//core::scoring::list_read_options_in_get_score_function( opts );
	//core::pack::task::PackerTask::list_options_read( opts );
	//core::pack::task::operation::ReadResfile::list_options_read( opts );
	//protocols::simple_moves::PackRotamersMover::list_options_read( opts );

	add_options( opts );

	/// Options defined in this app:
	//
	//add_option( basic::options::OptionKeys::minimize_sidechains );
	//add_option( basic::options::OptionKeys::min_pack );
	//add_option( basic::options::OptionKeys::off_rotamer_pack );
}

--class--::~--class--() {}

///@details
/// https://www.rosettacommons.org/docs/wiki/development_documentation/tutorials/jd3_derived_jq/creating_the_job_dag
///
/// Ask the base class to read the job-definition file or the command line;
/// next, we'll construct the Tags for each of the preliminary job nodes
/// and write down which Resources depend on which preliminary job nodes
///
protocols::jd3::JobDigraphOP
--class--::initial_job_dag()
{

    return StandardJobQueen::initial_job_dag();

}

///@details
/// https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/larval_jobs
///
std::list< jd3::LarvalJobOP >
--class--::determine_job_list (
        core::Size job_dag_node_index,
        core::Size max_njobs )
{
    return StandardJobQueen::determine_job_list(job_dag_node_index, max_njobs);
}

LarvalJobs
--class--::next_batch_of_larval_jobs_for_job_node( Size job_dag_node_index, Size max_njobs ){
    return StandardJobQueen::next_batch_of_larval_jobs_for_job_node( job_dag_node_index, max_njobs);
}

protocols::jd3::JobOP
--class--::complete_larval_job_maturation(
	protocols::jd3::LarvalJobCOP larval_job,
	utility::options::OptionCollectionCOP job_options,
	utility::vector1< protocols::jd3::JobResultCOP > const & )
{

	TR << "Completing larval job maturation" << std::endl;

	/// This is the core of your code.  Most likely, you will need the lines below.
	MoverJobOP mature_job( utility::pointer::make_shared< MoverJob >() );
	core::pose::PoseOP pose = StandardJobQueen::pose_for_job( larval_job, *job_options );
	mature_job->pose( pose );

	utility::tag::TagCOP job_tag;
	if ( larval_job->inner_job()->jobdef_tag() )  {
		job_tag = larval_job->inner_job()->jobdef_tag();
	}

	//Here, you create your mover(s) and call what you need to call. You need to give your protocol a LOCAL
	// Options collection from which to read from.
	//SomeMoverOP mover_protocol  = SomeMoverOP( new SomeMover());
	//mover_protocol->initialize_from_options(job_options); //You will need to write this in your protocol!
    //mature_job->mover( mover_protocol );


	// Setup your mover in whatever way you need.
	// SEE the FixBB JD3 example in src/apps/pilot/andrew/fixbb_jd3.cc
	// for more complex setup options, such as combining XML and command-line syntax for specifying job-level options.



	//Set any SimpleMetrics to run after the job and included in the JobSummary.
	//mature_job->metrics( simple_metrics_for_job, "job_node_1" );
	return mature_job;
}

///@details
/// https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/completed_job_summary
///
void
--class--::completed_job_summary(
	jd3::LarvalJobCOP job,
	core::Size result_index,
	jd3::JobSummaryOP summary)
{

	return StandardJobQueen::completed_job_summary( job, result_index, summary);
}

///@details
/// https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/discarding_job_results
///
std::list< jd3::JobResultID >
--class--::job_results_that_should_be_discarded() {

	return StandardJobQueen::job_results_that_should_be_discarded();
}

///@details
/// https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/outputting_results
///
std::list< jd3::output::OutputSpecificationOP >
--class--::jobs_that_should_be_output(){

    return StandardJobQueen::jobs_that_should_be_output();
}

///@details Here, you specify which tags that are defined in the <common> element.
///
void
--class--::append_common_tag_subelements(
	utility::tag::XMLSchemaDefinition & job_definition_xsd,
	utility::tag::XMLSchemaComplexTypeGenerator & ct_gen) const
{
	using namespace utility::tag;

    StandardJobQueen::append_common_tag_subelements( job_definition_xsd, ct_gen);

	/*
	ScoreFunctionLoader::provide_xml_schema( job_definition_xsd );
	TaskOperationLoader::provide_xml_schema( job_definition_xsd );

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_already_defined_subelement( ScoreFunctionLoader::loader_name(),     & ScoreFunctionLoader::score_function_loader_ct_namer )
		.add_already_defined_subelement( TaskOperationLoader::loader_name(),     & TaskOperationLoader::task_op_loader_ct_namer );
	ct_gen.add_ordered_subelement_set_as_repeatable( subelements );
	*/

}

///@details Here you can specify Schema's to use when creating your XSD for a specific JOB.
///
void
--class--::append_job_tag_subelements(
	utility::tag::XMLSchemaDefinition & job_definition_xsd,
	utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen) const
{
	using namespace utility::tag;

    StandardJobQueen::append_job_tag_subelements( job_definition_xsd, job_ct_gen);


	// Example for FixBB JD3 is below:

	/*
	PackRotamersMover::provide_xml_schema(   job_definition_xsd );
	MinMover::provide_xml_schema(            job_definition_xsd );
	ScoreFunctionLoader::provide_xml_schema( job_definition_xsd );
	TaskOperationLoader::provide_xml_schema( job_definition_xsd );

	//task operations and sfxns -- as many as you want
	XMLSchemaSimpleSubelementList task_and_sfxn_subelements;
	task_and_sfxn_subelements
		.add_already_defined_subelement( ScoreFunctionLoader::loader_name(), & ScoreFunctionLoader::score_function_loader_ct_namer )
		.add_already_defined_subelement( TaskOperationLoader::loader_name(), & TaskOperationLoader::task_op_loader_ct_namer );
	job_ct_gen.add_ordered_subelement_set_as_repeatable( task_and_sfxn_subelements );

	//pack -- at most one
	XMLSchemaSimpleSubelementList pack_subelement;
	pack_subelement.add_already_defined_subelement( PackRotamersMoverCreator::mover_name(), & complex_type_name_for_mover );
	job_ct_gen.add_ordered_subelement_set_as_optional( pack_subelement );

	//min -- at most one
	XMLSchemaSimpleSubelementList min_subelement;
	min_subelement.add_already_defined_subelement( MinMoverCreator::mover_name(), & complex_type_name_for_mover );
	job_ct_gen.add_ordered_subelement_set_as_optional( min_subelement );

	*/

}





--end_namespace--
