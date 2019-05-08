// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/StandardJobQueen.cc
/// @brief  StandardJobQueen class's method definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) Simplification/Refactor

//unit headers
#include <protocols/jd3/standard/StandardJobQueen.hh>

// package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/standard/PreliminaryLarvalJob.hh>
#include <protocols/jd3/standard/PreliminaryLarvalJobTracker.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/JobTracker.hh>
#include <protocols/jd3/jobs/MoverJob.hh>
#include <protocols/jd3/output/MultipleOutputSpecification.hh>
#include <protocols/jd3/output/MultipleOutputter.hh>
#include <protocols/jd3/pose_inputters/PoseInputSource.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.hh>
#include <protocols/jd3/pose_inputters/PoseInputterCreator.hh>
#include <protocols/jd3/pose_inputters/PoseInputterFactory.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.hh>
#include <protocols/jd3/deallocation/DeallocationMessage.hh>
#include <protocols/jd3/deallocation/InputPoseDeallocationMessage.hh>
#include <protocols/jd3/util.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose_options.hh>
// #include <basic/resource_manager/JobOptions.hh>

//utility headers
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
#include <utility/vector1.hh>

//basic headers
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>
#include <basic/options/option.cc.gen.hh>
#include <basic/datacache/ConstDataMap.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.jd3.standard.StandardJobQueen" );


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/list.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace standard {

using namespace utility::pointer;

StandardJobQueen::StandardJobQueen() :
	use_factory_provided_pose_inputters_( true ),
	use_factory_provided_pose_outputters_( true ),
	override_default_outputter_( false ),
	common_block_precedes_job_blocks_( true ),
	required_initialization_performed_( false ),
	preliminary_larval_jobs_determined_( false ),
	curr_inner_larval_job_index_( 0 ),
	njobs_made_for_curr_inner_larval_job_( 0 ),
	input_pose_counter_( 0 ),
	load_starting_poses_only_once_( false )
{
	// begin to populate the per-job options object
	pose_inputters::PoseInputterFactory::get_instance()->list_options_read( inputter_options_ );
	pose_outputters::PoseOutputterFactory::get_instance()->list_outputter_options_read( outputter_options_ );
	pose_outputters::PoseOutputterFactory::get_instance()->list_secondary_outputter_options_read( secondary_outputter_options_ );

	if ( basic::options::option[ basic::options::OptionKeys::jd3::load_input_poses_only_once ].value() ) {
		load_starting_poses_only_once_ = true;
	}
	prelim_job_tracker_ = make_shared< PreliminaryLarvalJobTracker >();
}

StandardJobQueen::~StandardJobQueen() = default;


std::string
StandardJobQueen::job_definition_xsd() const
{
	using namespace utility::tag;
	using namespace utility::options;
	using namespace basic::options;

	utility::options::OptionKeyList all_options = concatenate_all_options();

	XMLSchemaDefinition xsd;

	if ( use_factory_provided_pose_inputters_ ) {
		pose_inputters::PoseInputterFactory::get_instance()->define_pose_inputter_xml_schema( xsd );
	} else {
		utility::tag::define_xml_schema_group(
			inputter_creators_,
			pose_inputters::PoseInputterFactory::pose_inputter_xml_schema_group_name(),
			& pose_inputters::PoseInputterFactory::complex_type_name_for_pose_inputter,
			xsd );
	}

	XMLSchemaSimpleSubelementList input_subelements;
	input_subelements.add_group_subelement( & pose_inputters::PoseInputterFactory::pose_inputter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator input_ct;
	input_ct
		.element_name( "Input" )
		.description("XRW TO DO")
		.complex_type_naming_func( & job_def_complex_type_name )
		.set_subelements_pick_one( input_subelements )
		.write_complex_type_to_schema( xsd );

	if ( use_factory_provided_pose_outputters_ ) {
		pose_outputters::PoseOutputterFactory::get_instance()->define_pose_outputter_xml_schema( xsd );
	} else {
		pose_outputters::PoseOutputterFactory::get_instance()->define_secondary_pose_outputter_xml_schema( xsd );
		utility::tag::define_xml_schema_group(
			outputter_creators_,
			pose_outputters::PoseOutputterFactory::pose_outputter_xml_schema_group_name(),
			& pose_outputters::PoseOutputterFactory::complex_type_name_for_pose_outputter,
			xsd );
	}

	XMLSchemaSimpleSubelementList output_subelements;
	output_subelements.add_group_subelement( & pose_outputters::PoseOutputterFactory::pose_outputter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator output_ct;
	output_ct
		.element_name( "Output" )
		.description( "Give the (primary) outputter for this job. The primary outputter is responsible"
		"for giving a job its name and for writing the coordinates of the Pose to disk. Only a "
		"single (primary) outputter may be given for a job.")
		.complex_type_naming_func( & job_def_complex_type_name )
		.set_subelements_pick_one( output_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList secondary_output_subelements;
	secondary_output_subelements.add_group_subelement( & pose_outputters::PoseOutputterFactory::secondary_pose_outputter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator secondary_output_ct;
	secondary_output_ct
		.element_name( "SecondaryOutput" )
		.description( "List the secondary outputters that should be used for this job; as many secondary outputters as desired may be given for a job, and a particular secondary outputter may be given multiple times." )
		.complex_type_naming_func( & job_def_complex_type_name )
		.set_subelements_repeatable( secondary_output_subelements )
		.write_complex_type_to_schema( xsd );

	// write out option set -- this should be method-extracted into a place accessible to other job queens
	if  ( ! all_options.empty() ) {
		XMLSchemaComplexTypeGenerator option_generator;
		XMLSchemaSimpleSubelementList option_subelements;

		std::set< utility::keys::VariantKey< utility::options::OptionKey > > already_output_options;

		for ( auto const & iter : all_options ) {
			AttributeList attributes;
			utility::options::OptionKey const & opt_key( iter() );

			// only output each option once, even if it is read in more than
			// one context
			if ( already_output_options.count( opt_key ) ) continue;
			already_output_options.insert( opt_key );

			OptionTypes opt_type = option_type_from_key( opt_key );
			XMLSchemaType value_attribute_type = value_attribute_type_for_option( opt_type );
			if ( option[ opt_key ].has_default() ) {
				if ( opt_type == BOOLEAN_OPTION ) {
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, "XRW TO DO",  option[ opt_key ].raw_default_string( ) );
				} else {
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, "XRW TO DO",  option[ opt_key ].raw_default_string( ) );
				}
			} else { // no default; value is required, unless it's a boolean option
				if ( opt_type == BOOLEAN_OPTION ) {
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, "XRW TO DO",  "false"  );
				} else {
					attributes + XMLSchemaAttribute::required_attribute( "value", value_attribute_type , "XRW TO DO" );
				}
			}

			std::string decolonized_name = basic::options::replace_option_namespace_colons_with_underscores( iter );
			option_subelements.add_simple_subelement( decolonized_name, attributes , "");
		}
		option_generator.element_name( "Options" )
			.description( "XRW TO DO" )
			.complex_type_naming_func( & job_def_complex_type_name )
			.set_subelements_single_appearance_optional( option_subelements )
			.write_complex_type_to_schema( xsd );
	}

	XMLSchemaComplexTypeGenerator job_ct;

	// Job <Input> element is required
	XMLSchemaSimpleSubelementList job_input_subelement;
	job_input_subelement.add_already_defined_subelement( "Input", & job_def_complex_type_name );
	job_ct.add_ordered_subelement_set_as_optional( job_input_subelement );

	// Job <Output> element is optional
	XMLSchemaSimpleSubelementList job_output_subelement;
	job_output_subelement.add_already_defined_subelement( "Output", & job_def_complex_type_name );
	job_ct.add_ordered_subelement_set_as_optional( job_output_subelement );

	// Job <SecondaryOutput> element is optional
	XMLSchemaSimpleSubelementList job_secondary_output_subelement;
	job_secondary_output_subelement.add_already_defined_subelement( "SecondaryOutput", & job_def_complex_type_name );
	job_ct.add_ordered_subelement_set_as_optional( job_secondary_output_subelement );

	// Job <Options> element is optional
	XMLSchemaSimpleSubelementList job_options_subelement;
	job_options_subelement.add_already_defined_subelement( "Options", & job_def_complex_type_name );
	job_ct.add_ordered_subelement_set_as_optional( job_options_subelement );

	// Ask the derived class for what else belongs in the Job element.
	// NOTE: The derived class should only call the add_ordered_subelement_set_* functions
	// of the XMLSchemaComplexTypeGenerator or the <Input>, <Output>, and <Option>
	// subelements will be overwritten.
	append_job_tag_subelements( xsd, job_ct );

	// verify that the derived class did not call anything besides the add_ordered_subelement_set_*
	// functions.
	if ( job_ct.subelement_behavior() != se_ordered_sets ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Subclass of StandardJobQueen's append_job_tag_subelements"
			" method invokes a method of the XMLSchemaComplexTypeGenerator that overwrote the <Input>, <Output>, and"
			" <Options> elements.  It should only call methods named \"add_ordered_subelement_set_*\"" );
	}

	job_ct
		.element_name( "Job" )
		.description( "XRW TO DO" )
		.complex_type_naming_func( & job_def_complex_type_name )
		.add_attribute( XMLSchemaAttribute(  "nstruct", xsct_non_negative_integer, "The number of times to repeat this job -- i.e. the"
		" number of output structures to generate from this job. If not provided, then the command-line flag -nstruct will be read from." ))
		.write_complex_type_to_schema( xsd );

	XMLSchemaComplexTypeGenerator common_block_ct_gen;

	// Common block <Options> subelement
	XMLSchemaSimpleSubelementList common_block_option_subelement;
	common_block_option_subelement.add_already_defined_subelement( "Options", & job_def_complex_type_name );
	common_block_ct_gen.add_ordered_subelement_set_as_optional( common_block_option_subelement );

	// Ask the derived class for what else belongs in the Common element.
	// NOTE: The derived class should only call the add_ordered_subelement_set_* functions
	// of the XMLSchemaComplexTypeGenerator or the <Options> subelement will be
	// overwritten.
	append_common_tag_subelements( xsd, common_block_ct_gen );

	// verify that the derived class did not call anything besides the add_ordered_subelement_set_*
	// functions.
	if ( common_block_ct_gen.subelement_behavior() != se_ordered_sets ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Subclass of StandardJobQueen's append_job_tag_subelements"
			" method invokes a method of the XMLSchemaComplexTypeGenerator that overwrote the <Input>, <Output>, and"
			" <Options> elements.  It should only call methods named \"add_ordered_subelement_set_*\"" );
	}

	common_block_ct_gen
		.element_name( "Common" )
		.description( "XRW TO DO" )
		.complex_type_naming_func( & job_def_complex_type_name )
		.write_complex_type_to_schema( xsd );

	XMLSchemaComplexTypeGenerator job_def_file_ct;
	XMLSchemaSimpleSubelementList job_def_subelements;

	if ( common_block_precedes_job_blocks_ ) {
		job_def_subelements.add_already_defined_subelement( "Common", & job_def_complex_type_name, 0, 1 );
		job_def_subelements.add_already_defined_subelement( "Job", & job_def_complex_type_name, 1, xsminmax_unbounded );
	} else {
		job_def_subelements.add_already_defined_subelement( "Job", & job_def_complex_type_name, 1, xsminmax_unbounded );
		job_def_subelements.add_already_defined_subelement( "Common", & job_def_complex_type_name, 0, 1 );
	}
	job_def_file_ct.element_name( "JobDefinitionFile" )
		.complex_type_naming_func( & job_def_complex_type_name )
		.description( "XRW TO DO" )
		.set_subelements_single_appearance_required_and_ordered( job_def_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaElement root_element;
	root_element.name( "JobDefinitionFile" ).type_name( job_def_complex_type_name( "JobDefinitionFile" ));
	xsd.add_top_level_element( root_element );

	std::string xsd_string = xsd.full_definition();

	try {
		utility::tag::test_if_schema_is_valid( xsd_string );
	} catch ( utility::excn::Exception const & e ) {
		std::ostringstream oss;
		oss << "The XML Schema for the job definition file is invalid.  The error message is:\n" << e.msg()
			<< "\nAnd the whole schema is:\n" << xsd_string << "\nThis executable cannot be used in its"
			<< " current state.\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  oss.str() );
	}

	return xsd_string;

}

std::string
StandardJobQueen::resource_definition_xsd() const
{
	// TO DO!
	return "";
}

/// @details The base class provides a digraph with a single node -- that is, all the jobs
/// are independent of each other.  This is equivalent to the kind of jobs that could be
/// run in JD2.
JobDigraphOP
StandardJobQueen::create_initial_job_dag()
{
	determine_preliminary_job_list();

	// create a DAG with as many nodes in it as there are preliminary larval jobs
	JobDigraphOP job_graph = make_shared< JobDigraph >( preliminary_larval_jobs_.size() );
	return job_graph;
}

void StandardJobQueen::update_job_dag( JobDigraphUpdater & ) {}

/// @details The process begins by first constructing the job definition and resource definition
/// XSDs.  With these schemas, the %StandardJobQueen validates the input XML files (if present).
/// The %StandardJobQueen then populates preliminary versions of LarvalJob objects./ If the XSD
/// includes "command line options" (which may be specified either from the command line or in
/// the <options> sections of the Job XML file), the %StandardJobQueen loads the preliminary
/// LarvalJob objects with the options. These preliminary LarvalJob objects will not have been
/// nstruct expanded (i.e. if there are 100 nstruct for each of 5 different jobs, then there will
/// only be 5 preliminary larval jobs created). It then passes the preliminary LarvalJob list and
/// the TagOP objects for each preliminary LarvalJob to the derived class through the
/// refine_job_list method.
LarvalJobs
StandardJobQueen::determine_job_list( Size job_dag_node_index, Size max_njobs )
{
	// ok -- we're going to look for a job definition file, and failing that, fall back on
	// the PoseInputterFactory to determine where the input sources are coming from.

	if ( ! preliminary_larval_jobs_determined_ ) {
		determine_preliminary_job_list();
	}

	LarvalJobs larval_jobs;
	if ( job_dag_node_index <= preliminary_larval_jobs_.size() ) {

		// now that the PreliminaryLarvalJobs have been constructed, go ahead
		// and start delivering LarvalJobs to the JobDistributor.
		larval_jobs = next_batch_of_larval_jobs_from_prelim( job_dag_node_index, max_njobs );
	} else {
		larval_jobs = next_batch_of_larval_jobs_for_job_node( job_dag_node_index, max_njobs );
	}

	return larval_jobs;
}

LarvalJobs
StandardJobQueen::next_batch_of_larval_jobs_from_prelim( core::Size job_node_index, core::Size max_njobs )
{
	LarvalJobs jobs;
	if ( prelim_job_tracker_->get_job_node_assigned( job_node_index ) ) return jobs;

	while ( true ) {
		// Each iteration through this loop advances either  curr_inner_larval_job_, or
		// njobs_made_for_curr_inner_larval_job_ to ensure
		// that this loop does not continue indefinitely.

		if ( curr_inner_larval_job_index_ == 0 ) {

			PreliminaryLarvalJob curr_prelim_job = preliminary_larval_jobs_[ job_node_index ];
			inner_larval_jobs_for_curr_prelim_job_ = refine_preliminary_job( curr_prelim_job );
			curr_inner_larval_job_index_ = 1;
			njobs_made_for_curr_inner_larval_job_ = 0;

			// now initialize the inner larval jobs in the inner_larval_jobs_for_curr_prelim_job_ vector
			// setting the job tag.
			utility::tag::TagCOP output_tag;
			if ( curr_prelim_job.job_tag && curr_prelim_job.job_tag->hasTag( "Output" ) ) {
				output_tag = curr_prelim_job.job_tag->getTag( "Output" );
			}
			for ( InnerLarvalJobOP iijob : inner_larval_jobs_for_curr_prelim_job_ ) {
				pose_outputters::PoseOutputterOP ii_outputter = pose_outputter_for_job( *iijob );
				ii_outputter->determine_job_tag( output_tag, *curr_prelim_job.job_options, *iijob );
			}

		}

		if ( curr_inner_larval_job_index_ > inner_larval_jobs_for_curr_prelim_job_.size() ) {
			// prepare for the next time we call this function for a different job node
			curr_inner_larval_job_index_ = 0;
			prelim_job_tracker_->track_job_node_assigned(job_node_index, current_job_index() );
			return jobs;
		}

		if ( curr_inner_larval_job_index_ <= inner_larval_jobs_for_curr_prelim_job_.size() ) {
			// create LarvalJobs
			InnerLarvalJobOP curr_inner_larval_job = inner_larval_jobs_for_curr_prelim_job_[ curr_inner_larval_job_index_ ];
			core::Size max_to_make = max_njobs;

			if ( max_to_make > curr_inner_larval_job->nstruct_max() - njobs_made_for_curr_inner_larval_job_ ) {
				max_to_make = curr_inner_larval_job->nstruct_max() - njobs_made_for_curr_inner_larval_job_;
			}

			//Sets last pjn_index or first depending on whether this is the first nstruct being assigned.
			core::Size starting_index = 0; //Starting global ID we will be giving out
			core::Size ending_index   = 0; //Final global ID we will be giving out.


			starting_index = 1 + current_job_index();

			////////////////////////////////////////////////////////////
			///* Expand Job List AND increment the JobTracker's current job index during expansion! *///
			LarvalJobs curr_jobs = expand_job_list( curr_inner_larval_job, max_to_make );

			ending_index = current_job_index();

			core::Size n_made = curr_jobs.size();
			if ( max_njobs >= n_made ) {
				max_njobs -= n_made;
			} else {
				max_njobs = 0;
				// this should never happen!
				throw CREATE_EXCEPTION(utility::excn::Exception,  "expand_job_list returned " + utility::to_string( n_made ) + " jobs when it was given a maximum number of " + utility::to_string( max_to_make ) + " to make (with max_njobs of " + utility::to_string( max_njobs ) + ")\n" );
			}

			jobs.splice( jobs.end(), curr_jobs );

			if ( n_made + njobs_made_for_curr_inner_larval_job_ < curr_inner_larval_job->nstruct_max() ) {
				njobs_made_for_curr_inner_larval_job_ += n_made;
				// update the end index for this node; it is perhaps not the final job that will be submitted
				// for this node, but we want to be able to definitively answer the question for a job
				// that has been submitted as to whether it came from this node.

				//Now we know we are actually giving out jobs, so lets track them.
				prelim_job_tracker_->track_job_node_being_assigned( job_node_index, starting_index, ending_index );
				return jobs;
			} else {
				prelim_job_tracker_->track_job_node_being_assigned( job_node_index, starting_index, ending_index );
				++curr_inner_larval_job_index_;
				njobs_made_for_curr_inner_larval_job_ = 0;
			}
		}
	}
}

LarvalJobs
StandardJobQueen::next_batch_of_larval_jobs_for_job_node(core::Size /*job_node*/, core::Size /*max_jobs*/){
	LarvalJobs jobs;
	return jobs;
}

bool
StandardJobQueen::has_job_previously_been_output( protocols::jd3::LarvalJobCOP job )
{
	utility::options::OptionCollectionCOP job_options = options_for_job( *job->inner_job() );
	return pose_outputter_for_job( *job->inner_job() )->job_has_already_completed( *job, *job_options );
}

JobOP
StandardJobQueen::mature_larval_job(
	protocols::jd3::LarvalJobCOP larval_job,
	utility::vector1< JobResultCOP > const & input_results
)
{
	using namespace utility::options;
	using namespace utility::tag;
	using namespace basic::datacache;

	// There are a few pieces of initialization that need to occur on "remote nodes" upon the first
	// time they are asked to mature larval jobs.
	// E.g., the common_block_tags_ must be loaded prior to the call of complete_larval_job_maturation
	// and the PoseInputters require the Input TagCOP that corresponds to the preliminary job node
	// that this job came from.
	if ( ! required_initialization_performed_ ) {
		determine_preliminary_job_list();
	}

	// initialize the options collection for this job.
	utility::options::OptionCollectionCOP job_options = options_for_job( *larval_job->inner_job() );

	return complete_larval_job_maturation( larval_job, job_options, input_results );
}

/// @details Prepare this job for output by building an OutputSpecification for it and
/// storing this specification in the list of recent successes.
/// Note that tracking of completion of jobs from preliminary job nodes can only occur if
/// the DerivedJobQueen invoked the SJQ's version of next_batch_of_larval_jobs_from_prelim
/// in its determine_job_list method.
void StandardJobQueen::note_job_completed( LarvalJobCOP job, JobStatus status, Size nresults )
{
	prelim_job_tracker_->track_job_completed(job);

	Size pjn_index = prelim_job_tracker_->get_job_node_for_job_index( job->job_index() );
	if ( pjn_index != 0 ) {
		if ( prelim_job_tracker_->get_job_node_complete( pjn_index ) ) {
			// Let the derived JobQueen know that all of the jobs for a PJN have completed
			note_preliminary_job_node_is_complete( pjn_index );
		}
	}

	if ( status == jd3_job_status_success ) {
		utility::options::OptionCollectionOP job_options = options_for_job( *job->inner_job() );
		for ( Size ii = 1; ii <= nresults; ++ii ) {
			create_and_store_output_specification_for_job_result( job, *job_options, ii, nresults );
		}
	}

}

void StandardJobQueen::completed_job_summary( LarvalJobCOP, core::Size, JobSummaryOP ) {}

std::list< output::OutputSpecificationOP >
StandardJobQueen::jobs_that_should_be_output()
{
	std::list< output::OutputSpecificationOP > return_list;
	recent_successes_.swap( return_list );
	return return_list;
}

/// @details Default implementation does not discard any job results.
std::list< JobResultID >
StandardJobQueen::job_results_that_should_be_discarded() {
	return std::list< JobResultID >();
}

output::ResultOutputterOP
StandardJobQueen::result_outputter(
	output::OutputSpecification const & spec
)
{
	using MOS = output::MultipleOutputSpecification;
	using POS = pose_outputters::PoseOutputSpecification;
	using MO = output::MultipleOutputter;
	using MOOP = output::MultipleOutputterOP;
	debug_assert( dynamic_cast< MOS const * > ( &spec ) );
	auto const & mo_spec( static_cast< MOS const & > (spec) );
	MOOP outputters = make_shared< MO > ();
	for ( Size ii = 1; ii <= mo_spec.output_specifications().size(); ++ii ) {
		output::OutputSpecification const & ii_spec( *mo_spec.output_specifications()[ ii ] );
		debug_assert( dynamic_cast< POS const * > (&ii_spec) );
		auto const & ii_pos( static_cast< POS const & > (ii_spec) );
		// Note assumption that there is always a primary pose outputter -- return here
		// when trying to implement the -no_output option for JD3
		pose_outputters::PoseOutputterOP ii_outputter;
		if ( ii == 1 ) {
			ii_outputter = pose_outputter_for_job( ii_pos );
		} else {
			ii_outputter = secondary_outputter_for_job( ii_pos );
		}
		outputters->append_outputter( ii_outputter );
	}
	return outputters;
}


//void StandardJobQueen::completed_job_result(
// LarvalJobCOP job,
// core::Size result_index,
// JobResultOP job_result
//)
//{
// StandardInnerLarvalJobCOP inner_job = utility::pointer::dynamic_pointer_cast< StandardInnerLarvalJob const > ( job->inner_job() );
// if ( ! inner_job ) { throw bad_inner_job_exception(); }
// pose_outputters::PoseOutputterOP outputter = pose_outputter_for_job( *inner_job );
// PoseJobResultOP pose_result = utility::pointer::dynamic_pointer_cast< PoseJobResult >( job_result );
// if ( ! pose_result ) {
//  utility::excn::EXCN_Msg_Exception( "JobResult handed to StandardJobQueen::completed_job_result is not a PoseJobResult or derived from PoseJobResult" );
// }
// utility::options::OptionCollectionOP job_options = options_for_job( *inner_job );
// utility::tag::TagCOP outputter_tag;
// if ( job->inner_job()->jobdef_tag() ) {
//  utility::tag::TagCOP tag = job->inner_job()->jobdef_tag();
//  if ( tag->hasTag( "Output" ) ) {
//   outputter_tag = tag->getTag( "Output" )->getTags()[ 0 ];
//  }
// }
//
// Size const n_results_for_job = results_processed_for_job_[ job->job_index() ].n_results;
// JobOutputIndex output_index;
// output_index.primary_output_index   = job->nstruct_index();
// output_index.n_primary_outputs      = job->nstruct_max();
// output_index.secondary_output_index = result_index;
// output_index.n_secondary_outputs    = n_results_for_job;
//
// assign_output_index( job, result_index, n_results_for_job, output_index );
//
// outputter->write_output_pose( *job, output_index, *job_options, outputter_tag, *pose_result->pose() );
//
// std::list< pose_outputters::SecondaryPoseOutputterOP > secondary_outputters = secondary_outputters_for_job( *inner_job, *job_options );
// for ( std::list< pose_outputters::SecondaryPoseOutputterOP >::const_iterator
//   iter = secondary_outputters.begin(), iter_end = secondary_outputters.end();
//   iter != iter_end; ++iter ) {
//  (*iter)->write_output_pose( *job, output_index, *job_options, *pose_result->pose() );
// }
//
// note_job_result_output_or_discarded( job, result_index );
//}

/// @details Construct the XSD and then invoke the (private) determine_preliminary_job_list_from_xml_file method,
/// which is also invoked by determine_preliminary_job_list.
void
StandardJobQueen::determine_preliminary_job_list_from_xml_file( std::string const & job_def_string )
{
	std::string job_def_schema = job_definition_xsd();
	determine_preliminary_job_list_from_xml_file( job_def_string, job_def_schema );
}

void
StandardJobQueen::flush()
{
	for ( PoseOutputterMap::const_iterator iter = pose_outputter_map_.begin(); iter != pose_outputter_map_.end(); ++iter ) {
		iter->second->flush();
	}
	for ( SecondaryOutputterMap::const_iterator iter = secondary_outputter_map_.begin(); iter != secondary_outputter_map_.end(); ++iter ) {
		iter->second->flush();
	}
}

std::list< deallocation::DeallocationMessageOP >
StandardJobQueen::deallocation_messages()
{

	std::list< deallocation::DeallocationMessageOP > messages;
	// It will be safe to deallocate an input pose if all of the preliminary job nodes
	// that relied on that input pose have completed.

	utility::vector1< core::Size > poses_to_deallocate = prelim_job_tracker_->get_input_poses_to_deallocate();
	for ( core::Size const & pose_id : poses_to_deallocate  ) {

		deallocation::DeallocationMessageOP msg = make_shared< deallocation::InputPoseDeallocationMessage >( pose_id );
		messages.push_back( msg );

		input_pose_index_map_.erase( pose_id );
		prelim_job_tracker_->track_deallocated_input_pose(pose_id);
	}
	return messages;
}

void
StandardJobQueen::process_deallocation_message(
	deallocation::DeallocationMessageOP message
)
{
	using namespace deallocation;
	if ( message->deallocation_type() == input_pose_deallocation_msg ) {
		InputPoseDeallocationMessageOP pose_dealloc_msg =
			utility::pointer::dynamic_pointer_cast< InputPoseDeallocationMessage > ( message );
		if ( input_pose_index_map_.count( pose_dealloc_msg->pose_id() ) ) {
			TR << "Erasing previously-stored Pose with index " << pose_dealloc_msg->pose_id() << std::endl;
			input_pose_index_map_.erase( pose_dealloc_msg->pose_id() );
		}
	} else {
		// any other type of deallocation message, for now, will have to be
		// handled by the derived class.
		derived_process_deallocation_message( message );
	}
}


std::string
StandardJobQueen::job_def_complex_type_name( std::string const & type )
{ return "job_def_" + type + "_type"; }


void
StandardJobQueen::do_not_accept_all_pose_inputters_from_factory()
{
	if ( use_factory_provided_pose_inputters_ ) {
		use_factory_provided_pose_inputters_ = false;
		inputter_options_.clear();
	}
}

void
StandardJobQueen::allow_pose_inputter( pose_inputters::PoseInputterCreatorOP creator )
{
	if ( use_factory_provided_pose_inputters_ ) {
		use_factory_provided_pose_inputters_ = false;
		// the user has not told us they want to drop the original set of creators
		// so interpret this as them wanting to allow the user to provide another
		// inputter.
		inputter_creator_list_ = pose_inputters::PoseInputterFactory::get_instance()->pose_inputter_creators();
		for ( auto const & creator_ii : inputter_creator_list_ ) {
			inputter_creators_[ creator_ii->keyname() ] = creator_ii;
		}
	}
	if ( ! inputter_creators_.count( creator->keyname() ) ) {
		inputter_creator_list_.push_back( creator );
		inputter_creators_[ creator->keyname() ] = creator;
		creator->list_options_read( inputter_options_ );
	}
}

void
StandardJobQueen::do_not_accept_all_pose_outputters_from_factory()
{
	if ( use_factory_provided_pose_outputters_ ) {
		use_factory_provided_pose_outputters_ = false;
		override_default_outputter_ = true;
		outputter_options_.clear();
	}
}

void
StandardJobQueen::allow_pose_outputter( pose_outputters::PoseOutputterCreatorOP creator )
{
	if ( use_factory_provided_pose_outputters_ ) {
		use_factory_provided_pose_outputters_ = false;
		// the user has not told us they want to drop the original set of creators
		// so interpret this as them wanting to allow the user to provide another
		// outputter in addition to the original set.
		outputter_creator_list_ = pose_outputters::PoseOutputterFactory::get_instance()->pose_outputter_creators();
		for ( auto const & creator_ii : outputter_creator_list_ ) {
			outputter_creators_[ creator_ii->keyname() ] = creator_ii;
		}
	}
	if ( ! outputter_creators_.count( creator->keyname() ) ) {
		outputter_creator_list_.push_back( creator );
		outputter_creators_[ creator->keyname() ] = creator;
		creator->list_options_read( outputter_options_ );
		if ( override_default_outputter_ && ! default_outputter_creator_ ) {
			default_outputter_creator_ = creator;
		}
	}
}

void
StandardJobQueen::set_default_outputter( pose_outputters::PoseOutputterCreatorOP creator )
{
	override_default_outputter_ = true;
	default_outputter_creator_ = creator;
	if ( use_factory_provided_pose_outputters_ ) {
		// the user has not told us they want to drop the original set of creators
		// but to provide the desired functionality, the SJQ has to handle outputter
		// instantiation herself
		use_factory_provided_pose_outputters_ = false;
		outputter_creator_list_ = pose_outputters::PoseOutputterFactory::get_instance()->pose_outputter_creators();
		for ( auto const & creator_ii : outputter_creator_list_ ) {
			outputter_creators_[ creator_ii->keyname() ] = creator_ii;
		}
	}
}


/// @details Derived classes may choose to not override this method as a way to indicate that
/// they have no additional subtags of the <Job> tag they wish to add. This is a no-op
/// implementation.
void StandardJobQueen::append_job_tag_subelements(
	utility::tag::XMLSchemaDefinition &,
	utility::tag::XMLSchemaComplexTypeGenerator &
) const
{}

/// @details Derived classes may choose to not override this method as a way to indicate that
/// they have no common-block data that needs to be defined.
void
StandardJobQueen::append_common_tag_subelements(
	utility::tag::XMLSchemaDefinition &,
	utility::tag::XMLSchemaComplexTypeGenerator &
) const
{}

void
StandardJobQueen::determine_preliminary_job_list()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::tag;

	// only should be called once
	if ( preliminary_larval_jobs_determined_ ) return;
	preliminary_larval_jobs_determined_ = true;
	required_initialization_performed_ = true;

	// Always generate a job definition schema and in the process, validate the schema.
	// If derived JobQueens have defined an invalid schema, execution must stop.
	std::string job_def_schema = job_definition_xsd();
	// TO DO: validate the XSD

	if ( option[ in::file::job_definition_file ].user() ) {
		std::string job_def_string = utility::file_contents( option[ in::file::job_definition_file ] );
		determine_preliminary_job_list_from_xml_file( job_def_string, job_def_schema );
	} else {
		determine_preliminary_job_list_from_command_line();
	}

	parse_job_definition_tags( common_block_tags_, get_preliminary_larval_jobs() );
}

PreliminaryLarvalJobTracker &
StandardJobQueen::get_prelim_larval_job_tracker(){
	return *prelim_job_tracker_;
}

PreliminaryLarvalJobTracker const &
StandardJobQueen::get_prelim_larval_job_tracker() const {
	return *prelim_job_tracker_;
}

void
StandardJobQueen::note_preliminary_job_node_is_complete( core::Size )
{}

/// @details This base class implementation merely returns a one-element list containing the
/// input inner_job.  Derived classes have the flexibility to create several preliminary
/// jobs from this input job
InnerLarvalJobs
StandardJobQueen::refine_preliminary_job( PreliminaryLarvalJob const & prelim_job )
{
	InnerLarvalJobs one_job( 1, prelim_job.inner_job );
	return one_job;
}

LarvalJobs
StandardJobQueen::expand_job_list( InnerLarvalJobOP inner_job, core::Size max_larval_jobs_to_create )
{
	core::Size nstruct = inner_job->nstruct_max();
	LarvalJobs jobs;
	core::Size n_to_make = std::min( nstruct, max_larval_jobs_to_create );
	for ( core::Size jj = 1; jj <= n_to_make; ++jj ) {
		get_job_tracker().increment_current_job_index();
		LarvalJobOP job = make_shared< LarvalJob >( inner_job, njobs_made_for_curr_inner_larval_job_ + jj, current_job_index() );
		//TR << "Expand larval job " << current_job_index() << std::endl;
		jobs.push_back( job );
	}
	return jobs;
}

InnerLarvalJobOP
StandardJobQueen::create_and_init_inner_larval_job_from_preliminary( core::Size nstruct, core::Size prelim_job_node ) const
{
	InnerLarvalJobOP inner_job = make_shared< InnerLarvalJob> ( nstruct, prelim_job_node );

	PreliminaryLarvalJob const & prelim_job = preliminary_larval_jobs_[ prelim_job_node ];
	InnerLarvalJobCOP src = prelim_job.inner_job;
	inner_job->job_tag( src->job_tag() );
	inner_job->input_source( src->input_source_cop() );
	inner_job->jobdef_tag( src->jobdef_tag() );
	inner_job->outputter( src->outputter() );

	return inner_job;
}

void
StandardJobQueen::create_and_store_output_specification_for_job_result(
	LarvalJobCOP job,
	core::Size result_index,
	core::Size nresults
)
{

	utility::options::OptionCollectionOP job_options = options_for_job( *job->inner_job() );
	create_and_store_output_specification_for_job_result(
		job, *job_options, result_index, nresults );
}

void
StandardJobQueen::create_and_store_output_specification_for_job_result(
	LarvalJobCOP job,
	utility::options::OptionCollection const & job_options,
	core::Size result_index,
	core::Size nresults
)
{
	recent_successes_.push_back( create_output_specification_for_job_result( job, job_options, result_index, nresults ) );
}



/// @details Creates a MultipleOutputSpecification and packs it with one
/// OutputSpecification per PoseOutputter / SecondaryPoseOutputter. The order
/// of the OutputSpecifications that are given here needs to be recapitulated
/// for the PoseMultipleOutputter that will be created by the call to
/// result_outputter, as the OutputSpecifiations are iterated across in the
/// same order as the PoseOutputters in the PoseMultipleOutputter.
output::OutputSpecificationOP
StandardJobQueen::create_output_specification_for_job_result(
	LarvalJobCOP job,
	utility::options::OptionCollection const & job_options,
	core::Size result_index,
	core::Size nresults
)
{
	InnerLarvalJobCOP inner_job = job->inner_job() ;

	utility::tag::TagCOP outputter_tag;
	utility::tag::TagCOP secondary_output_tag;
	if ( inner_job->jobdef_tag() ) {
		utility::tag::TagCOP job_tag = inner_job->jobdef_tag();
		if ( job_tag->hasTag( "Output" ) ) {
			utility::tag::TagCOP output_tag = job_tag->getTag( "Output" );
			if ( output_tag->getTags().size() != 0 ) {
				outputter_tag = output_tag->getTags()[ 0 ];
			}
		}

		if ( job_tag->hasTag( "SecondaryOutput" ) ) {
			secondary_output_tag = job_tag->getTag( "SecondaryOutput" );
		}
	}

	JobOutputIndex output_index = build_output_index( job, result_index, nresults );

	pose_outputters::PoseOutputterOP outputter = pose_outputter_for_job( *inner_job );

	// Note assumption that there *is* a primary outputter for this job.
	// What does -no_output look like in JD3?
	pose_outputters::PoseOutputSpecificationOP primary_spec = outputter->create_output_specification(
		*job, output_index, job_options, outputter_tag );

	output::MultipleOutputSpecificationOP specs = make_shared< output::MultipleOutputSpecification >();
	specs->append_specification( primary_spec );

	// SecondaryOutputters are ordered:
	// 1st the set of (possibly-repeated) SecondaryOutputters that are given in the Job tag,
	// in the same order as they appear beneath the <SecondaryOutput> subtag.
	// 2nd the set of SecondaryOutputters that are specified on the command line and that are
	// not present in the Job tag.
	utility::vector1< pose_outputters::SecondaryPoseOutputterOP > secondary_outputters =
		secondary_outputters_for_job( *inner_job, job_options );
	for ( core::Size ii = 1; ii <= secondary_outputters.size(); ++ii ) {
		utility::tag::TagCOP secondary_outputter_tag;
		if ( secondary_output_tag ) {
			if ( secondary_output_tag->getTags().size() >= ii ) {
				secondary_outputter_tag = secondary_output_tag->getTags()[ ii-1 ];
			}
		}
		pose_outputters::PoseOutputSpecificationOP spec =
			secondary_outputters[ ii ]->create_output_specification(
			*job, output_index, job_options, secondary_outputter_tag );
		specs->append_specification( spec );
	}
	specs->output_index( output_index );
	specs->result_id( JobResultID( job->job_index(), result_index ) );
	return specs;
}

JobOutputIndex
StandardJobQueen::build_output_index(
	protocols::jd3::LarvalJobCOP job,
	Size result_index,
	Size n_results_for_job
)
{
	JobOutputIndex output_index;
	output_index.primary_output_index   = job->nstruct_index();
	output_index.n_primary_outputs      = std::max( job->nstruct_max(), (Size) 1000 );
	output_index.secondary_output_index = result_index;
	output_index.n_secondary_outputs    = n_results_for_job;

	assign_output_index( job, result_index, n_results_for_job, output_index );
	return output_index;
}

/// @brief No-op implementation -- leave the indices the way they were initialized
void
StandardJobQueen::assign_output_index(
	protocols::jd3::LarvalJobCOP,
	Size,
	Size,
	JobOutputIndex &
)
{}

void
StandardJobQueen::derived_process_deallocation_message(
	deallocation::DeallocationMessageOP
)
{}



void StandardJobQueen::add_options( utility::options::OptionKeyList const & opts )
{
	using namespace utility::options;
	for ( auto const & opt : opts ) {
		options_.push_back( opt );
	}
}

void StandardJobQueen::add_option( utility::options::OptionKey const & key )
{
	options_.emplace_back( key );
}

/// @details TO DO
void StandardJobQueen::remove_default_input_element() {}

utility::tag::TagCOP
StandardJobQueen::common_block_tags() const
{
	return common_block_tags_;
}

utility::options::OptionCollectionOP
StandardJobQueen::options_for_job( InnerLarvalJob const & inner_job ) const {
	return protocols::jd3::options_for_job( concatenate_all_options(), inner_job, common_block_tags_ );
}

utility::options::OptionCollectionOP
StandardJobQueen::options_from_tag(utility::tag::TagCOP job_options_tags) const {
	return protocols::jd3::options_from_tag( concatenate_all_options(), job_options_tags, common_block_tags_);
}

core::pose::PoseOP
StandardJobQueen::pose_for_job(
	LarvalJobCOP job,
	utility::options::OptionCollection const & options
)
{
	return pose_for_inner_job( job->inner_job(), options );
}

core::pose::PoseOP
StandardJobQueen::pose_for_inner_job(
	InnerLarvalJobCOP inner_job,
	utility::options::OptionCollection const & options
)
{
	// either read the Pose in using the pose_inputter (and then keep a copy
	// in the resource manager), or retrieve the Pose from the resource manager
	// initial version: just read the pose in for each job.
	auto const & input_source = dynamic_cast< pose_inputters::PoseInputSource const & >( inner_job->input_source() );
	TR << "Looking for input source " << input_source.input_tag() << " with pose_id " << input_source.source_id() << std::endl;

	if ( input_pose_index_map_.count( input_source.source_id() ) ) {

		core::pose::PoseOP pose = make_shared< core::pose::Pose > ();

		// Why does the standard job queen use detached_copy? Because it is important in multithreaded
		// contexts that no two Poses pointing to that same data end up in two separate threads,
		// and then try to modify that data at the same time.  It turns out there are places
		// in Pose where it covertly modifies data in some other Pose and this could lead to
		// race conditions.
		pose->detached_copy( *input_pose_index_map_[ input_source.source_id() ] );
		return pose;
	}

	utility::tag::TagCOP inputter_tag;
	if ( inner_job->jobdef_tag() && inner_job->jobdef_tag()->hasTag( "Input" ) ) {
		utility::tag::TagCOP input_tag = inner_job->jobdef_tag()->getTag( "Input" );
		if ( input_tag->getTags().size() != 0 ) {
			inputter_tag = input_tag->getTags()[ 0 ];
		}
	}

	core::pose::PoseOP input_pose = pose_inputter_for_job( *inner_job )->pose_from_input_source( input_source, options, inputter_tag );
	TR << "Storing Pose for input source " << input_source.input_tag() << " with pose_id " << input_source.source_id() << std::endl;
	input_pose_index_map_[ input_source.source_id() ] = input_pose;

	core::pose::PoseOP pose = make_shared< core::pose::Pose > ();
	pose->detached_copy( *input_pose );
	return pose;
}

core::pose::PoseOP
StandardJobQueen::pose_for_inner_job(
	InnerLarvalJobCOP inner_job
) {
	runtime_assert( inner_job );
	utility::options::OptionCollectionOP options = options_for_job(* inner_job );
	runtime_assert( options );
	return pose_for_inner_job( inner_job, * options );
}

//ResourceManagerOP StandardJobQueen::resource_manager()
//{}

/// @brief Access the pose inputter
pose_inputters::PoseInputterOP
StandardJobQueen::pose_inputter_for_job( InnerLarvalJob const & inner_job ) const
{
	// find the preliminary job node for this job, if available
	// and return the already-created pose inputter
	if ( inner_job.job_node() > preliminary_larval_jobs_.size() ) {
		if ( use_factory_provided_pose_inputters_ ) {
			return pose_inputters::PoseInputterFactory::get_instance()->new_pose_inputter( inner_job.input_source().origin() );
		} else {
			runtime_assert( inputter_creators_.count( inner_job.input_source().origin() ) );
			auto iter = inputter_creators_.find( inner_job.input_source().origin() );
			return iter->second->create_inputter();
		}
	} else {
		debug_assert( preliminary_larval_jobs_[ inner_job.job_node() ].pose_inputter );
		return preliminary_larval_jobs_[ inner_job.job_node() ].pose_inputter;
	}
}

/// @brief Access the pose outputter for a particular job; perhaps this outputter has already been created, or
/// perhaps it needs to be created and stored
pose_outputters::PoseOutputterOP
StandardJobQueen::pose_outputter_for_job( InnerLarvalJob const & inner_job )
{
	utility::options::OptionCollectionOP job_options = options_for_job( inner_job );
	return pose_outputter_for_job( inner_job, *job_options );
}

pose_outputters::PoseOutputterOP
StandardJobQueen::pose_outputter_for_job(
	InnerLarvalJob const & inner_job,
	utility::options::OptionCollection const & job_options
)
{
	if ( representative_pose_outputter_map_.empty() ) {
		utility_exit_with_message("SJQ: determine_preliminary_job_list() must be called prior to pose_outputter_for_job");
	}

	pose_outputters::PoseOutputterOP representative_outputter;

	representative_outputter = representative_pose_outputter_map_[ inner_job.outputter() ];

	utility::tag::TagCOP output_tag;
	if ( inner_job.jobdef_tag() ) {
		utility::tag::TagCOP job_tag = inner_job.jobdef_tag();
		if ( job_tag->hasTag( "Output" ) ) {
			output_tag = job_tag->getTag( "Output" );
		}
	}

	///Job-specific pose output
	std::string which_outputter = representative_outputter->outputter_for_job( output_tag, job_options, inner_job );

	if ( which_outputter == "" ) {
		return representative_outputter;
	}

	PoseOutputterMap::const_iterator iter = pose_outputter_map_.find( std::make_pair( inner_job.outputter(), which_outputter ) );
	if ( iter != pose_outputter_map_.end() ) {
		return iter->second;
	}

	pose_outputters::PoseOutputterOP outputter;
	if ( use_factory_provided_pose_outputters_ ) {
		outputter = pose_outputters::PoseOutputterFactory::get_instance()->new_pose_outputter( inner_job.outputter() );
	} else {
		outputter = outputter_creators_[ inner_job.outputter() ]->create_outputter();
	}

	pose_outputter_map_[ std::make_pair( inner_job.outputter(), which_outputter ) ] = outputter;
	return outputter;
}

pose_outputters::PoseOutputterOP
StandardJobQueen::pose_outputter_for_job( pose_outputters::PoseOutputSpecification const & spec )
{
	std::string const & outputter_type = spec.outputter_type();
	pose_outputters::PoseOutputterOP representative;
	auto rep_iter = representative_pose_outputter_map_.find( outputter_type );
	if ( rep_iter == representative_pose_outputter_map_.end() ) {
		if ( use_factory_provided_pose_outputters_ ) {
			representative = pose_outputters::PoseOutputterFactory::get_instance()->new_pose_outputter( outputter_type );
		} else {
			runtime_assert( outputter_creators_.count( outputter_type ) != 0 );
			auto iter = outputter_creators_.find( outputter_type );
			representative = iter->second->create_outputter();
		}
		representative_pose_outputter_map_[ outputter_type ] = representative;
	} else {
		representative = rep_iter->second;
	}

	// This outputter name may have a job-distributor-assigned suffix, and so it may not yet
	// have been created, even if this JobQueen has been doing some amount of output.
	std::string actual_outputter_name = representative->outputter_for_job( spec );
	if ( actual_outputter_name == "" ) return representative;

	std::pair< std::string, std::string > outputter_key =
		std::make_pair( outputter_type, actual_outputter_name );
	auto iter = pose_outputter_map_.find( outputter_key );
	if ( iter != pose_outputter_map_.end() ) {
		return iter->second;
	}

	pose_outputters::PoseOutputterOP outputter;
	if ( use_factory_provided_pose_outputters_ ) {
		outputter = pose_outputters::PoseOutputterFactory::get_instance()->new_pose_outputter( outputter_type );
	} else {
		debug_assert( outputter_creators_.count( outputter_type ) );
		outputter = outputter_creators_[ outputter_type ]->create_outputter();
	}
	pose_outputter_map_[ outputter_key ] = outputter;
	return outputter;

}


utility::vector1< pose_outputters::SecondaryPoseOutputterOP >
StandardJobQueen::secondary_outputters_for_job(
	InnerLarvalJob const & inner_job,
	utility::options::OptionCollection const & job_options
)
{
	std::set< std::string > secondary_outputters_added;
	utility::vector1< pose_outputters::SecondaryPoseOutputterOP > secondary_outputters;
	if ( inner_job.jobdef_tag() ) {
		utility::tag::TagCOP job_tag = inner_job.jobdef_tag();
		if ( job_tag->hasTag( "SecondaryOutput" ) ) {
			utility::tag::TagCOP secondary_output_tags = job_tag->getTag( "SecondaryOutput" );
			utility::vector0< utility::tag::TagCOP > const & subtags = secondary_output_tags->getTags();
			for ( core::Size ii = 0; ii < subtags.size(); ++ii ) {
				utility::tag::TagCOP iitag = subtags[ ii ];
				secondary_outputters_added.insert( iitag->getName() );

				// returns 0 if the secondary outputter is repressed for a particular job
				pose_outputters::SecondaryPoseOutputterOP outputter = secondary_outputter_for_job( inner_job, job_options, iitag->getName(), iitag );
				if ( outputter ) {
					secondary_outputters.push_back( outputter );
				}
			}
		}
	}

	if ( cl_outputters_.empty() ) {
		cl_outputters_ = pose_outputters::PoseOutputterFactory::get_instance()->secondary_pose_outputters_from_command_line();
	}
	for ( std::list< pose_outputters::SecondaryPoseOutputterOP >::const_iterator
			iter = cl_outputters_.begin(), iter_end = cl_outputters_.end();
			iter != iter_end; ++iter ) {

		if ( representative_secondary_outputter_map_.count( (*iter)->class_key() ) == 0 ) {
			representative_secondary_outputter_map_[ (*iter)->class_key() ] = *iter;
		}
		// Do not add a command-line specified secondary outputter if that one
		// has already been specified within the <SecondaryOutput> tag.
		if ( secondary_outputters_added.count( (*iter)->class_key() ) ) continue;

		// returns 0 if the secondary outputter is repressed for a particular job
		pose_outputters::SecondaryPoseOutputterOP outputter = secondary_outputter_for_job(
			inner_job, job_options, (*iter)->class_key(), nullptr );
		if ( outputter ) {
			secondary_outputters.push_back( outputter );
		}
	}

	return secondary_outputters;
}

void
StandardJobQueen::determine_preliminary_job_list_from_xml_file(
	std::string const & job_def_string,
	std::string const & job_def_schema
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::tag;
	using namespace pose_inputters;

	preliminary_larval_jobs_determined_ = true;

	load_job_definition_file( job_def_string, job_def_schema );

	// now iterate across all tags, and for each Job subtag, create a PreliminaryLarvalJob and load it
	// with all of the options that are within the <Option> subtag, if present -- and reading any options
	// not present in the tag from the (global) options system.
	Tag::tags_t const & subtags = job_definition_file_tags_->getTags();
	preliminary_larval_jobs_.reserve( subtags.size() );
	core::Size count_prelim_nodes( 0 );

	bool read_from_command_line = false;
	pose_inputters::PoseInputSourcesAndInputters cl_input_poses_and_inputters;
	for ( auto subtag : subtags ) {
		if ( subtag->getName() != "Job" ) {
			debug_assert( subtag->getName() == "Common" );
			continue;
		}

		TagCOP job_options_tag;
		if ( subtag->hasTag( "Options" ) ) {
			job_options_tag = subtag->getTag( "Options" );
		}
		utility::options::OptionCollectionCOP job_options = options_from_tag( job_options_tag );

		// ok -- look at input tag
		pose_inputters::PoseInputSourcesAndInputters input_poses_and_inputters;
		if ( subtag->hasTag( "Input" ) ) {
			TagCOP input_tag = subtag->getTag( "Input" );
			debug_assert( input_tag ); // XML schema validation should ensure that there is an "Input" subelement
			debug_assert( input_tag->getTags().size() == 1 ); // schema validation should ensure there is exactly one subelement
			TagCOP input_tag_child = input_tag->getTags()[ 0 ];

			// Get the right inputter for this Job.
			pose_inputters::PoseInputterOP inputter;
			if ( use_factory_provided_pose_inputters_ ) {
				inputter = pose_inputters::PoseInputterFactory::get_instance()->new_pose_inputter( input_tag_child->getName() );
			} else {
				runtime_assert( inputter_creators_.count( input_tag_child->getName() ) != 0 );
				inputter = inputter_creators_[ input_tag_child->getName() ]->create_inputter();
			}

			PoseInputSources input_poses = inputter->pose_input_sources_from_tag( *job_options, input_tag_child );
			input_poses_and_inputters.reserve( input_poses.size() );
			for ( auto input_source : input_poses ) {
				input_source->source_id( ++input_pose_counter_ );
				TR << "Assigning input_source " << input_source->input_tag() << " pose_id " << input_pose_counter_ << " and stored as " << input_source->source_id() << std::endl;
				input_poses_and_inputters.push_back( { input_source, inputter } );
			}
		} else if ( ! read_from_command_line || ! load_starting_poses_only_once_ ) {
			// We don't have an input tag -- instead, take the input Poses from the command line.
			input_poses_and_inputters = read_input_poses_from_command_line( );
			cl_input_poses_and_inputters = input_poses_and_inputters;
			read_from_command_line = true;
		} else {
			input_poses_and_inputters = cl_input_poses_and_inputters;
		}

		// Get the right outputter for this job.
		pose_outputters::PoseOutputterOP outputter = get_outputter_from_job_tag( subtag );

		if ( representative_pose_outputter_map_.count( outputter->class_key() ) == 0 ) {
			representative_pose_outputter_map_[ outputter->class_key() ] = outputter;
		}

		// now iterate across the input sources for this job and create
		// a preliminary job for each
		core::Size nstruct = nstruct_for_job( subtag );
		for ( auto pose_source_and_inputter : input_poses_and_inputters ) {
			PreliminaryLarvalJob prelim_job;
			InnerLarvalJobOP inner_job = make_shared< InnerLarvalJob >( nstruct, ++count_prelim_nodes );
			inner_job->input_source( pose_source_and_inputter.first );
			inner_job->jobdef_tag( subtag );
			inner_job->outputter( outputter->class_key() );

			prelim_job.inner_job = inner_job;
			prelim_job.job_tag = subtag;
			prelim_job.job_options = job_options;
			prelim_job.pose_inputter = pose_source_and_inputter.second;
			preliminary_larval_jobs_.push_back( prelim_job );
		}
	}

	///Initialize Preliminary Job Counters.
	prelim_job_tracker_->initialize_tracker(preliminary_larval_jobs_);
}


void
StandardJobQueen::load_job_definition_file(
	std::string const & job_def_string,
	std::string const & job_def_schema
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::tag;

	TR << "Loading job definition file" << std::endl;

	XMLValidationOutput validator_output;
	try {
		validator_output = validate_xml_against_xsd( job_def_string, job_def_schema );
	} catch ( utility::excn::Exception const & e ) {
		std::ostringstream oss;
		if ( option[ in::file::job_definition_file ].user() ) {
			oss << "Job definition file \"" << option[ in::file::job_definition_file ] << "\" failed to validate against"
				" the schema for this application\nUse the option -jd3::job_definition_schema <output filename> to output"
				" the schema to a file.\n" << e.msg() << "\n";
		} else {
			oss << "Job definition file  failed to validate against"
				" the schema for this application\nUse the option -jd3::job_definition_schema <output filename> to output"
				" the schema to a file.\n" << e.msg() << "\n";
		}
		throw CREATE_EXCEPTION(utility::excn::Exception,  oss.str() );
	}

	if ( ! validator_output.valid() ) {
		std::ostringstream oss;
		if ( option[ in::file::job_definition_file ].user() ) {
			oss << "Job definition file \"" << option[ in::file::job_definition_file ] << "\" failed to validate against"
				" the schema for this application\nUse the option -jd3::job_definition_schema <output filename> to output"
				" the schema to a file.\n";
		} else {
			oss << "Job definition file failed to validate against"
				" the schema for this application\nUse the option -jd3::job_definition_schema <output filename> to output"
				" the schema to a file.\n";
		}
		oss << "Error messages were: " << validator_output.error_messages() << "\n";
		oss << "Warning messages were: " << validator_output.warning_messages() << "\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  oss.str() );
	}


	job_definition_file_tags_ = Tag::create( job_def_string );

	// look for the Common block, if there is one
	Tag::tags_t const & subtags = job_definition_file_tags_->getTags();
	for ( auto subtag : subtags ) {
		if ( subtag->getName() == "Common" ) {
			common_block_tags_ = subtag;
			return;
		}
	}
}


void
StandardJobQueen::determine_preliminary_job_list_from_command_line()
{
	using namespace utility::tag;
	using namespace pose_inputters;

	// read from the command line a list of all of the input jobs
	PoseInputSourcesAndInputters input_poses_and_inputters =
		read_input_poses_from_command_line();

	pose_outputters::PoseOutputterOP outputter = get_outputter_from_job_tag( utility::tag::TagCOP() );

	if ( representative_pose_outputter_map_.count( outputter->class_key() ) == 0 ) {
		representative_pose_outputter_map_[ outputter->class_key() ] = outputter;
	}

	// pass in a null-pointing TagCOP and construct the job options object from the command line.
	utility::options::OptionCollectionCOP job_options = options_from_tag( utility::tag::TagCOP() );

	// now iterate across the input sources for this job and create
	// a preliminary job for each
	preliminary_larval_jobs_.reserve( input_poses_and_inputters.size() );
	core::Size count_prelim_nodes( 0 );
	for ( auto const & pose_source_and_inputter : input_poses_and_inputters ) {
		PreliminaryLarvalJob prelim_job;
		Size nstruct = nstruct_for_job( nullptr );
		InnerLarvalJobOP inner_job = make_shared< InnerLarvalJob >( nstruct, ++count_prelim_nodes );
		inner_job->input_source( pose_source_and_inputter.first );
		inner_job->outputter( outputter->class_key() );

		prelim_job.inner_job = inner_job;
		prelim_job.job_tag = TagCOP(); // null ptr
		prelim_job.job_options = job_options;
		prelim_job.pose_inputter = pose_source_and_inputter.second;
		preliminary_larval_jobs_.push_back( prelim_job );
	}

	///Initialize Preliminary Job Counters.
	prelim_job_tracker_->initialize_tracker(preliminary_larval_jobs_);
}

/// @details Requests that set of input structures that have been indicated on the command
/// line from the set of allowed PoseInputters (either all inputters registered with the
/// PoseInputterFactory, or a subset as indicated by the derived JobQueen). Then each pose is
/// givn a separate pose_id. Note that if this is called multiple times, then each structure
/// that is listed on the command line will be given multiple pose_ids. As a consequence,
/// each structure may be loaded multiple times.
pose_inputters::PoseInputSourcesAndInputters
StandardJobQueen::read_input_poses_from_command_line()
{
	using namespace pose_inputters;
	PoseInputSourcesAndInputters input_poses;
	if ( use_factory_provided_pose_inputters_ ) {
		input_poses = PoseInputterFactory::get_instance()->pose_inputs_from_command_line();
	} else {
		for ( auto inputter_creator : inputter_creator_list_ ) {
			PoseInputterOP inputter = inputter_creator->create_inputter();
			if ( inputter->job_available_on_command_line() ) {
				PoseInputSources iter_sources = inputter->pose_input_sources_from_command_line();
				input_poses.reserve( input_poses.size() + iter_sources.size() );
				for ( core::Size ii = 1; ii <= iter_sources.size(); ++ii ) {
					input_poses.push_back( std::make_pair( iter_sources[ ii ], inputter ) );
				}
			}
		}
	}
	for ( auto input_source : input_poses ) {
		input_source.first->source_id( ++input_pose_counter_ );
		TR << "Assigning input_source " << input_source.first->input_tag() << " pose_id " << input_pose_counter_ << std::endl;
	}
	return input_poses;
}

/// @details Retrieve a particular secondary outputter for a job, pointing to the
/// possibly (likely) shared outputter that several jobs will write their output to.
/// If a representative outputter has not been created for this job (as is sometimes
/// the case), then this function will update the representative_secondary_outputter_map_
/// member variable.
pose_outputters::SecondaryPoseOutputterOP
StandardJobQueen::secondary_outputter_for_job(
	InnerLarvalJob const & inner_job,
	utility::options::OptionCollection const & job_options,
	std::string const & secondary_outputter_type,
	utility::tag::TagCOP outputter_tag
)
{
	pose_outputters::SecondaryPoseOutputterOP representative_outputter;
	auto rep_iter = representative_secondary_outputter_map_.find( secondary_outputter_type );
	if ( rep_iter == representative_secondary_outputter_map_.end() ) {
		representative_outputter = pose_outputters::PoseOutputterFactory::get_instance()->
			new_secondary_outputter( secondary_outputter_type );
		representative_secondary_outputter_map_[ secondary_outputter_type ] =
			representative_outputter;
	} else {
		representative_outputter = rep_iter->second;
	}
	debug_assert( representative_outputter );

	std::string which_outputter = representative_outputter->outputter_for_job( outputter_tag, job_options, inner_job );


	if ( which_outputter == "" ) {
		return representative_outputter;
	}
	SecondaryOutputterMap::const_iterator iter =
		secondary_outputter_map_.find( std::make_pair( secondary_outputter_type, which_outputter ) );
	if ( iter != secondary_outputter_map_.end() ) {
		return iter->second;
	}

	pose_outputters::SecondaryPoseOutputterOP outputter =
		pose_outputters::PoseOutputterFactory::get_instance()->new_secondary_outputter( secondary_outputter_type );
	secondary_outputter_map_[ std::make_pair( secondary_outputter_type, which_outputter ) ] = outputter;
	return outputter;

}

pose_outputters::SecondaryPoseOutputterOP
StandardJobQueen::secondary_outputter_for_job( pose_outputters::PoseOutputSpecification const & spec )
{
	std::string const & outputter_type = spec.outputter_type();
	pose_outputters::SecondaryPoseOutputterOP representative;
	auto rep_iter = representative_secondary_outputter_map_.find( outputter_type );
	if ( rep_iter == representative_secondary_outputter_map_.end() ) {
		representative = pose_outputters::PoseOutputterFactory::get_instance()->new_secondary_outputter( outputter_type );
		representative_secondary_outputter_map_[ outputter_type ] = representative;
	} else {
		representative = rep_iter->second;
	}

	// This outputter name may have a job-distributor-assigned suffix, and so it may not yet
	// have been created, even if this JobQueen has been doing some amount of output.
	std::string actual_outputter_name = representative->outputter_for_job( spec );
	if ( actual_outputter_name == "" ) return representative;

	std::pair< std::string, std::string > outputter_key =
		std::make_pair( outputter_type, actual_outputter_name );
	auto iter = secondary_outputter_map_.find( outputter_key );
	if ( iter != secondary_outputter_map_.end() ) {
		return iter->second;
	}

	pose_outputters::SecondaryPoseOutputterOP outputter;
	outputter = pose_outputters::PoseOutputterFactory::get_instance()->new_secondary_outputter( outputter_type );
	secondary_outputter_map_[ outputter_key ] = outputter;
	return outputter;
}

utility::options::OptionKeyList
StandardJobQueen::concatenate_all_options() const
{
	utility::options::OptionKeyList all_options( options_ );
	all_options.insert( all_options.end(), inputter_options_.begin(), inputter_options_.end() );
	all_options.insert( all_options.end(), outputter_options_.begin(), outputter_options_.end() );
	all_options.insert( all_options.end(), secondary_outputter_options_.begin(), secondary_outputter_options_.end() );
	return all_options;
}

pose_outputters::PoseOutputterOP
StandardJobQueen::get_outputter_from_job_tag( utility::tag::TagCOP tag ) const
{
	using utility::tag::TagCOP;
	pose_outputters::PoseOutputterOP outputter;
	if ( tag && tag->hasTag( "Output" ) ) {
		TagCOP output_tag = tag->getTag( "Output" );
		debug_assert( output_tag );
		debug_assert( output_tag->getTags().size() == 1 );
		TagCOP output_subtag = output_tag->getTags()[ 0 ];
		if ( use_factory_provided_pose_outputters_ ) {
			outputter = pose_outputters::PoseOutputterFactory::get_instance()->new_pose_outputter( output_subtag->getName() );
		} else {
			runtime_assert( outputter_creators_.count( output_subtag->getName() ) != 0 );
			auto iter = outputter_creators_.find( output_subtag->getName() );
			outputter = iter->second->create_outputter();
		}
	} else {
		if ( use_factory_provided_pose_outputters_ ) {
			outputter = pose_outputters::PoseOutputterFactory::get_instance()->pose_outputter_from_command_line();
		} else {
			for ( auto outputter_creator : outputter_creator_list_ ) {
				if ( outputter_creator->outputter_specified_by_command_line() ) {
					outputter = outputter_creator->create_outputter();
					break;
				}
			}
			if ( ! outputter ) {
				if ( override_default_outputter_ ) {
					outputter = default_outputter_creator_->create_outputter();
				} else {
					runtime_assert( outputter_creators_.count( pose_outputters::PDBPoseOutputter::keyname() ) );
					outputter = make_shared< pose_outputters::PDBPoseOutputter >();
				}
			}
		}
	}
	return outputter;
}

// For the first node in the JobDAG, the %StandardJobQueen will spool out LarvalJobs
// slowly to the JobDistributor (in increments of the max_njobs parameter in the call
// to job_dag_node_index). Since max_njobs may be smaller than the nstruct parameter,
// the %StandardJobQueen will need to be able to interrupt the spooling of jobs until
// the JobDistributor is ready for them. For this reason, it keeps what is effectively
// a set of indices into a while loop for LarvalJob construction.
bool
StandardJobQueen::get_preliminary_larval_jobs_determined() const {
	return preliminary_larval_jobs_determined_;
}

utility::vector1< PreliminaryLarvalJob > const &
StandardJobQueen::get_preliminary_larval_jobs() const {
	return preliminary_larval_jobs_;
}

core::Size
StandardJobQueen::current_job_index() const {
	return get_job_tracker().current_job_index();
}

} // namespace standard
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::standard::StandardJobQueen::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobQueen >( this ) );
	// EXEMPT options_ inputter_options_ outputter_options_ secondary_outputter_options_
	//arc( CEREAL_NVP( options_ ) ); // utility::options::OptionKeyList
	//arc( CEREAL_NVP( inputter_options_ ) ); // utility::options::OptionKeyList
	//arc( CEREAL_NVP( outputter_options_ ) ); // utility::options::OptionKeyList
	//arc( CEREAL_NVP( secondary_outputter_options_ ) ); // utility::options::OptionKeyList
	arc( CEREAL_NVP( use_factory_provided_pose_inputters_ ) ); // _Bool
	arc( CEREAL_NVP( inputter_creators_ ) ); // std::map<std::string, pose_inputters::PoseInputterCreatorCOP>
	arc( CEREAL_NVP( inputter_creator_list_ ) ); // std::list<pose_inputters::PoseInputterCreatorCOP>
	arc( CEREAL_NVP( use_factory_provided_pose_outputters_ ) ); // _Bool
	arc( CEREAL_NVP( outputter_creators_ ) ); // std::map<std::string, pose_outputters::PoseOutputterCreatorCOP>
	arc( CEREAL_NVP( outputter_creator_list_ ) ); // std::list<pose_outputters::PoseOutputterCreatorCOP>
	arc( CEREAL_NVP( override_default_outputter_ ) ); // _Bool
	arc( CEREAL_NVP( default_outputter_creator_ ) ); // pose_outputters::PoseOutputterCreatorOP
	arc( CEREAL_NVP( pose_outputters_ ) ); // std::map<std::string, pose_outputters::PoseOutputterOP>
	arc( CEREAL_NVP( common_block_precedes_job_blocks_ ) ); // _Bool
	arc( CEREAL_NVP( required_initialization_performed_ ) ); // _Bool
	arc( CEREAL_NVP( job_definition_file_tags_ ) ); // utility::tag::TagCOP
	arc( CEREAL_NVP( common_block_tags_ ) ); // utility::tag::TagCOP
	arc( CEREAL_NVP( preliminary_larval_jobs_determined_ ) ); // _Bool
	arc( CEREAL_NVP( preliminary_larval_jobs_ ) ); // utility::vector1<PreliminaryLarvalJob>
	arc( CEREAL_NVP( prelim_job_tracker_ ) ); // PreliminaryLarvalJobTrackerOP
	arc( CEREAL_NVP( inner_larval_jobs_for_curr_prelim_job_ ) ); // InnerLarvalJobs
	arc( CEREAL_NVP( curr_inner_larval_job_index_ ) ); // Size
	arc( CEREAL_NVP( njobs_made_for_curr_inner_larval_job_ ) ); // Size
	arc( CEREAL_NVP( recent_successes_ ) ); // std::list<output::OutputSpecificationOP>
	arc( CEREAL_NVP( representative_pose_outputter_map_ ) ); // RepresentativeOutputterMap
	arc( CEREAL_NVP( representative_secondary_outputter_map_ ) ); // SecondaryRepresentativeOutputterMap
	arc( CEREAL_NVP( pose_outputter_map_ ) ); // PoseOutputterMap
	arc( CEREAL_NVP( secondary_outputter_map_ ) ); // SecondaryOutputterMap
	arc( CEREAL_NVP( cl_outputters_ ) ); // std::list<pose_outputters::SecondaryPoseOutputterOP>
	arc( CEREAL_NVP( input_pose_counter_ ) ); // core::Size
	arc( CEREAL_NVP( input_pose_index_map_ ) ); // std::map<core::Size, core::pose::PoseOP>
	arc( CEREAL_NVP( load_starting_poses_only_once_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::standard::StandardJobQueen::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobQueen >( this ) );
	// EXEMPT options_ inputter_options_ outputter_options_ secondary_outputter_options_
	//arc( options_ ); // utility::options::OptionKeyList
	//arc( inputter_options_ ); // utility::options::OptionKeyList
	//arc( outputter_options_ ); // utility::options::OptionKeyList
	//arc( secondary_outputter_options_ ); // utility::options::OptionKeyList
	arc( use_factory_provided_pose_inputters_ ); // _Bool
	arc( inputter_creators_ ); // std::map<std::string, pose_inputters::PoseInputterCreatorCOP>

	std::list< std::shared_ptr< protocols::jd3::pose_inputters::PoseInputterCreator > > local_inputter_creator_list;
	arc( local_inputter_creator_list ); // std::list<pose_inputters::PoseInputterCreatorCOP>
	//inputter_creator_list_ = local_inputter_creator_list; // copy the non-const pointer(s) into the const pointer(s) - THIS FAILS
	//Slow loop for now:
	for ( std::shared_ptr< protocols::jd3::pose_inputters::PoseInputterCreator > const & element : local_inputter_creator_list ) {
		inputter_creator_list_.push_back( element );
	}


	arc( use_factory_provided_pose_outputters_ ); // _Bool
	arc( outputter_creators_ ); // std::map<std::string, pose_outputters::PoseOutputterCreatorCOP>

	std::list< std::shared_ptr< protocols::jd3::pose_outputters::PoseOutputterCreator > > local_outputter_creator_list;
	arc( local_outputter_creator_list ); // std::list<pose_outputters::PoseOutputterCreatorCOP>
	//outputter_creator_list_ = local_outputter_creator_list; // copy the non-const pointer(s) into the const pointer(s)
	for ( std::shared_ptr< protocols::jd3::pose_outputters::PoseOutputterCreator > const & element : local_outputter_creator_list ) {
		outputter_creator_list_.push_back( element );
	}

	arc( override_default_outputter_ ); // _Bool
	arc( default_outputter_creator_ ); // pose_outputters::PoseOutputterCreatorOP
	arc( pose_outputters_ ); // std::map<std::string, pose_outputters::PoseOutputterOP>
	arc( common_block_precedes_job_blocks_ ); // _Bool
	arc( required_initialization_performed_ ); // _Bool
	std::shared_ptr< utility::tag::Tag > local_job_definition_file_tags;
	arc( local_job_definition_file_tags ); // utility::tag::TagCOP
	job_definition_file_tags_ = local_job_definition_file_tags; // copy the non-const pointer(s) into the const pointer(s)
	std::shared_ptr< utility::tag::Tag > local_common_block_tags;
	arc( local_common_block_tags ); // utility::tag::TagCOP
	common_block_tags_ = local_common_block_tags; // copy the non-const pointer(s) into the const pointer(s)
	arc( preliminary_larval_jobs_determined_ ); // _Bool
	arc( preliminary_larval_jobs_ ); // utility::vector1<PreliminaryLarvalJob>
	arc( prelim_job_tracker_ ); // PreliminaryLarvalJobTrackerOP
	arc( inner_larval_jobs_for_curr_prelim_job_ ); // InnerLarvalJobs
	arc( curr_inner_larval_job_index_ ); // Size
	arc( njobs_made_for_curr_inner_larval_job_ ); // Size
	arc( recent_successes_ ); // std::list<output::OutputSpecificationOP>
	arc( representative_pose_outputter_map_ ); // RepresentativeOutputterMap
	arc( representative_secondary_outputter_map_ ); // SecondaryRepresentativeOutputterMap
	arc( pose_outputter_map_ ); // PoseOutputterMap
	arc( secondary_outputter_map_ ); // SecondaryOutputterMap
	arc( cl_outputters_ ); // std::list<pose_outputters::SecondaryPoseOutputterOP>
	arc( input_pose_counter_ ); // core::Size
	arc( input_pose_index_map_ ); // std::map<core::Size, core::pose::PoseOP>
	arc( load_starting_poses_only_once_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::standard::StandardJobQueen );
CEREAL_REGISTER_TYPE( protocols::jd3::standard::StandardJobQueen )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_standard_StandardJobQueen )
#endif // SERIALIZATION
