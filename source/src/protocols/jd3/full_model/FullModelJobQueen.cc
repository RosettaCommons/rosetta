// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/full_model/FullModelJobQueen.cc
/// @brief  FullModelJobQueen class's method definitions
/// @author Andy Watkins (amw579@stanford.edu)

//unit headers
#include <protocols/jd3/full_model/FullModelJobQueen.hh>

// package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/full_model/FullModelInnerLarvalJob.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/full_model/MoverAndFullModelJob.hh>
#include <protocols/jd3/InputSource.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputSource.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputter.hh>
#include <protocols/jd3/full_model_inputters/PDBFullModelInputter.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputterCreator.hh>
#include <protocols/jd3/full_model_inputters/PDBFullModelInputterCreator.hh>
#include <protocols/jd3/full_model_inputters/SilentFileFullModelInputterCreator.hh>
#include <protocols/jd3/full_model_inputters/FullModelInputterFactory.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.hh>
#include <protocols/jd3/deallocation/DeallocationMessage.hh>
#include <protocols/jd3/deallocation/InputPoseDeallocationMessage.hh>

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
#include <basic/options/option.cc.gen.hh>
#include <basic/datacache/ConstDataMap.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.FullModelJobQueen" );


namespace protocols {
namespace jd3 {
namespace full_model {

PreliminaryLarvalJob::PreliminaryLarvalJob() {}
PreliminaryLarvalJob::~PreliminaryLarvalJob() = default;
PreliminaryLarvalJob::PreliminaryLarvalJob( PreliminaryLarvalJob const & src ) :
	inner_job( src.inner_job ),
	job_tag( src.job_tag ),
	job_options( src.job_options ),
	full_model_inputter( src.full_model_inputter )
{
}

PreliminaryLarvalJob &
PreliminaryLarvalJob::operator = ( PreliminaryLarvalJob const & rhs )
{
	if ( this != &rhs ) {
		inner_job = rhs.inner_job;
		job_tag   = rhs.job_tag;
		job_options = rhs.job_options;
		full_model_inputter = rhs.full_model_inputter;
	}
	return *this;
}



utility::excn::EXCN_Msg_Exception
bad_inner_job_exception()
{
	return utility::excn::EXCN_Msg_Exception( "The InnerLarvalJob provided to the FullModelJobQueen by the DerivedJobQueen must derive from FullModelInnerLarvalJob. Cannot proceed." );
}


FullModelJobQueen::FullModelJobQueen() :
	use_factory_provided_full_model_inputters_( true ),
	use_factory_provided_pose_outputters_( true ),
	override_default_outputter_( false ),
	required_initialization_performed_( false ),
	larval_job_counter_( 0 ),
	preliminary_larval_jobs_determined_( false ),
	curr_inner_larval_job_index_( 0 ),
	njobs_made_for_curr_inner_larval_job_( 0 ),
	input_pose_counter_( 0 )
{
	do_not_accept_all_full_model_inputters_from_factory();
	full_model_inputters::FullModelInputterCreatorOP full_model_inputter_creator( new full_model_inputters::PDBFullModelInputterCreator );
	full_model_inputters::FullModelInputterCreatorOP full_model_inputter_creator2( new full_model_inputters::SilentFileFullModelInputterCreator );
	allow_full_model_inputter( full_model_inputter_creator );
	allow_full_model_inputter( full_model_inputter_creator2 );

	// begin to populate the per-job options object
	full_model_inputters::FullModelInputterFactory::get_instance()->list_options_read( inputter_options_ );
	pose_outputters::PoseOutputterFactory::get_instance()->list_outputter_options_read( outputter_options_ );
	pose_outputters::PoseOutputterFactory::get_instance()->list_secondary_outputter_options_read( secondary_outputter_options_ );

	// add options?
	// In theory we'd use some function provided by the stepwise options classes here or something...
	add_option( basic::options::OptionKeys::in::file::fasta );
	/*utility::options::OptionKeyList opts;
	core::scoring::list_read_options_in_get_score_function( opts );
	core::pack::task::PackerTask::list_options_read( opts );
	core::pack::task::operation::ReadResfile::list_options_read( opts );
	protocols::simple_moves::PackRotamersMover::list_options_read( opts );
	add_options( opts );
	add_option( basic::options::OptionKeys::minimize_sidechains );
	add_option( basic::options::OptionKeys::min_pack );
	add_option( basic::options::OptionKeys::off_rotamer_pack );*/
}

FullModelJobQueen::~FullModelJobQueen() = default;

utility::options::OptionTypes
option_type_from_key(
	utility::options::OptionKey const & key
)
{
	using namespace utility::options;
	if ( dynamic_cast< StringOptionKey const * > ( &key ) ) {
		return STRING_OPTION;
	} else if ( dynamic_cast< StringVectorOptionKey const * > ( &key ) ) {
		return STRING_VECTOR_OPTION;
	} else if ( dynamic_cast< PathOptionKey const * > ( &key ) ) {
		return PATH_OPTION;
	} else if ( dynamic_cast< PathVectorOptionKey const * > ( &key ) ) {
		return PATH_VECTOR_OPTION;
	} else if ( dynamic_cast< FileOptionKey const * > ( &key ) ) {
		return FILE_OPTION;
	} else if ( dynamic_cast< FileVectorOptionKey const * > ( &key ) ) {
		return FILE_VECTOR_OPTION;
	} else if ( dynamic_cast< RealOptionKey const * > ( &key ) ) {
		return REAL_OPTION;
	} else if ( dynamic_cast< RealVectorOptionKey const * > ( &key ) ) {
		return REAL_VECTOR_OPTION;
	} else if ( dynamic_cast< IntegerOptionKey const * > ( &key ) ) {
		return INTEGER_OPTION;
	} else if ( dynamic_cast< IntegerVectorOptionKey const * > ( &key ) ) {
		return INTEGER_VECTOR_OPTION;
	} else if ( dynamic_cast< BooleanOptionKey const * > ( &key ) ) {
		return BOOLEAN_OPTION;
	} else if ( dynamic_cast< BooleanVectorOptionKey const * > ( &key ) ) {
		return BOOLEAN_VECTOR_OPTION;
	} else if ( dynamic_cast< ResidueChainVectorOptionKey const * > ( &key ) ) {
		return RESIDUE_CHAIN_VECTOR_OPTION;
	}
	return UNKNOWN_OPTION;
}

utility::tag::XMLSchemaType
value_attribute_type_for_option(
	utility::options::OptionTypes const & key
)
{
	using namespace utility::options;
	using namespace utility::tag;
	switch ( key ) {
	case STRING_OPTION :
	case STRING_VECTOR_OPTION :
	case PATH_OPTION :
	case PATH_VECTOR_OPTION :
	case FILE_OPTION :
	case FILE_VECTOR_OPTION :
	case RESIDUE_CHAIN_VECTOR_OPTION :
		return xs_string;
	case REAL_OPTION :
		return xsct_real;
	case REAL_VECTOR_OPTION :
		return xsct_real_wsslist;
	case INTEGER_OPTION :
		return xs_integer;
	case INTEGER_VECTOR_OPTION :
		return xsct_int_wsslist;
	case BOOLEAN_OPTION :
		return xsct_rosetta_bool;
	case BOOLEAN_VECTOR_OPTION :
		return xsct_bool_wsslist; // note: double check that options system uses utility/string_funcs.hh to cast from strings to bools.
	default :
		throw utility::excn::EXCN_Msg_Exception( "Unsupported option type hit in FullModelJobQueen.cc's value_attribute" );
	}
	return "ERROR";
}

std::string
FullModelJobQueen::job_definition_xsd() const
{
	using namespace utility::tag;
	using namespace utility::options;
	using namespace basic::options;

	utility::options::OptionKeyList all_options = concatenate_all_options();

	XMLSchemaDefinition xsd;

	if ( use_factory_provided_full_model_inputters_ ) {
		full_model_inputters::FullModelInputterFactory::get_instance()->define_full_model_inputter_xml_schema( xsd );
	} else {
		utility::tag::define_xml_schema_group(
			inputter_creators_,
			full_model_inputters::FullModelInputterFactory::full_model_inputter_xml_schema_group_name(),
			& full_model_inputters::FullModelInputterFactory::complex_type_name_for_full_model_inputter,
			xsd );
	}

	XMLSchemaSimpleSubelementList input_subelements;
	input_subelements.add_group_subelement( & full_model_inputters::FullModelInputterFactory::full_model_inputter_xml_schema_group_name );
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
		.description( "XRW TO DO" )
		.complex_type_naming_func( & job_def_complex_type_name )
		.set_subelements_pick_one( output_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList secondary_output_subelements;
	secondary_output_subelements.add_group_subelement( & pose_outputters::PoseOutputterFactory::secondary_pose_outputter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator secondary_output_ct;
	secondary_output_ct
		.element_name( "SecondaryOutput" )
		.description( "XRW TO DO" )
		.complex_type_naming_func( & job_def_complex_type_name )
		.set_subelements_single_appearance_optional( secondary_output_subelements )
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
	job_ct.add_ordered_subelement_set_as_required( job_input_subelement );

	// Job <Output> element is optional
	XMLSchemaSimpleSubelementList job_output_subelement;
	job_output_subelement.add_already_defined_subelement( "Output", & job_def_complex_type_name );
	job_ct.add_ordered_subelement_set_as_optional( job_output_subelement );

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
		throw utility::excn::EXCN_Msg_Exception( "Subclass of FullModelJobQueen's append_job_tag_subelements"
			" method invokes a method of the XMLSchemaComplexTypeGenerator that overwrote the <Input>, <Output>, and"
			" <Options> elements.  It should only call methods named \"add_ordered_subelement_set_*\"" );
	}

	job_ct
		.element_name( "Job" )
		.description( "XRW TO DO" )
		.complex_type_naming_func( & job_def_complex_type_name )
		.add_attribute( XMLSchemaAttribute::attribute_w_default(  "nstruct", xsct_non_negative_integer, "XRW TO DO",  "1"  ))
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
		throw utility::excn::EXCN_Msg_Exception( "Subclass of FullModelJobQueen's append_job_tag_subelements"
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
	job_def_subelements.add_already_defined_subelement( "Common", & job_def_complex_type_name, 0, 1 );
	job_def_subelements.add_already_defined_subelement( "Job", & job_def_complex_type_name, 1, xsminmax_unbounded );
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
	} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
		std::ostringstream oss;
		oss << "The XML Schema for the job definition file is invalid.  The error message is:\n" << e.msg()
			<< "\nAnd the whole schema is:\n" << xsd_string << "\nThis executable cannot be used in its"
			<< " current state.\n";
		throw utility::excn::EXCN_Msg_Exception( oss.str() );
	}

	return xsd_string;

}

std::string
FullModelJobQueen::resource_definition_xsd() const
{
	// TO DO!
	return "";
}

/// @details The base class provides a digraph with a single node -- that is, all the jobs
/// are independent of each other.  This is equivalent to the kind of jobs that could be
/// run in JD2.
JobDigraphOP
FullModelJobQueen::initial_job_dag()
{
	determine_preliminary_job_list();

	core::Size n_pjns = preliminary_larval_jobs_.size();
	preliminary_job_node_inds_.resize( n_pjns );
	preliminary_job_nodes_complete_.resize( n_pjns );
	pjn_job_ind_begin_.resize( n_pjns );
	pjn_job_ind_end_.resize( n_pjns );

	for ( core::Size ii = 1; ii <= n_pjns; ++ii ) {
		preliminary_job_node_inds_[ ii ] = ii;
		pjn_job_ind_begin_[ ii ] = pjn_job_ind_end_[ ii ] = 0;
		preliminary_job_nodes_complete_[ ii ] = 0;
	}
	// create a DAG with as many nodes in it as there are preliminary larval jobs
	job_graph_ = JobDigraphOP( new JobDigraph( preliminary_larval_jobs_.size() ) );
	return job_graph_;
}

void FullModelJobQueen::update_job_dag( JobDigraphUpdater & ) {}

/// @details The process begins by first constructing the job definition and resource definition
/// XSDs.  With these schemas, the %FullModelJobQueen validates the input XML files (if present).
/// The %FullModelJobQueen then populates preliminary versions of LarvalJob objects./ If the XSD
/// includes "command line options" (which may be specified either from the command line or in
/// the <options> sections of the Job XML file), the %FullModelJobQueen loads the preliminary
/// LarvalJob objects with the options. These preliminary LarvalJob objects will not have been
/// nstruct expanded (i.e. if there are 100 nstruct for each of 5 different jobs, then there will
/// only be 5 preliminary larval jobs created). It then passes the preliminary LarvalJob list and
/// the TagOP objects for each preliminary LarvalJob to the derived class through the
/// refine_job_list method.
LarvalJobs
FullModelJobQueen::determine_job_list( Size job_dag_node_index, Size max_njobs )
{
	// ok -- we're going to look for a job definition file, and failing that, fall back on
	// the FullModelInputterFactory to determine where the input sources are coming from.
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

	for ( auto job : larval_jobs ) {
		runtime_assert( job );
		core::Size pose_id = job->inner_job()->input_source().pose_id();
		if ( !last_job_for_input_pose_.count( pose_id ) ) {
			last_job_for_input_pose_[ pose_id ] = job->job_index();
		} else if ( last_job_for_input_pose_[ pose_id ] < job->job_index() ) {
			last_job_for_input_pose_[ pose_id ] = job->job_index();
		}
	}

	return larval_jobs;
}

bool
FullModelJobQueen::has_job_completed( protocols::jd3::LarvalJobCOP job )
{
	FullModelInnerLarvalJobCOP inner_job = utility::pointer::dynamic_pointer_cast< FullModelInnerLarvalJob const > ( job->inner_job() );
	if ( ! inner_job ) { throw bad_inner_job_exception(); }

	utility::options::OptionCollectionOP job_options = options_for_job( *inner_job );
	runtime_assert( job_options );
	return pose_outputter_for_job( *inner_job )->job_has_already_completed( *job, *job_options );
}

void
FullModelJobQueen::mark_job_as_having_begun( protocols::jd3::LarvalJobCOP job )
{
	FullModelInnerLarvalJobCOP inner_job = utility::pointer::dynamic_pointer_cast< FullModelInnerLarvalJob const > ( job->inner_job() );
	if ( ! inner_job ) { throw bad_inner_job_exception(); }
	utility::options::OptionCollectionOP job_options = options_for_job( *inner_job );
	return pose_outputter_for_job( *inner_job )->mark_job_as_having_started( *job, *job_options );
}

JobOP
FullModelJobQueen::mature_larval_job(
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
	// and the FullModelInputters require the Input TagCOP that corresponds to the preliminary job node
	// that this job came from.
	if ( ! required_initialization_performed_ ) {
		determine_preliminary_job_list();
	}

	// initialize the options collection for this job.
	FullModelInnerLarvalJobCOP inner_job = utility::pointer::dynamic_pointer_cast< FullModelInnerLarvalJob const > ( larval_job->inner_job() );
	if ( ! inner_job ) { throw bad_inner_job_exception(); }
	utility::options::OptionCollectionCOP job_options = options_for_job( *inner_job );

	return complete_larval_job_maturation( larval_job, job_options, input_results );
}

bool FullModelJobQueen::larval_job_needed_for_note_job_completed() const { return false; }

void FullModelJobQueen::note_job_completed( LarvalJobCOP job, JobStatus status, Size nresults)
{
	// pass the status to the other note_job_completed function
	note_job_completed( job->job_index(), status, nresults );
}

void FullModelJobQueen::note_job_completed( core::Size job_id, JobStatus status, Size nresults )
{
	completed_jobs_.insert( job_id );
	if ( status == jd3_job_status_success ) {
		for ( Size ii = 1; ii <= nresults; ++ii ) {
			recent_successes_.push_back( std::make_pair( job_id, ii ) );
		}
		successful_jobs_.insert( job_id );
		PartialOutputStatus output_status;
		output_status.n_results = nresults;
		output_status.n_output_or_discarded  = 0;
		results_processed_for_job_[ job_id ] = output_status;
	} else {
		failed_jobs_.insert( job_id );
	}
}

bool FullModelJobQueen::larval_job_needed_for_completed_job_summary() const { return false; }

void FullModelJobQueen::completed_job_summary( LarvalJobCOP, Size, JobSummaryOP ) {}

void FullModelJobQueen::completed_job_summary( core::Size, Size, JobSummaryOP ) {}

std::list< JobResultID > FullModelJobQueen::jobs_that_should_be_output()
{
	std::list< JobResultID > return_list = recent_successes_;
	recent_successes_.clear();
	return return_list;
}

/// @details Default implementation does not discard any job results.
std::list< JobResultID >
FullModelJobQueen::job_results_that_should_be_discarded() {
	return std::list< JobResultID >();
}

void FullModelJobQueen::completed_job_result( LarvalJobCOP job, core::Size result_index, JobResultOP job_result )
{
	FullModelInnerLarvalJobCOP inner_job = utility::pointer::dynamic_pointer_cast< FullModelInnerLarvalJob const > ( job->inner_job() );
	if ( ! inner_job ) { throw bad_inner_job_exception(); }
	pose_outputters::PoseOutputterOP outputter = pose_outputter_for_job( *inner_job );
	FullModelJobResultOP pose_result = utility::pointer::dynamic_pointer_cast< FullModelJobResult >( job_result );
	if ( ! pose_result ) {
		utility::excn::EXCN_Msg_Exception( "JobResult handed to StandardJobQueen::completed_job_result is not a PoseJobResult or derived from PoseJobResult" );
	}
	utility::options::OptionCollectionOP job_options = options_for_job( *inner_job );
	utility::tag::TagCOP outputter_tag;
	if ( job->inner_job()->jobdef_tag() ) {
		utility::tag::TagCOP tag = job->inner_job()->jobdef_tag();
		if ( tag->hasTag( "Output" ) ) {
			outputter_tag = tag->getTag( "Output" )->getTags()[ 0 ];
		}
	}
	std::pair< core::Size, core::Size > pose_of_total =
		{ result_index, results_processed_for_job_[ job->job_index() ].n_results };

	outputter->write_output_pose( *job, pose_of_total, *job_options, outputter_tag, *pose_result->pose() );

	std::list< pose_outputters::SecondaryPoseOutputterOP > secondary_outputters = secondary_outputters_for_job( *inner_job, *job_options );
	for ( auto const & secondary_outputter : secondary_outputters ) {
		secondary_outputter->write_output_pose( *job, *job_options, *pose_result->pose() );
	}

	note_job_result_output_or_discarded( job, result_index );
}

/// @details Construct the XSD and then invoke the (private) determine_preliminary_job_list_from_xml_file method,
/// which is also invoked by determine_preliminary_job_list.
void
FullModelJobQueen::determine_preliminary_job_list_from_xml_file( std::string const & job_def_string )
{
	std::string job_def_schema = job_definition_xsd();
	determine_preliminary_job_list_from_xml_file( job_def_string, job_def_schema );
}

void
FullModelJobQueen::flush()
{
	for ( auto const & elem : pose_outputter_map_ ) {
		elem.second->flush();
	}
	for ( auto const & elem : secondary_outputter_map_ ) {
		elem.second->flush();
	}
}

std::list< deallocation::DeallocationMessageOP >
FullModelJobQueen::deallocation_messages()
{
	// We will operate under the assumption that all of the times a Pose will
	// be used will be for sequential jobs, and that it is safe when this function
	// is called to deallocate all of the Poses whose last job index is less than
	// the number of jobs we have handed out so far - but that it would not be
	// safe to deallocate the Pose whose last job index is equal to the number of
	// jobs we've handed out so far because we will perhaps hand out more jobs
	// that use that Pose in just a moment.

	std::list< deallocation::DeallocationMessageOP > messages;
	for ( auto it = last_job_for_input_pose_.begin(); it != last_job_for_input_pose_.end(); /*no inc*/ ) {

		if ( it->second < larval_job_counter_ ) {
			deallocation::DeallocationMessageOP msg( new deallocation::InputPoseDeallocationMessage( it->first ));
			messages.push_back( msg );

			auto pose_iter = input_pose_index_map_.find( it->first );
			if ( pose_iter != input_pose_index_map_.end() ) {
				input_pose_index_map_.erase( pose_iter );
			}
			it = last_job_for_input_pose_.erase( it );
		} else {
			++it;
		}
	}
	return messages;
}

void
FullModelJobQueen::process_deallocation_message(
	deallocation::DeallocationMessageOP message
)
{
	using namespace deallocation;
	if ( message->deallocation_type() == input_pose_deallocation_msg ) {
		InputPoseDeallocationMessageOP pose_dealloc_msg =
			utility::pointer::dynamic_pointer_cast< InputPoseDeallocationMessage > ( message );
		if ( input_pose_index_map_.count( pose_dealloc_msg->pose_id() ) ) {
			input_pose_index_map_.erase( pose_dealloc_msg->pose_id() );
		}
	} else {
		// any other type of deallocation message, for now, will have to be
		// handled by the derived class.
		derived_process_deallocation_message( message );
	}
}


std::string
FullModelJobQueen::job_def_complex_type_name( std::string const & type )
{ return "job_def_" + type + "_type"; }


void
FullModelJobQueen::do_not_accept_all_full_model_inputters_from_factory()
{
	if ( use_factory_provided_full_model_inputters_ ) {
		use_factory_provided_full_model_inputters_ = false;
		inputter_options_.clear();
	}
}

void
FullModelJobQueen::allow_full_model_inputter( full_model_inputters::FullModelInputterCreatorOP creator )
{
	if ( use_factory_provided_full_model_inputters_ ) {
		use_factory_provided_full_model_inputters_ = false;
		// the user has not told us they want to drop the original set of creators
		// so interpret this as them wanting to allow the user to provide another
		// inputter.
		inputter_creator_list_ = full_model_inputters::FullModelInputterFactory::get_instance()->full_model_inputter_creators();
		for ( auto creator : inputter_creator_list_ ) {
			inputter_creators_[ creator->keyname() ] = creator;
		}
	}
	if ( ! inputter_creators_.count( creator->keyname() ) ) {
		inputter_creator_list_.push_back( creator );
		inputter_creators_[ creator->keyname() ] = creator;
		creator->list_options_read( inputter_options_ );
	}
}

void
FullModelJobQueen::do_not_accept_all_pose_outputters_from_factory()
{
	if ( use_factory_provided_pose_outputters_ ) {
		use_factory_provided_pose_outputters_ = false;
		override_default_outputter_ = true;
		outputter_options_.clear();
	}
}

void
FullModelJobQueen::allow_pose_outputter( pose_outputters::PoseOutputterCreatorOP creator )
{
	if ( use_factory_provided_pose_outputters_ ) {
		use_factory_provided_pose_outputters_ = false;
		// the user has not told us they want to drop the original set of creators
		// so interpret this as them wanting to allow the user to provide another
		// outputter in addition to the original set.
		outputter_creator_list_ = pose_outputters::PoseOutputterFactory::get_instance()->pose_outputter_creators();
		for ( auto creator : outputter_creator_list_ ) {
			outputter_creators_[ creator->keyname() ] = creator;
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
FullModelJobQueen::set_default_outputter( pose_outputters::PoseOutputterCreatorOP creator )
{
	override_default_outputter_ = true;
	default_outputter_creator_ = creator;
	if ( use_factory_provided_pose_outputters_ ) {
		// the user has not told us they want to drop the original set of creators
		// but to provide the desired functionality, the SJQ has to handle outputter
		// instantiation herself
		use_factory_provided_pose_outputters_ = false;
		outputter_creator_list_ = pose_outputters::PoseOutputterFactory::get_instance()->pose_outputter_creators();
		for ( auto creator : outputter_creator_list_ ) {
			outputter_creators_[ creator->keyname() ] = creator;
		}
	}
}


/// @details Derived classes may choose to not override this method as a way to indicate that
/// they have no additional subtags of the <Job> tag they wish to add. This is a no-op
/// implementation.
void FullModelJobQueen::append_job_tag_subelements(
	utility::tag::XMLSchemaDefinition &,
	utility::tag::XMLSchemaComplexTypeGenerator &
) const
{}

/// @details Derived classes may choose to not override this method as a way to indicate that
/// they have no common-block data that needs to be defined.
void
FullModelJobQueen::append_common_tag_subelements(
	utility::tag::XMLSchemaDefinition &,
	utility::tag::XMLSchemaComplexTypeGenerator &
) const
{}

/// @details In most JobQueens, this function will often generate a bunch of 'inner' jobs.
/// (Perhaps think of these as input states, or modeling problems.) For this one, the command
/// line can ONLY ever specify a single input state -- just because of how the CL is put
/// together -- and EVENTUALLY a job definition file can perhaps do more.
void
FullModelJobQueen::determine_preliminary_job_list()
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

}

/// @brief Read access for the subset of nodes in the job DAG which the %FullModelJobQueen
/// is responsible for producing the larval_jobs. They are called "preliminary" jobs because
/// they do not depend on outputs from any previous node in the graph. (The set of job nodes
/// that contain no incoming edges, though, could perhaps be different from the set of
/// preliminary job nodes, so the %FullModelJobQueen requires that the derived job queen
/// inform her of which nodes are the preliminary job nodes.)
utility::vector1< core::Size > const &
FullModelJobQueen::preliminary_job_nodes() const
{
	return preliminary_job_node_inds_;
}

bool
FullModelJobQueen::all_jobs_assigned_for_preliminary_job_node( core::Size node_id ) const
{
	return preliminary_job_nodes_complete_[ node_id ];
}


/// returns zero if no jobs for this node have yet been created
core::Size FullModelJobQueen::preliminary_job_node_begin_job_index( core::Size node_id ) const
{
	return pjn_job_ind_begin_[ node_id ];
}

/// @brief The index of the last job for a particular preliminary-job node; this function
/// returns zero if there are some jobs for this node that have not yet been created.
core::Size FullModelJobQueen::preliminary_job_node_end_job_index( core::Size node_id ) const
{
	return pjn_job_ind_end_[ node_id ];
}

/// @brief Read access to derived JobQueens to the preliminary job list.
/// This will return an empty list if  determine_preliminary_jobs has not yet
/// been called.
utility::vector1< PreliminaryLarvalJob > const &
FullModelJobQueen::preliminary_larval_jobs() const
{
	return preliminary_larval_jobs_;
}

numeric::DiscreteIntervalEncodingTree< core::Size > const &
FullModelJobQueen::completed_jobs() const
{
	return completed_jobs_;
}

numeric::DiscreteIntervalEncodingTree< core::Size > const &
FullModelJobQueen::successful_jobs() const
{
	return successful_jobs_;
}

numeric::DiscreteIntervalEncodingTree< core::Size > const &
FullModelJobQueen::failed_jobs() const
{
	return failed_jobs_;
}

numeric::DiscreteIntervalEncodingTree< core::Size > const &
FullModelJobQueen::processed_jobs() const {
	return processed_jobs_;
}

/// @details This base class implementation merely returns a one-element list containing the
/// input inner_job.  Derived classes have the flexibility to create several preliminary
/// jobs from this input job
FullModelInnerLarvalJobs
FullModelJobQueen::refine_preliminary_job( PreliminaryLarvalJob const & prelim_job )
{
	FullModelInnerLarvalJobs one_job( 1, prelim_job.inner_job );
	return one_job;
}

LarvalJobs
FullModelJobQueen::expand_job_list( FullModelInnerLarvalJobOP inner_job, core::Size max_larval_jobs_to_create )
{
	core::Size nstruct = inner_job->nstruct_max();
	LarvalJobs jobs;
	core::Size n_to_make = std::min( nstruct, max_larval_jobs_to_create );
	for ( core::Size jj = 1; jj <= n_to_make; ++jj ) {
		LarvalJobOP job = create_larval_job( inner_job, njobs_made_for_curr_inner_larval_job_ + jj, ++larval_job_counter_ );
		jobs.push_back( job );
	}
	return jobs;
}

FullModelInnerLarvalJobOP
FullModelJobQueen::create_inner_larval_job( core::Size nstruct, core::Size prelim_job_node ) const
{
	return FullModelInnerLarvalJobOP( new FullModelInnerLarvalJob( nstruct, prelim_job_node ) );
}

/// @details Factory method instantiates the base-class LarvalJob to start.
/// Non-const so that derived classes may do thing like keep track of each of the
/// larval jobs they instantiate.
LarvalJobOP
FullModelJobQueen::create_larval_job(
	FullModelInnerLarvalJobOP job,
	core::Size nstruct_index,
	core::Size larval_job_index
)
{
	return LarvalJobOP( new LarvalJob( job, nstruct_index, larval_job_index ));
}

JobOP
FullModelJobQueen::create_job( LarvalJobCOP ) const
{
	return JobOP( new MoverAndFullModelJob );
}

LarvalJobs
FullModelJobQueen::next_batch_of_larval_jobs_for_job_node( core::Size, core::Size )
{
	return LarvalJobs(); // default ctor to give an empty list
}

void
FullModelJobQueen::derived_process_deallocation_message(
	deallocation::DeallocationMessageOP
)
{}



void FullModelJobQueen::add_options( utility::options::OptionKeyList const & opts )
{
	using namespace utility::options;
	for ( auto const & opt : opts ) {
		options_.push_back( opt );
	}
}

void FullModelJobQueen::add_option( utility::options::OptionKey const & key )
{
	options_.emplace_back(key );
}

/// @details TO DO
void FullModelJobQueen::remove_default_input_element() {}

utility::tag::TagCOP
FullModelJobQueen::common_block_tags() const
{
	return common_block_tags_;
}

core::pose::PoseOP
FullModelJobQueen::pose_for_job(
	LarvalJobCOP job,
	utility::options::OptionCollection const & options
)
{
	using namespace full_model_inputters;
	// either read the Pose in using the full_model_inputter (and then keep a copy
	// in the resource manager), or retrieve the Pose from the resource manager
	// initial version: just read the pose in for each job.

	FullModelInnerLarvalJobCOP inner_job = utility::pointer::dynamic_pointer_cast< FullModelInnerLarvalJob const > ( job->inner_job() );
	if ( ! inner_job ) { throw bad_inner_job_exception(); }

	FullModelInputSource const & input_source = dynamic_cast< FullModelInputSource const & >( inner_job->input_source() );

	if ( input_pose_index_map_.count( input_source.pose_id() ) ) {
		core::pose::PoseOP pose( new core::pose::Pose );
		// Why does the full_model job queen use detached_copy? Because it is important in multithreaded
		// contexts that no two Poses pointing to that same data end up in two separate threads,
		// and then try to modify that data at the same time.  It turns out there are places
		// in Pose where it covertly modifies data in some other Pose and this could lead to
		// race conditions.
		pose->detached_copy( *input_pose_index_map_[ input_source.pose_id() ] );
		return pose;
	}

	utility::tag::TagCOP inputter_tag;
	if ( inner_job->jobdef_tag() && inner_job->jobdef_tag()->hasTag( "Input" ) ) {
		utility::tag::TagCOP input_tag = inner_job->jobdef_tag()->getTag( "Input" );
		if ( input_tag->getTags().size() != 0 ) {
			inputter_tag = input_tag->getTags()[ 0 ];
		}
	}

	core::pose::PoseOP input_pose = full_model_inputter_for_job( *inner_job )->full_model_from_input_source( input_source, options, inputter_tag );
	input_pose_index_map_[ input_source.pose_id() ] = input_pose;

	core::pose::PoseOP pose( new core::pose::Pose );
	pose->detached_copy( *input_pose );
	return pose;
}

//ResourceManagerOP FullModelJobQueen::resource_manager()
//{}

/// @brief Access the pose inputter
full_model_inputters::FullModelInputterOP
FullModelJobQueen::full_model_inputter_for_job( FullModelInnerLarvalJob const & inner_job ) const
{
	// find the preliminary job node for this job, if available
	// and return the already-created pose inputter
	if ( inner_job.prelim_job_node() == 0 ) {
		if ( use_factory_provided_full_model_inputters_ ) {
			return full_model_inputters::FullModelInputterFactory::get_instance()->new_full_model_inputter( inner_job.input_source().origin() );
		} else {
			runtime_assert( inputter_creators_.count( inner_job.input_source().origin() ) );
			auto iter = inputter_creators_.find( inner_job.input_source().origin() );
			return iter->second->create_inputter();
		}
	} else {
		debug_assert( preliminary_larval_jobs_[ inner_job.prelim_job_node() ].full_model_inputter );
		return preliminary_larval_jobs_[ inner_job.prelim_job_node() ].full_model_inputter;
	}
}

/// @brief Access the pose outputter for a particular job; perhaps this outputter has already been created, or
/// perhaps it needs to be created and stored
pose_outputters::PoseOutputterOP
FullModelJobQueen::pose_outputter_for_job( FullModelInnerLarvalJob const & inner_job )
{
	utility::options::OptionCollectionOP job_options = options_for_job( inner_job );
	return pose_outputter_for_job( inner_job, *job_options );
}

pose_outputters::PoseOutputterOP
FullModelJobQueen::pose_outputter_for_job(
	FullModelInnerLarvalJob const & inner_job,
	utility::options::OptionCollection const & job_options
)
{
	pose_outputters::PoseOutputterOP representative_outputter;
	representative_outputter = representative_pose_outputter_map_[ inner_job.outputter() ];

	utility::tag::TagCOP output_tag;
	if ( inner_job.jobdef_tag() ) {
		utility::tag::TagCOP job_tag = inner_job.jobdef_tag();
		if ( job_tag->hasTag( "Output" ) ) {
			output_tag = job_tag->getTag( "Output" );
		}
	}
	runtime_assert( representative_outputter );
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


std::list< pose_outputters::SecondaryPoseOutputterOP >
FullModelJobQueen::secondary_outputters_for_job(
	FullModelInnerLarvalJob const & inner_job,
	utility::options::OptionCollection const & job_options
)
{
	std::set< std::string > secondary_outputters_added;
	std::list< pose_outputters::SecondaryPoseOutputterOP > secondary_outputters;
	if ( inner_job.jobdef_tag() ) {
		utility::tag::TagCOP job_tag = inner_job.jobdef_tag();
		if ( job_tag->hasTag( "SecondaryOutput" ) ) {
			utility::tag::TagCOP secondary_output_tags = job_tag->getTag( "SecondaryOutput" );
			utility::vector0< utility::tag::TagCOP > const & subtags = secondary_output_tags->getTags();
			for ( core::Size ii = 0; ii < subtags.size(); ++ii ) {
				utility::tag::TagCOP iitag = subtags[ ii ];
				if ( secondary_outputters_added.count( iitag->getName() ) ) continue;
				secondary_outputters_added.insert( iitag->getName() );

				// returns 0 if the secondary outputter is repressed for a particular job
				pose_outputters::SecondaryPoseOutputterOP outputter = secondary_outputter_for_job( inner_job, job_options, iitag->getName() );
				if ( outputter ) {
					secondary_outputters.push_back( outputter );
				}
			}
		}
	}

	if ( cl_outputters_.empty() ) {
		cl_outputters_ = pose_outputters::PoseOutputterFactory::get_instance()->secondary_pose_outputters_from_command_line();
	}
	for ( pose_outputters::SecondaryPoseOutputterOP const & elem : cl_outputters_ ) {
		if ( representative_secondary_outputter_map_.count( elem->class_key() ) == 0 ) {
			representative_secondary_outputter_map_[ elem->class_key() ] = elem;
		}
		if ( secondary_outputters_added.count( elem->class_key() ) ) continue;

		// returns 0 if the secondary outputter is repressed for a particular job
		pose_outputters::SecondaryPoseOutputterOP outputter = secondary_outputter_for_job(
			inner_job, job_options, elem->class_key() );
		if ( outputter ) {
			secondary_outputters.push_back( outputter );
		}
	}

	return secondary_outputters;
}


/// @details the nstruct count is taken from either the Job tag, or the
/// command line. "nstruct" is not a valid option to be provided
/// in the <Options> element.
core::Size
FullModelJobQueen::nstruct_for_job( utility::tag::TagCOP job_tag ) const
{
	if ( job_tag && job_tag->hasOption( "nstruct" ) ) {
		return job_tag->getOption< core::Size >( "nstruct" );
	}
	return basic::options::option[ basic::options::OptionKeys::out::nstruct ];
}


utility::options::OptionCollectionOP
FullModelJobQueen::options_for_job( FullModelInnerLarvalJob const & inner_job ) const
{
	using namespace utility::tag;
	TagCOP job_options_tag;
	if ( inner_job.jobdef_tag() ) {
		TagCOP job_tags = inner_job.jobdef_tag();
		if ( job_tags && job_tags->hasTag( "Options" ) ) {
			job_options_tag = job_tags->getTag( "Options" );
		}
	}

	// now let's walk through all of the option keys and read their values
	// out of the global options system into the per-job options object
	return options_from_tag( job_options_tag );
}


/// @details Note that jobs specified on the command line but which the FullModelJobQueen does
/// not know about do not get set or added to the per-job options.  Functions trying to read
/// per-job options that the FullModelJobQueen does not know about will die. This is intentional.
utility::options::OptionCollectionOP
FullModelJobQueen::options_from_tag( utility::tag::TagCOP job_options_tag ) const
{
	using namespace utility::tag;
	using namespace utility::options;

	utility::options::OptionKeyList all_options = concatenate_all_options();
	OptionCollectionOP opts = basic::options::read_subset_of_global_option_collection( all_options );

	TagCOP common_options_tag;
	if ( common_block_tags_ && common_block_tags_->hasTag( "Options" ) ) {
		common_options_tag = common_block_tags_->getTag( "Options" );
	}

	for ( auto const & option : all_options ) {
		utility::options::OptionKey const & opt( option() );
		OptionTypes opt_type = option_type_from_key( opt );

		std::string opt_tag_name = basic::options::replace_option_namespace_colons_with_underscores( opt );
		if ( job_options_tag && job_options_tag->hasTag( opt.identifier() ) ) {
			TagCOP opt_tag = job_options_tag->getTag( opt.identifier() );
			if ( opt_type == BOOLEAN_OPTION ) {
				(*opts)[ opt ].set_value( opt_tag->getOption< std::string >( "value", "true" ) );
			} else {
				debug_assert( opt_tag->hasOption( "value" ) );
				(*opts)[ opt ].set_value( opt_tag->getOption< std::string >( "value" ) );
			}
		} else if ( common_options_tag && common_options_tag->hasTag( opt.identifier() ) ) {
			TagCOP opt_tag = common_options_tag->getTag( opt.identifier() );
			if ( opt_type == BOOLEAN_OPTION ) {
				(*opts)[ opt ].set_value( opt_tag->getOption< std::string >( "value", "true" ) );
			} else {
				debug_assert( opt_tag->hasOption( "value" ) );
				(*opts)[ opt ].set_value( opt_tag->getOption< std::string >( "value" ) );
			}
		}
	}

	return opts;
}

void
FullModelJobQueen::determine_preliminary_job_list_from_xml_file(
	std::string const & job_def_string,
	std::string const & job_def_schema
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::tag;
	using namespace full_model_inputters;

	preliminary_larval_jobs_determined_ = true;

	load_job_definition_file( job_def_string, job_def_schema );

	// now iterate across all tags, and for each Job subtag, create a PreliminaryLarvalJob and load it
	// with all of the options that are within the <Option> subtag, if present -- and reading any options
	// not present in the tag from the (global) options system.
	Tag::tags_t const & subtags = job_definition_file_tags_->getTags();
	preliminary_larval_jobs_.reserve( subtags.size() );
	core::Size count_prelim_nodes( 0 );
	for ( auto subtag : subtags ) {
		if ( subtag->getName() != "Job" ) {
			debug_assert( subtag->getName() == "Common" );
			continue;
		}

		utility::options::OptionCollectionCOP job_options = options_from_tag( subtag );

		// ok -- look at input tag
		TagCOP input_tag = subtag->getTag( "Input" );
		debug_assert( input_tag ); // XML schema validation should ensure that there is an "Input" subelement
		debug_assert( input_tag->getTags().size() == 1 ); // schema validation should ensure there is exactly one subelement
		TagCOP input_tag_child = input_tag->getTags()[ 0 ];

		// Get the right inputter for this Job.
		full_model_inputters::FullModelInputterOP inputter;
		if ( use_factory_provided_full_model_inputters_ ) {
			inputter = full_model_inputters::FullModelInputterFactory::get_instance()->new_full_model_inputter( input_tag_child->getName() );
		} else {
			runtime_assert( inputter_creators_.count( input_tag_child->getName() ) != 0 );
			inputter = inputter_creators_[ input_tag_child->getName() ]->create_inputter();
		}

		FullModelInputSources input_poses = inputter->full_model_input_sources_from_tag( *job_options, input_tag_child );
		for ( auto input_source : input_poses ) { input_source->pose_id( ++input_pose_counter_ ); }

		// Get the right outputter for this job.
		pose_outputters::PoseOutputterOP outputter = get_outputter_from_job_tag( subtag );

		if ( representative_pose_outputter_map_.count( outputter->class_key() ) == 0 ) {
			representative_pose_outputter_map_[ outputter->class_key() ] = outputter;
		}

		// now iterate across the input sources for this job and create
		// a preliminary job for each
		core::Size nstruct = nstruct_for_job( subtag );
		for ( auto const & input_pose : input_poses ) {
			PreliminaryLarvalJob prelim_job;
			FullModelInnerLarvalJobOP inner_job( create_inner_larval_job( nstruct, ++count_prelim_nodes ) );
			inner_job->input_source( input_pose );
			inner_job->jobdef_tag( subtag );
			inner_job->outputter( outputter->class_key() );

			prelim_job.inner_job = inner_job;
			prelim_job.job_tag = subtag;
			prelim_job.job_options = job_options;
			prelim_job.full_model_inputter = inputter;

			preliminary_larval_jobs_.push_back( prelim_job );
		}
	}
}


void
FullModelJobQueen::load_job_definition_file(
	std::string const & job_def_string,
	std::string const & job_def_schema
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::tag;

	XMLValidationOutput validator_output;
	try {
		validator_output = validate_xml_against_xsd( job_def_string, job_def_schema );
	} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
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
		throw utility::excn::EXCN_Msg_Exception( oss.str() );
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
		throw utility::excn::EXCN_Msg_Exception( oss.str() );
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
FullModelJobQueen::determine_preliminary_job_list_from_command_line()
{
	using namespace utility::tag;
	using namespace full_model_inputters;

	// read from the command line a list of all of the input jobs
	FullModelInputterFactory::FullModelInputSourcesAndInputters input_poses;
	if ( use_factory_provided_full_model_inputters_ ) {
		input_poses = FullModelInputterFactory::get_instance()->full_model_inputs_from_command_line();
	} else {
		for ( auto inputter_creator : inputter_creator_list_ ) {
			FullModelInputterOP inputter = inputter_creator->create_inputter();
			if ( inputter->job_available_on_command_line() ) {
				FullModelInputSources iter_sources = inputter->full_model_input_sources_from_command_line();
				input_poses.reserve( input_poses.size() + iter_sources.size() );
				for ( core::Size ii = 1; ii <= iter_sources.size(); ++ii ) {
					input_poses.push_back( std::make_pair( iter_sources[ ii ], inputter ) );
				}
			}
		}
	}
	for ( auto input_source : input_poses ) { input_source.first->pose_id( ++input_pose_counter_ ); }

	pose_outputters::PoseOutputterOP outputter = get_outputter_from_job_tag( utility::tag::TagCOP() );

	if ( representative_pose_outputter_map_.count( outputter->class_key() ) == 0 ) {
		representative_pose_outputter_map_[ outputter->class_key() ] = outputter;
	}

	// pass in a null-pointing TagCOP and construct the job options object from the command line.
	utility::options::OptionCollectionCOP job_options = options_from_tag( utility::tag::TagCOP() );

	// now iterate across the input sources for this job and create
	// a preliminary job for each
	preliminary_larval_jobs_.reserve( input_poses.size() );
	core::Size count_prelim_nodes( 0 );
	for ( auto const & input_pose : input_poses ) {
		PreliminaryLarvalJob prelim_job;
		Size nstruct = nstruct_for_job( 0 );
		FullModelInnerLarvalJobOP inner_job( create_inner_larval_job( nstruct, ++count_prelim_nodes ) );
		inner_job->input_source( input_pose.first );
		inner_job->outputter( outputter->class_key() );

		prelim_job.inner_job = inner_job;
		prelim_job.job_tag = TagCOP(); // null ptr
		prelim_job.job_options = job_options;
		prelim_job.full_model_inputter = input_pose.second;
		preliminary_larval_jobs_.push_back( prelim_job );
	}
}

LarvalJobs
FullModelJobQueen::next_batch_of_larval_jobs_from_prelim( core::Size job_node_index, core::Size max_njobs )
{
	LarvalJobs jobs;
	if ( preliminary_job_nodes_complete_[ job_node_index ] ) return jobs;

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
			for ( FullModelInnerLarvalJobOP iijob : inner_larval_jobs_for_curr_prelim_job_ ) {
				pose_outputters::PoseOutputterOP ii_outputter = pose_outputter_for_job( *iijob );
				ii_outputter->determine_job_tag( output_tag, *curr_prelim_job.job_options, *iijob );
			}

		}

		if ( curr_inner_larval_job_index_ > inner_larval_jobs_for_curr_prelim_job_.size() ) {
			// prepare for the next time we call this function for a different job node
			curr_inner_larval_job_index_ = 0;
			pjn_job_ind_end_[ job_node_index ] = larval_job_counter_;
			preliminary_job_nodes_complete_[ job_node_index ] = 1;
			return jobs;
		}

		if ( curr_inner_larval_job_index_ <= inner_larval_jobs_for_curr_prelim_job_.size() ) {
			// create LarvalJobs
			FullModelInnerLarvalJobOP curr_inner_larval_job = inner_larval_jobs_for_curr_prelim_job_[ curr_inner_larval_job_index_ ];
			core::Size max_to_make = max_njobs;

			if ( max_to_make > curr_inner_larval_job->nstruct_max() - njobs_made_for_curr_inner_larval_job_ ) {
				max_to_make = curr_inner_larval_job->nstruct_max() - njobs_made_for_curr_inner_larval_job_;
			}

			if ( pjn_job_ind_begin_[ job_node_index ] == 0 ) {
				pjn_job_ind_begin_[ job_node_index ] = 1 + larval_job_counter_;
			}
			LarvalJobs curr_jobs = expand_job_list( curr_inner_larval_job, max_to_make );

			core::Size n_made = curr_jobs.size();
			if ( max_njobs >= n_made ) {
				max_njobs -= n_made;
			} else {
				max_njobs = 0;
				// this should never happen!
				throw utility::excn::EXCN_Msg_Exception( "expand_job_list returned " + utility::to_string( n_made ) + " jobs when it was given a maximum number of " + utility::to_string( max_to_make ) + " to make (with max_njobs of " + utility::to_string( max_njobs ) + ")\n" );
			}

			jobs.splice( jobs.end(), curr_jobs );
			if ( n_made + njobs_made_for_curr_inner_larval_job_ < curr_inner_larval_job->nstruct_max() ) {
				njobs_made_for_curr_inner_larval_job_ += n_made;
				return jobs;
			} else {
				++curr_inner_larval_job_index_;
				njobs_made_for_curr_inner_larval_job_ = 0;
			}
		}
	}
}

/// @details helper function that should only be called by the above secondary_outputter_for_job function
/// because of its assumption that the representative_secondary_outputter_map_ map contains an entry for
/// the requested secondary_outputter_type
pose_outputters::SecondaryPoseOutputterOP
FullModelJobQueen::secondary_outputter_for_job(
	FullModelInnerLarvalJob const & inner_job,
	utility::options::OptionCollection const & job_options,
	std::string const & secondary_outputter_type
)
{
	pose_outputters::SecondaryPoseOutputterOP representative_outputter;
	representative_outputter = representative_secondary_outputter_map_[ secondary_outputter_type ];
	debug_assert( representative_outputter );

	utility::tag::TagCOP output_tag;
	if ( inner_job.jobdef_tag() ) {
		utility::tag::TagCOP job_tag = inner_job.jobdef_tag();
		if ( job_tag->hasTag( "SecondaryOutput" ) ) {
			output_tag = job_tag->getTag( "SecondaryOutput" );
		}
	}
	std::string which_outputter = representative_outputter->outputter_for_job( output_tag, job_options, inner_job );


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

void
FullModelJobQueen::note_job_result_output_or_discarded( LarvalJobCOP job, Size result_index )
{
	auto result_pos_iter = results_processed_for_job_.find( job->job_index() );
	if ( result_pos_iter == results_processed_for_job_.end() ) {
		std::ostringstream oss;
		oss << "From StandardJobQueen::note_job_result_output_or_discarded:\n";
		oss << "Tried to note that job " << job->job_index() << " result # " << result_index <<
			" was output; this job has already had all of its results output or the" <<
			" StandardJobQueen was unaware of its existence\n";
		throw utility::excn::EXCN_Msg_Exception( oss.str() );
	}

	if ( result_pos_iter->second.results_output_or_discarded.member( result_index ) ) {
		// we've already output this job; the derived class has messed up here
		std::ostringstream oss;
		oss << "From StandardJobQueen::note_job_result_output_or_discarded:\n";
		oss << "Tried to note that job " << job->job_index() << " result # " << result_index <<
			" was output; this result index for this job has already been output or discarded.\n";
		throw utility::excn::EXCN_Msg_Exception( oss.str() );
	}


	if ( result_pos_iter->second.n_results == ++result_pos_iter->second.n_output_or_discarded ) {
		results_processed_for_job_.erase( result_pos_iter );
		processed_jobs_.insert( job->job_index() );
	} else {
		result_pos_iter->second.results_output_or_discarded.insert( result_index );
	}
}

/*void
FullModelJobQueen::note_job_result_output( LarvalJobCOP job )
{
output_jobs_.insert( job->job_index() );
}*/

utility::options::OptionKeyList
FullModelJobQueen::concatenate_all_options() const
{
	utility::options::OptionKeyList all_options( options_ );
	all_options.insert( all_options.end(), inputter_options_.begin(), inputter_options_.end() );
	all_options.insert( all_options.end(), outputter_options_.begin(), outputter_options_.end() );
	all_options.insert( all_options.end(), secondary_outputter_options_.begin(), secondary_outputter_options_.end() );
	return all_options;
}

pose_outputters::PoseOutputterOP
FullModelJobQueen::get_outputter_from_job_tag( utility::tag::TagCOP tag ) const
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
					outputter = pose_outputters::PoseOutputterOP( new pose_outputters::PDBPoseOutputter );
				}
			}
		}
	}
	return outputter;
}

} // namespace full_model
} // namespace jd3
} // namespace protocols
