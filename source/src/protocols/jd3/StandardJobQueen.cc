// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/StandardJobQueen.cc
/// @brief  StandardJobQueen class's method definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/StandardJobQueen.hh>

// package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/MoverAndPoseJob.hh>
#include <protocols/jd3/PoseInputSource.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.hh>
#include <protocols/jd3/pose_inputters/PoseInputterFactory.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.hh>

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
#include <utility/vector1.hh>

//basic headers
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.cc.gen.hh>
#include <basic/datacache/ConstDataMap.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.StandardJobQueen" );


namespace protocols {
namespace jd3 {

PreliminaryLarvalJob::PreliminaryLarvalJob() {}
PreliminaryLarvalJob::~PreliminaryLarvalJob() = default;
PreliminaryLarvalJob::PreliminaryLarvalJob( PreliminaryLarvalJob const & src ) :
	inner_job( src.inner_job ),
	job_tag( src.job_tag ),
	job_options( src.job_options )
{}

PreliminaryLarvalJob &
PreliminaryLarvalJob::operator = ( PreliminaryLarvalJob const & rhs )
{
	if ( this != &rhs ) {
		inner_job = rhs.inner_job;
		job_tag   = rhs.job_tag;
		job_options = rhs.job_options;
	}
	return *this;
}


InputSourceAndImportPoseOptions::InputSourceAndImportPoseOptions() :
	input_source_( new PoseInputSource ),
	import_pose_options_( new core::import_pose::ImportPoseOptions )
{}

InputSourceAndImportPoseOptions::InputSourceAndImportPoseOptions(
	PoseInputSource const & input_source,
	core::import_pose::ImportPoseOptions const & options
) :
	input_source_( new PoseInputSource( input_source )),
	import_pose_options_( new core::import_pose::ImportPoseOptions( options ))
{}

InputSourceAndImportPoseOptions::InputSourceAndImportPoseOptions( InputSourceAndImportPoseOptions const & src ) :
	input_source_( src.input_source_ ),
	import_pose_options_( src.import_pose_options_ )
{}

InputSourceAndImportPoseOptions::~InputSourceAndImportPoseOptions()
{}

InputSourceAndImportPoseOptions &
InputSourceAndImportPoseOptions::operator = ( InputSourceAndImportPoseOptions const & rhs )
{
	if ( this != &rhs ) {
		// shallow copy
		input_source_ = rhs.input_source_;
		import_pose_options_ = rhs.import_pose_options_;
	}
	return *this;
}

bool InputSourceAndImportPoseOptions::operator == ( InputSourceAndImportPoseOptions const & rhs ) const
{
	if ( input_source_ != rhs.input_source_ && ! (*input_source_ == *rhs.input_source_ ) ) return false;
	if ( import_pose_options_ != rhs.import_pose_options_ && ! (*import_pose_options_ == *rhs.import_pose_options_ ) ) return false;
	return true;
}

bool InputSourceAndImportPoseOptions::operator <  ( InputSourceAndImportPoseOptions const & rhs ) const
{
	if ( input_source_ != rhs.input_source_ && *input_source_ <  *rhs.input_source_ ) return true;
	if ( input_source_ != rhs.input_source_ || *input_source_ == *rhs.input_source_ ) return false;

	if ( import_pose_options_ != rhs.import_pose_options_ && *import_pose_options_ < *rhs.import_pose_options_ ) return true;

	return false;
}

PoseInputSource const & InputSourceAndImportPoseOptions::input_source() const { return *input_source_; }
void InputSourceAndImportPoseOptions::input_source( PoseInputSource const & setting ) {
	input_source_ = PoseInputSourceOP( new PoseInputSource( setting ));
}

core::import_pose::ImportPoseOptions const & InputSourceAndImportPoseOptions::import_pose_options() const
{
	return *import_pose_options_;
}

void InputSourceAndImportPoseOptions::import_pose_options( core::import_pose::ImportPoseOptions const & setting )
{
	using namespace core::import_pose;
	import_pose_options_ = ImportPoseOptionsOP( new ImportPoseOptions( setting ));
}


StandardJobQueen::StandardJobQueen() :
	job_definition_file_read_( false ),
	larval_job_counter_( 0 ),
	preliminary_larval_jobs_determined_( false ),
	curr_inner_larval_job_index_( 0 ),
	njobs_made_for_curr_inner_larval_job_( 0 )
{
	// begin to populate the per-job options object
	pose_inputters::PoseInputterFactory::get_instance()->list_options_read( options_ );
	pose_outputters::PoseOutputterFactory::get_instance()->list_options_read( options_ );
}

StandardJobQueen::~StandardJobQueen() = default;

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
		return xs_decimal;
	case REAL_VECTOR_OPTION :
		return xsct_real_wsslist;
	case INTEGER_OPTION :
		return xs_integer;
	case INTEGER_VECTOR_OPTION :
		return xsct_int_wsslist;
	case BOOLEAN_OPTION :
		return xs_boolean;
	case BOOLEAN_VECTOR_OPTION :
		return xsct_bool_wsslist; // note: double check that options system uses utility/string_funcs.hh to cast from strings to bools.
	default :
		throw utility::excn::EXCN_Msg_Exception( "Unsupported option type hit in StandardJobQueen.cc's value_attribute" );
	}
	return "ERROR";
}

std::string job_def_group() { return "job_def_jobs"; }
std::string job_complex_type_name( std::string const & type ) { return "job_def_" + type + "_type"; }
//std::string job_input_complex_type_name( std::string const & type ) { return "job_def_input_" + type + "_type"; }
//std::string job_output_complex_type_name( std::string const & type ) { return "job_def_output_" + type + "_type"; }

std::string
StandardJobQueen::job_definition_xsd() const
{
	using namespace utility::tag;
	using namespace utility::options;
	using namespace basic::options;

	XMLSchemaDefinition xsd;

	pose_inputters::PoseInputterFactory::get_instance()->define_pose_inputter_xml_schema( xsd );

	XMLSchemaSimpleSubelementList input_subelements;
	input_subelements.add_group_subelement( & pose_inputters::PoseInputterFactory::pose_inputter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator input_ct;
	input_ct
		.element_name( "Input" )
		.complex_type_naming_func( & job_complex_type_name )
		.set_subelements_pick_one( input_subelements )
		.write_complex_type_to_schema( xsd );

	pose_outputters::PoseOutputterFactory::get_instance()->define_pose_outputter_xml_schema( xsd );

	XMLSchemaSimpleSubelementList output_subelements;
	output_subelements.add_group_subelement( & pose_outputters::PoseOutputterFactory::pose_outputter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator output_ct;
	output_ct
		.element_name( "Output" )
		.complex_type_naming_func( & job_complex_type_name )
		.set_subelements_pick_one( output_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList secondary_output_subelements;
	secondary_output_subelements.add_group_subelement( & pose_outputters::PoseOutputterFactory::secondary_pose_outputter_xml_schema_group_name );
	XMLSchemaComplexTypeGenerator secondary_output_ct;
	secondary_output_ct
		.element_name( "SecondaryOutput" )
		.complex_type_naming_func( & job_complex_type_name )
		.set_subelements_single_appearance_optional( secondary_output_subelements )
		.write_complex_type_to_schema( xsd );

	// write out option set -- this should be method-extracted into a place accessible to other job queens
	if  ( ! options_.empty() ) {
		XMLSchemaComplexTypeGenerator option_generator;
		XMLSchemaSimpleSubelementList option_subelements;

		std::set< utility::keys::VariantKey< utility::options::OptionKey > > already_output_options;
		for ( auto const & iter : options_ ) {
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
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, option[ opt_key ].raw_default_string() );
				} else {
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, option[ opt_key ].raw_default_string() );
				}
			} else { // no default; value is required, unless it's a boolean option
				if ( opt_type == BOOLEAN_OPTION ) {
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, "false" );
				} else {
					attributes + XMLSchemaAttribute::required_attribute( "value", value_attribute_type );
				}
			}

			std::string decolonized_name = basic::options::replace_option_namespace_colons_with_underscores( iter );
			option_subelements.add_simple_subelement( decolonized_name, attributes );
		}
		option_generator.element_name( "Options" )
			.complex_type_naming_func( & job_complex_type_name )
			.set_subelements_single_appearance_optional( option_subelements )
			.write_complex_type_to_schema( xsd );
	}

	XMLSchemaComplexTypeGenerator job_ct;

	// Job <Input> element is required
	XMLSchemaSimpleSubelementList job_input_subelement;
	job_input_subelement.add_already_defined_subelement( "Input", & job_complex_type_name );
	job_ct.add_ordered_subelement_set_as_required( job_input_subelement );

	// Job <Output> element is optional
	XMLSchemaSimpleSubelementList job_output_subelement;
	job_output_subelement.add_already_defined_subelement( "Output", & job_complex_type_name );
	job_ct.add_ordered_subelement_set_as_optional( job_output_subelement );

	// Job <Options> element is optional
	XMLSchemaSimpleSubelementList job_options_subelement;
	job_options_subelement.add_already_defined_subelement( "Options", & job_complex_type_name );
	job_ct.add_ordered_subelement_set_as_optional( job_options_subelement );

	// Ask the derived class for what else belongs in the Job element.
	// NOTE: The derived class should only call the add_ordered_subelement_set_* functions
	// of the XMLSchemaComplexTypeGenerator or the <Input>, <Output>, and <Option>
	// subelements will be overwritten.
	append_job_tag_subelements( xsd, job_ct );

	// verify that the derived class did not call anything besides the add_ordered_subelement_set_*
	// functions.
	if ( job_ct.subelement_behavior() != se_ordered_sets ) {
		throw utility::excn::EXCN_Msg_Exception( "Subclass of StandardJobQueen's append_job_tag_subelements"
			" method invokes a method of the XMLSchemaComplexTypeGenerator that overwrote the <Input>, <Output>, and"
			" <Options> elements.  It should only call methods named \"add_ordered_subelement_set_*\"" );
	}

	job_ct
		.element_name( "Job" )
		.complex_type_naming_func( & job_complex_type_name )
		.add_attribute( XMLSchemaAttribute::attribute_w_default(  "nstruct", xsct_non_negative_integer, "1" ))
		.write_complex_type_to_schema( xsd );

	XMLSchemaComplexTypeGenerator common_block_ct_gen;

	// Common block <Options> subelement
	XMLSchemaSimpleSubelementList common_block_option_subelement;
	common_block_option_subelement.add_already_defined_subelement( "Options", & job_complex_type_name );
	common_block_ct_gen.add_ordered_subelement_set_as_optional( common_block_option_subelement );

	// Ask the derived class for what else belongs in the Common element.
	// NOTE: The derived class should only call the add_ordered_subelement_set_* functions
	// of the XMLSchemaComplexTypeGenerator or the <Options> subelement will be
	// overwritten.
	append_common_tag_subelements( xsd, common_block_ct_gen );

	// verify that the derived class did not call anything besides the add_ordered_subelement_set_*
	// functions.
	if ( common_block_ct_gen.subelement_behavior() != se_ordered_sets ) {
		throw utility::excn::EXCN_Msg_Exception( "Subclass of StandardJobQueen's append_job_tag_subelements"
			" method invokes a method of the XMLSchemaComplexTypeGenerator that overwrote the <Input>, <Output>, and"
			" <Options> elements.  It should only call methods named \"add_ordered_subelement_set_*\"" );
	}

	common_block_ct_gen
		.element_name( "Common" )
		.complex_type_naming_func( & job_complex_type_name )
		.write_complex_type_to_schema( xsd );

	XMLSchemaComplexTypeGenerator job_def_file_ct;
	XMLSchemaSimpleSubelementList job_def_subelements;
	job_def_subelements.add_already_defined_subelement( "Common", & job_complex_type_name, 0, 1 );
	job_def_subelements.add_already_defined_subelement( "Job", & job_complex_type_name, 1, xsminmax_unbounded );
	job_def_file_ct.element_name( "JobDefinitionFile" )
		.complex_type_naming_func( & job_complex_type_name )
		.set_subelements_single_appearance_required_and_ordered( job_def_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaElement root_element;
	root_element.name( "JobDefinitionFile" ).type_name( job_complex_type_name( "JobDefinitionFile" ));
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
StandardJobQueen::resource_definition_xsd() const
{
	// TO DO!
	return "";
}

/// @details The base class provides a digraph with a single node -- that is, all the jobs
/// are independent of each other.  This is equivalent to the kind of jobs that could be
/// run in JD2.
JobDigraphOP
StandardJobQueen::initial_job_dag()
{
	determine_preliminary_job_list();
	preliminary_job_nodes_complete_.resize( preliminary_larval_jobs_.size() );
	for ( core::Size ii = 1; ii <= preliminary_larval_jobs_.size(); ++ii ) {
		preliminary_job_nodes_complete_[ ii ] = 0;
	}
	// create a DAG with as many nodes in it as there are preliminary larval jobs
	job_graph_ = JobDigraphOP( new JobDigraph( preliminary_larval_jobs_.size() ) );
	return job_graph_;
}

void StandardJobQueen::update_job_dag( JobDigraphUpdater & ) {}

/// @details The process begins by first constructing the job definition and resource definition
/// XSDs.  With these schemas, the %StandardJobQueen validates the input XML files (if present).
/// The %StandardJobQueen then populates preliminary versions of LarvalJob objects./ If the XSD
/// includes "command line options" (which may be specified either from the command line or in
/// the <options> section of the Job XML file), the %StandardJobQueen loads the preliminary
/// LarvalJob objects with the options. These preliminary LarvalJob objects will not have been
/// nstruct expanded (i.e. if there are 100 nstruct for each of 5 different jobs, then there will
/// only be 5 preliminary larval jobs created). It then passes the preliminary LarvalJob list and
/// the TagOP objects for each preliminary LarvalJob to the derived class through the
/// refine_job_list method.
///
/// @note This is only a temporary implementation of this function; future versions will
/// be compatible with the idea that the JobQueen defines not only a job for a single "round"
/// of work, but that she provides a DAG of work, recognizing nodes in the DAG, and providing
/// some but possibly not all of the jobs that correspond to a single node in the DAG.
/// This prototype is aimed solely at getting the XSD interface between the base class
/// %StandardJobQueen and the derived JobQueen hammerd out.
LarvalJobs
StandardJobQueen::determine_job_list( Size job_dag_node_index, Size max_njobs )
{
	// ok -- we're going to look for a job definition file, and failing that, fall back on
	// the PoseInputterFactory to determine where the input sources are coming from.

	if ( ! preliminary_larval_jobs_determined_ ) {
		determine_preliminary_job_list();
	}

	if ( job_dag_node_index <= preliminary_larval_jobs_.size() ) {

		// now that the PreliminaryLarvalJobs have been constructed, go ahead
		// and start delivering LarvalJobs to the JobDistributor.
		return next_batch_of_larval_jobs_from_prelim( job_dag_node_index, max_njobs );
	} else {
		return next_batch_of_larval_jobs_for_job_node( job_dag_node_index, max_njobs );
	}

}

bool
StandardJobQueen::has_job_completed( protocols::jd3::LarvalJobCOP job )
{
	return pose_outputter_for_job( *job->inner_job() )->job_has_already_completed( *job );
}

void
StandardJobQueen::mark_job_as_having_begun( protocols::jd3::LarvalJobCOP job )
{
	return pose_outputter_for_job( *job->inner_job() )->mark_job_as_having_started( *job );
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

	// The common_block_tags_ must be loaded prior to the call of complete_larval_job_maturation
	// and that may not have happened on this JobQueen.
	if ( ! job_definition_file_read_ ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		job_definition_file_read_ = true;
		if ( option[ in::file::job_definition_file ].user() ) {
			std::string job_def_schema = job_definition_xsd();
			std::string job_def_string = utility::file_contents( option[ in::file::job_definition_file ] );
			load_job_definition_file( job_def_string, job_def_schema );
		}
	}


	// initialize the options collection for this job.
	utility::options::OptionCollectionCOP job_options = options_for_job( * larval_job->inner_job() );

	return complete_larval_job_maturation( larval_job, job_options, input_results );
}

bool StandardJobQueen::larval_job_needed_for_note_job_completed() const { return false; }

void StandardJobQueen::note_job_completed( LarvalJobCOP job, JobStatus status )
{
	// pass the status to the other note_job_completed function
	note_job_completed( job->job_index(), status );
}

void StandardJobQueen::note_job_completed( core::Size job_id, JobStatus status )
{
	if ( status == jd3_job_status_success ) {
		recent_successes_.push_back( job_id );
		successful_jobs_.insert( job_id );
	} else {
		failed_jobs_.insert( job_id );
	}
}

bool StandardJobQueen::larval_job_needed_for_completed_job_summary() const { return false; }

void StandardJobQueen::completed_job_summary( LarvalJobCOP, JobSummaryOP ) {}

void StandardJobQueen::completed_job_summary( core::Size, JobSummaryOP ) {}

std::list< core::Size > StandardJobQueen::jobs_that_should_be_output()
{
	std::list< core::Size > return_list = recent_successes_;
	recent_successes_.clear();
	return return_list;
}

/// @details Default implementation does not discard any job results.
std::list< core::Size >
StandardJobQueen::job_results_that_should_be_discarded() {
	return std::list< core::Size >();
}

void StandardJobQueen::completed_job_result( LarvalJobCOP job, JobResultOP job_result )
{
	pose_outputters::PoseOutputterOP outputter = pose_outputter_for_job( *job->inner_job() );
	PoseJobResultOP pose_result = utility::pointer::dynamic_pointer_cast< PoseJobResult >( job_result );
	if ( ! pose_result ) {
		utility::excn::EXCN_Msg_Exception( "JobResult handed to StandardJobQueen::completed_job_result is not a PoseJobResult or derived from PoseJobResult" );
	}
	utility::options::OptionCollectionOP job_options = options_for_job( *job->inner_job() );
	outputter->write_output_pose( *job, *job_options, *pose_result->pose() );

	std::list< pose_outputters::SecondaryPoseOutputterOP > secondary_outputters = secondary_outputters_for_job( *job->inner_job(), *job_options );
	for ( std::list< pose_outputters::SecondaryPoseOutputterOP >::const_iterator
			iter = secondary_outputters.begin(), iter_end = secondary_outputters.end();
			iter != iter_end; ++iter ) {
		(*iter)->write_output_pose( *job, *job_options, *pose_result->pose() );
	}
}

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

	// Always generate a job definition schema and in the process, validate the schema.
	// If derived JobQueens have defined an invalid schema, execution must stop.
	std::string job_def_schema = job_definition_xsd();
	// TO DO: validate the XSD

	if ( option[ in::file::job_definition_file ].user() ) {
		std::string job_def_string = utility::file_contents( option[ in::file::job_definition_file ] );
		job_definition_file_read_ = true;
		determine_preliminary_job_list_from_xml_file( job_def_string, job_def_schema );
	} else {
		determine_preliminary_job_list_from_command_line();
	}

}

/// @brief Read access for the subset of nodes in the job DAG which the %StandardJobQueen
/// is responsible for producing the larval_jobs. They are called "preliminary" jobs because
/// they do not depend on outputs from any previous node in the graph. (The set of job nodes
/// that contain no incoming edges, though, could perhaps be different from the set of
/// preliminary job nodes, so the %StandardJobQueen requires that the deried job queen
/// inform her of which nodes are the preliminary job nodes.)
utility::vector1< core::Size > const &
StandardJobQueen::preliminary_job_nodes() const
{
	return preliminary_job_nodes_complete_;
}

/////// @brief Allow the derived job queen to specify a node in the JobDAG as being
/////// "preliminary" in the sense a) that the %StandardJobQueen is responsible for creating the
/////// list of larval jobs for this node, and b) there are no nodes that this node depends
/////// on having completed before it can run.
////void
////StandardJobQueen::declare_job_node_to_be_preliminary( core::Size job_node_index )
////{
//// if ( job_node_index > job_graph_.num_nodes() ) {
////  throw utility::excn::EXCN_Msg_Exception( "Could not declare job node " + utility::to_string( job_node_index )
////   + " as a preliminary job node as the job_graph_ has only " + utility::to_string( job_graph_.num_nodes() ) + "nodes." );
//// }
//// if ( job_graph_->get_node()->indegree() != 0 ) {
////  throw utility::excn::EXCN_Msg_Exception( "Could not declare job node " + utility::to_string( job_node_index )
////   + " as a preliminary job node as this node has an indegree of " + utility::to_string( job_graph_->get_node()->indegree() )
////   + " and it must instead have an indegree of 0." );
//// }
//// preliminary_job_nodes_.push_back( job_node_index );
////}

/// @brief Read access to derived JobQueens to the preliminary job list.
/// This will return an empty list if  determine_preliminary_jobs has not yet
/// been called.
utility::vector1< PreliminaryLarvalJob > const &
StandardJobQueen::preliminary_larval_jobs() const
{
	return preliminary_larval_jobs_;
}

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
		LarvalJobOP job = create_larval_job( inner_job, njobs_made_for_curr_inner_larval_job_ + jj, ++larval_job_counter_ );
		//TR << "Expand larval job " << larval_job_counter_ << std::endl;
		jobs.push_back( job );
	}
	return jobs;
}

InnerLarvalJobOP
StandardJobQueen::create_inner_larval_job( core::Size nstruct ) const
{
	return InnerLarvalJobOP( new InnerLarvalJob( nstruct ) );
}

/// @details Factory method instantiates the base-class LarvalJob to start.
/// Non-const so that derived classes may do thing like keep track of each of the
/// larval jobs they instantiate.
LarvalJobOP
StandardJobQueen::create_larval_job(
	InnerLarvalJobOP job,
	core::Size nstruct_index,
	core::Size larval_job_index
)
{
	return LarvalJobOP( new LarvalJob( job, nstruct_index, larval_job_index ));
}

JobOP
StandardJobQueen::create_job( LarvalJobCOP ) const
{
	return JobOP( new MoverAndPoseJob );
}

LarvalJobs
StandardJobQueen::next_batch_of_larval_jobs_for_job_node( core::Size, core::Size )
{
	return LarvalJobs(); // default ctor to give an empty list
}


void StandardJobQueen::add_options( utility::options::OptionKeyList const & opts )
{
	using namespace utility::options;
	for ( auto const & opt : opts ) {
		options_.push_back( opt );
	}
}

void StandardJobQueen::add_option( utility::options::OptionKey const & key )
{
	options_.emplace_back(key );
}

/// @details TO DO
void StandardJobQueen::remove_default_input_element() {}

utility::tag::TagCOP
StandardJobQueen::common_block_tags() const
{
	return common_block_tags_;
}

core::pose::PoseOP
StandardJobQueen::pose_for_job( LarvalJobCOP job, utility::options::OptionCollection const & options )
{
	// either read the Pose in using the pose_inputter (and then keep a copy
	// in the resource manager), or retrieve the Pose from the resource manager
	// initial version: just read the pose in for each job.

	typedef utility::pointer::shared_ptr< InputSourceAndImportPoseOptions > InputSourceAndImportPoseOptionsOP;

	InputSourceAndImportPoseOptionsOP settings_ptr;
	if ( job->inner_job()->input_source().origin() == "PDB" ) {
		core::import_pose::ImportPoseOptions import_opts( options );
		settings_ptr = InputSourceAndImportPoseOptionsOP( new InputSourceAndImportPoseOptions( job->inner_job()->input_source(), import_opts ));

		std::map< InputSourceAndImportPoseOptions, core::pose::PoseOP >::const_iterator iter =
			previously_read_in_poses_.find( *settings_ptr );
		if ( iter != previously_read_in_poses_.end() ) {
			core::pose::PoseOP pose( new core::pose::Pose );
			// Why does the standard job queen use detached_copy? Because it is important in multithreaded
			// contexts that no two Poses pointing to that same data end up in two separate threads,
			// and then try to modify that data at the same time.  It turns out there are places
			// in Pose where it covertly modifies data in some other Pose and this could lead to
			// race conditions.
			pose->detached_copy( *iter->second );
			return pose;
		}
	}

	core::pose::PoseOP pose = pose_inputter_for_job( *job->inner_job() )->pose_from_input_source( job->inner_job()->input_source(), options );
	if ( settings_ptr ) {
		previously_read_in_poses_[ *settings_ptr ] = pose;
	}
	return pose;
}

//ResourceManagerOP StandardJobQueen::resource_manager()
//{}

/// @brief Access the pose inputter
pose_inputters::PoseInputterOP
StandardJobQueen::pose_inputter_for_job( InnerLarvalJob const & inner_job ) const
{
	return pose_inputters::PoseInputterFactory::get_instance()->new_pose_inputter( inner_job.input_source().origin() );
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
	pose_outputters::PoseOutputterOP representative_outputter;
	representative_outputter = representative_pose_outputter_map_[ inner_job.outputter() ];

	utility::tag::TagCOP output_tag;
	if ( inner_job.jobdef_tag() ) {
		utility::tag::TagCOP job_tag = inner_job.jobdef_tag();
		if ( job_tag->hasTag( "Output" ) ) {
			output_tag = job_tag->getTag( "Output" );
		}
	}
	std::string which_outputter = representative_outputter->outputter_for_job( output_tag, job_options, inner_job );

	if ( which_outputter == "" ) {
		return representative_outputter;
	}
	PoseOutputterMap::const_iterator iter = pose_outputter_map_.find( std::make_pair( inner_job.outputter(), which_outputter ) );
	if ( iter != pose_outputter_map_.end() ) {
		return iter->second;
	}

	pose_outputters::PoseOutputterOP outputter = pose_outputters::PoseOutputterFactory::get_instance()->new_pose_outputter( inner_job.outputter() );
	pose_outputter_map_[ std::make_pair( inner_job.outputter(), which_outputter ) ] = outputter;
	return outputter;
}


std::list< pose_outputters::SecondaryPoseOutputterOP >
StandardJobQueen::secondary_outputters_for_job(
	InnerLarvalJob const & inner_job,
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
	for ( std::list< pose_outputters::SecondaryPoseOutputterOP >::const_iterator
			iter = cl_outputters_.begin(), iter_end = cl_outputters_.end();
			iter != iter_end; ++iter ) {

		if ( representative_secondary_outputter_map_.count( (*iter)->class_key() ) == 0 ) {
			representative_secondary_outputter_map_[ (*iter)->class_key() ] = *iter;
		}
		if ( secondary_outputters_added.count( (*iter)->class_key() ) ) continue;

		// returns 0 if the secondary outputter is repressed for a particular job
		pose_outputters::SecondaryPoseOutputterOP outputter = secondary_outputter_for_job(
			inner_job, job_options, (*iter)->class_key() );
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
StandardJobQueen::nstruct_for_job( utility::tag::TagCOP job_tag ) const
{
	if ( job_tag && job_tag->hasOption( "nstruct" ) ) {
		return job_tag->getOption< core::Size >( "nstruct" );
	}
	return basic::options::option[ basic::options::OptionKeys::out::nstruct ];
}


utility::options::OptionCollectionOP
StandardJobQueen::options_for_job( InnerLarvalJob const & inner_job ) const
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


/// @details Note that jobs specified on the command line but which the StandardJobQueen does
/// not know about do not get set or added to the per-job options.  Functions trying to read
/// per-job options that the StandardJobQueen does not know about will die. This is intentional.
utility::options::OptionCollectionOP
StandardJobQueen::options_from_tag( utility::tag::TagCOP job_options_tag ) const
{
	using namespace utility::tag;
	using namespace utility::options;

	OptionCollectionOP opts = basic::options::read_subset_of_global_option_collection( options_ );

	TagCOP common_options_tag;
	if ( common_block_tags_ && common_block_tags_->hasTag( "Options" ) ) {
		common_options_tag = common_block_tags_->getTag( "Options" );
	}

	for ( auto const & option : options_ ) {
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
StandardJobQueen::determine_preliminary_job_list_from_xml_file(
	std::string const & job_def_string,
	std::string const & job_def_schema
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::tag;

	load_job_definition_file( job_def_string, job_def_schema );

	// now iterate across all tags, and for each Job subtag, create a PreliminaryLarvalJob and load it
	// with all of the options that are within the <Option> subtag, if present -- and reading any options
	// not present in the tag from the (global) options system.
	Tag::tags_t const & subtags = job_definition_file_tags_->getTags();
	preliminary_larval_jobs_.reserve( subtags.size() );
	for ( auto subtag : subtags ) {
		if ( subtag->getName() != "Job" ) {
			debug_assert( subtag->getName() == "Common" );
			continue;
		}

		utility::options::OptionCollectionCOP job_options = options_from_tag( subtag );

		// ok -- look at input tag
		// which for now should be a PDB tag
		TagCOP input_tag = subtag->getTag( "Input" );
		debug_assert( input_tag ); // XML schema validation should ensure that there is an "Input" subelement
		debug_assert( input_tag->getTags().size() == 1 ); // schema validation should ensure there is exactly one subelement
		TagCOP input_tag_child = input_tag->getTags()[ 0 ];
		pose_inputters::PoseInputterOP inputter =
			pose_inputters::PoseInputterFactory::get_instance()->new_pose_inputter( input_tag_child->getName() );
		PoseInputSources input_poses = inputter->pose_input_sources_from_tag( input_tag_child );

		pose_outputters::PoseOutputterOP outputter;
		TagCOP output_tag;
		if ( subtag->hasTag( "Output" ) ) {
			output_tag = subtag->getTag( "Output" );
			debug_assert( output_tag );
			debug_assert( output_tag->getTags().size() == 1 );
			TagCOP output_subtag = output_tag->getTags()[ 0 ];
			outputter = pose_outputters::PoseOutputterFactory::get_instance()->new_pose_outputter( output_subtag->getName() );
		} else {
			outputter = pose_outputters::PoseOutputterFactory::get_instance()->pose_outputter_from_command_line();
		}

		if ( representative_pose_outputter_map_.count( outputter->class_key() ) == 0 ) {
			representative_pose_outputter_map_[ outputter->class_key() ] = outputter;
		}

		// now iterate across the input sources for this job and create
		// a preliminary job for each
		core::Size nstruct = nstruct_for_job( subtag );
		for ( PoseInputSources::const_iterator iter = input_poses.begin(); iter != input_poses.end(); ++iter ) {
			PreliminaryLarvalJob prelim_job;
			InnerLarvalJobOP inner_job( create_inner_larval_job( nstruct ) );
			inner_job->input_source( *iter );
			inner_job->jobdef_tag( subtag );
			inner_job->outputter( outputter->class_key() );

			prelim_job.inner_job = inner_job;
			prelim_job.job_tag = subtag;
			prelim_job.job_options = job_options;

			preliminary_larval_jobs_.push_back( prelim_job );
		}
	}
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
	for ( Tag::tags_t::const_iterator iter = subtags.begin(); iter != subtags.end(); ++iter ) {
		TagCOP subtag = *iter;
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

	// read from the command line a list of all of the input jobs
	PoseInputSources input_poses = pose_inputters::PoseInputterFactory::get_instance()->pose_inputs_from_command_line();

	pose_outputters::PoseOutputterOP outputter =
		pose_outputters::PoseOutputterFactory::get_instance()->pose_outputter_from_command_line();

	if ( representative_pose_outputter_map_.count( outputter->class_key() ) == 0 ) {
		representative_pose_outputter_map_[ outputter->class_key() ] = outputter;
	}

	// pass in a null-pointing TagCOP and construct the job options object from the command line.
	utility::options::OptionCollectionCOP job_options = options_from_tag( utility::tag::TagCOP() );

	// now iterate across the input sources for this job and create
	// a preliminary job for each
	preliminary_larval_jobs_.reserve( input_poses.size() );
	for ( PoseInputSources::const_iterator iter = input_poses.begin(); iter != input_poses.end(); ++iter ) {
		PreliminaryLarvalJob prelim_job;
		Size nstruct = nstruct_for_job( 0 );
		InnerLarvalJobOP inner_job( create_inner_larval_job( nstruct ) );
		inner_job->input_source( *iter );
		inner_job->outputter( outputter->class_key() );

		prelim_job.inner_job = inner_job;
		prelim_job.job_tag = TagCOP(); // null ptr
		prelim_job.job_options = job_options;
		preliminary_larval_jobs_.push_back( prelim_job );
	}
}

LarvalJobs
StandardJobQueen::next_batch_of_larval_jobs_from_prelim( core::Size job_node_index, core::Size max_njobs )
{
	LarvalJobs jobs;
	if ( preliminary_job_nodes_complete_[ job_node_index ] ) return jobs;


	core::Size njobs_already_made( 0 );
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
			return jobs;
		}

		if ( curr_inner_larval_job_index_ <= inner_larval_jobs_for_curr_prelim_job_.size() ) {
			// create LarvalJobs
			InnerLarvalJobOP curr_inner_larval_job = inner_larval_jobs_for_curr_prelim_job_[ curr_inner_larval_job_index_ ];
			core::Size max_to_make = max_njobs;
			if ( max_njobs > njobs_already_made ) {
				max_to_make = max_njobs - njobs_already_made;
			}
			LarvalJobs curr_jobs = expand_job_list( curr_inner_larval_job, max_to_make );
			core::Size n_made = curr_jobs.size();
			jobs.splice( jobs.end(), curr_jobs );
			if ( n_made + njobs_made_for_curr_inner_larval_job_ <= curr_inner_larval_job->nstruct_max() ) {
				njobs_made_for_curr_inner_larval_job_ += n_made;
				preliminary_job_nodes_complete_[ job_node_index ] = 1;
				return jobs;
			} else {
				njobs_already_made += n_made;
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
StandardJobQueen::secondary_outputter_for_job(
	InnerLarvalJob const & inner_job,
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


} // namespace jd3
} // namespace protocols
