// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/StandardJobQueen.cc
/// @brief  StandardJobQueen class's method definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/StandardJobQueen.hh>

// package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/MoverAndPoseJob.hh>
#include <protocols/jd3/PoseInputSource.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.hh>
#include <protocols/jd3/pose_inputters/PoseInputterFactory.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>

//project headers
#include <core/pose/Pose.hh>
// #include <basic/resource_manager/JobOptions.hh>

//utility headers
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/OptionCollection.hh>
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

namespace protocols {
namespace jd3 {

PreliminaryLarvalJob::PreliminaryLarvalJob() {}
PreliminaryLarvalJob::~PreliminaryLarvalJob() {}
PreliminaryLarvalJob::PreliminaryLarvalJob( PreliminaryLarvalJob const & src ) :
	inner_job( src.inner_job ),
	job_tag( src.job_tag )
{}

PreliminaryLarvalJob &
PreliminaryLarvalJob::operator = ( PreliminaryLarvalJob const & rhs )
{
	if ( this != &rhs ) {
		inner_job = rhs.inner_job;
		job_tag   = rhs.job_tag;
	}
	return *this;
}


StandardJobQueen::StandardJobQueen()
{
	// begin to populate the per-job options object
	pose_inputters::PoseInputterFactory::get_instance()->list_options_read( options_ );
	// TO DO: same thing for the PoseOutputterFactory.
}

StandardJobQueen::~StandardJobQueen() {}

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

	// write out option set
	if  ( ! options_.empty() ) {
		XMLSchemaComplexTypeGenerator option_generator;
		XMLSchemaSimpleSubelementList option_subelements;

		std::set< utility::keys::VariantKey< utility::options::OptionKey > > already_output_options;
		for ( OptionKeyList::const_iterator iter = options_.begin(); iter != options_.end(); ++iter ) {
			AttributeList attributes;
			utility::options::OptionKey const & opt_key( (*iter)() );

			// only output each option once, even if it is read in more than
			// one context
			if ( already_output_options.count( opt_key ) ) continue;
			already_output_options.insert( opt_key );

			OptionTypes opt_type = option_type_from_key( opt_key );
			XMLSchemaType value_attribute_type = value_attribute_type_for_option( opt_type );
			if ( option[ opt_key ].has_default() ) {
				if ( opt_type == BOOLEAN_OPTION ) {
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, option[ opt_key ].default_string() );
				} else {
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, option[ opt_key ].default_string() );
				}
			} else { // no default; value is required, unless it's a boolean option
				if ( opt_type == BOOLEAN_OPTION ) {
					attributes + XMLSchemaAttribute::attribute_w_default(  "value", value_attribute_type, "false" );
				} else {
					attributes + XMLSchemaAttribute::required_attribute( "value", value_attribute_type );
				}
			}
			option_subelements.add_simple_subelement( iter->identifier(), attributes );
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

/// @note This is only a temporary implementation of this function; future versions will
/// be compatible with the idea that the JobQueen defines not only a job for a single "round"
/// of work, but that she provides a DAG of work, recognizing nodes in the DAG, and providing
/// some but possibly not all of the jobs that correspond to a single node in the DAG.
/// This prototype is aimed solely at getting the XSD interface between the base class
/// %StandardJobQueen and the derived JobQueen hammerd out.
LarvalJobs
StandardJobQueen::determine_job_list()
{
	// ok -- we're going to look for a job definition file, and failing that, fall back on
	// the PoseInputterFactory to determine where the input sources are coming from.

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::tag;

	// Always generate a job definition schema and in the process, validate the schema.
	// If derived JobQueens have defined an invalid schema, execution must stop.
	std::string job_def_schema = job_definition_xsd();

	if ( option[ in::file::job_definition_file ].user() ) {
		std::string job_def_string = utility::file_contents( option[ in::file::job_definition_file ] );
		return determine_job_list_from_xml_file( job_def_string, job_def_schema );
	} else {
		return determine_job_list_from_command_line();
	}

}

bool
StandardJobQueen::has_job_completed( protocols::jd3::LarvalJobCOP job )
{
	return pose_outputter_for_job( *job->inner_job() )->job_has_already_completed( *job );
}

JobOP
StandardJobQueen::mature_larval_job( protocols::jd3::LarvalJobCOP larval_job )
{
	using namespace utility::options;
	using namespace utility::tag;
	using namespace basic::datacache;

	// initialize the options collection for this job.
	utility::options::OptionCollectionCOP job_options = options_for_job( * larval_job->inner_job() );

	return complete_larval_job_maturation( larval_job, job_options );
}

/// @details Construct the XSD and then invoke the (private) determine_job_list_from_xml_file method,
/// which is also invoked by determine_job_list.
LarvalJobs StandardJobQueen::determine_job_list_from_xml_file( std::string const & job_def_string )
{
	std::string job_def_schema = job_definition_xsd();
	return determine_job_list_from_xml_file( job_def_string, job_def_schema );
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

/// @details The base class implementation merely returns a one-element list containing the
/// input inner_job.  Derived classes have the flexibility to create several preliminary
/// jobs from this input job
InnerLarvalJobs
StandardJobQueen::refine_preliminary_job( PreliminaryLarvalJob const & prelim_job )
{
	InnerLarvalJobs one_job( 1, prelim_job.inner_job );
	return one_job;
}

LarvalJobs
StandardJobQueen::expand_job_list( InnerLarvalJobOP inner_job ) const
{
	core::Size nstruct = nstruct_for_job( *inner_job );
	LarvalJobs jobs;
	for ( core::Size jj = 1; jj <= nstruct; ++jj ) {
		LarvalJobOP job = create_larval_job( inner_job, jj );
		jobs.push_back( job );
	}
	return jobs;
}

InnerLarvalJobOP
StandardJobQueen::create_inner_larval_job() const
{
	return InnerLarvalJobOP( new InnerLarvalJob );
}

/// @details Factory method instantiates the base-class LarvalJob to start.
LarvalJobOP
StandardJobQueen::create_larval_job( InnerLarvalJobOP job, core::Size nstruct_index ) const
{
	return LarvalJobOP( new LarvalJob( job, nstruct_index ));
}

JobOP
StandardJobQueen::create_job( LarvalJobCOP ) const
{
	return JobOP( new MoverAndPoseJob );
}


void StandardJobQueen::add_options( utility::options::OptionKeyList const & opts )
{
	using namespace utility::options;
	for ( OptionKeyList::const_iterator iter = opts.begin(); iter != opts.end(); ++iter ) {
		options_.push_back( *iter );
	}
}

void StandardJobQueen::add_option( utility::options::OptionKey const & key )
{
	options_.push_back( key );
}

void StandardJobQueen::remove_default_input_element() {}

utility::tag::TagCOP
StandardJobQueen::common_block_tags() const
{
	return common_block_tags_;
}

core::pose::PoseOP
StandardJobQueen::pose_for_job( LarvalJobCOP job, utility::options::OptionCollection const & options ) const
{
	// either read the Pose in using the pose_inputter (and then keep a copy
	// in the resource manager), or retrieve the Pose from the resource manager
	// initial version: just read the pose in for each job.

	return pose_inputter_for_job( *job->inner_job() )->pose_from_input_source( job->inner_job()->input_source(), options );
}

//ResourceManagerOP StandardJobQueen::resource_manager()
//{}

/// @brief Access the pose inputter
pose_inputters::PoseInputterOP
StandardJobQueen::pose_inputter_for_job( InnerLarvalJob const & inner_job ) const
{
	return pose_inputters::PoseInputterFactory::get_instance()->new_pose_inputter( inner_job.input_source().origin() );
}

/// @brief Access the pose outputter
pose_outputters::PoseOutputterOP
StandardJobQueen::pose_outputter_for_job( InnerLarvalJob const & inner_job ) const
{
	return pose_outputters::PoseOutputterFactory::get_instance()->new_pose_outputter( inner_job.outputter() );
}

/// @details the nstruct count is taken from either the Job tag, or the
/// command line. "nstruct" is not a valid option to be provided
/// in the <Options> element.
core::Size
StandardJobQueen::nstruct_for_job( InnerLarvalJob const & inner_job ) const
{
	using namespace utility::tag;
	if ( inner_job.const_data_map().has( "tags", "job_tags" ) ) {
		TagCOP job_tags;
		job_tags = inner_job.const_data_map().get_ptr< Tag >( "tags", "job_tags" );
		if ( job_tags->hasOption( "nstruct" ) ) {
			return job_tags->getOption< core::Size >( "nstruct" );
		}
	}
	return basic::options::option[ basic::options::OptionKeys::out::nstruct ];
}


utility::options::OptionCollectionOP
StandardJobQueen::options_for_job( InnerLarvalJob const & inner_job ) const
{
	using namespace utility::tag;

	TagCOP job_options_tag;
	if ( inner_job.const_data_map().has( "tags", "job_tags" ) ) {
		TagCOP job_tags;
		job_tags = inner_job.const_data_map().get_ptr< Tag >( "tags", "job_tags" );
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

	for ( OptionKeyList::const_iterator iter = options_.begin(); iter != options_.end(); ++iter ) {
		utility::options::OptionKey const & opt( (*iter)() );
		OptionTypes opt_type = option_type_from_key( opt );

		if ( job_options_tag && job_options_tag->hasTag( opt.identifier() ) ) {
			TagCOP opt_tag = job_options_tag->getTag( opt.identifier() );
			if ( opt_type == BOOLEAN_OPTION ) {
				(*opts)[ opt ].set_cl_value( opt_tag->getOption< std::string >( "value", "true" ) );
			} else {
				debug_assert( opt_tag->hasOption( "value" ) );
				(*opts)[ opt ].set_cl_value( opt_tag->getOption< std::string >( "value" ) );
			}
		} else if ( common_options_tag && common_options_tag->hasTag( opt.identifier() ) ) {
			TagCOP opt_tag = common_options_tag->getTag( opt.identifier() );
			if ( opt_type == BOOLEAN_OPTION ) {
				(*opts)[ opt ].set_cl_value( opt_tag->getOption< std::string >( "value", "true" ) );
			} else {
				debug_assert( opt_tag->hasOption( "value" ) );
				(*opts)[ opt ].set_cl_value( opt_tag->getOption< std::string >( "value" ) );
			}
		}
	}

	return opts;
}

LarvalJobs
StandardJobQueen::determine_job_list_from_xml_file(
	std::string const & job_def_string,
	std::string const & job_def_schema
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::tag;

	try {
		validate_xml_against_xsd( job_def_string, job_def_schema );
	} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
		std::ostringstream oss;
		oss << "Job definition file \"" << option[ in::file::job_definition_file ]() << "\" failed to validate against"
			" the schema for this application\nUse the option -jd3::job_definition_schema <output filename> to output"
			" the schema to a file.\n" << e.msg() << "\n";
	}
	TagCOP job_def_tag = Tag::create( job_def_string );

	LarvalJobs jobs;
	// now iterate across all tags, and for each Job subtag, create a PreliminaryLarvalJob and load it
	// with all of the options that are within the <Option> subtag, if present -- and reading any options
	// not present in the tag from the (global) options system.
	Tag::tags_t const & subtags = job_def_tag->getTags();
	for ( Tag::tags_t::const_iterator iter = subtags.begin(); iter != subtags.end(); ++iter ) {
		TagCOP subtag = *iter;
		if ( subtag->getName() != "Job" ) {
			debug_assert( subtag->getName() == "Common" );
			common_block_tags_ = subtag;
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

		// now iterate across the input sources for this job and create
		// a preliminary job for each
		for ( PoseInputSources::const_iterator iter = input_poses.begin(); iter != input_poses.end(); ++iter ) {
			PreliminaryLarvalJob prelim_job;
			InnerLarvalJobOP inner_job( create_inner_larval_job() );
			inner_job->input_source( *iter );
			inner_job->const_data_map().add( "tags", "job_tags", subtag );
			inner_job->outputter( outputter->class_key() );

			prelim_job.inner_job = inner_job;
			prelim_job.job_tag = subtag;
			prelim_job.job_options = job_options;

			expand_preliminary_larval_job( prelim_job, outputter, job_options, output_tag, jobs );
		}
	}
	return jobs;
}

LarvalJobs
StandardJobQueen::determine_job_list_from_command_line()
{
	using namespace utility::tag;

	LarvalJobs jobs;
	// read from the command line a list of all of the input jobs
	PoseInputSources input_poses = pose_inputters::PoseInputterFactory::get_instance()->pose_inputs_from_command_line();

	pose_outputters::PoseOutputterOP outputter =
		pose_outputters::PoseOutputterFactory::get_instance()->pose_outputter_from_command_line();

	// pass in a null-pointing TagCOP and construct the job options object from the command line.
	utility::options::OptionCollectionCOP job_options = options_from_tag( utility::tag::TagCOP() );

	// now iterate across the input sources for this job and create
	// a preliminary job for each
	for ( PoseInputSources::const_iterator iter = input_poses.begin(); iter != input_poses.end(); ++iter ) {
		PreliminaryLarvalJob prelim_job;
		InnerLarvalJobOP inner_job( create_inner_larval_job() );
		inner_job->input_source( *iter );

		prelim_job.inner_job = inner_job;
		prelim_job.job_tag = TagCOP(); // null ptr
		prelim_job.job_options = job_options;

		expand_preliminary_larval_job( prelim_job, outputter, job_options, TagCOP(), jobs );
	}
	return jobs;
}

void
StandardJobQueen::expand_preliminary_larval_job(
	PreliminaryLarvalJob const & prelim_job,
	pose_outputters::PoseOutputterOP outputter,
	utility::options::OptionCollectionCOP job_options,
	utility::tag::TagCOP output_tag,
	LarvalJobs & jobs
)
{
	// ask the derived job queen to turn the expand or refine the list of inner larval jobs
	// for this preliminary job
	InnerLarvalJobs one_input_pose_inner_jobs = refine_preliminary_job( prelim_job );

	for ( InnerLarvalJobs::iterator inner_larval_job_iter = one_input_pose_inner_jobs.begin();
			inner_larval_job_iter != one_input_pose_inner_jobs.end(); ++inner_larval_job_iter ) {
		// now ask the job outputter to devise the job_tag for each inner job
		outputter->determine_job_tag( output_tag, *job_options, **inner_larval_job_iter );
		// and then expand the list of inner-larval jobs into a list of larval jobs, one for each nstruct
		LarvalJobs one_input_pose_jobs = expand_job_list( *inner_larval_job_iter );
		jobs.splice( jobs.end(), one_input_pose_jobs );
	}

}


} // namespace jd3
} // namespace protocols
