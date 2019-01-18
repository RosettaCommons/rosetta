// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/util.cc
/// @brief Utility functions for JD3.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd3/util.hh>

//Core headers
#include <core/select/residue_selector/util.hh>

//Protocol headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/parser/ResidueSelectorLoader.hh>
#include <protocols/parser/ScoreFunctionLoader.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>

// Basic headers
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

static basic::Tracer TR( "protocols.jd3.util" );

namespace protocols {
namespace jd3 {

using namespace utility::tag;
using namespace utility::options;
using namespace basic::options;

void
append_residue_selector_subelements(
	utility::tag::XMLSchemaDefinition & job_definition_xsd,
	utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen)
{
	parser::ResidueSelectorLoader::provide_xml_schema( job_definition_xsd );
	XMLSchemaSimpleSubelementList selector_subelements;
	selector_subelements.add_already_defined_subelement( parser::ResidueSelectorLoader::loader_name(), & parser::ResidueSelectorLoader::res_selector_loader_ct_namer );
	job_ct_gen.add_ordered_subelement_set_as_repeatable( selector_subelements);
}

void
append_single_scorefxn_subelement(
	utility::tag::XMLSchemaDefinition & job_definition_xsd,
	utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen )
{
	parser::ScoreFunctionLoader::provide_xml_schema( job_definition_xsd );
	XMLSchemaSimpleSubelementList sfxn_subelements;
	std::string scorefxn_name = "ScoreFunction";
	sfxn_subelements.add_already_defined_subelement( scorefxn_name, & parser::ScoreFunctionLoader::score_function_loader_ct_namer );
	job_ct_gen.add_ordered_subelement_set_as_optional( sfxn_subelements);
}

///@brief Add all residue selectors to a given datamap from a tag.
void
add_residue_selectors_to_datamap(
	core::pose::Pose const & pose,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{
	if ( ! tag ) return;

	using namespace utility::tag;
	using namespace protocols::parser;

	for ( auto iter = tag->getTags().begin(); iter != tag->getTags().end(); ++iter ) {
		if ( (*iter)->getName() == ResidueSelectorLoader::loader_name() ) {
			ResidueSelectorLoader selector_loader;
			selector_loader.load_data( pose, *iter, datamap );
		}
	}
}

///@brief Return a single residue selector from tag.
core::select::residue_selector::ResidueSelectorCOP
parse_single_residue_selector(
	core::pose::Pose const & pose,
	utility::tag::TagCOP tag,
	std::string const & selector_name )
{
	basic::datacache::DataMap datamap;
	add_residue_selectors_to_datamap(pose, tag, datamap);
	return core::select::residue_selector::get_residue_selector( selector_name, datamap );
}

/// @details the nstruct count is taken from either the Job tag, or the
/// command line. "nstruct" is not a valid option to be provided
/// in the <Options> element and so we do not have to worry about passing
/// in as a parameter a local OptionCollection object for this job.
core::Size
nstruct_for_job( utility::tag::TagCOP job_tag )
{
	if ( job_tag && job_tag->hasOption( "nstruct" ) ) {
		return job_tag->getOption< core::Size >( "nstruct" );
	}
	return basic::options::option[ basic::options::OptionKeys::out::nstruct ];
}

/// @details Note that jobs specified on the command line but which the StandardJobQueen does
/// not know about do not get set or added to the per-job options.  Functions trying to read
/// per-job options that the StandardJobQueen does not know about will die. This is intentional.
utility::options::OptionCollectionOP
options_from_tag(
	utility::options::OptionKeyList const & option_keys,
	utility::tag::TagCOP job_options_tag,
	utility::tag::TagCOP common_block_tags )
{
	using namespace utility::tag;
	using namespace utility::options;

	OptionCollectionOP opts = basic::options::read_subset_of_global_option_collection( option_keys );

	TagCOP common_options_tag;
	if ( common_block_tags && common_block_tags->hasTag( "Options" ) ) {
		common_options_tag = common_block_tags->getTag( "Options" );
	}

	for ( auto const & option : option_keys ) {
		utility::options::OptionKey const & opt( option() );
		OptionTypes opt_type = option_type_from_key( opt );

		std::string opt_tag_name = basic::options::replace_option_namespace_colons_with_underscores( opt );
		if ( job_options_tag && job_options_tag->hasTag( opt_tag_name ) ) {
			TagCOP opt_tag = job_options_tag->getTag( opt_tag_name );
			if ( opt_type == BOOLEAN_OPTION ) {
				(*opts)[ opt ].set_value( opt_tag->getOption< std::string >( "value", "true" ) );
			} else {
				debug_assert( opt_tag->hasOption( "value" ) );
				(*opts)[ opt ].set_value( opt_tag->getOption< std::string >( "value" ) );
			}
		} else if ( common_options_tag && common_options_tag->hasTag( opt_tag_name ) ) {
			TagCOP opt_tag = common_options_tag->getTag( opt_tag_name );
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

utility::options::OptionCollectionOP
options_for_job(
	utility::options::OptionKeyList const & option_keys,
	InnerLarvalJob const & inner_job,
	utility::tag::TagCOP common_block_tags )
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
	return options_from_tag( option_keys, job_options_tag, common_block_tags );
}




void
print_job_template(){



	TR << "The following is an empty (template) Job file:\n"
		<< "\n"
		<< "<JobDefinitionFile>\n"
		<< "\t<Job>\n"
		<< "\t\t<Input>\n"
		<< "\t\t</Input>\n"
		<< "\t\t<Output>\n"
		<< "\t\t</Output>\n"
		<< "\t\t<Options>\n"
		<< "\t\t</Options>\n"
		<< "\t</Job>\n"
		<< "</JobDefinitionFile>\n\n";

	TR << "Common Input Tags: \n "<<
		"<PDB filename=\"my_fname.pdb\"/>\n"<< std::endl;

	TR << "Common Ouptut Tags: \n" <<
		"<PDB filename_pattern=\"Some_job-specific_name_$\"/>\n" << std::endl;

	TR << std::endl;
	TR << "The rosetta_scripts application will now exit." << std::endl;
	TR.flush();

}

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
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Unsupported option type hit in JD3 value_attribute" );
	}
	return "ERROR";
}



} //protocols
} //jd3


