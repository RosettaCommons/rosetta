// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ScreeningFeatures.cc
/// @brief  report JSON object with information needed for vHTS postprocessing
/// @author Sam DeLuca

#include <protocols/features/ScreeningFeatures.hh>

#include <protocols/jd2/util.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/string_util.hh>
#include <utility/json_spirit/json_spirit_utils.h>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/tag/Tag.hh>
#include <algorithm>

#include <core/chemical/ResidueType.hh>
#include <utility/tools/make_vector.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/ScreeningFeaturesCreator.hh>

namespace protocols {
namespace features {

using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

ScreeningFeatures::ScreeningFeatures() = default;

ScreeningFeatures::ScreeningFeatures(ScreeningFeatures const & ) = default;


ScreeningFeatures::~ScreeningFeatures() = default;

// XRW TEMP std::string ScreeningFeatures::type_name() const
// XRW TEMP {
// XRW TEMP  return "ScreeningFeatures";
// XRW TEMP }

void ScreeningFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;
	Column struct_id("struct_id",DbDataTypeOP( new DbBigInt() ),false);
	Column chain_name("chain_id",DbDataTypeOP( new DbText(1) ),false);
	Column residue_number("residue_number",DbDataTypeOP( new DbInteger() ),false);
	Column name3("name3",DbDataTypeOP( new DbText() ),false);
	Column experiment_group("group_name",DbDataTypeOP( new DbText() ),false);
	Column descriptor_data("descriptor_data",DbDataTypeOP( new DbText() ),false);

	utility::vector1<Column> primary_keys;
	primary_keys.push_back(struct_id);
	primary_keys.push_back(residue_number);

	Schema screening_features("screening_features",PrimaryKey(primary_keys));
	screening_features.add_column(struct_id);
	screening_features.add_column(chain_name);
	screening_features.add_column(residue_number);
	screening_features.add_column(name3);
	screening_features.add_column(experiment_group);
	screening_features.add_column(descriptor_data);
	screening_features.add_foreign_key(
		ForeignKey(struct_id,"structures","struct_id")
	);
	screening_features.write(db_session);

}


utility::vector1<std::string> ScreeningFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	return dependencies;
}

core::Size
ScreeningFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & /*relevant_residues*/,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	InsertGenerator screening_insert("screening_features");
	screening_insert.add_column("struct_id");
	screening_insert.add_column("chain_id");
	screening_insert.add_column("residue_number");
	screening_insert.add_column("name3");
	screening_insert.add_column("group_name");
	screening_insert.add_column("descriptor_data");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id",struct_id) );
	RowDataBaseOP chain_id_data( new RowData<std::string>("chain_id",chain_) );

	std::vector<utility::json_spirit::Pair>  descriptor_json_data(get_desriptor_data());


	std::string descriptor_string(utility::json_spirit::write(utility::json_spirit::Value(descriptor_json_data)));
	RowDataBaseOP descriptor_data_column( new RowData<std::string>("descriptor_data",descriptor_string) );

	std::string group_name;

	std::map<  std::string, std::string > const & sspairmap( protocols::jd2::get_string_string_pairs_from_current_job() );
	if ( sspairmap.count( "input_group_name" ) ) {
		group_name = sspairmap.at( "input_group_name" );
	}

	RowDataBaseOP group_name_data( new RowData<std::string>("group_name",group_name) );

	core::Size chain_id = core::pose::get_chain_id_from_chain(chain_,pose);

	for ( core::Size resnum= pose.conformation().chain_begin(chain_id); resnum <= pose.conformation().chain_end(chain_id); ++resnum ) {
		RowDataBaseOP resnum_data( new RowData<core::Size>("residue_number",resnum) );
		RowDataBaseOP name3_data( new RowData<std::string>("name3",pose.residue_type(resnum).name3()) );

		screening_insert.add_row(utility::tools::make_vector(struct_id_data,chain_id_data,descriptor_data_column,group_name_data,resnum_data,name3_data));
	}
	screening_insert.write_to_database(db_session);

	return 0;


}


void
ScreeningFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/)
{

	if ( !basic::options::option[basic::options::OptionKeys::in::file::screening_job_file].user() ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "The ScreeningFeatures reporter is only usable with the ScreeningJobInputter. specify input using -in:file:screening_job_inputter");
	}

	if ( !tag->hasOption("chain") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScreeningFeatures requires the 'chain' tag");
	}

	chain_ = tag->getOption<std::string>("chain");

	core::Size descriptor_count = 0;
	for ( utility::tag::TagCOP const & sub_tag : tag->getTags() ) {
		if ( sub_tag->getName() != "descriptor" ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScreeningFeatures only supports subtags with the name 'descriptor");
		}

		if ( !sub_tag->hasOption("type") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScreeningFeatures descriptor subtags require a 'type' option");
		}

		std::string descriptor_type(sub_tag->getOption<std::string>("type"));
		descriptors_.push_back(descriptor_type);
		descriptor_count++;
	}
	if ( descriptor_count == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScreeningFeatures requires at least one 'descriptor' subtag");
	}

}

std::vector<utility::json_spirit::Pair>  ScreeningFeatures::get_desriptor_data() const
{
	std::vector<utility::json_spirit::Pair>  descriptor_data;

	for ( std::pair< std::string, std::string > const & sspair: protocols::jd2::get_string_string_pairs_from_current_job() ) {
		if ( std::find(descriptors_.begin(),descriptors_.end(),sspair.first) != descriptors_.end() ) {
			utility::json_spirit::Pair data_pair(sspair.first,sspair.second);
			descriptor_data.push_back(data_pair);
		}
	}

	for ( std::pair< std::string, core::Real > const & srpair: protocols::jd2::get_string_real_pairs_from_current_job()  ) {
		if ( std::find(descriptors_.begin(),descriptors_.end(),srpair.first) != descriptors_.end() ) {
			utility::json_spirit::Pair data_pair(srpair.first,srpair.second);
			descriptor_data.push_back(data_pair);
		}
	}

	return descriptor_data;



}

std::string ScreeningFeatures::type_name() const {
	return class_name();
}

std::string ScreeningFeatures::class_name() {
	return "ScreeningFeatures";
}

void ScreeningFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"chain", xs_string,
		"chain tag");

	AttributeList descriptor_attlist;
	descriptor_attlist + XMLSchemaAttribute::required_attribute(
		"type", xs_string,
		"descriptor type");

	XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement("descriptor", descriptor_attlist, "");

	protocols::features::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, class_name(),
		"report JSON object with information needed for vHTS postprocessing",
		attlist, ssl );
}

std::string ScreeningFeaturesCreator::type_name() const {
	return ScreeningFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ScreeningFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new ScreeningFeatures );
}

void ScreeningFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ScreeningFeatures::provide_xml_schema( xsd );
}



}
}
