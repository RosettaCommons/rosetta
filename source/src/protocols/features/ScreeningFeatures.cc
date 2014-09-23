// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ScreeningFeatures.cc
/// @brief  report JSON object with information needed for vHTS postprocessing
/// @author Sam DeLuca

#include <protocols/features/ScreeningFeatures.hh>


#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/string_util.hh>
#include <utility/json_spirit/json_spirit_utils.h>
#include <utility/json_spirit/json_spirit_writer.h>
#include <algorithm>

#include <core/chemical/ResidueType.hh>
#include <utility/tools/make_vector.hh>
#include <boost/foreach.hpp>

namespace protocols {
namespace features {

using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

ScreeningFeatures::ScreeningFeatures()
{

}

ScreeningFeatures::ScreeningFeatures(ScreeningFeatures const & src) : protocols::features::FeaturesReporter(src),
		chain_(src.chain_),
		descriptors_(src.descriptors_)
{
}


ScreeningFeatures::~ScreeningFeatures()
{

}

std::string ScreeningFeatures::type_name() const
{
	return "ScreeningFeatures";
}

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
	
	jd2::JobDistributor* job_distributor = jd2::JobDistributor::get_instance();
	jd2::JobCOP current_job = job_distributor->current_job();
	std::string group_name;

	jd2::Job::StringStringPairs::const_iterator string_string_begin = current_job->output_string_string_pairs_begin();
	jd2::Job::StringStringPairs::const_iterator string_string_end = current_job->output_string_string_pairs_end();
	
	for(jd2::Job::StringStringPairs::const_iterator it = string_string_begin; it != string_string_end;++it)
	{
		if(it->first == "input_group_name")
		{
			group_name = it->second;
		}
	}

	RowDataBaseOP group_name_data( new RowData<std::string>("group_name",group_name) );

	core::Size chain_id = core::pose::get_chain_id_from_chain(chain_,pose);

	for(core::Size resnum= pose.conformation().chain_begin(chain_id); resnum <= pose.conformation().chain_end(chain_id); ++resnum)
	{
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

	if(!basic::options::option[basic::options::OptionKeys::in::file::screening_job_file].user())
	{
		throw utility::excn::EXCN_BadInput("The ScreeningFeatures reporter is only usable with the ScreeningJobInputter. specify input using -in:file:screening_job_inputter");
	}

	if(!tag->hasOption("chain"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("ScreeningFeatures requires the 'chain' tag");
	}

	chain_ = tag->getOption<std::string>("chain");

	core::Size descriptor_count = 0;
	BOOST_FOREACH( utility::tag::TagCOP const sub_tag, tag->getTags() )
	{
		if(sub_tag->getName() != "descriptor")
		{
			throw utility::excn::EXCN_RosettaScriptsOption("ScreeningFeatures only supports subtags with the name 'descriptor");
		}

		if(!sub_tag->hasOption("type"))
		{
			throw utility::excn::EXCN_RosettaScriptsOption("ScreeningFeatures descriptor subtags require a 'type' option");
		}

		std::string descriptor_type(sub_tag->getOption<std::string>("type"));
		descriptors_.push_back(descriptor_type);
		descriptor_count++;
	}
	if(descriptor_count == 0)
	{
		throw utility::excn::EXCN_RosettaScriptsOption("ScreeningFeatures requires at least one 'descriptor' subtag");
	}

}

std::vector<utility::json_spirit::Pair>  ScreeningFeatures::get_desriptor_data() const
{
	using jd2::Job;
	
	jd2::JobDistributor* job_distributor = jd2::JobDistributor::get_instance();
	jd2::JobCOP current_job = job_distributor->current_job();

	Job::StringStringPairs::const_iterator string_string_begin = current_job->output_string_string_pairs_begin();
	Job::StringStringPairs::const_iterator string_string_end = current_job->output_string_string_pairs_end();
	
	Job::StringRealPairs::const_iterator string_real_begin = current_job->output_string_real_pairs_begin();
	Job::StringRealPairs::const_iterator string_real_end = current_job->output_string_real_pairs_end();
	
	std::vector<utility::json_spirit::Pair>  descriptor_data;
	
	for(Job::StringStringPairs::const_iterator it = string_string_begin; it != string_string_end;++it)
	{
		if(std::find(descriptors_.begin(),descriptors_.end(),it->first) != descriptors_.end())
		{
			utility::json_spirit::Pair data_pair(it->first,it->second);
			descriptor_data.push_back(data_pair);
		}
	}
	
	for(Job::StringRealPairs::const_iterator it = string_real_begin; it != string_real_end;++it)
	{
		if(std::find(descriptors_.begin(),descriptors_.end(),it->first) != descriptors_.end())
		{
			utility::json_spirit::Pair data_pair(it->first,it->second);
			descriptor_data.push_back(data_pair);
		}
	}
	
	return descriptor_data;
	
	
	
}


}
}
