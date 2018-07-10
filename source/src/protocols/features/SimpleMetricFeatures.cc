// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/SimpleMetricFeatures.cc
/// @brief  Report a set of simple metrics into a features database.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit Headers
#include <protocols/features/SimpleMetricFeatures.hh>
#include <protocols/features/util.hh>
#include <core/simple_metrics/SimpleMetricFactory.hh>
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/CompositeRealMetric.hh>
#include <core/simple_metrics/CompositeStringMetric.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/PerResidueStringMetric.hh>
#include <core/simple_metrics/util.hh>

// Project Headers
#include <core/pose/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/excn/Exceptions.hh>

// Basic Headers
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>


// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>
#include <map>
#include <list>
#include <sstream>

#include <utility/vector0.hh>
#include <utility/string_util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/SimpleMetricFeaturesCreator.hh>

//cppdb External.
#include <cppdb/errors.h>


namespace protocols {
namespace features {

using core::pose::Pose;
using core::pose::PoseCOP;
using core::pose::PoseOP;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using utility::vector1;
using utility::sql_database::sessionOP;
using utility::tag::TagCOP;
using cppdb::statement;
using cppdb::cppdb_error;

using namespace core::simple_metrics;

static basic::Tracer TR( "protocols.features.SimpleMetricFeatures" );

// XRW TEMP string
// XRW TEMP SimpleMetricFeatures::type_name() const { return "SimpleMetricFeatures"; }



SimpleMetricFeatures::SimpleMetricFeatures():
	FeaturesReporter()
{}

SimpleMetricFeatures::SimpleMetricFeatures( vector1< SimpleMetricCOP > metrics ):
	FeaturesReporter()
{
	set_simple_metrics(metrics);
}

SimpleMetricFeatures::SimpleMetricFeatures(
	SimpleMetricFeatures const & src ) :
	FeaturesReporter(),
	prefix_(src.prefix_),
	suffix_(src.suffix_),
	table_name_(src.table_name_)
{
	metrics_ = src.metrics_;
}


void
SimpleMetricFeatures::add_simple_metric( SimpleMetricCOP metric ){
	metrics_.push_back( metric );
}

void
SimpleMetricFeatures::set_simple_metrics( vector1< SimpleMetricCOP > metrics ){
	metrics_ = metrics;
}

void
SimpleMetricFeatures::set_table_name( std::string table_name ){
	table_name_ = table_name;
}

void
SimpleMetricFeatures::set_prefix_suffix( std::string prefix, std::string suffix){
	prefix_ = prefix;
	suffix_ = suffix;
}

utility::vector1<std::string>
SimpleMetricFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	//dependencies.push_back("StructureFeatures");
	return dependencies;
}

void
SimpleMetricFeatures::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const &
) {
	runtime_assert(tag->getName() == type_name());

	vector1< SimpleMetricCOP > metrics = get_metrics_from_datamap_and_subtags(tag, data);
	set_simple_metrics(metrics);

	set_table_name( tag->getOption< std::string >("table_name", table_name_));
	set_prefix_suffix(
		tag->getOption< std::string >("prefix", prefix_),
		tag->getOption< std::string >("suffix", suffix_));
}

void
SimpleMetricFeatures::write_schema_to_db(
	sessionOP db_session
) const {

	write_general_schema_to_db( db_session );
	write_per_residue_schema_to_db( db_session );

}

void
SimpleMetricFeatures::write_general_schema_to_db(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	std::map< std::string, std::string > name_value_type;

	core::Size n_metrics = 0;
	for ( SimpleMetricCOP metric : metrics_ ) {
		utility::vector1< std::string > names = metric->get_metric_names();
		std::string  metric_class = metric->metric();

		std::string metric_type="";

		std::string metric_name = metric->simple_metric_type();
		std::string custom_type = metric->get_custom_type();
		if ( (metric_name == "PerResidueRealMetric" )|| ( metric_name == "PerResidueStringMetric") ) {
			continue;
		} else if ( (metric_name == "StringMetric" ) || ( metric_name == "CompositeStringMetric" ) ) {

			metric_type = "string";
			n_metrics+=1;

		} else if ( ( metric_name == "RealMetric") || ( metric_name =="CompositeRealMetric") ) {

			metric_type = "real";
			n_metrics+=1;

		} else {
			utility_exit_with_message("Unknown SimpleMetric Type " + metric_name);
		}

		for ( std::string m_name : names ) {
			if ( (metric_name == "CompositeStringMetric" )|| (metric_name == "CompositeRealMetric") ) {
				name_value_type[m_name+"_"+custom_type + metric_class] = metric_type;
			} else {
				name_value_type[custom_type  + m_name] = metric_type;
			}
		}
	}
	if ( n_metrics == 0 ) return;

	//Create the table using the set SimpleMetrics
	utility::vector1< Column > new_columns;

	Column struct_id_column("struct_id", DbDataTypeOP( new DbBigInt() ));
	new_columns.push_back(struct_id_column);

	//Primary Key
	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id_column);
	PrimaryKey primary_key(primary_key_columns);

	//Foreign Key
	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id_column);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	ForeignKey foreign_key(foreign_key_columns, "structures", reference_columns, true);
	Schema table(table_name_, primary_key);
	table.add_foreign_key(foreign_key);

	//C++ 11 we have ordered maps.
	for ( auto name_data : name_value_type ) {
		if ( name_data.second == "real" ) {
			Column new_column(name_data.first, DbDataTypeOP( new DbReal() ));
			new_columns.push_back( new_column );
		} else if ( name_data.second == "string" ) {
			Column new_column(name_data.first, DbDataTypeOP( new DbText() ));
			new_columns.push_back( new_column );
		}
	}

	Column prefix("prefix", DbDataTypeOP( new DbText()));
	Column suffix("suffix", DbDataTypeOP( new DbText()));
	table.add_column(prefix);
	table.add_column(suffix);
	for ( auto col : new_columns ) {
		table.add_column( col );
	}
	table.write(db_session);
}

void
SimpleMetricFeatures::write_per_residue_schema_to_db(utility::sql_database::sessionOP db_session) const {
	using namespace basic::database::schema_generator;

	//Make sure we have data to write.
	std::map< std::string, std::string > name_value_type;
	for ( SimpleMetricCOP metric : metrics_ ) {
		utility::vector1< std::string > const names = metric->get_metric_names();
		std::string metric_type="";
		std::string specific_metric_name = ""; //Ex SASA, RMSD, etc.
		std::string const metric_name = metric->simple_metric_type();

		if ( metric_name == "PerResidueRealMetric" ) {
			assert( names.size() == 1);
			specific_metric_name = names[1];
			metric_type = "real";
		} else if ( metric_name == "PerResidueStringMetric" ) {
			assert( names.size() == 1);
			specific_metric_name = names[1];
			metric_type = "string";
		} else {
			continue;
		}
		name_value_type[specific_metric_name] = metric_type;
	}

	if ( name_value_type.size() == 0 ) return; //No PerResidueMetrics.

	//Create the table using the set SimpleMetrics
	utility::vector1< Column > new_columns;

	Column struct_id_column("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column res_id_column("resNum", DbDataTypeOP( new DbInteger() ) );
	new_columns.push_back(struct_id_column);
	new_columns.push_back(res_id_column);

	//Primary Key
	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id_column);
	primary_key_columns.push_back(res_id_column);

	PrimaryKey primary_key(primary_key_columns);

	//Foreign Key
	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id_column);
	foreign_key_columns.push_back(res_id_column);

	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey struct_foreign_key(struct_id_column, "structures", "struct_id", true);
	ForeignKey res_foreign_key(foreign_key_columns, "residues", reference_columns, true);

	std::string const table_name = "per_residue_"+table_name_;

	Schema table(table_name, primary_key);
	table.add_foreign_key(struct_foreign_key);
	table.add_foreign_key(res_foreign_key);

	//C++ 11 we have ordered maps.
	for ( auto name_data : name_value_type ) {
		if ( name_data.second == "real" ) {
			Column new_column(name_data.first, DbDataTypeOP( new DbReal() ));
			new_columns.push_back( new_column );
		} else if ( name_data.second == "string" ) {
			Column new_column(name_data.first, DbDataTypeOP( new DbText() ));
			new_columns.push_back( new_column );
		}
	}

	Column prefix("prefix", DbDataTypeOP( new DbText()));
	Column suffix("suffix", DbDataTypeOP( new DbText()));
	table.add_column(prefix);
	table.add_column(suffix);
	for ( auto col : new_columns ) {
		table.add_column( col );
	}
	table.write(db_session);
}

core::Size
SimpleMetricFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	using namespace core::scoring;
	using namespace basic::database::schema_generator;

	report_general_features(      pose, relevant_residues, struct_id, db_session );
	report_per_residue_features(  pose, relevant_residues, struct_id, db_session );

	return 0;
}


core::Size
SimpleMetricFeatures::report_general_features(
	Pose const & pose,
	vector1< bool > const & ,
	StructureID struct_id,
	sessionOP db_session
){
	using namespace core::scoring;
	using namespace basic::database::schema_generator;

	//Calculate the data, figure out the names of and types of the columns.
	// We want to be able to work on any simple metric, but I don't know how to do this elegantly.

	std::map< std::string, MetricData > name_data_map;

	//This is needed as there is a bug in the struct where the text will be completely garbled if put into the database.
	// I think its some kind of memory issue, but regardless, storing the string data in this map works.
	// I'm keeping the original code, but this note exists as to why its bit funky here.
	std::map< std::string, std::string > string_data;


	for ( SimpleMetricCOP metric: metrics_ ) {

		std::string custom_type = metric->get_custom_type();
		if ( custom_type != "" ) custom_type=custom_type+"_";

		std::string const metric_type = metric->simple_metric_type();


		//Can't use a switch for strings for some reason.
		if      ( metric_type == "RealMetric" ) {
			RealMetric const & r_metric = dynamic_cast<RealMetric const & >( *metric );
			MetricData new_data = MetricData(r_metric.calculate(pose));
			name_data_map[custom_type+r_metric.metric() ] = new_data;

		} else if ( metric_type == "StringMetric" ) {
			StringMetric const & r_metric = dynamic_cast<StringMetric const & >( *metric );
			std::string value = r_metric.calculate(pose);
			MetricData new_data( value );
			name_data_map[custom_type+r_metric.metric() ] = new_data;
			string_data[custom_type+r_metric.metric() ] = value;
		} else if ( metric_type == "CompositeRealMetric" ) {
			CompositeRealMetric const & r_metric = dynamic_cast<CompositeRealMetric const & >( *metric );
			std::map< std::string, core::Real > result = r_metric.calculate(pose);
			for ( auto result_pair : result ) {
				MetricData new_data( result_pair.second);
				name_data_map[ result_pair.first+"_"+custom_type+r_metric.metric() ] = new_data;
			}
		} else if ( metric_type == "CompositeStringMetric" ) {
			CompositeStringMetric const & r_metric = dynamic_cast<CompositeStringMetric const & >( *metric );
			std::map< std::string, std::string > result = r_metric.calculate(pose);
			for ( auto result_pair : result ) {
				MetricData new_data( result_pair.second);
				name_data_map[ result_pair.first+"_"+custom_type+r_metric.metric() ] = new_data;
				string_data[ result_pair.first+"_"+custom_type+r_metric.metric() ] = result_pair.second;
			}
		} else if ( (metric_type == "PerResidueRealMetric") || (metric_type == "PerResidueStringMetric" ) ) {
			continue;
		} else {
			utility_exit_with_message("Metric Type not understood! "+ metric_type);
		}

	}

	if ( name_data_map.size() == 0 ) return 0; //No data.

	//Make the statement string.
	std::string statement_string = "INSERT INTO "+table_name_+"(struct_id,prefix,suffix,";
	core::Size n = 3;
	for ( auto data_pair : name_data_map ) {
		n+=1;
		statement_string = statement_string + data_pair.first;
		if ( n < name_data_map.size() + 3 ) {
			statement_string = statement_string + ",";
		}
	}
	statement_string = statement_string + ") VALUES "+get_question_mark_string(n);

	///Better error handling - so that if someone tries to do something funky, we catch it.
	statement stmt;
	try {
		stmt = basic::database::prepare_statement_no_catch(statement_string,db_session);
	}
catch (cppdb_error const & error){

	TR << std::endl << std::endl;
	TR.Error << "Input If you are calling ReportToDB on the same database, make sure either the table names are different, or the same exact simple metrics are being used for both calls!" << std::endl;
	TR << std::endl << " Here is the full error message: " << std::endl ;

	TR << "   Failed to safely prepare the following statement: " << std::endl;
	TR << "    " << statement_string << std::endl;
	TR << "    " << error.what() << std::endl;
	TR << std::endl;
	throw CREATE_EXCEPTION(utility::excn::BadInput, "Columns for successive ReportToDb must match.");
}

	//Add the data to the database.
	stmt.bind(1,struct_id);
	stmt.bind(2,prefix_);
	stmt.bind(3,suffix_);

	n=3;
	for ( auto data_pair : name_data_map ) {
		n+=1;
		//TR << data_pair.first << " " << data_pair.second.data_type << std::endl;
		if ( data_pair.second.data_type == "real" ) {
			core::Real value = data_pair.second.real_value;
			stmt.bind(n, value);
			//TR << data_pair.second.real_value << std::endl;
		} else if ( data_pair.second.data_type == "string" ) {
			stmt.bind(n, string_data[data_pair.first]);
			//TR << data_pair.second.string_value << std::endl;
		}

	}

	basic::database::safely_write_to_database(stmt);
	return 0;
}

core::Size
SimpleMetricFeatures::report_per_residue_features(
	Pose const & pose,
	vector1< bool > const & ,
	StructureID struct_id,
	sessionOP db_session
){
	using namespace core::scoring;
	using namespace basic::database::schema_generator;

	utility::vector1< std::string > data_names;
	std::map< std::string, std::map< core::Size, std::string > > string_data;
	std::map< std::string, std::map< core::Size, core::Real > > real_data;

	///Clone the metrics, wipeout their ResidueSelectors so that they are calculated on ALL residues.
	/// This is so that data is valid instead of having a "NA" for a string, which could totally be something real.
	for ( SimpleMetricCOP metric_original: metrics_ ) {

		SimpleMetricOP metric = metric_original->clone();

		std::string custom_type = metric->get_custom_type();
		if ( custom_type != "" ) custom_type=custom_type+"_";

		std::string const metric_type = metric->simple_metric_type();

		if ( metric_type == "PerResidueStringMetric" ) {
			PerResidueStringMetric & r_metric = dynamic_cast<PerResidueStringMetric & >( *metric );
			r_metric.set_residue_selector( nullptr ); //reset residue selector.

			//name_data_map[custom_type+r_metric.metric() ] = new_data;
			//string_data[ ] = value;
			std::map< core::Size, std::string > const data = r_metric.calculate(pose);
			std::string data_name = custom_type+r_metric.metric();
			data_names.push_back( data_name );
			string_data[ data_name ] = data;
		} else if ( metric_type == "PerResidueRealMetric" ) {
			PerResidueRealMetric & r_metric = dynamic_cast<PerResidueRealMetric & >( *metric );
			r_metric.set_residue_selector( nullptr );
			//MetricData new_data = MetricData(r_metric.calculate(pose));
			//name_data_map[custom_type+r_metric.metric() ] = new_data;

			std::map< core::Size, core::Real > const data = r_metric.calculate(pose);
			std::string const data_name = custom_type+r_metric.metric();
			data_names.push_back( data_name );
			real_data[ data_name ] = data;
		} else {
			continue;
		}
	}

	if ( data_names.size() == 0 ) return 0; //No PerResidue Data

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		write_per_residue_data_row( struct_id, i, data_names, string_data, real_data, db_session );
	}
	return 0;

}

void
SimpleMetricFeatures::write_per_residue_data_row(
	StructureID struct_id,
	core::Size resnum,
	utility::vector1< std::string > const & data_names,
	std::map< std::string, std::map< core::Size, std::string > > const & string_data,
	std::map< std::string, std::map< core::Size, core::Real > > const & real_data,
	sessionOP db_session )
{

	//Make the statement string.
	std::string statement_string = "INSERT INTO "+table_name_+"(struct_id,resNum,prefix,suffix,";
	core::Size current_row = 4;
	for ( auto name : data_names ) {
		current_row+=1;
		statement_string = statement_string + name;
		if ( current_row < data_names.size() + 4 ) {
			statement_string = statement_string + ",";
		}
	}
	statement_string = statement_string + ") VALUES "+get_question_mark_string( current_row );

	///Better error handling - so that if someone tries to do something funky, we catch it.
	statement stmt;
	try {
		stmt = basic::database::prepare_statement_no_catch(statement_string,db_session);
	}
catch (cppdb_error const & error){

	TR << std::endl << std::endl;
	TR.Error << "Input If you are calling ReportToDB on the same database, make sure either the table names are different, or the same exact simple metrics are being used for both calls!" << std::endl;
	TR << std::endl << " Here is the full error message: " << std::endl ;

	TR << "   Failed to safely prepare the following statement: " << std::endl;
	TR << "    " << statement_string << std::endl;
	TR << "    " << error.what() << std::endl;
	TR << std::endl;
	throw CREATE_EXCEPTION(utility::excn::BadInput, "Columns for successive ReportToDb must match.");
}

	//Add the data to the database.
	stmt.bind(1,struct_id);
	stmt.bind(2, resnum);
	stmt.bind(3,prefix_);
	stmt.bind(4,suffix_);

	current_row=4;
	for ( auto name : data_names ) {
		current_row+=1;
		if ( real_data.count(name) ) {
			core::Real data = real_data.at(name).at(resnum);
			stmt.bind(current_row, data );

		} else if ( string_data.count(name) ) {
			std::string data = string_data.at(name).at(resnum);
			stmt.bind( current_row, data );
		}
	}
	basic::database::safely_write_to_database(stmt);
}



std::string SimpleMetricFeatures::type_name() const {
	return class_name();
}

std::string SimpleMetricFeatures::class_name() {
	return "SimpleMetricFeatures";
}

void SimpleMetricFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	SimpleMetricFactory::get_instance()->define_simple_metric_xml_schema( xsd );

	attlist
		+ XMLSchemaAttribute( "metrics", xs_string, "Comma-separated list of previously defined simple_metrics to be added." )
		+ XMLSchemaAttribute( "prefix", xs_string, "Prefix tag for the data.  Added as an extra column in the data. " )
		+ XMLSchemaAttribute( "suffix", xs_string, "suffix tag for the data.  Added as an extra column in the data. " )
		+ XMLSchemaAttribute::attribute_w_default( "table_name", xs_string, "The table to add metrics to.  Should match the same exact data that you had before if the table has already been created!", "simple_metrics" );


	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & SimpleMetricFactory::get_instance()->simple_metric_xml_schema_group_name );

	protocols::features::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "Run a set of SimpleMetrics on the pose, and put the data into a features table ", attlist, subelements );

}

std::string SimpleMetricFeaturesCreator::type_name() const {
	return SimpleMetricFeatures::class_name();
}

protocols::features::FeaturesReporterOP
SimpleMetricFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new SimpleMetricFeatures );
}

void SimpleMetricFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SimpleMetricFeatures::provide_xml_schema( xsd );
}


} //namesapce
} //namespace
