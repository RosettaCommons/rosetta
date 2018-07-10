// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/SimpleMetricFeatures.hh
/// @brief  Report a set of simple metrics into a features database.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_features_SimpleMetricFeatures_hh
#define INCLUDED_protocols_features_SimpleMetricFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/SimpleMetricFeatures.fwd.hh>
#include <core/simple_metrics/SimpleMetric.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace features {

struct MetricData {

	core::Real  real_value;
	std::string string_value;

	std::string data_type;

	MetricData(){}

	MetricData(core::Real r_value):
		real_value(r_value),
		data_type("real"){}

	MetricData(std::string r_value):
		string_value(r_value),
		data_type("string"){}

};

class SimpleMetricFeatures : public FeaturesReporter {
public:

	SimpleMetricFeatures();

	SimpleMetricFeatures(
		utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics
	);

	SimpleMetricFeatures( SimpleMetricFeatures const & src );

	~SimpleMetricFeatures() override= default;

public:

	/// @brief return string with class name

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const override;

public:

	///@brief Add a simple_metric
	void
	add_simple_metric( core::simple_metrics::SimpleMetricCOP metric );

	///@brief Set the simple metrics we will run.
	void
	set_simple_metrics( utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics );

	///@brief Set extra columns, prefix and suffix, for each data added to the database.  Used for multiple runs of the features.
	void
	set_prefix_suffix( std::string prefix="", std::string suffix="");

	void
	set_table_name( std::string = "simple_metrics");

public:

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose) override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;



	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	void
	write_general_schema_to_db(
		utility::sql_database::sessionOP db_session) const ;

	void
	write_per_residue_schema_to_db(
		utility::sql_database::sessionOP db_session) const ;

	/// @brief collect all the feature data for the pose
	core::Size
	report_general_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) ;

	core::Size
	report_per_residue_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) ;

	void
	write_per_residue_data_row(
		StructureID struct_id,
		core::Size resnum,
		utility::vector1< std::string > const & data_names,
		std::map< std::string, std::map< core::Size, std::string > > const & string_data,
		std::map< std::string, std::map< core::Size, core::Real > > const & real_data,
		utility::sql_database::sessionOP db_session );

private:

	utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics_;
	std::string prefix_="";
	std::string suffix_="";
	std::string table_name_="simple_metrics";

};

} // features namespace
} // protocols namespace

#endif // include guard
