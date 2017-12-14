// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/TaskOperationFeatures.hh
/// @brief  report
/// @author Kevin Houlihan (khouli@unc.edu)

#ifndef INCLUDED_protocols_features_TaskOperationFeatures_hh
#define INCLUDED_protocols_features_TaskOperationFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/TaskOperationFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pack/task/TaskFactory.hh>

// Utility Headers
#include <utility>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace features {

class TaskOperationFeatures : public protocols::features::FeaturesReporter {
public:
	TaskOperationFeatures();

	TaskOperationFeatures( TaskOperationFeatures const & src );

	~TaskOperationFeatures() override;

	/// @brief return string with class name

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const override;

private:
	void
	write_task_operations_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_task_operation_residue_effects_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/) override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

	void
	insert_task_operations_row(
		core::Size const taskop_id,
		std::string const & taskop_name,
		utility::sql_database::sessionOP db_session);

	void
	insert_task_operation_residue_effects_row(
		StructureID const struct_id,
		core::Size const resNum,
		core::Size const taskop_id,
		bool pack,
		bool design,
		utility::sql_database::sessionOP db_session);

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool run_once_;

	//std::map<std::string, core::pack::task::TaskFactoryCOP> taskop_keys_factories_;

	struct Taskop_id_name_factory_ {
		Taskop_id_name_factory_(
			Size i,
			std::string const & n,
			core::pack::task::TaskFactoryCOP t
		) : id(i), name(n), tf(t) {}
		core::Size id;
		std::string name;
		core::pack::task::TaskFactoryCOP tf;
	};

	utility::vector1<Taskop_id_name_factory_> taskops_;

};

} // namespace
} // namespace

#endif
