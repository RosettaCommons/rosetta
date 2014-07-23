// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/TaskOperationFeatures.hh
/// @brief  report
/// @author Kevin Houlihan (khouli@unc.edu)

#ifndef INCLUDED_protocols_features_TaskOperationFeatures_hh
#define INCLUDED_protocols_features_TaskOperationFeatures_hh

//External

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pack/task/TaskFactory.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>
#include <vector>

#include <boost/unordered_map.hpp>

#include <devel/metawidget/MetaWidget.hh>
//#include <WidgetTypes.hh>

namespace protocols{
namespace features{

//class TaskOperationFeatures : public protocols::features::FeaturesReporter {
class TaskOperationFeatures : public MetaWidget::MetaWidget<MetaWidget::FeaturesReporter, TaskOperationFeatures> {
public:
#include <devel/metawidget/WidgetTypes.hh>
	typedef MetaWidget<FeaturesReporter, TaskOperationFeatures> Base;
	typedef Base::RegType RegType;
	static const char* name;

	TaskOperationFeatures();
	~TaskOperationFeatures();

	TaskOperationFeatures( TaskOperationFeatures const & src );

	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	void
	write_task_operations_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_task_operation_residue_effects_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	core::Size
	report_task_operations(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	core::Size
	report_task_operation_effects(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_task_operations_row(
		core::Size const taskop_id,
		std::string const& type_name,
		std::string const& user_name,
		utility::sql_database::sessionOP db_session);

	void
	insert_task_operation_residue_effects_row(
		StructureID const struct_id,
		core::Size const resNum,
		core::Size const taskop_id,
		bool pack,
		bool design,
		utility::sql_database::sessionOP db_session);

private:
	typedef core::pack::task::operation::TaskOperationCOP TaskOperationCOP;
	struct named_taskop {
		named_taskop(const std::string& type_name, const std::string& user_name, TaskOperation* to) : type_name(type_name), user_name(user_name), taskopOP(to) {}
		std::string type_name;
		std::string user_name;
		TaskOperationCOP taskopOP;
	};
	std::vector<named_taskop> named_taskops;

};

} // namespace
} // namespace

#endif
