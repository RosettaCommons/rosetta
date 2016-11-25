// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/TaskOperationFeatures.cc
/// @brief  report
/// @author Kevin Houlihan (khouli@unc.edu)

// Unit Headers
#include <protocols/features/TaskOperationFeatures.hh>

//External
#include <boost/assign/list_of.hpp>

// Project Headers
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/assert.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector0.hh>
#include <utility/string_util.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/assign/list_of.hpp>

// C++ Headers
#include <cmath>
#include <utility/excn/Exceptions.hh>
#include <algorithm>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/TaskOperationFeaturesCreator.hh>


//Auto Headers


namespace protocols {
namespace features {

using std::string;
using std::stringstream;
using std::endl;
using core::Size;
using core::Real;
using core::Vector;
using core::pack::task::PackerTaskCOP;
using core::pose::Pose;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using utility::sql_database::sessionOP;
using core::pack::task::TaskFactoryOP;
using core::pack::task::TaskFactoryCOP;
using core::pack::task::TaskFactory;
using core::pack::task::operation::TaskOperation;
using utility::tag::TagCOP;
using cppdb::statement;
using utility::vector1;
using basic::Tracer;
using basic::database::insert_or_ignore;
using utility::tag::TagPtr;
using boost::assign::list_of;

static Tracer TR("protocols.features.TaskOperationFeatures");

TaskOperationFeatures::TaskOperationFeatures()
// : run_once_(true)
{}

TaskOperationFeatures::TaskOperationFeatures(TaskOperationFeatures const &) :
	FeaturesReporter()
	// run_once_(src.run_once_)
{}

TaskOperationFeatures::~TaskOperationFeatures() = default;

// XRW TEMP string
// XRW TEMP TaskOperationFeatures::type_name() const { return "TaskOperationFeatures"; }

void
TaskOperationFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	TR << "Writing schemas" << std::endl;
	write_task_operations_table_schema(db_session);
	write_task_operation_residue_effects_table_schema(db_session);
}

void
TaskOperationFeatures::write_task_operations_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	TR << "Writing task_operations schema" << std::endl;
	Column taskop_id("taskop_id", DbDataTypeOP( new DbInteger() ));
	Column taskop_name("taskop_name", DbDataTypeOP( new DbText() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(taskop_id);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("task_operations", primary_key);

	table.add_column(taskop_name);

	table.write(db_session);
}

void
TaskOperationFeatures::write_task_operation_residue_effects_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	TR << "Writing task_operation_residue_effects schema" << std::endl;
	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum("resNum", DbDataTypeOP( new DbInteger() ));
	Column taskop_id("taskop_id", DbDataTypeOP( new DbInteger() ));
	Column pack("pack", DbDataTypeOP( new DbInteger() ));
	Column design("design", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	primary_key_columns.push_back(taskop_id);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(resNum);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(taskop_id);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("taskop_id");
	ForeignKey foreign_key2(foreign_key_columns2, "task_operations", reference_columns2, true);

	Schema table("task_operation_residue_effects", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);

	table.add_column(pack);
	table.add_column(design);

	table.write(db_session);
}

utility::vector1<std::string>
TaskOperationFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

void
TaskOperationFeatures::parse_my_tag(
	TagCOP tag,
	DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	if ( ! tag->hasOption("task_operations") ) {
		stringstream error_msg;
		error_msg
			<< "The " << type_name() << " reporter requires a 'task_operations' tag:" << endl
			<< endl
			<< "    <feature name=" << type_name() <<" task_operations=comma_separated_task_operation_names />" << endl;
		throw utility::excn::EXCN_RosettaScriptsOption(error_msg.str());
	}

	// Fill task_factories_, one TaskFactoryCOP per TaskOperator
	Size taskop_id(0);
	std::string taskop_list = tag->getOption< std::string >("task_operations");
	utility::vector0< std::string > const taskop_keys( utility::string_split( taskop_list, ',' ) );
	for ( auto const & taskop_key : taskop_keys ) {
		taskop_id++;
		if ( data.has( "task_operations", taskop_key ) ) {
			TaskFactoryOP new_task_factory( new TaskFactory );
			new_task_factory->push_back( data.get_ptr< TaskOperation >
				( "task_operations", taskop_key ) );
			taskops_.push_back(Taskop_id_name_factory_(taskop_id, taskop_key, new_task_factory) );
			//   taskop_keys_factories_.insert(
			//    std::pair<std::string const, TaskFactoryCOP>(*taskop_key, new_task_factory));
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + taskop_key + " not found in DataMap.");
		}
	}
}

Size
TaskOperationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	// assert pose.update_residue_neighbors() has been called:
	runtime_assert(
		!pose.conformation().structure_moved() &&
		pose.energies().residue_neighbors_updated());

	// if(taskop_keys_factories_.empty()) utility_exit_with_message("List of TaskOperators"
	//   " given to TaskFeatures is empty.");
	if ( taskops_.empty() ) {
		utility_exit_with_message("List of TaskOperators"
			" given to TaskFeatures is empty.");
	}

	// task_operations table rows, only run once
	// if (run_once_) {
	TR << "Inserting TaskOperations rows" << std::endl;
	for ( vector1<Taskop_id_name_factory_>::const_iterator
			top = taskops_.begin(),
			top_end = taskops_.end(); top != top_end; ++top ) {
		insert_task_operations_row(top->id, top->name, db_session);
	}
	//  run_once_ = false;
	// } else {
	//  TR << "Skipping task_operations rows" << std::endl;
	// }

	TR << "Generating tasks from saved task factories" << std::endl;
	// gets tasks
	std::map<core::Size, PackerTaskCOP> tasks;
	for ( vector1<Taskop_id_name_factory_>::const_iterator
			top = taskops_.begin(),
			top_end = taskops_.end(); top != top_end; ++top ) {
		//Size taskop_id = top->id;
		tasks[top->id] = (top->tf)->create_task_and_apply_taskoperations( pose );
	}

	TR << "Inserting task_operation_residue_effects rows" << std::endl;
	// task_operation_residue_effects rows
	for ( Size resNum = 1; resNum <= pose.size(); ++resNum ) {
		if ( !check_relevant_residues( relevant_residues, resNum ) ) continue;

		for ( auto & task : tasks ) {
			bool const pack = (task.second)->pack_residue(resNum);
			bool const design = (task.second)->design_residue(resNum);
			Size taskop_id = task.first;
			insert_task_operation_residue_effects_row(struct_id, resNum,
				taskop_id, pack, design, db_session);
		}
	}

	return 0;
}

void
TaskOperationFeatures::insert_task_operations_row(
	Size const taskop_id,
	std::string const & taskop_name,
	sessionOP db_session
){
	TR << "writing task_operation_residue_effects row" << std::endl;

	string statement_string;
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3 :
		statement_string = "INSERT OR IGNORE INTO task_operations"
			" (taskop_id, taskop_name) VALUES (?,?);";
		break;
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres :
		statement_string = "INSERT IGNORE INTO task_operations"
			" (taskop_id, taskop_name) VALUES (?,?);";
		break;
	default :
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}

	statement stmt(basic::database::safely_prepare_statement(statement_string, db_session));
	stmt.bind(1,taskop_id);
	stmt.bind(2,taskop_name);
	basic::database::safely_write_to_database(stmt);

}

void
TaskOperationFeatures::insert_task_operation_residue_effects_row(
	StructureID const struct_id,
	Size const resNum,
	Size const taskop_id,
	bool pack,
	bool design,
	sessionOP db_session
){
	TR << "writing task_operation_residue_effects row" << std::endl;

	std::string statement_string = "INSERT INTO task_operation_residue_effects"
		" (struct_id, resNum, taskop_id, pack, design) VALUES (?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string, db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,resNum);
	stmt.bind(3,taskop_id);
	stmt.bind(4,pack);
	stmt.bind(5,design);
	basic::database::safely_write_to_database(stmt);
}

std::string TaskOperationFeatures::type_name() const {
	return class_name();
}

std::string TaskOperationFeatures::class_name() {
	return "TaskOperationFeatures";
}

void TaskOperationFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute(
		"task_operation", xs_string,
		"comma separated TaskOperation list");

	protocols::features::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"TaskOperation feature reporter",
		attlist );
}

std::string TaskOperationFeaturesCreator::type_name() const {
	return TaskOperationFeatures::class_name();
}

protocols::features::FeaturesReporterOP
TaskOperationFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new TaskOperationFeatures );
}

void TaskOperationFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TaskOperationFeatures::provide_xml_schema( xsd );
}


} // namesapce
} // namespace
