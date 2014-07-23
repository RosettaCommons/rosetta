// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/pack/task/TaskFactory.hh>
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
#include <boost/foreach.hpp>
//#include <boost/unordred_set.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <cmath>
#include <utility/excn/Exceptions.hh>
#include <algorithm>
#include <iostream>


namespace protocols{
namespace features{
static basic::Tracer TR("protocols.features.TaskOperationFeatures");
const char* TaskOperationFeatures::name = "TaskOperationFeatures";
namespace{TaskOperationFeatures::RegType SelfRegistration;}

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
using cppdb::statement;
using utility::vector1;
using basic::Tracer;
using basic::database::insert_or_ignore;
using boost::assign::list_of;


TaskOperationFeatures::TaskOperationFeatures()
{}
//	: run_once_(true)

TaskOperationFeatures::TaskOperationFeatures(TaskOperationFeatures const &)
{}
//	run_once_(src.run_once_)

TaskOperationFeatures::~TaskOperationFeatures() {}

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
	Column taskop_id("taskop_id", new DbInteger());
	Column taskop_typename("taskop_typename", new DbText());
	Column taskop_username("taskop_username", new DbText());

	Columns primary_key_columns;
	primary_key_columns.push_back(taskop_id);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("task_operations", primary_key);

	table.add_column(taskop_typename);
	table.add_column(taskop_username);

	table.write(db_session);
}

void
TaskOperationFeatures::write_task_operation_residue_effects_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	TR << "Writing task_operation_residue_effects schema" << std::endl;
	Column struct_id("struct_id", new DbBigInt());
	Column resNum("resNum", new DbInteger());
	Column taskop_id("taskop_id", new DbInteger());
	Column pack("packable", new DbInteger());
	Column design("designable", new DbInteger());

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
TaskOperationFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

void
TaskOperationFeatures::parse_my_tag(
    utility::tag::TagCOP tag,
    basic::datacache::DataMap& data,
    const protocols::filters::Filters_map&,
    const protocols::moves::Movers_map&,
    const core::pose::Pose&
) {
    // build map of TaskOperation type name to user name given in the Rosetta script
    TagCOP root = tag->getParent();
    std::string tagName = root->getName();
    // move up the tag tree until we get to the root
    while (tagName != "ROSETTASCRIPTS" && tagName != "dock_design") {
        root = root->getParent();
        tagName = root->getName();
    }

    TR << "Got root tag, name = " << root->getName() << std::endl;
    utility::vector0<TagCOP> root_subtags = root->getTags();
    TR << "The root tag has " << root_subtags.size() << " subtags, diving in" << std::endl;
    // from the root of the tag tree, find the TASKOPERATIONS tag
    boost::unordered_map<std::string, std::string> named_taskop_types;
    for (utility::vector0<TagCOP>::const_iterator it = root_subtags.begin(), root_subtags_end = root_subtags.end();
            it != root_subtags_end; ++it) {
        TR << "Looking for TASKOPERATIONS, inspecting tag: " << (*it)->getName() << std::endl;
        if ((*it)->getName() == "TASKOPERATIONS") {
            TR << "Found TASKOPERATIONS tag: " << (*it)->getName() << std::endl;
            utility::vector0<TagCOP> taskop_subtags = (*it)->getTags();
            for (utility::vector0<TagCOP>::const_iterator jt = taskop_subtags.begin(), taskop_subtags_end = taskop_subtags.end();
                    jt != taskop_subtags_end; ++jt) {
                TagCOP taskopTag = *jt;
                const std::string taskoptype = taskopTag->getName();
                const std::string uname = taskopTag->getOption<std::string>("name");
                TR << "1) taskop type: " << taskoptype << std::endl;
                TR << "1) taskop user name: " << uname << std::endl;
                // map user provided name to type name
                named_taskop_types[uname] = taskoptype;
            }
        }
    }

    const std::string task_op_opt = "task_operations";

    // if an explicit list of task ops was given, find the given task ops and only the given task ops
    if (tag->hasOption(task_op_opt)) {
        utility::vector1<std::string> requested_taskops = utility::string_split(tag->getOption<std::string>(task_op_opt), ',');
        for (utility::vector1<std::string>::const_iterator it = requested_taskops.begin(), end = requested_taskops.end(); it != end; ++it) {
            const std::string user_name = *it;
            const std::string type_name = named_taskop_types[user_name];
            TaskOperation* top = data.get<core::pack::task::operation::TaskOperation*>( "task_operations", user_name );
            named_taskops.push_back(named_taskop(type_name, user_name, top));
        }
    }
    // if not explicit list of task ops was given, find them all
    else {
        TR << "iterating over DataMap" << std::endl;
        typedef basic::datacache::DataMap DataMap;
        typedef std::map<std::string, utility::pointer::ReferenceCountOP> SubDataMap;
        for (DataMap::const_iterator it = data.begin(); it != data.end(); ++it) {
            if (it->first == "task_operations") {
                const SubDataMap& submap = it->second;
                for (SubDataMap::const_iterator jt = submap.begin(); jt != submap.end(); ++jt) {
                    const std::string user_name = jt->first;
                    const std::string type_name = named_taskop_types[user_name];
                    // TR << "\tKey: " << user_name << '\n';
                    named_taskops.push_back( named_taskop(type_name, user_name, dynamic_cast<TaskOperation*>(jt->second.get())) );
                }
                break;
            }
        }
    }

    // get these from functions
    // std::vector<std::string> requested_taskops = requested_taskops();
    // std::vector<std::string> available_taskops = available_taskops();

    // Debugger
    /*
    typedef utility::tag::Tag::options_t options_t;
    const options_t options_map = tag->getOptions();
    TR << "iterating over options in my tag\n";
    for (options_t::const_iterator it = options_map.begin(); it != options_map.end(); ++it) {
        TR << "option: first = " << it->first << ", second = " << it->second << std::endl;
        // swallow extra options, cast to string will always be accepted
        tag->getOption<std::string>(it->first);
    }
     */
    // Debugger

    // The above code qcquires a mapping of user name to TaskOperation* but the original
    // derived types of the those TaskOperation classes are unknown.

    // Work back up the tag system to find the original name, map it


    for (std::vector<named_taskop>::const_iterator it = named_taskops.begin(), end = named_taskops.end();
            it != end; ++it) {
        TR << "Named taskop" << std::endl;
        TR << '\t' << "user provided name: " << it->user_name << std::endl;
        TR << '\t' << "Rosetta type: " << it->type_name << std::endl;
    }

    //data_taskops = data.get< core::pack::task::operation::TaskOperation * >( "task_operations", *t_o_key );

    // TR << "Adding the following task operations\n";
    // for ( vec1_str::const_iterator to_it = data_taskops.begin(), end = t_o_keys.end();
    //         to_it != end; ++to_it ) {
    //     if ( data.has( "task_operations", *t_o_key ) ) {
    //         TR << *t_o_key << ' ';
    //         named_taskops.push_back( named_taskop(*t_o_key, data.get<TaskOperation*>("task_operations", *t_o_key)) );
    //     }
    //     else {
    //         throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + *t_o_key + " not found in basic::datacache::DataMap.");
    //     }
    // }
    // TR << std::endl;
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

	if (named_taskops.empty()) {
        TR << "No TaskOperations found, none will be reported" << std::endl;
        return 0;
    }
    report_task_operations(pose, relevant_residues, struct_id, db_session);
    report_task_operation_effects(pose, relevant_residues, struct_id, db_session);

    // task_operations table describes the taskoperations themselves
    // taskop_ids the indices of whatever order the taskops happen to be in

	return 0;
}

Size
TaskOperationFeatures::report_task_operations(
	Pose const&,
	vector1< bool > const&,
	StructureID,
	sessionOP db_session
){
    core::Size task_op_id = 1;
    for (std::vector<named_taskop>::const_iterator it = named_taskops.begin(), end = named_taskops.end(); it != end; ++it) {
        insert_task_operations_row(task_op_id, it->type_name, it->user_name, db_session);
        ++task_op_id;
    }
    return 0;
}

Size
TaskOperationFeatures::report_task_operation_effects(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
    // how the task operations affect individual residues
    core::Size taskop_id = 1;
    for (std::vector<named_taskop>::const_iterator it = named_taskops.begin(), end = named_taskops.end(); it != end; ++it) {
        TaskFactory tf;
        tf.push_back(it->taskopOP);
        PackerTaskCOP task = tf.create_task_and_apply_taskoperations(pose);
        for (core::Size resi = 1, end = pose.n_residue() + 1; resi != end; ++resi) {
            if (!relevant_residues[resi])
                continue;
            bool pack = task->pack_residue(resi);
            bool design = task->design_residue(resi);
            insert_task_operation_residue_effects_row(struct_id, resi, taskop_id, pack, design, db_session);
        }
        ++taskop_id;
    }
    return 0;
}

void
TaskOperationFeatures::insert_task_operations_row(
	Size const taskop_id,
	std::string const& type_name,
	std::string const& user_name,
	sessionOP db_session
){
	string statement_string;
	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3:
		statement_string = "INSERT OR IGNORE INTO task_operations"
				" (taskop_id, taskop_typename, taskop_username) VALUES (?,?,?);";
		break;
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres:
		statement_string = "INSERT IGNORE INTO task_operations"
				" (taskop_id, taskop_typename, taskop_username) VALUES (?,?,?);";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}

	statement stmt(basic::database::safely_prepare_statement(statement_string, db_session));
	stmt.bind(1,taskop_id);
	stmt.bind(2,type_name);
	stmt.bind(3,user_name);
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
	std::string statement_string = "INSERT INTO task_operation_residue_effects"
		" (struct_id, resNum, taskop_id, packable, designable) VALUES (?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string, db_session));
	stmt.bind(1,struct_id);
	stmt.bind(2,resNum);
	stmt.bind(3,taskop_id);
	stmt.bind(4,pack);
	stmt.bind(5,design);
	basic::database::safely_write_to_database(stmt);
}

} // namesapce
} // namespace
