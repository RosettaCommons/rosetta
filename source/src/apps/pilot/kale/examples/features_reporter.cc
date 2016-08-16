// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers {{{1
// Core headers
#include <core/pose/Pose.hh>

// Protocol Headers
#include <protocols/moves/Mover.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ReportToDB.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>
#include <basic/database/sql_utils.hh>
#include <boost/noncopyable.hpp>
#include <devel/init.hh>

// C++ headers
#include <string>

// Global Names {{{1
using namespace std;
using namespace basic::options;

using core::Size;
using core::pose::Pose;
using core::import_pose::pose_from_file;
using protocols::moves::Mover;
using protocols::features::ReportToDB;
using protocols::features::ReportToDBOP;
using protocols::features::StructureID;
using utility::vector1;
using utility::tools::make_vector;
using utility::tools::make_vector1;
using utility::pointer::owning_ptr;
// }}}1

class MyFeaturesReporter; // {{{1
typedef owning_ptr<MyFeaturesReporter> MyFeaturesReporterOP;
typedef owning_ptr<MyFeaturesReporter const> MyFeaturesReporterCOP;

class MyFeaturesReporter
	: public protocols::features::FeaturesReporter,
	  private boost::noncopyable {

public:

	string type_name() const {
		return "MyFeaturesReporter";
	}

	vector1<string> features_reporter_dependencies() const {
		vector1<string> dependencies;
		return dependencies;
	}

	virtual void write_schema_to_db(
			utility::sql_database::sessionOP db_session) const;

	Size report_features(
			Pose const & pose,
			utility::vector1<bool> const & relevant_residues,
			StructureID struct_id,
			utility::sql_database::sessionOP db_session);

};

void MyFeaturesReporter::write_schema_to_db( // {{{1
		utility::sql_database::sessionOP db_session) const {

	using namespace basic::database::schema_generator;

	Column id("id", new DbBigInt(), false);
	Column my_number("my_number", new DbBigInt(), false);
	Column my_string("my_string", new DbText(), false);
	Columns key_columns = make_vector1(id);

	Schema trajectories("my_table", PrimaryKey(key_columns));
	trajectories.add_column(my_number);
	trajectories.add_column(my_string);
	trajectories.write(db_session);
}

Size MyFeaturesReporter::report_features( // {{{1
		Pose const & /*pose_orig*/,
		vector1<bool> const & /*relevant_residues*/,
		StructureID id_value,
		utility::sql_database::sessionOP db_session) {

	using namespace basic::database::insert_statement_generator;

	InsertGenerator trajectory_insert("my_table");
	trajectory_insert.add_column("id");
	trajectory_insert.add_column("my_number");
	trajectory_insert.add_column("my_string");

	RowDataBaseOP id = new RowData<StructureID>("id", id_value);
	RowDataBaseOP my_number = new RowData<Size>("my_number", 42);
	RowDataBaseOP my_string = new RowData<string>("my_string", "hello world");

	trajectory_insert.add_row(make_vector(id, my_number, my_string));
	trajectory_insert.write_to_database(db_session);
	cout << "report_features()" << endl;
	return 0;
}

// }}}1

class MyApp; // {{{1

typedef utility::pointer::owning_ptr<MyApp> MyAppOP;
typedef utility::pointer::owning_ptr<MyApp const> MyAppCOP;

class MyApp : public Mover {

public:
	void apply(Pose & pose);
	string get_name() const { return "My Application"; }
	MoverOP fresh_instance() const { return new MyApp(); }

};

void MyApp::apply(Pose & pose) { // {{{1
	ReportToDBOP report_to_db = new ReportToDB;
	report_to_db->set_batch_name("MyExampleDatabase");
	report_to_db->set_batch_description("Just an example program.");
	report_to_db->add_features_reporter(new MyFeaturesReporter);
	report_to_db->apply(pose);
}
// }}}1

// To run this program, a few options need to be specified:
// $ rr && rx -- -out:use_database -inout:dbms:database_name traj.db
// $ sqlite3 -header -column traj.db 'select * from table my_table'

int main(int argc, char** argv) {
	devel::init(argc, argv);
	protocols::jd2::JobDistributor::get_instance()->go(new MyApp);
}
