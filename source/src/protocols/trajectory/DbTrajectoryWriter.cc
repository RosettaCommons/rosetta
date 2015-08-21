// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Headers {{{1
#include <protocols/trajectory/DbTrajectoryWriter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>

// }}}1

namespace protocols {
namespace trajectory {

// Global Names {{{1
using namespace std;

using core::Real;
using core::Size;
using core::pose::Pose;

using utility::vector1;
using utility::tools::make_vector;
using utility::tools::make_vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

static thread_local basic::Tracer TR( "protocols.trajectory.DbTrajectoryWriter" );
// }}}1

DbTrajectoryWriter::DbTrajectoryWriter( // {{{1
	Size job_id, Pose const & pose, Size frequency, Size cache_limit) {

	job_id_ = job_id;
	iteration_ = 0;
	frequency_ = frequency;
	cache_limit_ = cache_limit;

	write_schema_to_db();
	update(pose);
}

void DbTrajectoryWriter::set_frequency(Size setting) { // {{{1
	frequency_ = setting;
}

void DbTrajectoryWriter::set_cache_limit(Size setting) { // {{{1
	cache_limit_ = setting;
}

// {{{1
/// @details This method should be called on every iteration.  If you don't
/// want to actually store a frame every iteration, the optional frequency
/// argument to the constructor can be used to specify how often frames should
/// be recorded.
void DbTrajectoryWriter::update(Pose const & pose) {
	if ( iteration_ % frequency_ == 0 ) {
		Frame frame;
		frame.iteration = iteration_;
		frame.pose = pose;
		frame_cache_.push_back(frame);
	}

	iteration_ += 1;

	if ( frame_cache_.size() >= cache_limit_ ) {
		write_cache_to_db();
	}

}

void DbTrajectoryWriter::finalize() const { // {{{1
	write_cache_to_db();
}

void DbTrajectoryWriter::write_schema_to_db() const { // {{{1
	using utility::sql_database::sessionOP;
	using namespace basic::database::schema_generator;

	sessionOP db_session = basic::database::get_db_session();

	Column job_id("job_id", DbDataTypeOP( new DbBigInt() ), false);
	Column iteration("iteration", DbDataTypeOP( new DbBigInt() ), false);
	Column score("score", DbDataTypeOP( new DbReal() ), false);
	Column silent_pose("silent_pose", DbDataTypeOP( new DbText() ), false);

	PrimaryKey composite_key(make_vector1(job_id, iteration));
	Schema trajectories("trajectories", composite_key);

	trajectories.add_column(score);
	trajectories.add_column(silent_pose);
	trajectories.write(db_session);
}

void DbTrajectoryWriter::write_cache_to_db() const { // {{{1
	using utility::sql_database::sessionOP;
	using namespace core::io::silent;
	using namespace basic::database::insert_statement_generator;

	sessionOP db_session = basic::database::get_db_session();

	InsertGenerator trajectory_insert("trajectories");
	trajectory_insert.add_column("job_id");
	trajectory_insert.add_column("iteration");
	trajectory_insert.add_column("score");
	trajectory_insert.add_column("silent_pose");

	RowDataBaseOP job( new RowData<Size>("job_id", job_id_) );

	BOOST_FOREACH ( Frame frame, frame_cache_ ) {
		stringstream string_stream;
		SilentFileData silent_file;
		SilentStructOP silent_data( new BinarySilentStruct(frame.pose, "db") );
		silent_file._write_silent_struct(*silent_data, string_stream);

		RowDataBaseOP iteration( new RowData<Size>(
			"iteration", frame.iteration) );
		RowDataBaseOP score( new RowData<Real>(
			"score", frame.pose.energies().total_energy()) );
		RowDataBaseOP silent_pose( new RowData<string>(
			"silent_pose", string_stream.str()) );

		trajectory_insert.add_row(make_vector(job, iteration, score, silent_pose));
	}

	frame_cache_.clear();
	trajectory_insert.write_to_database(db_session);
}
// }}}1

} // namespace trajectory
} // namespace protocols
