// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Kale Kundert (kale.kundert@ucsf.edu)

// Headers {{{1
#include <protocols/canonical_sampling/DbTrajectoryRecorder.hh>
#include <protocols/canonical_sampling/DbTrajectoryRecorderCreator.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

// Protocol headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>  // required for Windows build

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

// External headers
#include <cppdb/frontend.h>

// MPI headers
#ifdef USEMPI
#include <mpi.h>
#include <protocols/jd2/util.hh>
#endif

// C++ headers
#include <sstream>

// }}}1

namespace protocols {
namespace canonical_sampling {

// Global Names {{{1
using namespace std;
using namespace core;
using core::pose::Pose;
using protocols::moves::MoverOP;
using utility::vector1;
using utility::tools::make_vector;
using utility::tools::make_vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.DbTrajectoryRecorder" );
// }}}1

string DbTrajectoryRecorderCreator::keyname() const { // {{{1
	return DbTrajectoryRecorderCreator::mover_name();
}

MoverOP DbTrajectoryRecorderCreator::create_mover() const { // {{{1
	return MoverOP( new DbTrajectoryRecorder );
}

string DbTrajectoryRecorderCreator::mover_name() { // {{{1
	return "DbTrajectoryRecorder";
}
// }}}1

DbTrajectoryRecorder::DbTrajectoryRecorder() // {{{1
: TrajectoryRecorder(),
	job_id_(-1) {}

DbTrajectoryRecorder::DbTrajectoryRecorder(Size job_id) // {{{1
: TrajectoryRecorder(),
	job_id_(job_id) {}

DbTrajectoryRecorder::DbTrajectoryRecorder( DbTrajectoryRecorder const & ) = default;

MoverOP DbTrajectoryRecorder::clone() const { // {{{1
	return MoverOP( new protocols::canonical_sampling::DbTrajectoryRecorder( *this ) );
}

MoverOP DbTrajectoryRecorder::fresh_instance() const { // {{{1
	return MoverOP( new DbTrajectoryRecorder );
}

string DbTrajectoryRecorder::get_name() const { // {{{1
	return "DbTrajectoryRecorder";
}

void DbTrajectoryRecorder::initialize_simulation( // {{{1
	core::pose::Pose & pose,
	MetropolisHastingsMover const & mover,
	core::Size cycle) {

	TrajectoryRecorder::initialize_simulation(pose, mover, cycle);
	write_schema_to_db();
	write_model(pose);
}

void DbTrajectoryRecorder::finalize_simulation( // {{{1
	core::pose::Pose & pose,
	MetropolisHastingsMover const & mover) {

	TrajectoryRecorder::finalize_simulation(pose, mover);
	write_cache_to_db();
}

bool DbTrajectoryRecorder::restart_simulation( // {{{1
	core::pose::Pose &,
	MetropolisHastingsMover&,
	core::Size&,
	core::Size&,
	core::Real&) {

	utility_exit_with_message("DbTrajectoryRecorder does not support restarting trajectories.");

	return false;
}

void DbTrajectoryRecorder::write_schema_to_db() const { // {{{1
	using utility::sql_database::sessionOP;
	using namespace basic::database::schema_generator;

#ifdef USEMPI
	// Avoid writing the schema to the database from a bunch of different
	// processes at the same time.  The problem isn't that the schema would
	// somehow get written twice.  Rather, it's that sqlite can't handle any
	// parallelism at all.

	int rank; MPI_Comm_rank(jd2::current_mpi_comm(), &rank);
	if (rank > 0) return;
#endif

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

void DbTrajectoryRecorder::write_cache_to_db() const { // {{{1
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

	for ( Frame const & frame : frame_cache_ ) {
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

void DbTrajectoryRecorder::write_model( // {{{1
	core::pose::Pose const & pose,
	MetropolisHastingsMover const *) {

	Frame frame;
	frame.iteration = step_count();
	frame.pose = pose;
	frame_cache_.push_back(frame);

	if ( frame_cache_.size() >= cache_limit() ) {
		write_cache_to_db();
	}
}
// }}}1

} // namespace canonical_sampling
} // namespace protocols
