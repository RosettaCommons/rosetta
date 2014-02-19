// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_canonical_sampling_DbTrajectoryRecorder_hh
#define INCLUDED_protocols_canonical_sampling_DbTrajectoryRecorder_hh

// Project forward headers
#include <protocols/canonical_sampling/DbTrajectoryRecorder.fwd.hh>
#include <protocols/canonical_sampling/TrajectoryRecorder.hh>

// Project headers
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <string>

namespace protocols {
namespace canonical_sampling {

/// @brief Record a trajectory to an SQL database.
///
/// @details This class builds upon Rosetta's database framework, which means 
/// that there is support for SQLite3, MySQL and PostgreSQL.  Database options 
/// must be specified on the command line (i.e. there's no API for this):
///
/// @code{.sh}
/// # SQLite3
/// -out:use_database
/// -inout:dbms:database_name traj.sqlite
///
/// # MySQL
/// -out:use_database
/// -inout:dbms:mode mysql
/// -inout:dbms:database_name ...
/// -inout:dbms:user ...
/// -inout:dbms:port ...
/// -inout:dbms:password ...
/// @endcode

class DbTrajectoryRecorder : public TrajectoryRecorder {

public:
	/// @brief Default constructor.
	DbTrajectoryRecorder();

	/// @brief Constructor with a @a job_id parameter to use as a foreign key.
	DbTrajectoryRecorder( core::Size job_id );

	/// @brief Constructor with @a job_id and @a temp_level parameters.
	DbTrajectoryRecorder( core::Size job_id, core::Size temp_level );

	/// @brief Copy constructor.
	DbTrajectoryRecorder( DbTrajectoryRecorder const & );

public:
	protocols::moves::MoverOP	clone() const;

	protocols::moves::MoverOP	fresh_instance() const;

	std::string	get_name() const;

	/// @brief Return the job id that will be used as a foreign key in the 
	/// trajectory table that gets generated.
	core::Size job_id() const { return job_id_; }

	/// @brief Set the job id that will be used as a foreign key in the 
	/// trajectory table that gets generated.
	void set_job_id( core::Size id ) { job_id_ = id; }

	/// @brief Return the temperature level that this observer is focusing on.
	core::Size temp_level() const { return temp_level_; }

	/// @brief Set the temperature level that this observer should focus on.
	/// @details Temperature levels start counting from one.  If zero is given, 
	/// it will be taken to mean that every level should be recorded.
	void set_temp_level( core::Size level ) { temp_level_ = level; }

	void initialize_simulation(
			core::pose::Pose & pose,
			MetropolisHastingsMover const & mover,
			core::Size cycle
	);

	void finalize_simulation(
			core::pose::Pose & pose,
			MetropolisHastingsMover const & mover
	);

	/// @brief Not implemented, except to complain if accidentally used.
	bool restart_simulation(
			core::pose::Pose & pose,
			MetropolisHastingsMover& mover,
			core::Size& cycle,
			core::Size& temp_level,
			core::Real& temperature
	);

private:

	/// @brief Generate the table schemas and write them to the database.
	void write_schema_to_db() const;

	/// @brief Write any cached poses into the database, then clear the cache.
	void write_cache_to_db() const;

	/// @brief Append the given model to the trajectory being written.
	/// @details The pose may be cached and not immediately written to the 
	/// database.  
	void 	write_model(
		core::pose::Pose const & pose,
		MetropolisHastingsMoverCAP mover=0
	);

	/// @brief Write the starting structure to the database.
	/// @details This case is handled specially to guarantee that the frame 
	/// labeled as "iteration 0" always represents the starting structure.  In 
	/// the future, this method will make sure only the root node writes this 
	/// frame.
	void write_first_model(
		core::pose::Pose const & pose,
		MetropolisHastingsMoverCAP mover
	);

private:

	core::Size job_id_;
	core::Size temp_level_;

	/// @brief Helper struct used store cached poses.
	struct Frame {
		core::Size temp_level;
		core::Size iteration;
		core::pose::Pose pose;
	};
	mutable utility::vector1<Frame> frame_cache_;

}; // DbTrajectoryRecorder


} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_DbTrajectoryRecorder_HH
