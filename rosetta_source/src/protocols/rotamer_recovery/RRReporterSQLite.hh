// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRReporterSQLite.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_rotamer_recovery_RRReporterSQLite_hh
#define INCLUDED_protocols_rotamer_recovery_RRReporterSQLite_hh

// Unit Headers
#include <protocols/rotamer_recovery/RRReporter.hh>
#include <protocols/rotamer_recovery/RRReporterSQLite.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

// C++ Headers
#include <ostream>

namespace protocols {
namespace rotamer_recovery {

class RRReporterSQLite : public RRReporter {

public: // constructors destructors

	struct OutputLevel{
		enum e {
			full = 1,
			features,
			none
		};
	};

RRReporterSQLite();

	RRReporterSQLite(
		std::string const & database_fname,
		OutputLevel::e const output_level = OutputLevel::full
	);

	RRReporterSQLite(
		utility::sql_database::sessionOP db_session,
		OutputLevel::e const output_level = OutputLevel::full
	);

	~RRReporterSQLite();

	RRReporterSQLite( RRReporterSQLite const & );

public: // public interface

	static
	std::string
	schema(
		OutputLevel::e output_level = OutputLevel::full);

	virtual
	void
	set_comparer_info(
		std::string const & comparer_name,
		std::string const & comparer_params);

	void
	set_output_level(
		OutputLevel::e const output_level );

	OutputLevel::e
	get_output_level() const;

	void
	set_struct_id1(
		core::Size const struct_id1);

	core::Size
	get_struct_id1() const;

	void
	set_struct_id2(
		Size const struct_id1);

	core::Size
	get_struct_id2() const;

	virtual
	void
	reset_recovery();

	virtual
	void
	report_rotamer_recovery(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real score,
		bool recovered
	);


	virtual
	void
	report_rotamer_recovery_full(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real score,
		bool recovered
	);

	virtual
	void
	report_rotamer_recovery_features(
		core::Size struct_id1,
		core::conformation::Residue const & res1,
		core::Real score
	);

	virtual
	core::Real
	recovery_rate() const;

	virtual
	void
	show(std::ostream & out ) const;

	virtual
	void
	show() const;

private: // data members

	OutputLevel::e output_level_;

	core::Size struct_id1_;
	core::Size struct_id2_;

	std::string comparer_name_;
	std::string comparer_params_;

	core::Size residues_considered_;
	core::Size rotamers_recovered_;

	utility::sql_database::sessionOP db_session_;

};

} // namespace rotamer_recovery
} // namespace protocols

#endif // include guard
