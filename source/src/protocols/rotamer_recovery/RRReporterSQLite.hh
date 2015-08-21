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
#include <protocols/features/FeaturesReporter.fwd.hh>
#include <protocols/features/ReportToDB.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <string>
// C++ Headers


namespace protocols {
namespace rotamer_recovery {

enum OutputLevel {
	OL_full = 1,
	OL_features,
	OL_none
};


class RRReporterSQLite : public RRReporter {

public: // constructors destructors

	RRReporterSQLite();

	RRReporterSQLite(
		std::string const & database_name,
		std::string const & database_pq_schema = "",
		OutputLevel output_level = protocols::rotamer_recovery::OL_full
	);

	RRReporterSQLite(
		utility::sql_database::sessionOP db_session,
		OutputLevel output_level = protocols::rotamer_recovery::OL_full
	);

	~RRReporterSQLite();

	RRReporterSQLite( RRReporterSQLite const & );

public: // public interface

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the nchi table schema
	void
	write_nchi_table_schema(
		utility::sql_database::sessionOP db_session) const;

	/// @brief generate the full rotamer_recovery table schema
	void
	write_rotamer_recovery_full_table_schema(
		utility::sql_database::sessionOP db_session) const;

	/// @brief generate the features rotamer_recovery table schema
	void
	write_rotamer_recovery_features_table_schema(
		utility::sql_database::sessionOP db_session) const;

	/// @brief generate the predicted_features table schema
	void
	write_predicted_features_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	void
	set_protocol_info(
		std::string const & protocol_name,
		std::string const & protocol_params);

	void
	set_comparer_info(
		std::string const & comparer_name,
		std::string const & comparer_params);

	void
	db_session(
		utility::sql_database::sessionOP db_session);

	utility::sql_database::sessionOP
	db_session();

	void
	set_output_level(OutputLevel output_level );

	OutputLevel get_output_level() const;

	void
	set_struct_id1(
		protocols::features::StructureID const struct_id1);

	protocols::features::StructureID
	get_struct_id1() const;

	void
	set_predicted_report_to_db(
		features::ReportToDBOP report_to_db);

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
		protocols::features::StructureID struct_id1,
		core::conformation::Residue const & res1,
		core::Real score,
		bool recovered
	);

	virtual
	void
	report_predicted_features(
		features::StructureID struct_id1,
		core::conformation::Residue const & res1,
		core::pose::Pose const & predicted_pose,
		core::conformation::Residue const & predicted_res
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

	OutputLevel output_level_;

	protocols::features::StructureID struct_id1_;

	//Additional features can be reported for the predicted conformations
	protocols::features::ReportToDBOP report_to_db_;

	std::string protocol_name_;
	std::string protocol_params_;

	std::string comparer_name_;
	std::string comparer_params_;

	core::Size residues_considered_;
	core::Size rotamers_recovered_;

	std::string database_name_;
	std::string database_pq_schema_;
	utility::sql_database::sessionOP db_session_;

};

} // namespace rotamer_recovery
} // namespace protocols

#endif // include guard
