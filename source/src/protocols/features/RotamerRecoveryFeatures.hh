// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/RotamerRecoveryFeatures.hh
/// @brief  report Rotamer Recovery Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_RotamerRecoveryFeatures_hh
#define INCLUDED_protocols_features_RotamerRecoveryFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/RotamerRecoveryFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/rotamer_recovery/RRReporterSQLite.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRProtocol.hh>


namespace protocols{
namespace features{

class RotamerRecoveryFeatures : public protocols::features::FeaturesReporter {
public:
	RotamerRecoveryFeatures();

	RotamerRecoveryFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	RotamerRecoveryFeatures(
		RotamerRecoveryFeatures const & src);

	virtual ~RotamerRecoveryFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	void
	initialize_task_factory(
		core::pack::task::TaskFactoryOP task_factory);

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & /*pose*/);

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

private:
	core::scoring::ScoreFunctionOP scfxn_;
	protocols::rotamer_recovery::RRReporterSQLiteOP reporter_;
	protocols::rotamer_recovery::RRProtocolOP protocol_;
	protocols::rotamer_recovery::RRComparerOP comparer_;
	core::pack::task::TaskFactoryOP task_factory_;

};

} // namespace
} // namespace

#endif // include guard
