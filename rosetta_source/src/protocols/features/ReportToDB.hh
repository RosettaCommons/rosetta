// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ReportToDB.hh
/// @brief  report feature data to sqlite database
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ReportToDB_hh
#define INCLUDED_protocols_features_ReportToDB_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/features/ReportToDB.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/features/FeaturesReporter.fwd.hh>
#include <protocols/features/FeaturesReporterFactory.fwd.hh>
#include <protocols/features/ProtocolFeatures.fwd.hh>
#include <protocols/features/BatchFeatures.fwd.hh>
#include <protocols/features/StructureFeatures.fwd.hh>


// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// Boost Headers
#include <boost/uuid/uuid.hpp>

// C++ Headers
#include <string>


namespace protocols{
namespace features{

class ReportToDB : public protocols::moves::Mover {

public:
	ReportToDB();

	ReportToDB(std::string const & name);

	ReportToDB(
		std::string const & name,
		std::string const & database_fname,
		std::string const & sample_source,
		core::scoring::ScoreFunctionOP scfxn,
		bool use_transactions=true,
		core::Size cache_size=2000);

	ReportToDB(ReportToDB const & src);

	virtual ~ReportToDB();

	virtual
	void
	register_options() const;

	virtual moves::MoverOP fresh_instance() const;

	virtual moves::MoverOP clone() const;

	virtual std::string get_name() const { return "ReportToDB"; }


	void
	parse_db_tag_item(
		utility::tag::TagPtr const tag);

	void
	parse_sample_source_tag_item(
		utility::tag::TagPtr const tag);

	void
	parse_protocol_id_tag_item(
		utility::tag::TagPtr const tag);

	/* Undefined, commenting out to fix PyRosetta build  void
	parse_struct_id_type_tag_item(
		utility::tag::TagPtr const tag);
	*/

	/* Undefined, commenting out to fix PyRosetta build  void
	parse_first_struct_id_tag_item(
		utility::tag::TagPtr const tag);
	*/

	void
	parse_db_mode_tag_item(
		utility::tag::TagPtr const tag);

	void
	parse_separate_db_per_mpi_process_tag_item(
		utility::tag::TagPtr const tag);

	void
	parse_use_transactions_tag_item(
		utility::tag::TagPtr const tag);

	void
	parse_cache_size_tag_item(
		utility::tag::TagPtr const tag);

	void
	parse_name_tag_item(TagPtr const tag);

	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose );

	void
	check_features_reporter_dependencies(
		FeaturesReporterOP test_feature_reporter
	) const;

	void
	initialize_reporters();

	utility::sql_database::sessionOP
	initialize_database();

	///@brief initialize the pose and return the relevant residues
	utility::vector1< bool >
	initialize_pose(
		Pose & pose
	) const;

	///@brief Add the defined features reporters to the
	///'features_reporters' table in the database
	void
	write_features_reporters_table(
		utility::sql_database::sessionOP db_session) const;

	///@brief Link the defined features reporters to the batch of
	///structures extracted with this invocation of the ReportToDB mover
	void
	write_batch_reports_table(
		utility::sql_database::sessionOP db_session) const;

	///@brief write tables linking the batches table with the features
	///datababase
	void
	write_linking_tables(
		utility::sql_database::sessionOP db_session) const;

	void
	apply(
		Pose& pose);

	boost::uuids::uuid
	report_structure_features(
		utility::vector1<bool> const & relevant_residues,
		utility::sql_database::sessionOP db_session) const;

	void
	report_features(
		core::pose::Pose const & pose,
		boost::uuids::uuid struct_id,
		utility::vector1<bool> const & relevant_residues,
		utility::sql_database::sessionOP db_session) const;

private:
	std::string database_fname_;
	std::string database_mode_;
	bool separate_db_per_mpi_process_;
	std::string sample_source_;
	std::string name_;

	core::scoring::ScoreFunctionOP scfxn_;
	bool use_transactions_;

	core::Size cache_size_;

	core::Size protocol_id_;
	core::Size batch_id_;

	// initialized in parse_my_tag
	core::pack::task::TaskFactoryOP task_factory_;

	protocols::features::FeaturesReporterFactory * features_reporter_factory_;
	protocols::features::ProtocolFeaturesOP protocol_features_;
	protocols::features::BatchFeaturesOP batch_features_;
	protocols::features::StructureFeaturesOP structure_features_;
	utility::vector1< protocols::features::FeaturesReporterOP > features_reporters_;
	bool initialized;
};

} // namespace
} // namespace

#endif //include guard
