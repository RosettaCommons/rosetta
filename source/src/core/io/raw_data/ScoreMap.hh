// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ScoreMap.hh
/// @brief  A place to put some common functions for score-file output
/// @author Monica Berrondo

#ifndef INCLUDED_core_io_raw_data_ScoreMap_hh
#define INCLUDED_core_io_raw_data_ScoreMap_hh

#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <map>


namespace core {
namespace io {
namespace raw_data {

class ScoreMap : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ScoreMap() override;
	/// @brief full atom energies for output

	///@brief Scores and grabs (nonzero, weighted) energies map.
	static
	void
	score_and_add_energies_to_map(
		std::map< std::string, core::Real> & score_map,
		core::scoring::ScoreFunctionCOP score_fxn,
		core::pose::Pose & pose
	);

	/// @brief generates a (nonzero) (weighted) scoremap assuming the pose is already scored (note const w.r.t. pose)
	static
	std::map< std::string, core::Real>
	get_energies_map_from_scored_pose( core::pose::Pose const & pose );

	///@brief Return by value version of extra string data.
	///
	///@details
	/// Get data set as PoseExtraScores (ARBITRARY_STRING_DATA) and SimpleMetrics (SIMPLE_METRIC_DATA).
	static
	std::map< std::string, std::string >
	get_arbitrary_string_data_from_pose( pose::Pose const & pose );

	///@brief Return by value version of extra score data.
	///
	///@details
	/// Get data set as PoseExtraScores (ARBITRARY_FLOAT_DATA) and SimpleMetrics (SIMPLE_METRIC_DATA).
	static
	std::map< std::string, core::Real >
	get_arbitrary_score_data_from_pose( pose::Pose const & pose );


	/////////////////// Add Data ///////////////////////

	///@brief Add energies to scores map.
	static
	void
	add_energies_data_from_scored_pose(
		core::pose::Pose const & pose,
		std::map< std::string, core::Real> & scores
	);

	///@brief Add data set as PoseExtraScores (ARBITRARY_STRING_DATA) and SimpleMetrics (SIMPLE_METRIC_DATA).
	static
	void
	add_arbitrary_string_data_from_pose(
		core::pose::Pose const & pose,
		std::map < std::string, std::string > & string_map
	);

	///@brief Add data set as PoseExtraScores (ARBITRARY_FLOAT_DATA) and SimpleMetrics (SIMPLE_METRIC_DATA).
	static
	void
	add_arbitrary_score_data_from_pose(
		core::pose::Pose const & pose,
		std::map < std::string, core::Real > & score_map
	);


	/////////// Etc //////////////

	/// @brief print out the values in the scoremap
	static
	void
	print(
		std::map < std::string, core::Real > const & score_map,
		std::ostream & out
	);


};

} // raw_data
} // io
} // core

#endif
