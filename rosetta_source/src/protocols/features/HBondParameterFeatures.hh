// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/HBondParameterFeatures.hh
/// @brief  report HBond Parameter features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_HBondParameterFeatures_hh
#define INCLUDED_protocols_features_HBondParameterFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/HBondParameterFeatures.fwd.hh>

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <utility/sql_database/DatabaseSessionManager.hh>

// C++ Headers
#include <string>

#include <core/scoring/ScoreFunction.fwd.hh>


namespace protocols{
namespace features{

class HBondParameterFeatures : public protocols::features::FeaturesReporter {
public:
	HBondParameterFeatures();

	HBondParameterFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	HBondParameterFeatures( HBondParameterFeatures const & src );

	virtual ~HBondParameterFeatures();

	///@breif return string with class name
	std::string
	type_name() const;

	///@breif return sql statements that setup the right tables
	std::string
	schema() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	///@breif collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

private:
	core::scoring::ScoreFunctionOP scfxn_;
};

} // namespace
} // namespace

#endif // include guard
