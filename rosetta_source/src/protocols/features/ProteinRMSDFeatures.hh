// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinRMSDFeatures.hh
/// @brief  report ProteinRMSD similarity of structure against supplied reference structure
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ProteinRMSDFeatures_hh
#define INCLUDED_protocols_features_ProteinRMSDFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ProteinRMSDFeatures.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

namespace protocols{
namespace features{

class ProteinRMSDFeatures : public FeaturesReporter {
public:

	/// @details If you use this constructore, before applying, make
	/// sure the reference pose is specified, either through the
	/// set_reference_pose method or parse_my_tag.
	ProteinRMSDFeatures() {}

	ProteinRMSDFeatures(
		core::pose::PoseCOP reference_pose
	);

	ProteinRMSDFeatures(
		ProteinRMSDFeatures const & ) :
		FeaturesReporter()
	{}

	virtual ~ProteinRMSDFeatures(){}

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief return sql statements that setup the right tables
	std::string
	schema() const;

	core::pose::PoseCOP reference_pose() const;
	void reference_pose(core::pose::PoseCOP);

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose);

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session
	);

private:
	core::pose::PoseCOP reference_pose_;

};

} // features namespace
} // protocols namespace

#endif // include guard
