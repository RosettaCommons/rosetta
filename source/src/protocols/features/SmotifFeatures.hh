// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SmotifFeatures.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_SmotifFeature_hh
#define INCLUDED_protocols_features_SmotifFeature_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/SmotifFeatures.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.fwd.hh>

#include <numeric/xyzVector.fwd.hh>

// C++ Headers
#include <string>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace features {

struct SecondaryStructureSegment{
	core::Size segment_id;
	core::Size residue_begin;
	core::Size residue_end;
	std::string dssp;
};

class SmotifFeatures : public protocols::features::FeaturesReporter {
public:
	SmotifFeatures();

	SmotifFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	SmotifFeatures( SmotifFeatures const & src );

	~SmotifFeatures() override;

	/// @brief return string with class name
	// XRW TEMP  std::string
	// XRW TEMP  type_name() const override;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/) override;

	void
	calculate_angles(
		utility::vector1< numeric::xyzVector<core::Real> > const & ss1_coords,
		utility::vector1< numeric::xyzVector<core::Real> > const & ss2_coords,
		core::Real & distance, /*output*/
		core::Real & hoist_angle_degrees, /*output*/
		core::Real & packing_angle_degrees, /*output*/
		core::Real & meridian_angle_degrees /*output*/
	);

	utility::vector1<SecondaryStructureSegment>
	get_ss_segments(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	/* Undefined, commenting out to fix PyRosetta build  numeric::xyzVector<core::Real>
	pca(utility::vector1< numeric::xyzVector< core::Real > > coords); */

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

} // namespace
} // namespace

#endif // include guard
