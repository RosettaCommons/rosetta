// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/UnrecognizedAtomFeatures.hh
/// @brief  report residue to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_UnrecognizedAtomFeatures_hh
#define INCLUDED_protocols_features_UnrecognizedAtomFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/UnrecognizedAtomFeatures.fwd.hh>

//External

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>


namespace protocols{
namespace features{

class UnrecognizedAtomFeatures : public protocols::features::FeaturesReporter {
public:
	UnrecognizedAtomFeatures();

	UnrecognizedAtomFeatures(core::Real neighbor_distance_cutoff);

	UnrecognizedAtomFeatures(UnrecognizedAtomFeatures const & src);

	virtual ~UnrecognizedAtomFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the unrecognized_atoms table schema
	virtual void
	write_unrecognized_atoms_table_schema(
		utility::sql_database::sessionOP db_session
	) const;

	/// @brief generate the unrecognized_residues table schema
	virtual void
	write_unrecognized_residues_table_schema(
		utility::sql_database::sessionOP db_session
	) const;

	/// @brief generate the unrecognized_neighbors table schema
	virtual void
	write_unrecognized_neighbors_table_schema(
		utility::sql_database::sessionOP db_session
	) const;

public:
	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

private:
	void
	insert_unrecognized_residues_rows(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_unrecognized_atoms_rows(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_unrecognized_neighbors_rows(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

public:
	void
	delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_sesion);

private:

	core::Real neighbor_distance_cutoff_;

};

} // namespace
} // namespace

#endif // include guard
