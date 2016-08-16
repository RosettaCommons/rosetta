// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProteinBackboneAtomAtomPairFeatures.hh
/// @brief  report atom-atom pair distances between atoms in protein backbones to features Statistics Scientific Benchmark
/// @author Matthew O'Meara (mattjoomeara@gmail.com)

#ifndef INCLUDED_protocols_features_ProteinBackboneAtomAtomPairFeatures_hh
#define INCLUDED_protocols_features_ProteinBackboneAtomAtomPairFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ProteinBackboneAtomAtomPairFeatures.fwd.hh>

//External

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class ProteinBackboneAtomAtomPairFeatures : public protocols::features::FeaturesReporter {
public:
	ProteinBackboneAtomAtomPairFeatures();

	ProteinBackboneAtomAtomPairFeatures( ProteinBackboneAtomAtomPairFeatures const & src );

	virtual ~ProteinBackboneAtomAtomPairFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the protein_backbone_atom_atom_pairs table schema
	void
	write_protein_backbone_atom_atom_pairs_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

};

} // namespace
} // namespace

#endif // include guard
