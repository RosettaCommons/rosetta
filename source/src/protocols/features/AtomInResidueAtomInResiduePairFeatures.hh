// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/AtomInResidueAtomInResiduePairFeatures.hh
/// @brief  report atom-atom pair geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_AtomInResidueAtomInResiduePairFeatures_hh
#define INCLUDED_protocols_features_AtomInResidueAtomInResiduePairFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/AtomInResidueAtomInResiduePairFeatures.fwd.hh>

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

class AtomInResidueAtomInResiduePairFeatures : public protocols::features::FeaturesReporter {
public:
	AtomInResidueAtomInResiduePairFeatures();

	AtomInResidueAtomInResiduePairFeatures( AtomInResidueAtomInResiduePairFeatures const & src );

	virtual ~AtomInResidueAtomInResiduePairFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the atom_in_residue_pairs table schema
	void
	write_atom_in_residue_pairs_table_schema(
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

	/// @brief Similar to radial density distributions defined in
	///Lu H, Skolnick J. A distance-dependent atomic knowledge-based
	///potential for improved protein structure
	///selection. Proteins. 2001;44(3):223-32. Available at:
	///http://www.ncbi.nlm.nih.gov/pubmed/11455595.
	void
	report_atom_pairs(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

};

} // namespace
} // namespace

#endif // include guard
