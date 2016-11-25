// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.hh
/// @brief --brief--
/// @author --name-- (--email--)


#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

#include <--path--/--class--.fwd.hh>
#include <protocols/features/FeaturesReporter.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

--namespace--

///@brief --brief--
class --class-- : public protocols::features::FeaturesReporter {

public:

	--class--();
	--class--(--class-- const & src);

	~--class--() override;

	--class--OP
	clone() const override;

	/// @brief return string with class name
	std::string
	type_name() const;

public: // Feature-Reporter Specific Methods

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session
	) const override;

	/// @brief return the set of features reporters that are required to
	/// also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session ) override;


private:

};

//The Creator class will need these functions:
//protocols::features::FeaturesReporterOP create_features_reporter() const override;
//std::string type_name() const override;
//void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

--end_namespace--


#endif //INCLUDED_--path_underscore--_--class--_hh
