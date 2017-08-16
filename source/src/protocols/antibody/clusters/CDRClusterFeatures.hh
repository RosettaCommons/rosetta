// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/cluster/CDRClusterFeatures.hh
/// @brief FeaturesReporter for North/Dunbrack Cannonical CDRClusters
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_clusters_CDRCLUSTERFEATURES_HH
#define INCLUDED_protocols_antibody_clusters_CDRCLUSTERFEATURES_HH

#include <protocols/antibody/clusters/CDRClusterFeatures.fwd.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace antibody {
namespace clusters {

class CDRClusterFeatures : public features::FeaturesReporter {

public:
	CDRClusterFeatures();

	virtual ~CDRClusterFeatures();


	//Required
	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	utility::vector1< std::string >
	features_reporter_dependencies() const override;

	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & residues,
		features::StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

	/// @brief Specify specific CDRs to load and analyze, with cdr_definition and scheme
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/) override;

public:

	/// @brief Limit CDRs being analyzed.
	void
	set_cdrs_to_use(utility::vector1< CDRNameEnum > cdrs );

	/// @brief set the numbering scheme used by the pose.
	void
	set_numbering_scheme(AntibodyNumberingSchemeEnum const & numbering_scheme);

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	utility::vector1< CDRNameEnum > cdrs_;
	AntibodyNumberingSchemeEnum numbering_scheme_;

};




}
}
}


#endif //#ifndef INCLUDED_protocols/antibody_clusters_CDRCLUSTERFEATURES_HH

