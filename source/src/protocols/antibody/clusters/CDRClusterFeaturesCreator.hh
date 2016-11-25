// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/clusters/CDRClusterFeaturesCreator.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_clusters_CDRCLUSTERFEATURESCREATOR_HH
#define INCLUDED_protocols_antibody_clusters_CDRCLUSTERFEATURESCREATOR_HH

#include <protocols/features/FeaturesReporterCreator.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
namespace clusters {

class CDRClusterFeaturesCreator: public features::FeaturesReporterCreator {


public:
	// XRW TEMP  CDRClusterFeaturesCreator();

	// XRW TEMP  virtual ~CDRClusterFeaturesCreator();

	// XRW TEMP  virtual features::FeaturesReporterOP create_features_reporter() const;

	// XRW TEMP  virtual std::string type_name() const;
	protocols::features::FeaturesReporterOP create_features_reporter() const override;
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


} //clusters
} //antibody
} //protocols


#endif //#ifndef INCLUDED_protocols/antibody_clusters_CDRCLUSTERFEATURESCREATOR_HH

