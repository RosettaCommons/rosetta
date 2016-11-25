// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/utilities/PrepareForCentroid.hh>
#include <protocols/loop_modeling/utilities/PrepareForCentroidCreator.hh>

// Core headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

// XRW TEMP moves::MoverOP PrepareForCentroidCreator::create_mover() const {
// XRW TEMP  return moves::MoverOP( new PrepareForCentroid );
// XRW TEMP }

// XRW TEMP string PrepareForCentroidCreator::keyname() const {
// XRW TEMP  return "PrepareForCentroid";
// XRW TEMP }

PrepareForCentroid::PrepareForCentroid() {}

bool PrepareForCentroid::do_apply(Pose & pose) {
	using core::util::switch_to_residue_type_set;
	using core::chemical::CENTROID_t;

	if ( ! pose.is_centroid() ) {
		switch_to_residue_type_set(pose, CENTROID_t);
	}

	return true;
}

std::string PrepareForCentroid::get_name() const {
	return mover_name();
}

std::string PrepareForCentroid::mover_name() {
	return "PrepareForCentroid";
}

void PrepareForCentroid::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Switch to centroid mode", attlist );
}

std::string PrepareForCentroidCreator::keyname() const {
	return PrepareForCentroid::mover_name();
}

protocols::moves::MoverOP
PrepareForCentroidCreator::create_mover() const {
	return protocols::moves::MoverOP( new PrepareForCentroid );
}

void PrepareForCentroidCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PrepareForCentroid::provide_xml_schema( xsd );
}


}
}
}

