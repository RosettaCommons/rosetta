// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/symmetry/SetupForSequenceSymmetryMover.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/symmetry/SetupForSequenceSymmetryMover.hh>
#include <protocols/symmetry/SetupForSequenceSymmetryMoverCreator.hh>

#include <protocols/moves/mover_schemas.hh>
#include <protocols/residue_selectors/StoreResidueSubsetMover.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>

// Utility Headers
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace symmetry {

static basic::Tracer TR( "protocols.simple_moves.symmetry.SetupForSequenceSymmetryMover" );

SetupForSequenceSymmetryMover::SetupForSequenceSymmetryMover()
{}

SetupForSequenceSymmetryMover::~SetupForSequenceSymmetryMover() = default;


void
SetupForSequenceSymmetryMover::apply( core::pose::Pose & pose ) {
	TR << "Adding " << independent_region_selectors_.size() << " selectors" << std::endl;
	for ( core::Size ii = 1; ii <= independent_region_selectors_.size(); ++ii ) {
		std::string const magic_selector_name =
			"SequenceSymmetricAnnealer_" + std::to_string( ii );
		residue_selectors::StoreResidueSubsetMover add_subset( independent_region_selectors_[ ii ], magic_selector_name, true );
		add_subset.apply( pose );
#ifndef NDEBUG
		TR << "Subset " << magic_selector_name << " is being added" << std::endl;
#endif
	}
}

void SetupForSequenceSymmetryMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) {

	using namespace basic::options;

	std::string const selectors = tag->getOption< std::string >( "independent_regions", "" );
	if ( selectors.size() == 0 ) return;

	utility::vector1< std::string > const elements =
		utility::string_split( selectors, ',' );

	independent_region_selectors_.clear();
	for ( std::string const & selector_name : elements ) {
		independent_region_selectors_.emplace_back(
			core::select::residue_selector::get_residue_selector( selector_name, datamap )
		);
	}
}

void SetupForSequenceSymmetryMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "independent_regions", xs_string , "Comma-separated list of residue selectors. Each one defines a region of sequence symmetry. For example, say you had a dimer and a trimer in the same system. You would pass one residue selector that selects all of the residues in the dimer and a second residue selector that selects all of the residues in the trimer. independent_regions=\"dimer_selector,trimer_selector\". If the user provides residue selectors, this will not enforce sequence symmetry on regions not covered by any selection. Otherwise the entire protein will be treated as one single region of sequence symmetry." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TODO", attlist );
}

std::string SetupForSequenceSymmetryMoverCreator::keyname() const {
	return SetupForSequenceSymmetryMover::mover_name();
}

protocols::moves::MoverOP
SetupForSequenceSymmetryMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SetupForSequenceSymmetryMover >();
}

void SetupForSequenceSymmetryMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetupForSequenceSymmetryMover::provide_xml_schema( xsd );
}



} //symmetry
} // protocols
