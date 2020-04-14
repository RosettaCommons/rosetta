// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/IdealizeAnomericHydrogens.cc
/// @brief Mover to idealize carbohydrate anomeric hydrogens
/// @author Brandon Frenz


// Unit headers
#include <protocols/carbohydrates/IdealizeAnomericHydrogens.hh>
#include <protocols/carbohydrates/IdealizeAnomericHydrogensCreator.hh>

// Package headers
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>

// Project headers
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/kinematics/util.hh>
#include <core/conformation/Conformation.hh>

#include <protocols/moves/MoverStatus.hh>
#include <protocols/rosetta_scripts/util.hh>


// Basic headers
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>

// Utility header
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.carbohydrates.IdealizeAnomericHydrogens" );


namespace protocols {
namespace carbohydrates {

IdealizeAnomericHydrogens::~IdealizeAnomericHydrogens()= default;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

protocols::moves::MoverOP
IdealizeAnomericHydrogens::clone() const{
	return utility::pointer::make_shared< IdealizeAnomericHydrogens >( *this );
}

moves::MoverOP
IdealizeAnomericHydrogens::fresh_instance() const
{
	return utility::pointer::make_shared< IdealizeAnomericHydrogens >();
}


void
IdealizeAnomericHydrogens::apply( core::pose::Pose & pose )
{
	using namespace core::pose::carbohydrates;
	using namespace core::chemical::carbohydrates;
	using namespace core::chemical;

	for ( core::Size i=1; i<=pose.size(); i++ ) {
		if ( !pose.residue_type(i).is_carbohydrate() ) continue;

		CarbohydrateInfoCOP carbo_info( pose.residue(i).carbohydrate_info() );
		core::Size anomeric_sidechain_index = carbo_info->anomeric_sidechain_index();
		if ( anomeric_sidechain_index == 0 ) {
			TR.Warning << "Cannot idealize residue " << i << " [" << pose.residue_type(i).name() << "] as it is missing an anomeric designation." << std::endl;
			continue;
		}
		core::chemical::AtomICoor icoor = pose.residue_type(i).icoor(anomeric_sidechain_index);
		core::id::AtomID id(anomeric_sidechain_index,i);
		numeric::xyzVector<core::Real> icoor_xyz = icoor.build( pose.residue(i), pose.conformation() );
		pose.conformation().set_xyz( id, icoor_xyz );
	}
}

/////////////// Creator ///////////////




std::string IdealizeAnomericHydrogens::get_name() const {
	return mover_name();
}

std::string IdealizeAnomericHydrogens::mover_name() {
	return "IdealizeAnomericHydrogens";
}

void IdealizeAnomericHydrogens::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_residue_selector( attlist );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Mover to idealize carbohydrate anomeric hydrogens", attlist );
}

void IdealizeAnomericHydrogens::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &
) {} // No configuration needed.

std::string IdealizeAnomericHydrogensCreator::keyname() const {
	return IdealizeAnomericHydrogens::mover_name();
}

protocols::moves::MoverOP
IdealizeAnomericHydrogensCreator::create_mover() const {
	return utility::pointer::make_shared< IdealizeAnomericHydrogens >();
}

void IdealizeAnomericHydrogensCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	IdealizeAnomericHydrogens::provide_xml_schema( xsd );
}

//Citation Manager

bool IdealizeAnomericHydrogens::mover_provides_citation_info() const {
	return true;
}

/// @brief Provide the citation.
/// @returns A vector of citation collections.  This allows the mover to provide citations for
/// itself and for any modules that it invokes.
/// @details Also provides citations for movers called by the BundleGridSampler.
/// @author Brandon Frenz (brandon.frenz@brandon.frenz@gmail.com)
utility::vector1< basic::citation_manager::CitationCollectionCOP >
IdealizeAnomericHydrogens::provide_citation_info() const {
	basic::citation_manager::CitationCollectionOP cc(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
		get_name(), basic::citation_manager::CitedModuleType::Mover
		)
	);

	cc->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi( "10.1016/j.str.2018.09.006" ) );
	utility::vector1< basic::citation_manager::CitationCollectionCOP > returnvec{ cc };
	return returnvec;
}


} //carbohydrates
} //protocols
