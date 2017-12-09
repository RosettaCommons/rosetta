// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/IdealizeAnomericHydrogens.cc
/// @brief This code changes all of the dihedrals of a particular glycosidic linkage based on database info,
///   esentially sampling carbohydrate dihedral conformers of two residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>

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

// Utility header
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.carbohydrates.IdealizeAnomericHydrogens" );


namespace protocols {
namespace carbohydrates {

IdealizeAnomericHydrogens::~IdealizeAnomericHydrogens()= default;

//void
//IdealizeAnomericHydrogens::parse_my_tag(
// utility::tag::TagCOP tag,
// basic::datacache::DataMap& datamap,
// protocols::filters::Filters_map const & ,
// protocols::moves::Movers_map const & ,
// core::pose::Pose const & )
//{
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

protocols::moves::MoverOP
IdealizeAnomericHydrogens::clone() const{
	return protocols::moves::MoverOP( new IdealizeAnomericHydrogens( *this ) );
}

moves::MoverOP
IdealizeAnomericHydrogens::fresh_instance() const
{
	return protocols::moves::MoverOP( new IdealizeAnomericHydrogens );
}

// XRW TEMP std::string
// XRW TEMP IdealizeAnomericHydrogens::get_name() const {
// XRW TEMP  return "IdealizeAnomericHydrogens";
// XRW TEMP }

void
IdealizeAnomericHydrogens::apply( core::pose::Pose & pose )
{
	using namespace core::pose::carbohydrates;
	using namespace core::chemical::carbohydrates;
	using namespace core::chemical;

	for ( core::Size i=1; i<=pose.size(); i++ ) {
		if ( !pose.residue_type(i).is_carbohydrate() ) continue;
		core::Size anomeric_h = pose.residue_type(i).atom_index("H1");
		core::chemical::AtomICoor icoor = pose.residue_type(i).icoor(anomeric_h);
		core::id::AtomID id(anomeric_h,i);
		numeric::xyzVector<core::Real> icoor_xyz = icoor.build(pose.residue(i), pose.conformation());
		pose.conformation().set_xyz( id, icoor_xyz);
	}
}

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP IdealizeAnomericHydrogensCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new IdealizeAnomericHydrogens );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP IdealizeAnomericHydrogensCreator::keyname() const {
// XRW TEMP  return IdealizeAnomericHydrogens::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP IdealizeAnomericHydrogens::mover_name(){
// XRW TEMP  return "IdealizeAnomericHydrogens";
// XRW TEMP }

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

std::string IdealizeAnomericHydrogensCreator::keyname() const {
	return IdealizeAnomericHydrogens::mover_name();
}

protocols::moves::MoverOP
IdealizeAnomericHydrogensCreator::create_mover() const {
	return protocols::moves::MoverOP( new IdealizeAnomericHydrogens );
}

void IdealizeAnomericHydrogensCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	IdealizeAnomericHydrogens::provide_xml_schema( xsd );
}


} //carbohydrates
} //protocols
