// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/DeclareStructureDataCovalentBondMover.cc
/// @brief declares covalent bond
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/movers/DeclareStructureDataCovalentBondMover.hh>
#include <protocols/denovo_design/movers/DeclareStructureDataCovalentBondMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.denovo_design.movers.DeclareStructureDataCovalentBondMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
DeclareStructureDataCovalentBondMover::DeclareStructureDataCovalentBondMover():
	protocols::moves::Mover( DeclareStructureDataCovalentBondMover::mover_name() ),
	atom1_( "" ),
	atom2_( "" )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
DeclareStructureDataCovalentBondMover::~DeclareStructureDataCovalentBondMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
DeclareStructureDataCovalentBondMover::apply( core::pose::Pose & pose )
{
	if ( atom1_.find(',') == std::string::npos ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "Atom1 description " + atom1_ + " is malformed. It should be of the form segment#residue,atom." );
	}
	if ( atom2_.find(',') == std::string::npos ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "Atom2 description " + atom2_ + " is malformed. It should be of the form segment#residue,atom." );
	}

	utility::vector1< std::string > res_plus_atom1 = utility::string_split( atom1_, ',' );
	if ( res_plus_atom1.size() != 2 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "Atom1 description " + atom1_ + " is malformed. It should be of the form segment#residue,atom." );
	}

	utility::vector1< std::string > res1_strs = utility::string_split( res_plus_atom1[1], '.' );
	if ( res1_strs.size() == 1 ) {
		res1_strs.push_back( "0" );
	}
	if ( res1_strs.size() != 2 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "Atom1 description " + atom1_ + " is malformed. It should be of the form segment#residue,atom." );
	}

	std::string const seg1 = res1_strs[1];
	core::Size const res1 = boost::lexical_cast< core::Size >(res1_strs[2]);
	std::string const a1 = res_plus_atom1[2];

	utility::vector1< std::string > res_plus_atom2 = utility::string_split( atom2_, ',' );
	if ( res_plus_atom2.size() != 2 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "Atom2 description " + atom2_ + " is malformed. It should be of the form segment#residue,atom." );
	}
	utility::vector1< std::string > res2_strs = utility::string_split( res_plus_atom2[1], '.' );
	if ( res2_strs.size() == 1 ) {
		res2_strs.push_back( "0" );
	}
	if ( res2_strs.size() != 2 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "Atom2 description " + atom2_ + " is malformed. It should be of the form segment#residue,atom." );
	}

	std::string const seg2 = res2_strs[1];
	core::Size const res2 = boost::lexical_cast< core::Size >(res2_strs[2]);
	std::string const a2 = res_plus_atom2[2];

	components::StructureData perm = components::StructureDataFactory::get_instance()->get_from_pose( pose );
	core::Size const resid1 = perm.pose_residue( seg1, res1 );
	std::string const newname = pose.residue(resid1).type().name() + ":connect" + a1;
	protocols::simple_moves::MutateResidue mut1( resid1, newname );
	mut1.set_preserve_atom_coords( true );
	mut1.apply( pose );

	core::Size const resid2 = perm.pose_residue( seg2, res2 );
	std::string const newname2 = pose.residue(resid2).type().name() + ":connect" + a2;
	protocols::simple_moves::MutateResidue mut2( resid2, newname2 );
	mut2.set_preserve_atom_coords( true );
	mut2.apply( pose );

	perm.declare_covalent_bond( seg1, res1, a1, seg2, res2, a2 );
	components::StructureDataFactory::get_instance()->save_into_pose( pose, perm );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
DeclareStructureDataCovalentBondMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &
operator<<( std::ostream & os, DeclareStructureDataCovalentBondMover const & mover )
{
	mover.show(os);
	return os;
}


////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
DeclareStructureDataCovalentBondMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	atom1_ = tag->getOption< std::string >( "atom1", atom1_ );
	atom2_ = tag->getOption< std::string >( "atom2", atom2_ );
	if ( ! atom1_.size() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Atom1 option must be specified to DeclareCovalentBond." );
	}
	if ( ! atom2_.size() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Atom2 option must be specified to DeclareCovalentBond." );
	}
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
DeclareStructureDataCovalentBondMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new DeclareStructureDataCovalentBondMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
DeclareStructureDataCovalentBondMover::clone() const
{
	return protocols::moves::MoverOP( new DeclareStructureDataCovalentBondMover( *this ) );
}

// XRW TEMP std::string
// XRW TEMP DeclareStructureDataCovalentBondMover::get_name() const
// XRW TEMP {
// XRW TEMP  return DeclareStructureDataCovalentBondMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DeclareStructureDataCovalentBondMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "DeclareStructureDataCovalentBondMover";
// XRW TEMP }

////////////////////////////////////////////////////////////////////////////////
/// Creator ///
///////////////

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DeclareStructureDataCovalentBondMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new DeclareStructureDataCovalentBondMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DeclareStructureDataCovalentBondMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DeclareStructureDataCovalentBondMover::mover_name();
// XRW TEMP }

std::string DeclareStructureDataCovalentBondMover::get_name() const {
	return mover_name();
}

std::string DeclareStructureDataCovalentBondMover::mover_name() {
	return "DeclareStructureDataCovalentBondMover";
}

void DeclareStructureDataCovalentBondMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "atom1", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "atom2", xs_string, "XRW TO DO" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string DeclareStructureDataCovalentBondMoverCreator::keyname() const {
	return DeclareStructureDataCovalentBondMover::mover_name();
}

protocols::moves::MoverOP
DeclareStructureDataCovalentBondMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DeclareStructureDataCovalentBondMover );
}

void DeclareStructureDataCovalentBondMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DeclareStructureDataCovalentBondMover::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

} //protocols
} //denovo_design
} //movers

