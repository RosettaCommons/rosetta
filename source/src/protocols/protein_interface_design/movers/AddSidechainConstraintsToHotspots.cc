// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspots.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspots.hh>
#include <protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspotsCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility/tag/Tag.hh>
#include <protocols/protein_interface_design/movers/PlaceUtils.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/kinematics/FoldTree.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using core::pose::Pose;
using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.AddSidechainConstraintsToHotspots" );
static THREAD_LOCAL basic::Tracer TR_cst( "protocols.protein_interface_design.movers.AddSidechainConstraintsToHotspots_csts" );

// XRW TEMP std::string
// XRW TEMP AddSidechainConstraintsToHotspotsCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AddSidechainConstraintsToHotspots::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddSidechainConstraintsToHotspotsCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AddSidechainConstraintsToHotspots );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddSidechainConstraintsToHotspots::mover_name()
// XRW TEMP {
// XRW TEMP  return "AddSidechainConstraintsToHotspots";
// XRW TEMP }


AddSidechainConstraintsToHotspots::AddSidechainConstraintsToHotspots() :
	protocols::moves::Mover( AddSidechainConstraintsToHotspots::mover_name() ),
	chain_( 2 ),
	coord_sdev_( 1.0 )
{ }

AddSidechainConstraintsToHotspots::~AddSidechainConstraintsToHotspots() {}

void
AddSidechainConstraintsToHotspots::apply( Pose & pose )
{
	core::Size const begin( pose.conformation().chain_begin( chain() ) );
	core::Size const end( pose.conformation().chain_end( chain() ) );

	Size const num_cutpoints( pose.fold_tree().num_cutpoint() );
	if ( num_cutpoints <= 2 && residues().size() == 0 ) {
		TR<<"Not enough cutpoints in pose and no residues defined by user. Doing nothing"<<std::endl;
		return;
	}
	for ( Size i=2; i<=pose.fold_tree().num_cutpoint(); ++i ) {
		core::Size const cutpoint = pose.fold_tree().cutpoint( i );
		core::Size const cutpoint_i_1 = pose.fold_tree().cutpoint( i - 1 );
		if ( cutpoint - 1 != cutpoint_i_1 ) continue;//only mark residues that are cut on both ends
		if ( cutpoint <= end && cutpoint >= begin ) {
			add_residue( i );
		}
	}
	for ( core::Size const residue : residues() ) {
		using namespace core::scoring::constraints;

		core::scoring::func::HarmonicFuncOP dummy_cst;
		ConstraintCOPs constraint;
		constraint = add_coordinate_constraints( pose, pose.conformation().residue( residue ), chain(), residue, coord_sdev(), dummy_cst );
		for ( ConstraintCOP cst : constraint ) {
			cst->show_def( TR_cst, pose );
		}
	}
	TR.flush();
}

// XRW TEMP std::string
// XRW TEMP AddSidechainConstraintsToHotspots::get_name() const {
// XRW TEMP  return AddSidechainConstraintsToHotspots::mover_name();
// XRW TEMP }

void
AddSidechainConstraintsToHotspots::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	Movers_map const &,
	Pose const & pose )
{
	chain( tag->getOption< core::Size >( "chain", 2 ) );
	coord_sdev( tag->getOption< core::Real >( "coord_sdev", 1.0 ) );
	utility::vector1< core::Size > v1 = core::pose::get_resnum_list( tag, "resnums", pose );
	for ( core::Size const r : v1 ) { add_residue( r ); }
}

core::Size
AddSidechainConstraintsToHotspots::chain() const{
	return chain_;
}

void
AddSidechainConstraintsToHotspots::chain( core::Size const c ){
	chain_ = c;
}

core::Real
AddSidechainConstraintsToHotspots::coord_sdev() const{
	return coord_sdev_;
}

void
AddSidechainConstraintsToHotspots::coord_sdev( core::Real const c ){
	coord_sdev_ = c;
}

void
AddSidechainConstraintsToHotspots::add_residue( core::Size const res ){
	residues_.insert( res );
}

std::string AddSidechainConstraintsToHotspots::get_name() const {
	return mover_name();
}

std::string AddSidechainConstraintsToHotspots::mover_name() {
	return "AddSidechainConstraintsToHotspots";
}

void AddSidechainConstraintsToHotspots::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "chain", xsct_non_negative_integer, "Chain containing the relevant hotspots, numbered sequentially from 1", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "coord_sdev", xsct_real, "Size of coordinate standard deviations for constraints", "1.0" )
		+ XMLSchemaAttribute( "resnums", xsct_residue_number_cslist, "List of residues in residue numbering" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string AddSidechainConstraintsToHotspotsCreator::keyname() const {
	return AddSidechainConstraintsToHotspots::mover_name();
}

protocols::moves::MoverOP
AddSidechainConstraintsToHotspotsCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddSidechainConstraintsToHotspots );
}

void AddSidechainConstraintsToHotspotsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddSidechainConstraintsToHotspots::provide_xml_schema( xsd );
}


std::set< core::Size > const &
AddSidechainConstraintsToHotspots::residues() const{
	return residues_;
}

} //movers
} //protein_interface_design
} //protocols
