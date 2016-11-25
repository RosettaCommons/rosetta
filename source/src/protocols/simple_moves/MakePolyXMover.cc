// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./protocols/simple_moves/MakePolyXMover.cc
/// @brief  convert pose to poly XXX: any amino acid, default poly Ala
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit headers
#include <protocols/simple_moves/MakePolyXMover.hh>
#include <protocols/simple_moves/MakePolyXMoverCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.MakePolyXMover" );

namespace protocols {
namespace simple_moves {

// XRW TEMP std::string
// XRW TEMP MakePolyXMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return MakePolyXMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MakePolyXMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new MakePolyXMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP MakePolyXMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "MakePolyX";
// XRW TEMP }

MakePolyXMover::MakePolyXMover():
	protocols::moves::Mover( MakePolyXMover::mover_name() ),
	aa_( "ALA" ), //SML Jul 10 2016, changed AA to ALA.  Unclear if a typo or deliberate non amino acid?
	keep_pro_( false ),
	keep_gly_( true ),
	keep_disulfide_cys_( false ),
	selector_( new core::select::residue_selector::TrueResidueSelector )
{}

MakePolyXMover::MakePolyXMover( std::string aa, bool keep_pro, bool keep_gly, bool keep_disulfide_cys ):
	protocols::moves::Mover( MakePolyXMover::mover_name() ),
	aa_(std::move( aa )),
	keep_pro_( keep_pro ),
	keep_gly_( keep_gly ),
	keep_disulfide_cys_( keep_disulfide_cys ),
	selector_( new core::select::residue_selector::TrueResidueSelector )
{}

MakePolyXMover::~MakePolyXMover() = default;

/// @brief clone this object
protocols::moves::MoverOP MakePolyXMover::clone() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::MakePolyXMover( *this ) );
}

/// @brief create this type of object
protocols::moves::MoverOP MakePolyXMover::fresh_instance() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::MakePolyXMover() );
}

/// @details virtual main
void MakePolyXMover::apply( Pose & pose )
{
	using core::select::residue_selector::ResidueSubset;
	ResidueSubset const subset = selector_->apply( pose );
	// flip to poly-ala-gly-pro-disulf pose
	utility::vector1< Size > protein_residues;
	for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
		if ( subset[ i ] && pose.residue( i ).is_protein() ) {
			protein_residues.push_back( i );
		}
	}
	using protocols::toolbox::pose_manipulation::construct_poly_XXX_pose;
	TR << "Pose is converted to poly " << aa_ << std::endl;
	construct_poly_XXX_pose( aa_, pose, protein_residues, keep_pro_, keep_gly_, keep_disulfide_cys_ );
}

/// @brief parse xml
void
MakePolyXMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{
	aa_ = tag->getOption<std::string>( "aa", "ALA" );
	keep_pro_  = tag->getOption<bool>( "keep_pro", 0 );
	keep_gly_  = tag->getOption<bool>( "keep_gly", 1 );
	keep_disulfide_cys_  = tag->getOption<bool>( "keep_disulfide_cys", 0 );

	ResidueSelectorCOP selector = core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) selector_ = selector;

	TR << "MakePolyXMover was loaded" << std::endl;

	if ( keep_pro_ || keep_gly_ || keep_disulfide_cys_ ) {
		TR << "but keep AA types of ";
		if ( keep_pro_ ) TR << "Pro ";
		if ( keep_gly_ ) TR << "Gly  ";
		if ( keep_disulfide_cys_ ) TR << "Disulfide Cys";
		TR << std::endl;
	}

}

std::string MakePolyXMover::get_name() const {
	return mover_name();
}

std::string MakePolyXMover::mover_name() {
	return "MakePolyX";
}

void MakePolyXMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"aa", xs_string,
		"using amino acid type for converting",
		"ALA" )
		+ XMLSchemaAttribute::attribute_w_default(
		"keep_pro", xsct_rosetta_bool,
		"Pro is not converted to XXX",
		"false" )
		+ XMLSchemaAttribute::attribute_w_default(
		"keep_gly", xsct_rosetta_bool,
		"Gly is not converted to XXX",
		"true" )
		+ XMLSchemaAttribute::attribute_w_default(
		"keep_disulfide_cys", xsct_rosetta_bool,
		"disulfide CYS is not converted to XXX",
		"false" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist );
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Convert pose into poly XXX ( XXX can be any amino acid )",
		attlist );
}

std::string MakePolyXMoverCreator::keyname() const {
	return MakePolyXMover::mover_name();
}

protocols::moves::MoverOP
MakePolyXMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakePolyXMover );
}

void MakePolyXMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakePolyXMover::provide_xml_schema( xsd );
}


}  // namespace simple_moves
}  // namespace protocols
