// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/protein_interface_design/movers/ShoveResidueMover.cc
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/ShoveResidueMover.hh>
#include <protocols/protein_interface_design/movers/ShoveResidueMoverCreator.hh>

// Project headers
#include <utility/graph/Graph.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/selection.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/protein_interface_design/util.hh>
#include <protocols/rosetta_scripts/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.ShoveResidueMover" );

// XRW TEMP std::string
// XRW TEMP ShoveResidueMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ShoveResidueMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ShoveResidueMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ShoveResidueMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ShoveResidueMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "ShoveResidueMover";
// XRW TEMP }


ShoveResidueMover::ShoveResidueMover() :
	protocols::moves::Mover( ShoveResidueMover::mover_name() )
{}

ShoveResidueMover::ShoveResidueMover( Size resnum ) :
	protocols::moves::Mover( ShoveResidueMover::mover_name() ),
	shove_residues_( new core::select::residue_selector::ResidueIndexSelector( resnum ) )
{}

void
ShoveResidueMover::apply ( pose::Pose & pose )
{
	//using namespace rotamer_set;
	using namespace core::scoring;
	using namespace core::pack::task;
	using namespace core::pack::rotamer_set;
	runtime_assert( shove_residues_ != nullptr );

	for ( core::Size const resid : core::select::get_residues_from_subset( shove_residues_->apply(pose) ) ) {
		if ( pose.residue(resid).name3() == "GLY" ) {
			//The SHOVE_BB patch does not properly work with Glycine.  Looking at it, it's not clear that it even
			//Makes sense to use this mover with Glycine.  Given this, it's probably best to exit here.
			//If this isn't true, a special case should be added to the patch.  I added this if statement
			//Because running this mover on a glycine causes select_atoms_to_orient() fails during repacking
			// - Sam DeLuca

			utility_exit_with_message("ERROR: ShoveResidueMover does not currently work properly with glycines, and residue " + utility::to_string(resid) + " is a glycine.");

		}
		if ( remove_shove_variant_ ) {
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::SHOVE_BB, resid );
		} else {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::SHOVE_BB, resid );
		}
	}
}

// XRW TEMP std::string
// XRW TEMP ShoveResidueMover::get_name() const {
// XRW TEMP  return ShoveResidueMover::mover_name();
// XRW TEMP }

void
ShoveResidueMover::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	remove_shove_variant_ = tag->getOption<bool>( "remove_shove_variant", false );
	if ( tag->hasOption( "shove" ) ) {
		shove_residues_ = core::pose::get_resnum_selector( tag, "shove" );
	} else {
		using namespace core::select::residue_selector;
		std::string resnum = core::pose::get_resnum_string( tag );
		shove_residues_ = ResidueSelectorOP( new ResidueIndexSelector( resnum ) );
	}
}

std::string ShoveResidueMover::get_name() const {
	return mover_name();
}

std::string ShoveResidueMover::mover_name() {
	return "ShoveResidueMover";
}

void ShoveResidueMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	core::pose::attributes_for_get_resnum_string( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "remove_shove_variant", xsct_rosetta_bool, "Remove the shove variant from the residue in question?", "false" )
		+ XMLSchemaAttribute( "shove", xsct_refpose_enabled_residue_number_cslist, "Residues to which to add shove variant, especially if resnum isn't provided" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string ShoveResidueMoverCreator::keyname() const {
	return ShoveResidueMover::mover_name();
}

protocols::moves::MoverOP
ShoveResidueMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ShoveResidueMover );
}

void ShoveResidueMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ShoveResidueMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols

