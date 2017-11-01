// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/AddChainBreak.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#include <boost/range/algorithm.hpp>

// Unit headers
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>
#include <protocols/protein_interface_design/movers/AddChainBreakCreator.hh>
#include <core/pose/variant_util.hh>
// Package headers
#include <core/chemical/VariantType.hh>
// Project headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
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

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.AddChainBreak" );

// XRW TEMP std::string
// XRW TEMP AddChainBreakCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AddChainBreak::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddChainBreakCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AddChainBreak );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddChainBreak::mover_name()
// XRW TEMP {
// XRW TEMP  return "AddChainBreak";
// XRW TEMP }

AddChainBreak::AddChainBreak() :
	protocols::moves::Mover( AddChainBreak::mover_name() ),
	resnum_( "" ),
	change_foldtree_( true ),
	change_conformation_( false ),
	find_automatically_( false ),
	automatic_distance_cutoff_( 2.5 ),
	remove_( false )
{}

AddChainBreak::~AddChainBreak() {}

protocols::moves::MoverOP
AddChainBreak::clone() const {
	return (protocols::moves::MoverOP( new AddChainBreak( *this ) ) );
}

void
AddChainBreak::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	/// resnum & pdb_num are now equivalent
	if ( tag->hasOption( "resnum" ) ) {
		resnum_ = tag->getOption< std::string > ("resnum" );
	} else if ( tag->hasOption( "pdb_num" ) ) {
		resnum_ = tag->getOption< std::string > ("pdb_num" );
	}
	if ( tag->hasOption( "find_automatically" ) ) {
		find_automatically( tag->getOption< bool >( "find_automatically" ) );
		automatic_distance_cutoff( tag->getOption< core::Real >( "distance_cutoff", 2.5 ));
	}
	remove( tag->getOption< bool >( "remove", false ) );
	change_foldtree( tag->getOption< bool >( "change_foldtree", true ) );
	change_conformation( tag->getOption< bool >( "change_conformation", false ) );
	TR<<"resnum: "<<resnum_<<" change foldtree "<<change_foldtree()<<" find cutpoints automatically "<<find_automatically()<<" remove: "<<remove()<<std::endl;
}//end parse my tag

void
AddChainBreak::apply( core::pose::Pose & pose )
{
	using namespace core::chemical;
	using namespace pose;
	core::kinematics::FoldTree f( pose.fold_tree() );

	if ( resnum() != "" ) {
		core::Size const resn( core::pose::parse_resnum( resnum(), pose ) );

		if ( change_foldtree() ) {
			f.new_jump( resn, resn+1, resn );
		}

		if ( change_conformation() ) {
			pose.conformation().insert_chain_ending(resn);
		}

		if ( pose.residue(resn  ).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(pose,resn);
		if ( pose.residue(resn+1).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(pose,resn+1);
		if ( remove() ) {
			remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, resn );
			remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, resn +1);
		} else {
			add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, resn );
			add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, resn +1);
		}
	}
	utility::vector1< core::Size > cuts;
	cuts.clear();
	if ( find_automatically() ) {
		for ( core::Size i = 1; i < pose.size(); ++i ) {
			core::Real const distance( pose.residue( i ).xyz( "C" ).distance( pose.residue( i + 1 ).xyz( "N" ) ) );
			if ( distance >= automatic_distance_cutoff() ) {
				cuts.push_back( i );
				TR<<"Detecting cut at "<<i<<" with distance "<<distance<<std::endl;
			}
		}
	}
	for ( core::Size const res : cuts ) {
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, res );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, res +1);
		if ( change_foldtree() &&  !f.is_cutpoint(res) ) {
			f.new_jump( res, res+1, res );
		}

		if ( change_conformation() && boost::range::count(pose.conformation().chain_endings(), res) == 0 ) {
			pose.conformation().insert_chain_ending(res);
		}
	}
	if ( change_foldtree() ) {
		pose.fold_tree( f );
		TR<<"New fold tree: "<<pose.fold_tree()<<std::endl;
	}
}

// XRW TEMP std::string
// XRW TEMP AddChainBreak::get_name() const {
// XRW TEMP  return AddChainBreak::mover_name();
// XRW TEMP }

std::string AddChainBreak::get_name() const {
	return mover_name();
}

std::string AddChainBreak::mover_name() {
	return "AddChainBreak";
}

void AddChainBreak::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"remove", xsct_rosetta_bool, "Residue numbering (pdb); provide this or resnum",
		"false")
		+ XMLSchemaAttribute::attribute_w_default(
		"change_foldtree", xsct_rosetta_bool, "Add chainbreak as new jump in fold tree.",
		"true")
		+ XMLSchemaAttribute::attribute_w_default(
		"change_conformation", xsct_rosetta_bool, "Add chainbreak as chain ending in pose conformation.",
		"true")
		+ XMLSchemaAttribute( "resnum", xs_string, "Residue numbering (pdb); provide this or resnum" )
		+ XMLSchemaAttribute( "find_automatically", xsct_rosetta_bool, "Automatically find chain breaks" )
		+ XMLSchemaAttribute::attribute_w_default( "distance_cutoff", xsct_real, "Distance cutoff past which chain breaks will be automatically assigned", "2.5" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string AddChainBreakCreator::keyname() const {
	return AddChainBreak::mover_name();
}

protocols::moves::MoverOP
AddChainBreakCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddChainBreak );
}

void AddChainBreakCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddChainBreak::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
