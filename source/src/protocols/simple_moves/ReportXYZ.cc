// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ReportXYZ.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/ReportXYZ.hh>
#include <protocols/simple_moves/ReportXYZCreator.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.simple_moves.ReportXYZ" );

#include <utility/tag/Tag.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/exit.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

ReportXYZ::ReportXYZ()
: moves::Mover("ReportXYZ")
{
}

void
ReportXYZ::apply( Pose & pose )
{
	using namespace core::pose;
	numeric::xyzVector<core::Real> selected_pos(0.,0.,0.);
	if ( resnum_ !=999999 ) {
		selected_pos = pose.residue(resnum_).xyz("CA");
	} else {
		core::Size chain_id = core::pose::get_chain_id_from_chain(chain_,pose);
		core::conformation::Conformation const & conf = pose.conformation();
		core::Size const chain_begin( conf.chain_begin( chain_id ) );
		core::Size const chain_end(   conf.chain_end(   chain_id ) );
		core::Size selected_res=9999;
		if ( term_ =="n_term" ) {
			selected_res = chain_begin;
		}
		if ( term_=="c_term" ) {
			selected_res = chain_end;
		}
		selected_pos = pose.residue(selected_res).xyz("CA");
	}
	setPoseExtraScore( pose, "X", selected_pos.x());
	setPoseExtraScore( pose, "Y", selected_pos.y());
	setPoseExtraScore( pose, "Z", selected_pos.z());
}

moves::MoverOP
ReportXYZ::clone() const
{
	return moves::MoverOP( new ReportXYZ( *this ) );
}

moves::MoverOP
ReportXYZ::fresh_instance() const
{
	return moves::MoverOP( new ReportXYZ );
}

void
ReportXYZ::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	resnum_ =999999;
	term_ = "none";
	chain_ = tag->getOption<std::string>( "chain", "A" );
	if ( !tag->hasOption("resnum")&& !tag->hasOption("term") ) {
		utility_exit_with_message("User input of either resnum or term must be entered");
	}
	if ( tag->hasOption("term") ) {
		term_ = tag->getOption<std::string>( "term");
	}
	if ( tag->hasOption("resnum") ) {
		resnum_  = tag->getOption<core::Size>( "resnum" );
	}
}

std::string ReportXYZ::get_name() const {
	return mover_name();
}

std::string ReportXYZ::mover_name() {
	return "ReportXYZ";
}

void ReportXYZ::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "term" , xs_string , "Which terminus to get the xyz postions options n_term or c_term" )
		+ XMLSchemaAttribute( "resnum" , xsct_positive_integer , "Which residue to get the xyz of" )
		+ XMLSchemaAttribute::attribute_w_default( "chain" , xs_string , "which chain if there are multiple chains" , "A" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Adds the X,Y,Z position to the scorefile", attlist );

}


std::string ReportXYZCreator::keyname() const {
	return ReportXYZ::mover_name();
}

protocols::moves::MoverOP
ReportXYZCreator::create_mover() const {
	return protocols::moves::MoverOP( new ReportXYZ );
}

void ReportXYZCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReportXYZ::provide_xml_schema( xsd );
}


} // simple_moves
} // protocols
