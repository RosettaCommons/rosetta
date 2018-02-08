// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/WriteSSEMover.cc
/// @brief Writes SSE assignation from DSSP or prediction from PSIPRED as REMARK.
/// @author Jaume Bonet (jaume.bonet@gmail.com)


//Unit Headers
#include <protocols/simple_moves/WriteSSEMover.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/simple_moves/WriteSSEMoverCreator.hh>
//Project Headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>
#include <core/io/external/PsiPredInterface.hh>
#include <core/scoring/dssp/Dssp.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
#include <core/pose/extra_pose_info_util.hh>

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace core::scoring;


WriteSSEMover::WriteSSEMover():
	psipred_interface_(),
	cmd_( std::string() ),
	dssp_( false )
{}

WriteSSEMover::~WriteSSEMover() = default;


void WriteSSEMover::apply( core::pose::Pose &pose ) {

	std::string dssp_ss(pose.size(), 'L');
	if ( dssp_ ) {
		core::scoring::dssp::Dssp dssp( pose );
		dssp_ss = dssp.get_dssp_secstruct();
		core::pose::add_comment( pose, "DSSP", dssp_ss );
	}
	if ( cmd_ != "" ) {
		core::io::external::PsiPredResult const psipred_result = psipred_interface_->run_psipred( pose, dssp_ss );
		core::pose::add_comment( pose, "PSIPRED", psipred_result.pred_ss );
	}
}

void WriteSSEMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	cmd( tag->getOption< std::string >( "cmd", "" ) );
	dssp( tag->getOption< bool >( "dssp", false ) );
}

void WriteSSEMover::cmd( std::string const & cmd ) {
	cmd_ = cmd;
	if ( cmd_ != "" ) {
		psipred_interface_ = core::io::external::PsiPredInterfaceOP( new core::io::external::PsiPredInterface( cmd_ ) );
	}
}

std::string WriteSSEMover::get_name() const{
	return "WriteSSEMover";
}

protocols::moves::MoverOP
WriteSSEMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new WriteSSEMover );
}
std::string WriteSSEMoverCreator::keyname() const {
	return WriteSSEMover::mover_name();
}

std::string WriteSSEMover::mover_name() {
	return "WriteSSEMover";
}

void WriteSSEMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "cmd", xs_string, "Path to PSIPRED. If not provided, psipred is not executed", "")
		+ XMLSchemaAttribute::attribute_w_default( "dssp", xsct_rosetta_bool, "Set to true to get DSSP prediction (default false)", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Write PSIPRED and/or DSSP to silent file as REMARK or to the PDB with "
		"-pdb_comments option set to true.", attlist );
}

void WriteSSEMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	WriteSSEMover::provide_xml_schema( xsd );
}

}
}
