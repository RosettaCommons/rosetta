// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/AddPDBInfoMover.cc
/// @brief add PDBInfo to a pose
/// @author Jonathan Weinstein jonathan.weinstein@weizmann.ac.il

// Unit headers
#include <protocols/simple_moves/AddPDBInfoMover.hh>
#include <protocols/simple_moves/AddPDBInfoMoverCreator.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <core/id/AtomID.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <utility/io/ozstream.hh>

#include <utility/sys_util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#ifdef USEMPI
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#endif

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.simple_moves.AddPDBInfoMover" );

// XRW TEMP std::string
// XRW TEMP AddPDBInfoMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AddPDBInfoMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddPDBInfoMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AddPDBInfoMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddPDBInfoMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "AddPDBInfoMover";
// XRW TEMP }

AddPDBInfoMover::AddPDBInfoMover():
	protocols::moves::Mover( AddPDBInfoMover::mover_name() )
{}

AddPDBInfoMover::~AddPDBInfoMover() = default;

void AddPDBInfoMover::apply( core::pose::Pose & pose ) {
	//core::pose::PDBInfoOP pdb_info = new core::pose::PDBInfo(pose, true);
	//pose.pdb_info( pdb_info );
	pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );
	//TR << "set PDBInfo to " << std::endl << pose.pdb_info()->show( std::ostream & output=std::cout ) << std::endl;
}

void
AddPDBInfoMover::parse_my_tag( TagCOP const, basic::datacache::DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
}

protocols::moves::MoverOP
AddPDBInfoMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AddPDBInfoMover );
}

protocols::moves::MoverOP
AddPDBInfoMover::clone() const {
	return( protocols::moves::MoverOP( new AddPDBInfoMover( *this ) ) );
}
// XRW TEMP std::string
// XRW TEMP AddPDBInfoMover::get_name() const {
// XRW TEMP  return AddPDBInfoMover::mover_name();
// XRW TEMP }

std::string AddPDBInfoMover::get_name() const {
	return mover_name();
}

std::string AddPDBInfoMover::mover_name() {
	return "AddPDBInfoMover";
}

void AddPDBInfoMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//attlist + XMLSchemaAttribute::attribute_w_default( "fname" , xs_string , "Filename of dumped PDB." , "dump.pdb" )
	//  + XMLSchemaAttribute::attribute_w_default( "tag_time" , xsct_rosetta_bool , "If true, adds timestamp to name of pdb file." , "false" ) ;

	//protocols::rosetta_scripts::attributes_for_parse_score_function( attlist) ;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "create a PDBInfo for pose, from pose. useful when working with denovo protocls", attlist );
}

std::string AddPDBInfoMoverCreator::keyname() const {
	return AddPDBInfoMover::mover_name();
}

protocols::moves::MoverOP
AddPDBInfoMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddPDBInfoMover );
}

void AddPDBInfoMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddPDBInfoMover::provide_xml_schema( xsd );
}


} //simple_moves
} //protocols
