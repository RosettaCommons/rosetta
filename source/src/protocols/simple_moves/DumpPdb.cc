// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/DumpPdb.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/simple_moves/DumpPdb.hh>
#include <protocols/simple_moves/DumpPdbCreator.hh>


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

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.DumpPdb" );

// XRW TEMP std::string
// XRW TEMP DumpPdbCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DumpPdb::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DumpPdbCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new DumpPdb );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DumpPdb::mover_name()
// XRW TEMP {
// XRW TEMP  return "DumpPdb";
// XRW TEMP }

DumpPdb::DumpPdb():
	protocols::moves::Mover( DumpPdb::mover_name() ),
	fname_("dump.pdb"),
	scorefxn_(/* 0 */),
	addtime_(false)
{}

DumpPdb::DumpPdb( std::string  fname ) :
	protocols::moves::Mover( DumpPdb::mover_name() ),
	fname_(std::move(fname)),
	scorefxn_(/* 0 */),
	addtime_(false)
{}

DumpPdb::~DumpPdb() = default;

void DumpPdb::apply( core::pose::Pose & pose ) {
	std::string name( fname_ );
	if ( addtime_ ) {
#ifdef USEMPI
		name += "_" + utility::to_string(utility::mpi_rank());
#endif
		name += "_" + utility::timestamp_short() + ".pdb";
		TR << "Dumping PDB " << name << std::endl;
	}
	if ( scorefxn_ ) {
		pose.dump_scored_pdb( name, *scorefxn_ );
	} else {
		pose.dump_pdb( name );
	}
}

void DumpPdb::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn) {
	scorefxn_ = scorefxn;
}
void
DumpPdb::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	fname_ = tag->getOption<std::string>( "fname", "dump.pdb" );

	//JAB - XRW - bfactors are automatically output now.
	if ( tag->hasOption("bfactor") ) {
		TR.Warning << "Bfactor output by default now.  Tag not needed.  Continueing." << std::endl;
	}

	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}

	tag_time( tag->getOption<bool>( "tag_time", false ) );
	TR.Warning << "DEFINED DUMP_PDB MOVER. THIS IS USUALLY ONLY GOOD FOR DEBUGGING." << std::endl;
	TR << "dump_pdb mover with filename "<<fname_<<std::endl;
}

// XRW TEMP std::string
// XRW TEMP DumpPdb::get_name() const {
// XRW TEMP  return DumpPdb::mover_name();
// XRW TEMP }

std::string DumpPdb::get_name() const {
	return mover_name();
}

std::string DumpPdb::mover_name() {
	return "DumpPdb";
}

void DumpPdb::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "fname" , xs_string , "Filename of dumped PDB." , "dump.pdb" )
		+ XMLSchemaAttribute::attribute_w_default( "tag_time" , xsct_rosetta_bool , "If true, adds timestamp to name of pdb file." , "false" ) ;

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist) ;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Dumps a pdb. Recommended ONLY for debugging as you can't change the name of the file during a run, although if tag_time is true a timestamp with second resolution will be added to the filename, allowing for a limited amount of multi-dumping. If scorefxn is specified, a scored pdb will be dumped.", attlist );
}

std::string DumpPdbCreator::keyname() const {
	return DumpPdb::mover_name();
}

protocols::moves::MoverOP
DumpPdbCreator::create_mover() const {
	return protocols::moves::MoverOP( new DumpPdb );
}

void DumpPdbCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DumpPdb::provide_xml_schema( xsd );
}


} //simple_moves
} //protocols
