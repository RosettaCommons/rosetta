// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.simple_moves.DumpPdb" );

std::string
DumpPdbCreator::keyname() const
{
	return DumpPdbCreator::mover_name();
}

protocols::moves::MoverOP
DumpPdbCreator::create_mover() const {
	return new DumpPdb;
}

std::string
DumpPdbCreator::mover_name()
{
	return "DumpPdb";
}

DumpPdb::DumpPdb():
	protocols::moves::Mover( DumpPdbCreator::mover_name() ),
	fname_("dump.pdb"),
	scorefxn_(0)
{}

DumpPdb::DumpPdb( std::string const fname ) :
	protocols::moves::Mover( DumpPdbCreator::mover_name() ),
  fname_(fname),
	scorefxn_(0)
{}

DumpPdb::~DumpPdb() {}

void DumpPdb::apply( core::pose::Pose & pose ) {
	if ( scorefxn_ ) {
		pose.dump_scored_pdb( fname_, *scorefxn_ );
	}	else {
		pose.dump_pdb( fname_ );
	}
}

void DumpPdb::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn) {
	scorefxn_ = scorefxn;
}
void
DumpPdb::parse_my_tag( TagPtr const tag, DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	fname_ = tag->getOption<std::string>( "fname", "dump.pdb" );
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}
	TR<<"dump pdb\n";
	TR<<"WARNING: DEFINED DUMP_PDB MOVER. THIS IS USUALLY ONLY GOOD FOR DEBUGGING.\n";
	TR<<"with filename "<<fname_<<std::endl;
}

std::string
DumpPdb::get_name() const {
	return DumpPdbCreator::mover_name();
}

} //simple_moves
} //protocols
