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
#include <core/id/AtomID.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <utility/io/ozstream.hh>

#include <utility/basic_sys_util.hh>

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
	scorefxn_(0),
	addtime_(false)
{}

DumpPdb::DumpPdb( std::string const fname ) :
	protocols::moves::Mover( DumpPdbCreator::mover_name() ),
  fname_(fname),
	scorefxn_(0),
	addtime_(false)
{}

DumpPdb::~DumpPdb() {}

void DumpPdb::apply( core::pose::Pose & pose ) {
	std::string name( fname_ );
	if (bfactor_){ /* Added by gideonla june13, prints out the the bfactor column */
		 utility::io::ozstream out(name);
		core::id::AtomID_Map< Real > bfactors;
		core::pose::initialize_atomid_map( bfactors, pose, 0.0 );
		for (core::Size seqpos=1;seqpos<=pose.total_residue();++seqpos) { // go over all pose residues
			for( core::Size idx = 1; idx <= pose.residue( seqpos ).natoms(); ++idx ){
		 	  core::id::AtomID atom_id(idx,seqpos);
		 	 bfactors[atom_id] = pose.pdb_info()->temperature (  seqpos, idx);
			}
		} //end for loop over residues
		core::io::pdb::dump_bfactor_pdb( pose, bfactors, out );
		return;
	}
	if ( addtime_ ) {
#ifdef USEMPI
		name += "_" + utility::to_string(utility::mpi_rank()); 
#endif
		name += "_" + utility::timestamp_short() + ".pdb";
		TR << "Dumping PDB " << name << std::endl;
	}
	if ( scorefxn_ ) {
		pose.dump_scored_pdb( name, *scorefxn_ );
	}	else {
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
	bfactor_ = ( tag->getOption<bool>( "bfactor", false ) );
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}

	tag_time( tag->getOption<bool>( "tag_time", false ) );
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
