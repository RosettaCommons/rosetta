// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/StructureRestrictor.cc
///
/// @brief  lookup relevant chains for a structure in a table.
/// @author Matthew O'Meara

/// This should probably be a pilot app, but the way Rosetta Scripts
/// is set up, it can't be in the pilot apps

#include <boost/algorithm/string.hpp>
#include <core/conformation/Conformation.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/StructureRestrictor.hh>


#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>


#include <string>

#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using namespace std;
using namespace boost;
using namespace core;
using namespace pose;
using namespace protocols::jd2;

static thread_local basic::Tracer TR_SR( "protocols.moves.StructureRestrictor" );

namespace protocols {
namespace moves {

StructureRestrictor::StructureRestrictor():
	Mover("StructureRestrictor"),
	//  relevant_chains_fname( basic::options::option[ basic::options::OptionKeys::StructureRestrictor::relevant_chains].value() ),
	initialized( false )
{}

StructureRestrictor::StructureRestrictor( string const & name):
	Mover(name),
	//  relevant_chains_fname( basic::options::option[ basic::options::OptionKeys::StructureRestrictor::relevant_chains].value() ),
	initialized( false )
{}

StructureRestrictor::StructureRestrictor( StructureRestrictor const & src):
	//utility::pointer::ReferenceCount(),
	Mover(src)
{
	chain_map = std::map< std::string const, std::string const >( src.chain_map );
	relevant_chains_fname = src.relevant_chains_fname;
	initialized = src.initialized;
}

StructureRestrictor::~StructureRestrictor(){}

MoverOP StructureRestrictor::fresh_instance() const { return MoverOP( new StructureRestrictor ); }

MoverOP StructureRestrictor::clone() const
{
	return MoverOP( new StructureRestrictor( *this ) );
}


// So this this can be called from RosettaScripts
void
StructureRestrictor::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/ )
{
	if ( tag->hasOption("relevant_chains") ) {
		relevant_chains_fname = tag->getOption<string>("relevant_chains");
	}
}

void
StructureRestrictor::setup_relevant_chains(
	string const & relevant_chains_fname,
	map<string const, string const> & chain_map
){
	if ( relevant_chains_fname.length() == 0 ) {
		TR_SR.Error << " Cannot open relevant_chains_file '"<< relevant_chains_fname << "'" << endl;
		return;
	}
	utility::io::izstream relevant_chains_file( relevant_chains_fname );
	if ( !relevant_chains_file ) {
		TR_SR.Error << " Cannot open relevant_chains_file '"<< relevant_chains_fname << "'" << endl;
		return;
	} else {
		TR_SR << "Reading in relevant chains from file '"<< relevant_chains_fname << "'" << endl;
	}
	string line;
	//string chains = "*";
	vector<string> tokens;
	getline(relevant_chains_file,line); // header
	while ( getline( relevant_chains_file, line ) ) {
		string tab("\t");
		split(tokens, line, is_any_of(tab) );
		chain_map.insert(std::pair<string const, string const>(tokens[0], tokens[1]));
	}
	initialized = true;
}

// this is a hack because poses do not have canonical names!
string
StructureRestrictor::pose_name(Pose const & pose){
	//silent files and pdbs set the name of the pose differently
	string name = "No_Name_Found";
	if ( pose.pdb_info() ) {
		name = pose.pdb_info()->name();
	} else if ( pose.data().has( datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		name = static_cast< basic::datacache::CacheableString const & >
			( pose.data().get( datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	} else {
		name = JobDistributor::get_instance()->current_job()->input_tag();
	}
	return name;
}


void
StructureRestrictor::apply( Pose& pose ){
	if ( !initialized ) {
		setup_relevant_chains(relevant_chains_fname, chain_map);
	}


	string const & name = pose_name(pose);
	map<string const, string const>::iterator i(chain_map.find(name));
	if ( i == chain_map.end() ) {
		TR_SR << "No chain information found for structure " << name << "." << endl;
		return;
	}
	string chains = i->second;
	TR_SR << "Restricting structure " << name << " to chains " << chains << "." << endl;
	Size res_begin_delete = 1;
	for ( Size i=1; i <= pose.total_residue(); ++i ) {
		//INVARIANT: if we're in a stretch to delete then res_begin_delete
		//indicates the first residue in this stretch to delete
		if ( chains.find( pose.pdb_info()->chain(i), 0) != string::npos ) {
			//keep this position
			if ( res_begin_delete != i ) {
				pose.conformation().delete_residue_range_slow(res_begin_delete, i-1);
			}
			res_begin_delete = i+1;
		}
	}
	// don't for get the last section to delete
	if ( res_begin_delete <= pose.total_residue() ) {
		pose.conformation().delete_residue_range_slow(res_begin_delete, pose.total_residue());
	}
}


} // moves namespace
} // protocols namespace
