// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/SheetDB.cc
/// @brief Heirarchical Database of beta sheets
/// @detailed
/// @author Tom Linsky

// Unit Headers
#include <protocols/denovo_design/components/SheetDB.hh>

// Protocol Headers
#include <protocols/denovo_design/util.hh>
#include <protocols/fldsgn/topology/Sheet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/util.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

// Core Headers
#include <core/conformation/ResidueFactory.hh>
#include <core/graph/Graph.hh>
#include <core/graph/graph_util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <boost/lexical_cast.hpp>
#include <map>
#include <stack>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.SheetDB" );

namespace protocols {
namespace denovo_design {
namespace components {

///////////////////////////////////////////////////////////////////////////////
/// SheetDB Implementations
///////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
SheetDB::SheetDB():
	idealize_( true ),
	db_path_( "" )
{
	sheet_db_.clear();
}

/// @brief adds a sheetlist to the database with given descriptors
void
SheetDB::set_sheetlist(
	SheetList const & list,
	core::Size const nstrands,
	std::string const & in_orientations,
	std::string const & in_lengths )
{
	std::pair< std::string, std::string > canonical = canonicalize( in_orientations, in_lengths );

	set_sheetlist_as_is( list, nstrands, canonical.first, canonical.second );

	// now insert the opposite of the canonical
	std::string const rev_len = reverse_lengths( canonical.first, canonical.second );
	std::string const rev_or = reverse_orientations( canonical.first );
	SheetList newlist;
	for ( SheetList::const_iterator s=list.begin(), end=list.end(); s != end; ++s ) {
		newlist.push_back( reverse_chains( **s ) );
	}
	set_sheetlist_as_is( newlist, nstrands, rev_or, rev_len );
}

/// @brief adds a sheet list as it is presented
void
SheetDB::set_sheetlist_as_is(
	SheetList const & list,
	core::Size const nstrands,
	std::string const & orientations,
	std::string const & lengths )
{
	MapByStrands::iterator m1 = sheet_db_.find( nstrands );
	if ( m1 == sheet_db_.end() ) {
		MapByOrientations newmap;
		sheet_db_[nstrands] = newmap;
		m1 = sheet_db_.find( nstrands );
	}

	assert( m1 != sheet_db_.end() );
	MapByOrientations & omap = m1->second;
	MapByOrientations::iterator m2 = omap.find( orientations );
	if ( m2 == omap.end() ) {
		MapByLengths newmap;
		omap[orientations] = newmap;
		m2 = omap.find( orientations );
	}

	assert( m2 != omap.end() );
	MapByLengths & lmap = m2->second;
	lmap[lengths] = list;

}

/// @brief tells whether the given pose matches up with the given descriptor
bool
pose_matches_description( core::pose::Pose const & pose, std::string const & in_lengths_shifts )
{
	std::pair< utility::vector1< core::Size >, utility::vector1< int > > strands = parse_lengths( in_lengths_shifts );
	assert( strands.first.size() == strands.second.size() );
	bool matches = true;
	if ( strands.first.size() != pose.conformation().num_chains() ) {
		matches = false;
	} else {
		// check for proper terminal variants and chain endings
		for ( core::Size c=1; c<=pose.conformation().num_chains(); ++c ) {
			core::Size const chain_start = pose.conformation().chain_begin(c);
			core::Size const chain_end = pose.conformation().chain_end(c);
			core::Size const chain_len = chain_end - chain_start + 1;
			if ( chain_len != strands.first[c] + 2 ) {
				matches = false;
				TR << " length of chain " << c << " is wrong -- this doesn't match up" << std::endl;
				break;
			}

			core::Size const predicted_chain_end = chain_start + strands.first[c] + 1;
			if ( predicted_chain_end != chain_end ) {
				matches = false;
				TR << "Predicted chain end for chain " << c << " is wrong " << predicted_chain_end << " vs actual " << chain_end << std::endl;
				break;
			}

			if ( !pose.residue(chain_start).is_lower_terminus() ) {
				matches = false;
				TR << " residue " << chain_start << " does not have lower terminal variant" << std::endl;
				break;
			}

			if ( !pose.residue(predicted_chain_end).is_upper_terminus() ) {
				matches = false;
				TR << " residue " << predicted_chain_end << " does not have upper terminal variant" << std::endl;
				break;
			}

			// all other residues in chain should not be terminal
			for ( core::Size r=chain_start+1; r<chain_end; ++r ) {
				if ( pose.residue(r).is_terminus() ) {
					matches = false;
					TR << " residue " << r << " has a terminus variant and it shouldn't." << std::endl;
					break;
				}
			}
			if ( !matches ) {
				break;
			}
		}
	}
	return matches;
}

/// @brief add a sheet to the database based on descriptors
void
SheetDB::add_sheet(
	core::pose::PoseOP in_pose,
	core::Size const nstrands,
	std::string const & in_orientations,
	std::string const & in_lengths_shifts,
	bool const check_descriptor_against_pose = true )
{
	assert( in_pose );

	//TODO: check pose vs. descriptors
	std::string matching_orient_str = in_orientations;
	std::string matching_len_str = in_lengths_shifts;
	if ( check_descriptor_against_pose ) {
		if ( !pose_matches_description( *in_pose, in_lengths_shifts ) ) {
			matching_orient_str = reverse_orientations( in_orientations );
			matching_len_str = reverse_lengths( in_orientations, in_lengths_shifts );
		}
		if ( !pose_matches_description( *in_pose, matching_len_str ) ) {
			in_pose->dump_pdb( "bad_description.pdb" );
			TR << " ADD Sheet FAILED DUE TO BAD description: " << in_lengths_shifts << std::endl;
			assert( false );
			runtime_assert( false );
		}
		assert( pose_matches_description( *in_pose, matching_len_str ) );
	}

	// make canonical and add only canonical variant to DB
	std::pair< std::string, std::string > canonical = canonicalize( matching_orient_str, matching_len_str );
	std::string const & orientations = canonical.first;
	std::string const & lengths_shifts = canonical.second;

	// if canonicalization flipped the strings, we need to flip the pose
	core::pose::PoseOP pose;
	if ( ( orientations == matching_orient_str ) && ( lengths_shifts == matching_len_str ) ) {
		pose = in_pose;
	} else {
		pose = reverse_chains( *in_pose );
	}

	core::util::switch_to_residue_type_set( *pose, "centroid" );

	MapByStrands::iterator m1 = sheet_db_.find( nstrands );
	if ( m1 == sheet_db_.end() ) {
		MapByOrientations newmap;
		sheet_db_[nstrands] = newmap;
		m1 = sheet_db_.find( nstrands );
	}
	assert( m1 != sheet_db_.end() );
	MapByOrientations & omap = m1->second;
	MapByOrientations::iterator m2 = omap.find( orientations );
	if ( m2 == omap.end() ) {
		MapByLengths newmap;
		omap[orientations] = newmap;
		m2 = omap.find( orientations );
	}
	assert( m2 != omap.end() );

	MapByLengths & lmap = m2->second;
	MapByLengths::iterator m3 = lmap.find( lengths_shifts );
	if ( m3 == lmap.end() ) {
		SheetList newlist;
		lmap[lengths_shifts] = newlist;
		m3 = lmap.find( lengths_shifts );
	}
	assert( m3 != lmap.end() );

	m3->second.push_back( pose );
	TR << "Added to DB: N=" << nstrands << " orientations=" << orientations << " lengths=" << lengths_shifts << std::endl;
}

std::string
orientations_str( StrandOrientations const & orientations )
{
	std::stringstream output;
	for ( StrandOrientations::const_iterator o=orientations.begin(); o!=orientations.end(); ++o ) {
		if ( *o == UP ) output << "1";
		else if ( *o == DOWN ) output << "0";
		else {
			std::stringstream msg;
			msg << "Invalid orientation type " << *o << std::endl;
			utility_exit_with_message( msg.str() );
		}
	}
	return output.str();
}

/// @brief query the database for a set of orientations and shifts
bool
SheetDB::has_data( Lengths const & lengths, StrandOrientations const & orientations, RegisterShifts const & shifts ) const
{
	if ( ( lengths.size() != orientations.size() ) || ( lengths.size() != shifts.size() ) ) {
		std::stringstream msg;
		msg << "SheetDB: lengths, orientations, and register shift vectors given to has_data() must be the same size.  You provided lengths "
			<< lengths << ", orientations " << orientations << ", and shifts " << shifts << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return has_data( lengths.size(), orientations_str( orientations ), make_strand_info_str( lengths, shifts ) );
}

/// @brief query the database for a set of orientations and shifts
SheetList
SheetDB::sheet_list( Lengths const & lengths, StrandOrientations const & orientations, RegisterShifts const & shifts )
{
	if ( ( lengths.size() != orientations.size() ) || ( lengths.size() != shifts.size() ) ) {
		std::stringstream msg;
		msg << "SheetDB: lengths, orientations, and register shift vectors given to sheet_list() must be the same size.  You provided lengths "
			<< lengths << ", orientations " << orientations << ", and shifts " << shifts << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return sheet_list( lengths.size(), orientations_str( orientations ), make_strand_info_str( lengths, shifts ) );
}

/// @brief query the database for a set of orientations and shifts
SheetList
SheetDB::sheet_list_const( Lengths const & lengths, StrandOrientations const & orientations, RegisterShifts const & shifts ) const
{
	if ( ( lengths.size() != orientations.size() ) || ( lengths.size() != shifts.size() ) ) {
		std::stringstream msg;
		msg << "SheetDB: lengths, orientations, and register shift vectors given to sheet_list() must be the same size.  You provided lengths "
			<< lengths << ", orientations " << orientations << ", and shifts " << shifts << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return sheet_list_const( lengths.size(), orientations_str( orientations ), make_strand_info_str( lengths, shifts ) );
}

/// @brief query the database for a sheet and retrieve the list
SheetList
SheetDB::sheet_list( core::Size const nstrands, std::string const & orientations ) const
{
	MapByStrands::const_iterator it1 = sheet_db_.find( nstrands );
	if ( it1 == sheet_db_.end() ) {
		return SheetList();
	}
	assert( it1 != sheet_db_.end() );
	MapByOrientations const & omap = it1->second;
	MapByOrientations::const_iterator it2 = omap.find( orientations );
	if ( it2 == omap.end() ) {
		return SheetList();
	}
	assert( it2 != omap.end() );
	MapByLengths const & lmap = it2->second;
	SheetList retval;
	for ( MapByLengths::const_iterator it3 = lmap.begin(); it3 != lmap.end(); ++it3 ) {
		SheetList const & sheetlist = it3->second;
		for ( SheetList::const_iterator s = sheetlist.begin(); s != sheetlist.end(); ++s ) {
			retval.push_back( *s );
		}
	}
	return retval;
}

/// @brief reads sheetlist from the file database, and saves it to this database -- should always be called with canonical descriptors
void
SheetDB::load_sheetlist_from_file( core::Size const nstrands, std::string const & orientations, std::string const & strandinfo )
{
	std::pair< std::string, std::string > canonical = canonicalize( orientations, strandinfo );
	TR << "Requested: " << strandinfo << " Canonical: " << canonical.second << std::endl;

	if ( db_path_.empty() ) return;
	std::string const filename = db_path_ + "/" + boost::lexical_cast< std::string >(nstrands) + "/" +
		canonical.first + "/" + canonical.second + "/clusters_firstmember.out.gz";
	TR << "loading " << filename << std::endl;
	utility::io::izstream infile;
	infile.open( filename );
	if ( !infile.good() ) {
		return;
	}
	infile.close();

	core::io::silent::SilentFileData sfd;

	// try to read file, if there is an error, return empty list
	if ( ! sfd.read_file( filename ) ) {
		return;
	}

	// create sheetlist from file
	SheetList retval;
	for ( core::io::silent::SilentFileData::iterator s = sfd.begin(), end = sfd.end(); s != end; ++s ) {
		core::pose::PoseOP pose = core::pose::PoseOP( new core::pose::Pose() );
		s->fill_pose( *pose );
		retval.push_back( pose );
	}

	// save sheetlist
	set_sheetlist( retval, nstrands, canonical.first, canonical.second );
}

/// @brief queries the database to see if the requested data exists
bool
SheetDB::has_data( core::Size const nstrands, std::string const & orient, std::string const & strandinfo ) const
{
	MapByStrands::const_iterator it1 = sheet_db_.find( nstrands );
	if ( it1 == sheet_db_.end() ) {
		return false;
	}
	MapByOrientations::const_iterator it2 = it1->second.find( orient );
	if ( it2 == it1->second.end() ) {
		return false;
	}
	MapByLengths::const_iterator it3 = it2->second.find( strandinfo );
	if ( it3 == it2->second.end() ) {
		return false;
	}
	return true;
}

/// @brief query the database for a sheet and retrieve the list
SheetList
SheetDB::sheet_list( core::Size const nstrands, std::string const & in_orientations, std::string const & in_strandinfo )
{
	if ( has_data( nstrands, in_orientations, in_strandinfo ) ) {
		TR << "Returning existing data for " << in_strandinfo << std::endl;
		return sheet_list_const( nstrands, in_orientations, in_strandinfo );
	}

	MapByStrands::iterator it1 = sheet_db_.find( nstrands );
	if ( it1 == sheet_db_.end() ) {
		load_sheetlist_from_file( nstrands, in_orientations, in_strandinfo );
		it1 = sheet_db_.find( nstrands );
		if ( it1 == sheet_db_.end() ) {
			return SheetList();
		}
	}
	assert( it1 != sheet_db_.end() );
	MapByOrientations & omap = it1->second;
	MapByOrientations::iterator it2 = omap.find( in_orientations );
	if ( it2 == omap.end() ) {
		// try to load/save the list from file
		load_sheetlist_from_file( nstrands, in_orientations, in_strandinfo );
		it2 = omap.find( in_orientations );
		if ( it2 == omap.end() ) {
			return SheetList();
		}
	}
	assert( it2 != omap.end() );
	MapByLengths & lmap = it2->second;
	MapByLengths::iterator it3 = lmap.find( in_strandinfo );
	if ( it3 == lmap.end() ) {
		// if not found in memory, try to load the list from file
		load_sheetlist_from_file( nstrands, in_orientations, in_strandinfo );
		if ( it3 == lmap.end() ) {
			return SheetList();
		}
	}
	assert( has_data( nstrands, in_orientations, in_strandinfo  ) );
	return sheet_list_const( nstrands, in_orientations, in_strandinfo );
}

/// @brief query the database for a sheet and retrieve the list
SheetList
SheetDB::sheet_list_const( core::Size const nstrands, std::string const & in_orientations, std::string const & in_strandinfo ) const
{
	MapByStrands::const_iterator it1 = sheet_db_.find( nstrands );
	if ( it1 == sheet_db_.end() ) {
		return SheetList();
	}
	assert( it1 != sheet_db_.end() );
	MapByOrientations const & omap = it1->second;
	MapByOrientations::const_iterator it2 = omap.find( in_orientations );
	if ( it2 == omap.end() ) {
		return SheetList();
	}
	assert( it2 != omap.end() );
	MapByLengths const & lmap = it2->second;
	MapByLengths::const_iterator it3 = lmap.find( in_strandinfo );
	if ( it3 == lmap.end() ) {
		return SheetList();
	}
	assert( it3 != lmap.end() );
	return it3->second;
}

/// @brief queries the database for all lengths/shifts with nstrands strands and given orientations
utility::vector1< std::string >
SheetDB::lengths_shifts( core::Size const nstrands, std::string const & orientations ) const
{
	utility::vector1< std::string > retval;
	MapByStrands::const_iterator it1 = sheet_db_.find( nstrands );
	if ( it1 == sheet_db_.end() ) {
		return retval;
	}

	assert( it1 != sheet_db_.end() );
	MapByOrientations::const_iterator it2 = it1->second.find( orientations );
	if ( it2 == it1->second.end() ) {
		return retval;
	}

	assert( it2 != it1->second.end() );
	MapByLengths::const_iterator it3_end = it2->second.end();
	for ( MapByLengths::const_iterator it3 = it2->second.begin(); it3 != it3_end; ++it3 ) {
		retval.push_back( it3->first );
	}
	return retval;
}

/// @brief queries the database for all orientations with nstrands strands
utility::vector1< std::string >
SheetDB::orientations( core::Size const nstrands ) const
{
	utility::vector1< std::string > retval;
	MapByStrands::const_iterator it1 = sheet_db_.find( nstrands );
	if ( it1 == sheet_db_.end() ) {
		return retval;
	}
	assert( it1 != sheet_db_.end() );
	MapByOrientations::const_iterator it2_end = it1->second.end();
	for ( MapByOrientations::const_iterator it2 = it1->second.begin(); it2 != it2_end; ++it2 ) {
		retval.push_back( it2->first );
	}
	return retval;
}

/// @brief queries the database for all strand lengths
utility::vector1< core::Size >
SheetDB::nstrand_values() const
{
	utility::vector1< core::Size > retval;
	for ( MapByStrands::const_iterator it = sheet_db_.begin(); it != sheet_db_.end(); ++it ) {
		retval.push_back( it->first );
	}
	return retval;
}

/// @brief adds sheet(s) based on a strand pairing set and pose
int
SheetDB::add_sheets_from_pose( core::pose::Pose & pose )
{
	using namespace protocols::fldsgn::topology;
	if ( pose.empty() ) {
		return 0;
	}
	utility::vector1< core::Size > positions;
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		positions.push_back(i);
	}
	core::chemical::ResidueTypeSetCOP restype_set =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	protocols::toolbox::pose_manipulation::construct_poly_XXX_pose( "ALA", pose, positions, restype_set, false, false, false );
	protocols::moves::DsspMover dssp;
	dssp.apply( pose );
	SS_Info2_OP ssinfo = SS_Info2_OP( new SS_Info2( pose ) );
	StrandPairingSet spairs = calc_strand_pairing_set( pose, ssinfo );
	TR << spairs << std::endl;
	return add_sheets_from_spairs( pose, spairs, *ssinfo );
}

core::graph::Graph
graph_from_strandpairings(
	protocols::fldsgn::topology::StrandPairings const & spairs,
	std::map< core::Size, core::Size > & node_to_strand )
{
	std::map< core::Size, core::Size > strand_to_node;
	utility::vector1< std::pair< core::Size, core::Size > > g_edges;
	for ( core::Size i=1; i<=spairs.size(); ++i ) {
		if ( spairs[i]->size1() < SheetDB::MIN_STRAND_LEN ) {
			TR << "Skipping strand " << spairs[i]->s1() << " becuase it is too short: " << spairs[i]->size1() << std::endl;
			continue;
		}
		if ( spairs[i]->size2() < SheetDB::MIN_STRAND_LEN ) {
			TR << "Skipping strand " << spairs[i]->s2() << " becuase it is too short: " << spairs[i]->size2() << std::endl;
			continue;
		}
		if ( strand_to_node.find(spairs[i]->s1()) == strand_to_node.end() ) {
			core::Size const nodeidx = strand_to_node.size()+1;
			strand_to_node[spairs[i]->s1()] = nodeidx;
			assert( node_to_strand.find(nodeidx) == node_to_strand.end() );
			node_to_strand[nodeidx] = spairs[i]->s1();
		}
		if ( strand_to_node.find(spairs[i]->s2()) == strand_to_node.end() ) {
			core::Size const nodeidx = strand_to_node.size()+1;
			strand_to_node[spairs[i]->s2()] = nodeidx;
			assert( node_to_strand.find(nodeidx) == node_to_strand.end() );
			node_to_strand[nodeidx] = spairs[i]->s2();
		}
		assert( strand_to_node.find(spairs[i]->s1()) != strand_to_node.end() );
		assert( strand_to_node.find(spairs[i]->s2()) != strand_to_node.end() );
		if ( spairs[i]->rgstr_shift() != 99 ) {
			g_edges.push_back( std::pair< core::Size, core::Size >( strand_to_node[spairs[i]->s1()], strand_to_node[spairs[i]->s2()] ) );
		}
	}

	TR.Debug << node_to_strand << strand_to_node << " number of nodes=" << strand_to_node.size() << std::endl;
	core::graph::Graph g( strand_to_node.size() );
	for ( core::Size i=1; i<=g_edges.size(); ++i ) {
		TR.Debug << "Adding edge " << g_edges[i].first << " <--> " << g_edges[i].second << std::endl;
		g.add_edge( g_edges[i].first, g_edges[i].second );
	}

	return g;
}

int
add_to_pose(
	core::pose::PoseOP newpose,
	core::pose::Pose const & pose,
	core::Size s_start,
	core::Size s_stop )
{
	assert( newpose );
	if ( s_start > s_stop ) {
		core::Size const tmp = s_stop;
		s_stop = s_start;
		s_start = tmp;
	}
	assert( s_start <= s_stop );
	core::Size const len = s_stop - s_start + 1;
	assert( len > 0 );
	core::Size start = s_start;
	if ( ( start > 1 ) && // not off end of pose
			( pose.chain(start) == pose.chain(start-1) ) // same chain
			) {
		start = start - 1;
	}
	core::Size stop = s_stop;
	if ( ( stop < pose.size() ) && // not off end of pose
			( pose.chain(stop) == pose.chain(stop+1) ) // same chain
			) {
		stop = stop + 1;
	}

	core::pose::Pose tmppose( pose, start, stop );
	if ( tmppose.fold_tree().num_jump() ) {
		TR.Error << "pose segment being copied has a jump... skipping it." << std::endl;
		return -1;
	}
	assert( tmppose.fold_tree().num_jump() == 0 );
	core::kinematics::FoldTree ft = tmppose.fold_tree();
	ft.reorder( tmppose.size() );
	tmppose.fold_tree( ft );
	core::Size local_res = 1;
	// add terminal residues if necessary
	if ( s_start == start ) {
		// prepend residue here
		core::chemical::ResidueType const restype =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA");
		core::conformation::ResidueOP new_rsd = core::conformation::ResidueFactory::create_residue( restype );
		assert( new_rsd );
		if ( tmppose.residue(1).is_lower_terminus() ) {
			core::pose::remove_lower_terminus_type_from_pose_residue( tmppose, 1 );
		}
		tmppose.prepend_polymer_residue_before_seqpos( *new_rsd, 1, true );
		core::pose::add_lower_terminus_type_to_pose_residue( tmppose, 1 );
		start = start - 1;
	}
	for ( core::Size r=start; r<s_start; ++r ) {
		tmppose.set_phi( local_res, 180.0 );
		tmppose.set_psi( local_res, 180.0 );
		tmppose.set_omega( local_res, 180.0 );
		++local_res;
	}
	ft = tmppose.fold_tree();
	ft.reorder( 1 );
	tmppose.fold_tree( ft );
	if ( s_stop == stop ) {
		// append residue here
		core::chemical::ResidueType const restype =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA");
		core::conformation::ResidueOP new_rsd = core::conformation::ResidueFactory::create_residue( restype );
		assert( new_rsd );
		if ( tmppose.residue(tmppose.size()).is_upper_terminus() ) {
			core::pose::remove_upper_terminus_type_from_pose_residue( tmppose, tmppose.size() );
		}
		tmppose.append_polymer_residue_after_seqpos( *new_rsd, tmppose.size(), true );
		core::pose::add_upper_terminus_type_to_pose_residue( tmppose, tmppose.size() );
		stop = stop + 1;
	}
	local_res = tmppose.size();
	for ( core::Size r=stop; r>s_stop; --r ) {
		tmppose.set_phi( local_res, 180.0 );
		tmppose.set_psi( local_res, 180.0 );
		--local_res;
	}

	TR.Debug << "Copied residues " << start << " to " << stop << " s_start=" << s_start << " s_stop=" << s_stop << " len=" << len << std::endl;
	if ( newpose->size() > 0 ) {
		newpose->conformation().insert_conformation_by_jump(
			tmppose.conformation(),
			newpose->size()+1,
			newpose->num_jump(),
			1 );
	} else {
		*newpose = tmppose;
	}
	return 0;
}

utility::vector1< core::pose::PoseOP >
extract_sheets_from_strandlist(
	utility::vector1< core::pose::PoseOP > const & strands,
	core::Size const nstrands )
{
	utility::vector1< core::pose::PoseOP > poses;
	for ( core::Size stop_strand = nstrands; stop_strand <= strands.size(); ++stop_strand ) {
		core::pose::PoseOP pose = core::pose::PoseOP( new core::pose::Pose() );
		for ( core::Size i=nstrands; i>=1; --i ) {
			if ( pose->size() ) {
				pose->conformation().insert_conformation_by_jump( strands[stop_strand-i+1]->conformation(),
					pose->size()+1,
					pose->num_jump()+2,
					1,
					pose->num_jump()+1 );
			} else {
				pose = strands[stop_strand-i+1]->clone();
			}
		}
		TR.Debug << "Added pose with ss=" << pose->secstruct() << std::endl;
		poses.push_back( pose );
	}
	assert( poses.size() + nstrands - 1 == strands.size() );
	return poses;
}

core::pose::PoseOP
reorder_chains( core::pose::Pose & pose )
{
	using namespace protocols::fldsgn::topology;
	protocols::moves::DsspMover dssp;
	dssp.apply( pose );
	SS_Info2_OP ss_info = SS_Info2_OP( new SS_Info2( pose, pose.secstruct() ) );
	StrandPairingSet spairset = calc_strand_pairing_set( pose, ss_info );
	TR << "After ssinfo " << spairset << std::endl;
	core::Size const nstrands = pose.conformation().num_chains();

	// don't even bother if some chains aren't detected as paired strands
	if ( nstrands != spairset.num_strands() ) {
		TR.Debug << "Fewer paired strands than chains... skipping reordering." << std::endl;
		return core::pose::PoseOP();
	}

	core::Size start_strand = 0;
	// find strand with only one other paired strand
	for ( core::Size i=1; i<=nstrands; ++i ) {
		core::Size const n_paired = spairset.neighbor_strands( i ).size();
		if ( n_paired == 1 ) {
			start_strand = i;
			break;
		} else if ( !start_strand && ( n_paired == 2 ) ) {
			// in the case of barrels, all strands have two paired partners. We choose the first one here
			start_strand = i;
		}
	}
	assert( start_strand );
	std::set< core::Size > visited;
	std::stack< core::Size > nodes;
	utility::vector1< core::Size > neworder;
	nodes.push( start_strand );
	while ( nodes.size() ) {
		core::Size const strand = nodes.top();
		nodes.pop();
		if ( visited.find( strand ) != visited.end() ) {
			continue;
		}
		visited.insert( strand );
		neworder.push_back( strand );
		protocols::fldsgn::topology::StrandPairingSet::VecSize neighbors = spairset.neighbor_strands( strand );
		for ( core::Size i=1; i<=neighbors.size(); ++i ) {
			nodes.push( neighbors[i] );
		}
	}
	TR << "New strand order = " << neworder << std::endl;
	if ( visited.size() != nstrands ) {
		return core::pose::PoseOP();
	}
	assert( visited.size() == nstrands );
	assert( neworder.size() == nstrands );

	utility::vector1< core::pose::PoseOP > strands = pose.split_by_chain();
	core::pose::PoseOP newpose;
	for ( core::Size i=1; i<=neworder.size(); ++i ) {
		if ( !newpose ) {
			newpose = strands[neworder[i]]->clone();
		} else {
			newpose->conformation().insert_conformation_by_jump(
				strands[neworder[i]]->conformation(),
				newpose->size()+1,
				newpose->fold_tree().num_jump(),
				1 );
		}
	}
	assert( newpose->size() == pose.size() );
	assert( newpose->conformation().num_chains() == pose.conformation().num_chains() );
	return newpose;
}

utility::vector1< core::pose::PoseOP >
extract_sheets_from_pose(
	core::pose::Pose const & pose,
	protocols::fldsgn::topology::StrandPairings const & spairs,
	protocols::fldsgn::topology::SS_Info2 const & ss_info,
	bool const idealize )
{
	using namespace protocols::fldsgn::topology;

	// mover for idealizing bonds
	protocols::idealize::IdealizeMover idealizemover;

	utility::vector1< core::pose::PoseOP > sheets;
	std::map< core::Size, core::Size > node_to_strand;
	core::graph::Graph g = graph_from_strandpairings( spairs, node_to_strand );
	g.print_vertices();
	utility::vector1< std::pair< core::Size, core::Size > > connected_sheets = core::graph::find_connected_components( g );
	for ( core::Size i=1; i<=connected_sheets.size(); ++i ) {
		TR << "Tomponent Group " << connected_sheets[i].first << " " << connected_sheets[i].second << std::endl;
		bool success = true;
		std::stack< core::Size > nodes;
		std::set< core::Size > visited;
		nodes.push( connected_sheets[i].first );
		core::pose::PoseOP newpose = core::pose::PoseOP( new core::pose::Pose() );
		while ( nodes.size() ) {
			core::Size const nodenum = nodes.top();
			nodes.pop();
			if ( visited.find( nodenum ) != visited.end() ) {
				continue;
			}
			visited.insert( nodenum );
			core::Size const strand = node_to_strand[nodenum];
			core::Size start = ss_info.strand(strand)->begin();
			core::Size stop = ss_info.strand(strand)->end();
			if ( add_to_pose( newpose, pose, start, stop ) ) {
				success = false;
				break;
			}
			core::graph::Node const * n = g.get_node(nodenum);
			assert( n );
			for ( core::graph::Graph::EdgeListConstIter e=n->const_edge_list_begin(); e != n->const_edge_list_end(); ++e ) {
				core::Size const othernode = (*e)->get_other_ind(nodenum);
				nodes.push(othernode);
			}
		}
		if ( success ) {
			TR << "Reordering chains" << std::endl;
			newpose = reorder_chains( *newpose );
			if ( newpose ) {
				if ( idealize ) {
					idealizemover.apply( *newpose );
				}
				sheets.push_back( newpose );
			}
		}
	}
	return sheets;
}

/// @brief creates a list of whether a pose residue is beta-sheet-type paired to anything
utility::vector1< bool >
calc_paired_residues( core::pose::Pose & pose, protocols::fldsgn::topology::StrandPairingSet const & spairs )
{
	TR.Debug << "new set= " << spairs << std::endl;
	utility::vector1< bool > paired( pose.size(), false );
	for ( core::Size c=1; c<=pose.conformation().num_chains(); ++c ) {
		core::Size s1 = c;
		core::Size s2 = c+1;
		if ( s2 > pose.conformation().num_chains() ) {
			s2 = s1 - 1;
		}
		// if shift comes back as 99, there probably aren't valid pairings
		// these sheets are currently thrown out until I can learn more about shift=99
		if ( spairs.strand_pairing( s1, s2 )->rgstr_shift() == 99 ) {
			return utility::vector1< bool >( pose.size(), false );
		}
		// move from beginning to end looking for paired residues
		for ( core::Size i=pose.conformation().chain_begin(c); i<=pose.conformation().chain_end(c); ++i ) {
			TR.Debug << "Trying res " << i << "...";
			if ( spairs.strand_pairing( s1, s2 )->has_paired_residue( i ) ) {
				core::Size const pair_res = spairs.strand_pairing( s1, s2 )->residue_pair( i );
				paired[pair_res] = true;
				paired[i] = true;
				TR.Debug << " paired to " << pair_res << std::endl;
			} else {
				TR.Debug << " not paired" << std::endl;
			}
		}
	}
	return paired;
}

std::string
make_strand_info_str( Lengths const & lengths, RegisterShifts const & shifts )
{
	std::string retval;
	debug_assert( lengths.size() == shifts.size() );
	for ( core::Size i=1; i<=lengths.size(); ++i ) {
		if ( retval.size() ) {
			retval += ',';
		}
		retval += boost::lexical_cast< std::string >( lengths[i] );
		retval += ':';
		retval += boost::lexical_cast< std::string >( shifts[i] );
	}
	return retval;
}

std::string
make_lengths_str(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & paired,
	Lengths const & lengths,
	StrandOrientations const & orients,
	protocols::fldsgn::topology::StrandPairingSet const & spairs )
{
	RegisterShifts shifts;
	core::Size prev_start = 0;
	core::Size prev_stop = 0;

	assert( pose.conformation().num_chains() == orients.size() );
	for ( core::Size c=1; c<=pose.conformation().num_chains(); ++c ) {
		core::Size start = 0;
		core::Size stop = 0;
		for ( core::Size i=pose.conformation().chain_begin(c); i<=pose.conformation().chain_end(c); ++i ) {
			if ( paired[i] ) {
				start = i;
				break;
			}
		}
		for ( core::Size i=pose.conformation().chain_end(c); i>=pose.conformation().chain_begin(c); --i ) {
			if ( paired[i] ) {
				stop = i;
				break;
			}
		}
		if ( !start || !stop ) {
			TR << " no strand start/stop found. Chain start=" << pose.conformation().chain_begin(c) << " end=" << pose.conformation().chain_end(c) << " paired=" << paired << std::endl;
			return "";
		}
		TR << "Chain start=" << pose.conformation().chain_begin(c) << " end=" << pose.conformation().chain_end(c) << " start=" << start << " stop=" << stop << " paired=" << paired << std::endl;
		assert( start );
		assert( stop );
		int shift = 0;
		if ( c == 1 ) {
			shift = 0;
		} else {
			assert( prev_start );
			assert( prev_stop );
			// find a residue paired to the previous strand
			core::Size prev_res = prev_start;
			int offset = 0;
			if ( ( orients[c-1] == UP ) && ( orients[c] == UP ) ) {
				// find residue with pairing to this strand
				for ( prev_res=prev_start; prev_res<=prev_stop; ++prev_res ) {
					if ( spairs.strand_pairing( c-1, c )->has_paired_residue(prev_res) ) {
						TR << prev_res << " " << offset << " +paired" << std::endl;
						break;
					}
					++offset;
				}
				assert( spairs.strand_pairing( c-1, c )->has_paired_residue(prev_res) );
				if ( offset > 0 ) {
					shift = offset;
				} else {
					assert( prev_res == prev_start );
					int paired_with_prev = spairs.strand_pairing( c-1, c )->residue_pair(prev_res);
					shift = static_cast<int>(start) - paired_with_prev;
				}
			} else if ( ( orients[c-1] == DOWN ) && ( orients[c] == UP ) ) {
				for ( prev_res=prev_stop; prev_res>=prev_start; --prev_res ) {
					if ( spairs.strand_pairing( c-1, c )->has_paired_residue(prev_res) ) {
						TR << prev_res << " " << offset << " -paired" << std::endl;
						break;
					}
					++offset;
				}
				assert( spairs.strand_pairing( c-1, c )->has_paired_residue(prev_res) );
				if ( offset > 0 ) {
					shift = offset;
				} else {
					assert( prev_res == prev_stop );
					int paired_with_prev = spairs.strand_pairing( c-1, c )->residue_pair(prev_res);
					shift = static_cast<int>(start) - paired_with_prev;
				}
			} else if ( ( orients[c-1] == UP ) && ( orients[c] == DOWN ) ) {
				// find residue with pairing to this strand
				for ( prev_res=prev_start; prev_res<=prev_stop; ++prev_res ) {
					if ( spairs.strand_pairing( c-1, c )->has_paired_residue(prev_res) ) {
						TR << prev_res << " " << offset << " +paired" << std::endl;
						break;
					}
					++offset;
				}
				assert( spairs.strand_pairing( c-1, c )->has_paired_residue(prev_res) );
				if ( offset > 0 ) {
					shift = offset;
				} else {
					assert( prev_res == prev_start );
					int paired_with_prev = spairs.strand_pairing( c-1, c )->residue_pair(prev_res);
					shift = paired_with_prev - static_cast<int>(stop);
				}
			} else if ( ( orients[c-1] == DOWN ) && ( orients[c] == DOWN ) ) {
				for ( prev_res=prev_stop; prev_res>=prev_start; --prev_res ) {
					if ( spairs.strand_pairing( c-1, c )->has_paired_residue(prev_res) ) {
						TR << prev_res << " " << offset << " -paired" << std::endl;
						break;
					}
					++offset;
				}
				assert( spairs.strand_pairing( c-1, c )->has_paired_residue(prev_res) );
				if ( offset > 0 ) {
					shift = offset;
				} else {
					assert( prev_res == prev_stop );
					int paired_with_prev = spairs.strand_pairing( c-1, c )->residue_pair(prev_res);
					shift = static_cast<int>(stop) - paired_with_prev;
				}
			}
		}
		shifts.push_back( shift );
		TR << "Strand " << c << " prevstart=" << prev_start << " prevstop=" << prev_stop << " start=" << start << " stop=" << stop << " len=" << lengths[c] << " shift=" << shift << std::endl;
		prev_start = start;
		prev_stop = stop;
	}
	return make_strand_info_str( lengths, shifts );
}

/// @brief using the pose secondary structure, remove residues that aren't paired. returns a list of paired strand lengths, which do not necessarily correspond to pose lengths due to beta-bulges and other weird things.
std::pair< std::string, std::string >
prune_unpaired_residues( core::pose::Pose & pose )
{
	using namespace protocols::fldsgn::topology;
	protocols::moves::DsspMover dssp;
	dssp.apply( pose );
	SS_Info2_OP ss_info = SS_Info2_OP( new SS_Info2( pose, pose.secstruct() ) );
	StrandPairingSet spairs = calc_strand_pairing_set( pose, ss_info );

	utility::vector1< bool > paired = calc_paired_residues( pose, spairs );
	Lengths lengths( pose.conformation().num_chains(), 0 );
	StrandOrientations orients( pose.conformation().num_chains(), UP );
	core::Size prevchain_start = 0;
	// determine which residues to copy from each chain
	core::pose::PoseOP newpose = core::pose::PoseOP( new core::pose::Pose() );
	for ( core::Size c=1; c<=pose.conformation().num_chains(); ++c ) {
		if ( prevchain_start >= pose.conformation().chain_begin(c) ) {
			TR.Error << " chain " << c << " starts at res " << pose.conformation().chain_begin(c) << " but previous chain began at res " << prevchain_start << std::endl;
			assert( pose.conformation().chain_begin(c) > prevchain_start );
			runtime_assert( pose.conformation().chain_begin(c) > prevchain_start );
		}
		prevchain_start = pose.conformation().chain_begin(c);
		core::Size start = 0;
		core::Size stop = 0;
		core::Size paired_res_count = 0;
		// go forward looking for paired
		for ( core::Size i=pose.conformation().chain_begin(c); i<=pose.conformation().chain_end(c); ++i ) {
			if ( paired[i] ) {
				if ( !start ) {
					start = i;
				}
				++paired_res_count;
			}
		}

		// go backward looking for paired
		for ( core::Size i=pose.conformation().chain_end(c); i>=pose.conformation().chain_begin(c); --i ) {
			if ( paired[i] ) {
				stop = i;
				break;
			}
		}

		if ( stop == 0 || start == 0 ) {
			pose.clear();
			return std::make_pair( "", "" );
		}

		// check to make sure paired residues are contiguous
		assert( start <= stop );
		bool broken = false;
		for ( core::Size i=start; i<=stop; ++i ) {
			if ( !paired[i] ) {
				broken = true;
				TR.Debug << " strand is broken since residue " << i << " is not paired but is inside a strand." << std::endl;
				break;
			}
		}
		if ( broken ) {
			pose.clear();
			return std::make_pair( "", "" );
		}

		TR.Debug << "for strand " << c << ", start=" << start << " and stop=" << stop << std::endl;
		if ( add_to_pose( newpose, pose, start, stop ) ) {
			return std::make_pair( "", "" );
		}

		// if too short, skip
		if ( paired_res_count < SheetDB::MIN_STRAND_LEN ) {
			return std::make_pair( "", "" );
		}
		lengths[c] = paired_res_count;

		if ( c > 1 ) {
			char orient = spairs.strand_pairing( c-1, c )->orient();
			if ( orient == 'P' ) {
				orients[c] = orients[c-1];
			} else if ( orient == 'A' ) {
				if ( orients[c-1] == UP ) orients[c] = DOWN;
				else if ( orients[c-1] == DOWN ) orients[c] = UP;
				else {
					std::stringstream msg;
					msg << "Invalid orientation for strand " << c-1 << ": " << orients[ c-1 ] << std::endl;
					utility_exit_with_message( msg.str() );
				}
			} else {
				TR.Error << "Invalid orientation between strands " << c-1 << " and " << c << " : " << orient << std::endl;
			}
		}
	}
	std::string const lengths_str = make_lengths_str( pose, paired, lengths, orients, spairs );
	std::string const orient_str = orientations_str( orients );
	pose = *newpose;
	TR.Debug << "orientations=" << orient_str << " lengths=" << lengths_str << std::endl;
	return std::make_pair( orient_str, lengths_str );
}

/// @brief adds sheet(s) based on a strand pairing set and pose
int
SheetDB::add_sheets_from_spairs(
	core::pose::Pose const & pose,
	protocols::fldsgn::topology::StrandPairingSet & spairs, // should be const, but the api says end() is nonconst
	protocols::fldsgn::topology::SS_Info2 const & ss_info )
{
	using namespace protocols::fldsgn::topology;
	utility::vector1< core::pose::PoseOP > sheets = extract_sheets_from_pose( pose, spairs.strand_pairings(), ss_info, idealize_ );
	int added_count = 0;
	for ( core::Size i=1; i<=sheets.size(); ++i ) {
		TR << "Started sheet " << i << std::endl;
		utility::vector1< core::pose::PoseOP > strandlist = sheets[i]->split_by_chain();
		for ( core::Size nstrands=2; nstrands <= strandlist.size(); ++nstrands ) {
			utility::vector1< core::pose::PoseOP > sheets_n_strands =
				extract_sheets_from_strandlist( strandlist, nstrands );
			// for each sheet generated
			for ( core::Size s=1; s<=sheets_n_strands.size(); ++s ) {
				// redo dssp and eliminate unpaired residues
				TR.Debug << "Pruning ss= " << sheets_n_strands[s]->secstruct() << " nstrands = " << num_strands( *(sheets_n_strands[s]) ) << std::endl;
				std::pair< std::string, std::string > or_len = prune_unpaired_residues( *(sheets_n_strands[s]) );
				if ( ( or_len.first == "" ) || ( sheets_n_strands[s]->empty() ) ) {
					continue;
				}
				assert( nstrands == num_strands( *(sheets_n_strands[s]) ) );
				add_sheet( sheets_n_strands[s], nstrands, or_len.first, or_len.second );
				++added_count;
			}
		}
		TR << "Done with sheet " << i << std::endl;
	}
	return added_count;
}

/// @brief counts number of strands based on number of jumps. Pose must be a disconnected sheet
core::Size num_strands( core::pose::Pose const & pose )
{
	return pose.fold_tree().num_jump() + 1;
}

/// @brief gets a pair of strings, which refer to orientations, and lengths/shifts. Pose secstruct MUST be set
std::pair< std::string, std::string >
find_orientations_and_lengths( core::pose::Pose const & pose )
{
	// gets a vector of three strings, which refer to orientations, lengths, and shifts
	std::string const ss = pose.secstruct();
	assert( ss.size() && ( ss[0] != ' ' ) );
	protocols::fldsgn::topology::SS_Info2_OP sheet_ss_info =
		protocols::fldsgn::topology::SS_Info2_OP( new protocols::fldsgn::topology::SS_Info2( pose, ss ) );
	protocols::fldsgn::topology::StrandPairingSet pairs = protocols::fldsgn::topology::calc_strand_pairing_set( pose, sheet_ss_info );
	if ( pairs.size() == 0 ) {
		return std::make_pair( "", "" );
	}
	core::Size const strand_count = num_strands( pose );
	TR << "Nstrands = " << strand_count << std::endl;
	TR << "PairSet  = " << pairs << std::endl;
	assert( strand_count > 1 );
	std::map< core::Size, std::string > orientations_str;
	utility::vector1< bool > orientations( strand_count, true );
	utility::vector1< core::Size > lengths( strand_count, 0 );
	int const NOT_SET = -99999;
	utility::vector1< int > shifts( strand_count, NOT_SET );
	std::set< core::Size > parsed;
	parsed.insert( 1 ); // strand 1 is always "true" (aka up)
	shifts[1] = 0;
	lengths[1] = sheet_ss_info->strand(1)->length();
	for ( protocols::fldsgn::topology::StrandPairings::const_iterator p = pairs.begin(); p != pairs.end(); ++p ) {
		core::Size const s1 = (*p)->s1();
		core::Size const s2 = (*p)->s2();
		char const orient = (*p)->orient();
		core::Size const len1 = sheet_ss_info->strand(s1)->length();
		core::Size const len2 = sheet_ss_info->strand(s2)->length();
		int const shift = (*p)->rgstr_shift();
		TR << "s1=" << s1 << " s2=" << s2 << " len1=" << len1 << " len2=" << len2 << " orient=" << orient << " shift=" << shift << std::endl;
		assert( (*p)->begin1() < (*p)->end1() );
		assert( s1 > 0 );
		assert( s2 > 0 );
		assert( s1 <= strand_count );
		assert( s2 <= strand_count );
		assert( parsed.find(s1) != parsed.end() );

		// set lengths
		assert( lengths[s1] == len1 );
		assert( !lengths[s2] );
		TR << "Setting lengths for " << s2 << " to " << len2 << std::endl;
		lengths[s2] = len2;

		// set register shift
		assert( shifts[s1] != NOT_SET );
		assert( shifts[s2] == NOT_SET );
		shifts[s2] = shift;

		// set orientations
		if ( orient == 'A' ) {
			orientations[s2] = ! orientations[s1];
		} else if ( orient == 'P' ) {
			orientations[s2] = orientations[s1];
		} else {
			TR.Error << "Invalid orientation thrown: " << orient << std::endl;
			utility_exit();
		}
		parsed.insert(s2);
	}

	std::string orient_str;
	TR << "orientations:" << orientations << std::endl;
	TR << "lengths:" << lengths << std::endl;
	TR << "shifts:" << shifts << std::endl;
	for ( core::Size i=1; i<=orientations.size(); ++i ) {
		if ( orientations[i] ) {
			orient_str += '1';
		} else {
			orient_str += '0';
		}
	}
	assert( shifts.size() == lengths.size() );

	std::string len_str;
	for ( core::Size i=1; i<=lengths.size(); ++i ) {
		if ( len_str.size() ) {
			len_str += ',';
		}
		len_str += boost::lexical_cast< std::string >( lengths[i] ) + '_' + boost::lexical_cast< std::string >( shifts[i] );
	}
	return std::make_pair( orient_str, len_str );
}

std::string reverse_orientations( std::string const & orient )
{
	std::string retval = "";
	bool flip = ( orient[orient.size()-1] == '0' );
	for ( core::Size i=orient.size(); i>=1; --i ) {
		if ( flip ) {
			if ( orient[i-1] == '0' ) {
				retval += '1';
			} else if ( orient[i-1] == '1' ) {
				retval += '0';
			} else {
				TR << "Something other than 0 or 1 found in orientations string :" << orient << std::endl;
				assert( false );
			}
		} else {
			retval += orient[i-1];
		}
	}
	assert( retval.size() == orient.size() );
	return retval;
}

utility::vector1< bool >
parse_orientations( std::string const & orientations )
{
	utility::vector1< bool > retval;
	for ( core::Size i=0; i<orientations.size(); ++i ) {
		if ( orientations[i] == '0' ) {
			retval.push_back( false );
		} else if ( orientations[i] == '1' ) {
			retval.push_back( true );
		} else {
			TR << "Something other than 0 or 1 found in orientations string :" << orientations << std::endl;
			assert( false );
		}
	}
	assert( orientations.size() == retval.size() );
	return retval;
}

std::pair< utility::vector1< core::Size >, utility::vector1< int > >
parse_lengths( std::string const & lengths )
{
	std::pair< utility::vector1< core::Size >, utility::vector1< int > > retval;
	utility::vector1< std::string > strands = utility::string_split( lengths, ',' );
	for ( core::Size i=1; i<=strands.size(); ++i ) {
		utility::vector1< std::string > info = utility::string_split( strands[i], ':' );
		assert( info.size() == 2 );
		retval.first.push_back( boost::lexical_cast< core::Size >( info[1] ) );
		retval.second.push_back( boost::lexical_cast< int >( info[2] ) );
	}
	assert( retval.first.size() == retval.second.size() );
	assert( retval.first.size() == strands.size() );
	return retval;
}

std::string reverse_lengths( std::string const & orientations, std::string const & lengths )
{
	utility::vector1< bool > orients = parse_orientations( orientations );
	std::pair< utility::vector1< core::Size >, utility::vector1< int > > lens = parse_lengths( lengths );
	assert( orients.size() == lens.first.size() );
	assert( orients.size() == lens.second.size() );
	utility::vector1< core::Size > rlengths;
	utility::vector1< int > rshifts;
	bool flip = !orients[orients.size()];
	std::string retval;
	for ( core::Size i=lens.first.size(); i>=1; --i ) {
		rlengths.push_back( lens.first[i] );
		if ( i == lens.first.size() ) {
			rshifts.push_back( 0 );
			continue;
		}
		if ( !flip ) {
			rshifts.push_back( -lens.second[i+1] );
		} else {
			// flip the sheet
			int const top = lens.first[i];
			int const top_prev = lens.first[i+1] + lens.second[i+1];
			rshifts.push_back( top_prev - top );
		}
	}
	return protocols::denovo_design::components::make_strand_info_str( rlengths, rshifts );
}

std::pair< std::string, std::string >
canonicalize( std::string const & orientations, std::string const & lengths )
{
	std::string const rev_orientations = reverse_orientations( orientations );
	std::string const rev_lengths = reverse_lengths( orientations, lengths );

	std::pair< std::string, std::string > retval;
	retval.second = choose_canonical( lengths, rev_lengths );
	if ( retval.second == lengths ) {
		retval.first = orientations;
	} else if ( retval.second == rev_lengths ) {
		retval.first = rev_orientations;
	} else {
		assert( false );
	}
	return retval;
}

std::string const &
choose_canonical( std::string const & l1, std::string const & l2 )
{
	if ( l1 <= l2 ) {
		return l1;
	} else {
		return l2;
	}
}

core::pose::PoseOP
reverse_chains( core::pose::Pose const & pose )
{
	core::pose::PoseOP newpose;
	utility::vector1< core::pose::PoseOP > chains = pose.split_by_chain();
	for ( core::Size i=chains.size(); i>=1; --i ) {
		if ( !newpose ) {
			newpose = chains[i]->clone();
		} else {
			newpose->conformation().insert_conformation_by_jump(
				chains[i]->conformation(),
				newpose->size()+1,
				newpose->num_jump(),
				1 );
		}
	}
	assert( newpose->size() == pose.size() );
	assert( newpose->conformation().num_chains() == pose.conformation().num_chains() );
	return newpose;
}


} // components
} // denovo_design
} // protocols
