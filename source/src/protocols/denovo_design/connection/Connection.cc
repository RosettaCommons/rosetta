// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/Connection.cc
/// @brief Connection objects for piece-wise assembly of structures
/// @detailed
/// @author Tom Linsky


//Unit Headers
#include <protocols/denovo_design/connection/Connection.hh>
#include <protocols/denovo_design/connection/ConnectionCreators.hh>

//Project Headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/denovo_design/util.hh>

//Protocol Headers
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>

//Core Headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/sequence/ABEGOManager.hh>

//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers

// Boost headers
#include <boost/algorithm/string/predicate.hpp>
#include <boost/assign.hpp>

//C++ Headers
static thread_local basic::Tracer TR( "protocols.denovo_design.Connection" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace connection {
////////////////////////////////////////////////////////////////////////////////////////////////////

/// Motif code:

std::ostream & operator<<( std::ostream & os, Connection::Motif const & motif ) {
	os << motif.len << " " << motif.ss << " " << motif.abego;
	return os;
}

///  ---------------------------------------------------------------------------------
///  Connection main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
Connection::Connection() :
	components::NamedMover(),
	chain1_( 0 ),
	chain2_( 0 ),
	trials_( 1 ),
	overlap_( 0 ),
	do_remodel_( true ),
	allow_cyclic_( false ),
	connecting_bond_dist_( 1.5 ),
	idealized_abego_( false ),
	performs_orientation_( false ),
	check_abego_( true )
{
	motifs_.clear();
	motifs_.push_back( Motif( 0, 'L', "X" ) );
	cut_resis_.clear();
	comp1_ids_.clear();
	comp2_ids_.clear();
	disallowed_pairs_.clear();
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
Connection::~Connection() {}

/// @brief setup the parameters via an xml tag
void
Connection::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	components::NamedMover::parse_my_tag( tag, data, filters, movers, pose );
	comp1_ids_.clear();
	comp2_ids_.clear();
	if ( tag->hasOption( "segment1" ) ) {
		set_comp1_ids( tag->getOption< std::string >( "segment1" ) );
	}
	if ( tag->hasOption( "segment2" ) ) {
		set_comp2_ids( tag->getOption< std::string >( "segment2" ) );
	}

	if ( tag->hasOption( "chain1" ) && tag->hasOption( "segment1" ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "chain and comp can't both be specified in " + id() );
	}
	if ( tag->hasOption( "chain2" ) && tag->hasOption( "segment2" ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "chain and comp can't both be specified in " + id() );
	}
	if ( tag->hasOption( "chain1" ) ) {
		chain1_ = tag->getOption< core::Size >( "chain1" );
	}
	if ( tag->hasOption( "chain2" ) ) {
		chain2_ = tag->getOption< core::Size >( "chain2" );
	}
	set_check_abego( tag->getOption< core::Size >( "check_abego", check_abego_ ) );
	set_allow_cyclic( tag->getOption< bool >( "allow_cyclic", allow_cyclic_ ) );
	set_trials( tag->getOption< core::Size >( "trials", trials_ ) );
	set_overlap( tag->getOption< core::Size >( "overlap", overlap_ ) );
	set_idealized_abego( tag->getOption< bool >( "idealized_abego", idealized_abego_ ) );

	if ( tag->hasOption( "length" ) && tag->hasOption( "motif" ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Both lengths and motifs cannot be simultaneously specified in " + id());
	}

	if ( tag->hasOption( "idealized_abego" ) && tag->hasOption( "motif" ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Specific motifs cannot be specified when building idealized abegos " + id());
	}

	motifs_.clear();
	if ( tag->hasOption( "length" ) ) {
		std::string const length_str( tag->getOption< std::string >( "length" ) );
		set_lengths( length_str );
	}

	if ( tag->hasOption( "motif" ) ) {
		set_motifs( tag->getOption< std::string >( "motif" ) );
	}

	if ( !motifs_.size() ) {
		motifs_.push_back( Motif( 0, 'L', "X" ) );
	}

	cut_resis_.clear();
	std::string const cut_resi_str( tag->getOption< std::string >( "cut_resi", "" ) );
	if ( cut_resi_str != "" ) {
		set_cut_resis( cut_resi_str );
	}
}

/// @brief sets overlap
void
Connection::set_overlap( core::Size const overlap_val )
{
	overlap_ = overlap_val;
}

/// @brief sets comp1 ids from a string
void
Connection::set_comp1_ids( std::string const & id_str )
{
	if ( !id_str.empty() ) {
		set_comp1_ids( utility::string_split( id_str, ',' ) );
	} else {
		comp1_ids_.clear();
	}
}

/// @brief sets comp2 ids from a string
void
Connection::set_comp2_ids( std::string const & id_str )
{
	if ( !id_str.empty() ) {
		set_comp2_ids( utility::string_split( id_str, ',' ) );
	} else {
		comp2_ids_.clear();
	}
}

/// @brief sets the list of options for component 1
void
Connection::set_comp1_ids( utility::vector1< std::string > const & id_list )
{
	comp1_ids_.clear();
	if ( ! id_list.empty() ) {
		comp1_ids_ = id_list;
	}
}

/// @brief sets the list of options for component 1
void
Connection::set_comp2_ids( utility::vector1< std::string > const & id_list )
{
	comp2_ids_.clear();
	if ( ! id_list.empty() ) {
		comp2_ids_ = id_list;
	}
}

/// @brief sets the list of acceptable lengths by parsing a string
void
Connection::set_lengths( std::string const & length_str )
{
	set_lengths( parse_length_string( length_str ) );
}

void note_fwd_element_extensions(
	std::string & ss,
	std::string const & abego,
	char const sschar,
	char const element_abego )
{
	assert( ss.size() == abego.size() );
	// starting at beginning, locate extending elements
	core::Size i=0;
	core::Size const len = ss.size();
	while ( (i < len) && (abego[i] == element_abego) ) {
		ss[i] = sschar;
		++i;
	}
}

void note_rev_element_extensions(
	std::string & ss,
	std::string const & abego,
	char const sschar,
	char const element_abego )
{
	assert( ss.size() == abego.size() );
	// starting at end, locate extending elements
	int i=ss.size()-1;
	while ( (i >= 0) && (abego[i] == element_abego) ) {
		ss[i] = sschar;
		--i;
	}
}

/// @brief Given desired lengths, compute a set of idealized loop motifs via Nobu/Rie/YuRu rules
Connection::MotifList
Connection::calc_idealized_motifs(
	std::string const & abego1,
	std::string const & abego2,
	std::set< core::Size > const & lenset ) const
{
	static std::string const ab[] = {
		"GB" ,    "GBA" ,    "BAA" , // Basic loops
		"GBB",    "GBAB",    "BAAB", // Extend strand by one residue
		"AGB" ,   "AGBA" ,   "ABAA" , // Extend helix by one residue
		"AGBB",   "AGBAB",   "ABAAB",
		"AAGB" ,  "AAGBA" ,  "AABAA" , // Extend helix by two residues
		"AAGBB",  "AAGBAB",  "AABAAB",
		"AAAGB" , "AAAGBA" , "AAABAA" , // Extend helix by three residues
		"AAAGBB", "AAAGBAB", "AAABAAB" };
	static core::Size const numAB = 24;

	static std::string const ba[] = {
		"AB"   , "GBB"   , "BAB"   , "BGBB"   ,  // Basic loops plus extend strand by one residue
		"ABA"  , "GBBA"  , "BABA"  , "BGBBA"  ,  // Extend helix by one residue
		"ABAA" , "GBBAA" , "BABAA" , "BGBBAA" ,  // Extend helix by two residues
		"ABAAA", "BGGAAA", "BABAAA", "BGBBAAA" };
	static core::Size const numBA = 16;

	static std::string const bb[] = {
		"GG" ,  "EA" ,  "AA" ,  "AAG" ,   // Basic loops
		"BGG" , "BEA" , "BAA" , "BAAG" ,   // Extend first strand by one
		"GGB",  "EAB",  "AAB",  "AAGB",   // Extend second strand by one
		"BGGB", "BEAB", "BAAB", "BAAGB" }; // Extend both strands by one
	static core::Size const numBB = 16;

	static std::string const alpha_alpha[] = {
		"B", "AB", "AAB", "AAAB", "BA", "ABA", "AABA", "AAABA", "BAA", "ABAA",   // one-residue B and extensions
		"AABAA", "AAABAA", "BAAA", "ABAAA", "AABAAA", "AAABAAA",
		"G", "AG", "AAG", "AAAG", "GA", "AGA", "AAGA", "AAAGA", "GAA", "AGAA",   // one-residue G and extensions
		"AAGAA", "AAAGAA", "GAAA", "AGAAA", "AAGAAA", "AAAGAAA",
		"BB"   ,    "BG"   ,    "GB"   ,    "BAB"   ,    "BAG"   ,    "GAB"   ,  //loops not containing E or start/ending with A
		"BBA"  ,    "BGA"  ,    "GBA"  ,    "BABA"  ,    "BAGA"  ,    "GABA"  ,  // Extend h2 by one residue
		"BBAA" ,    "BGAA" ,    "GBAA" ,    "BABAA" ,    "BAGAA" ,    "GABAA" ,  // Extend h2 by two residues
		"BBAAA",    "BGAAA",    "GBAAA",    "BABAAA",    "BAGAAA",    "GABAAA",  // Extend h2 by three residues
		"ABB"   ,   "ABG"   ,   "AGB"   ,   "ABAB"   ,   "ABAG"   ,   "AGAB"   ,  // Extend h1 by one residue
		"ABBA"  ,   "ABGA"  ,   "AGBA"  ,   "ABABA"  ,   "ABAGA"  ,   "AGABA"  ,  //
		"ABBAA" ,   "ABGAA" ,   "AGBAA" ,   "ABABAA" ,   "ABAGAA" ,   "AGABAA" ,  //
		"ABBAAA",   "ABGAAA",   "AGBAAA",   "ABABAAA",   "ABAGAAA",   "AGABAAA",  //
		"AABB"   ,  "AABG"   ,  "AAGB"   ,  "AABAB"   ,  "AABAG"   ,  "AAGAB"   ,  // Extend h1 by two residue
		"AABBA"  ,  "AABGA"  ,  "AAGBA"  ,  "AABABA"  ,  "AABAGA"  ,  "AAGABA"  ,  //
		"AABBAA" ,  "AABGAA" ,  "AAGBAA" ,  "AABABAA" ,  "AABAGAA" ,  "AAGABAA" ,  //
		"AABBAAA",  "AABGAAA",  "AAGBAAA",  "AABABAAA",  "AABAGAAA",  "AAGABAAA",  //
		"AAABB"   , "AAABG"   , "AAAGB"   , "AAABAB"   , "AAABAG"   , "AAAGAB"   ,  // Extend h1 by three residue
		"AAABBA"  , "AAABGA"  , "AAAGBA"  , "AAABABA"  , "AAABAGA"  , "AAAGABA"  ,  //
		"AAABBAA" , "AAABGAA" , "AAAGBAA" , "AAABABAA" , "AAABAGAA" , "AAAGABAA" ,  //
		"AAABBAAA", "AAABGAAA", "AAAGBAAA", "AAABABAAA", "AAABAGAAA", "AAAGABAAA" };
	static core::Size const numAA = 96;

	MotifList retval;
	// determine previous residue by abego type
	if ( ( abego1 == "A" ) && ( abego2 == "A" ) ) {
		for ( core::Size i=0; i<numAA; ++i ) {
			if ( lenset.find( alpha_alpha[i].size() ) != lenset.end() ) {
				Motif m;
				for ( std::string::const_iterator c=alpha_alpha[i].begin(), endc=alpha_alpha[i].end(); c != endc; ++c ) {
					m.add( 1, 'L', *c );
				}
				note_fwd_element_extensions( m.ss, m.abego, 'H', 'A' );
				note_rev_element_extensions( m.ss, m.abego, 'H', 'A' );
				retval.push_back(m);
			}
		}
	} else if ( ( abego1 == "A" ) && ( abego2 == "B" ) ) {
		for ( core::Size i=0; i<numAB; ++i ) {
			if ( lenset.find( ab[i].size() ) != lenset.end() ) {
				Motif m;
				for ( std::string::const_iterator c=ab[i].begin(), endc=ab[i].end(); c != endc; ++c ) {
					m.add( 1, 'L', *c );
				}
				note_fwd_element_extensions( m.ss, m.abego, 'H', 'A' );
				retval.push_back(m);
			}
		}
	} else if ( ( abego1 == "B" ) && ( abego2 == "A" ) ) {
		for ( core::Size i=0; i<numBA; ++i ) {
			if ( lenset.find( ba[i].size() ) != lenset.end() ) {
				Motif m;
				for ( std::string::const_iterator c=ba[i].begin(), endc=ba[i].end(); c != endc; ++c ) {
					m.add( 1, 'L', *c );
				}
				note_rev_element_extensions( m.ss, m.abego, 'H', 'A' );
				retval.push_back(m);
			}
		}
	} else if ( ( abego1 == "B" ) && ( abego2 == "B" ) ) {
		for ( core::Size i=0; i<numBB; ++i ) {
			if ( lenset.find( bb[i].size() ) != lenset.end() ) {
				Motif m;
				for ( std::string::const_iterator c=bb[i].begin(), endc=bb[i].end(); c != endc; ++c ) {
					m.add( 1, 'L', *c );
				}
				retval.push_back(m);
			}
		}
	} else {
		// by default, throw 'X' abegos
		TR << "Warning: using idealized_abego, but secondary structure combination we are connecting (" << abego1 << " --> " << abego2 << " is not AA, AB, BA, or BB. Loop will be all LX" << std::endl;
		for ( std::set< core::Size >::const_iterator l=lenset.begin(), endl=lenset.end(); l != endl; ++l ) {
			retval.push_back( Motif( *l, 'L', "X" ) );
		}
	}
	return retval;
}

/// @brief set lengths based on a vector of possible lengths
void
Connection::set_lengths( utility::vector1< core::Size > const & lengths_val )
{
	motifs_.clear();
	for ( core::Size i=1; i<=lengths_val.size(); ++i ) {
		motifs_.push_back( Motif( lengths_val[i], 'L', "X" ) );
	}
}

/// @brief set possible motifs from a comma-separated string
void
Connection::set_motifs( std::string const & motif_str )
{
	set_motifs( utility::string_split( motif_str, ',' ) );
}

/// @brief set possible motifs from a list of motif strings
void
Connection::set_motifs( utility::vector1< std::string > const & motif_strs )
{
	motifs_.clear();
	for ( core::Size i=1; i<=motif_strs.size(); ++i ) {
		motifs_.push_back( parse_motif( motif_strs[i] ) );
	}
}

/// @brief set possible lengths of the connection with a string
void
Connection::set_cut_resis( std::string const & cut_str )
{
	set_cut_resis( parse_length_string(cut_str) );
}

/// @brief set possible lengths of the connection with a vector
void
Connection::set_cut_resis( utility::vector1< core::Size > const & cuts_val )
{
	cut_resis_ = cuts_val;
}

/// @brief parse a motif string, return a list of paired lengths and ss+abegos
Connection::Motif
Connection::parse_motif( std::string const & motif_str ) const
{
	Motif retval;
	utility::vector1< std::string > motifs = utility::string_split( motif_str, '-' );
	for ( core::Size i=1; i<=motifs.size(); ++i ) {
		// here, we can accept "3LX" or "3:LX"
		std::string motif( "" );
		for ( core::Size j=0; j<motifs[i].size(); ++j ) {
			if ( motifs[i][j] != ':' ) {
				motif += motifs[i][j];
			}
		}
		char const ss_type( motif[motif.size()-2] );
		if ( (ss_type != 'H') && (ss_type != 'L') && (ss_type != 'E') ) {
			TR.Error << "Invalid SS type in motif " << motif << std::endl;
			assert( (ss_type == 'H') || (ss_type == 'E') || (ss_type == 'L') );
			utility_exit();
		}
		std::string const abego_type( motif.substr( motif.size()-1, 1 ) );
		assert( abego_type.size() == 1 );
		if ( abego_type[0] > 'Z' || abego_type[0] < 'A' ) {
			TR.Error << "Invalid abego type in motif " << motif << std::endl;
			utility_exit();
		}
		int const len( utility::string2int( motif.substr( 0, motif.size()-2 ) ) );
		assert( len > 0 );
		retval.add( static_cast< core::Size >(len), ss_type, abego_type );
	}
	return retval;
}

/// @brief length of the connection to be built
core::Size
Connection::build_len( components::StructureData const & perm ) const
{
	return perm.get_data_int( id(), "length" );
}

/// @brief SS of the connection to be built
std::string
Connection::build_ss( components::StructureData const & perm ) const
{
	return perm.get_data_str( id(), "ss" );
}

/// @brief abego of the connection to be built
std::string
Connection::build_abego( components::StructureData const & perm ) const
{
	return perm.get_data_str( id(), "abego" );
}

/// @brief location of the cut to be placed (within the loop)
core::Size
Connection::cut_resi( components::StructureData const & perm ) const
{
	return perm.get_data_int( id(), "cut_resi" );
}

/// @brief stores cut point within the loop into the permutation
void
Connection::set_cut_resi( components::StructureData & perm, core::Size const cut_val ) const
{
	perm.set_data_int( id(), "cut_resi", cut_val );
}

/// @brief stores build length in the permutation
void
Connection::set_build_len( components::StructureData & perm, core::Size const len_val ) const
{
	perm.set_data_int( id(), "length", len_val );
}

/// @brief stores build secondary structure in the permutation
void
Connection::set_build_ss( components::StructureData & perm, std::string const & ss_val ) const
{
	perm.set_data_str( id(), "ss" , ss_val );
}

/// @brief stores build abego in the permutation
void
Connection::set_build_abego( components::StructureData & perm, std::string const & abego_val ) const
{
	perm.set_data_str( id(), "abego", abego_val );
}

/// @brief get id of the component 1 to be built
std::string const &
Connection::lower_segment_id( components::StructureData const & perm ) const
{
	return perm.get_data_str( id(), "segment1" );
}

/// @brief get id of the component 2 to be built
std::string const &
Connection::upper_segment_id( components::StructureData const & perm ) const
{
	return perm.get_data_str( id(), "segment2" );
}

/// @brief get id of the component 2 to be built
std::string const &
Connection::comp1_lower( components::StructureData const & perm ) const
{
	return perm.get_data_str( id(), "segment1_lower" );
}

/// @brief get id of the component 2 to be built
std::string const &
Connection::comp2_upper( components::StructureData const & perm ) const
{
	return perm.get_data_str( id(), "segment2_upper" );
}

/// @brief get id of the lower component actually being connected
std::string const &
Connection::loop_lower( components::StructureData const & perm ) const
{
	return perm.get_data_str( id(), "loop_lower" );
}

/// @brief get id of the upper component actually being connected
std::string const &
Connection::loop_upper( components::StructureData const & perm ) const
{
	return perm.get_data_str( id(), "loop_upper" );
}

/// @brief set id of the component 1 to be built
void
Connection::set_lower_segment_id( components::StructureData & perm, std::string const & comp ) const
{
	//perm.set_data_str( data_name( "segment1" ), add_parent_prefix(comp) );
	perm.set_data_str( id(), "segment1", comp );
}

/// @brief set id of the component 2 to be built
void
Connection::set_upper_segment_id( components::StructureData & perm, std::string const & comp ) const
{
	//perm.set_data_str( data_name( "segment2" ), add_parent_prefix(comp) );
	perm.set_data_str( id(), "segment2", comp );
}

/// @brief set id of the component 1 to be built
void
Connection::set_comp1_lower( components::StructureData & perm, std::string const & comp ) const
{
	//perm.set_data_str( data_name( "segment1_lower" ), add_parent_prefix(comp) );
	perm.set_data_str( id(), "segment1_lower", comp );
}

/// @brief set id of the component 2 to be built
void
Connection::set_comp2_upper( components::StructureData & perm, std::string const & comp ) const
{
	//perm.set_data_str( id(), "segment2_upper" ), add_parent_prefix(comp) );
	perm.set_data_str( id(), "segment2_upper", comp );
}

/// @brief set id of the lower component actually being connected
void
Connection::set_loop_lower( components::StructureData & perm, std::string const & comp ) const
{
	//perm.set_data_str( data_name( "loop_lower" ), add_parent_prefix(comp) );
	perm.set_data_str( id(), "loop_lower", comp );
}

/// @brief set id of the upper component actually being connected
void
Connection::set_loop_upper( components::StructureData & perm, std::string const & comp ) const
{
	//perm.set_data_str( data_name( "loop_upper" ), add_parent_prefix(comp) );
	perm.set_data_str( id(), "loop_upper", comp );
}

/// @brief performs setup and applies the connection
void
Connection::apply( core::pose::Pose & pose )
{
	// must have an ID
	if ( id().empty() ) {
		throw utility::excn::EXCN_Msg_Exception( "No name is specified to connection -- connections must be named!" );
	}
	bool success = false;
	components::StructureDataOP orig = components::StructureData::create_from_pose( pose, "" );
	components::StructureDataOP perm;
	for ( core::Size i = 1; i <= trials_; ++i ) {
		TR.Debug << "Starting trial " << i << " for connection " << id() << std::endl;
		perm = orig->clone();
		debug_assert( perm );
		set_last_move_status( setup_permutation( *perm ) );
		if ( get_last_move_status() != protocols::moves::MS_SUCCESS ) {
			continue;
		}
		apply_permutation( *perm );
		if ( get_last_move_status() == protocols::moves::MS_SUCCESS ) {
			success = true;
			break;
		}
	}
	if ( !success ) {
		TR.Error << "Unable to connect, connection " << id() << " failed." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}
	if ( perm->pose() ) {
		pose = *(perm->pose());
		perm->save_into_pose( pose );
	} else {
		TR.Error << "we should never be here." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "We should never be here." );
	}
}

/// @brief generates permutation information for length -- assemblies like to know length in advance for checking purposes
protocols::moves::MoverStatus
Connection::setup_permutation( components::StructureData & perm ) const
{
	assert( id() != "" );
	core::Real random = numeric::random::rg().uniform();
	return setup_from_random( perm, random );
}

/// @brief applies the connetion to a permutation
void
Connection::apply_permutation( components::StructureData & perm )
{
	set_last_move_status( protocols::moves::MS_SUCCESS );

	components::StructureDataOP orig = perm.clone();

	//perform setup and add necessary residues
	setup( perm );

	// do the work
	apply_connection( perm );

	// only proceed if the move succeeded.
	if ( get_last_move_status() != protocols::moves::MS_SUCCESS ) {
		TR << "Mover reported something other than success -- giving back unmodified pose." << std::endl;
		// reset perm
		perm = *(orig->clone());
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}
	std::string const & c1 = loop_lower(perm);
	std::string const & c2 = loop_upper(perm);
	TR << id() << " : " << "abego=" << perm.segment(c1).abego() << perm.segment(c2).abego()
		<< " ss=" << perm.segment(c1).ss() << " " << perm.segment(c2).ss() << std::endl;
	perm.connect_segments( c1, c2 );
	if ( ( c2 == c1 + "_1" ) &&
			( perm.segment(c1).cterm_resi()+1 == perm.segment(c2).nterm_resi() ) ) {
		perm.merge_segments( c1, c2, c1 );
	}

	if ( build_len(perm) > 0 ) {
		perm.set_cutpoint( c1, cut_resi(perm) );
		if ( do_remodel() && perm.pose() ) {
			protocols::moves::DsspMover dssp;
			perm.apply_mover(dssp);
			bool const checkval = check(perm);
			if ( !checkval ) {
				TR << "Failed ss/abego check." << std::endl;
				set_last_move_status( protocols::moves::FAIL_RETRY );
				return;
			}
		}
	}

	if ( perm.pose() ) {
		TR << "Final fold tree after " << id() << " = " << perm.pose()->fold_tree() << std::endl;
	}
}

/// @brief easier way to call setup() below
void
Connection::setup( components::StructureData & perm )
{
	setup( perm, comp1_lower(perm), lower_segment_id(perm), upper_segment_id(perm), comp2_upper(perm) );
}

/// @brief does jump setup work and creates loop residues
void
Connection::setup(
	components::StructureData & perm,
	std::string const & comp1_n,
	std::string const & comp1_c,
	std::string const & comp2_n,
	std::string const & comp2_c )
{
	// primary root is component 1, secondary root is component 2
	utility::vector1< std::string > roots;
	roots.push_back( lower_segment_id(perm) );
	roots.push_back( upper_segment_id(perm) );
	TR.Debug << "Consolidating with roots = " << roots << std::endl;
	perm.consolidate_movable_groups( roots );

	// make the loop
	create_loop( perm, comp1_n, comp1_c, comp2_n, comp2_c );
	// only do variants/bond declaration if the pose exists
	if ( perm.pose() ) {
		// remove terminal variants
		if ( perm.pose()->residue( perm.segment(loop_lower(perm)).cterm_resi() ).is_upper_terminus() ) {
			perm.remove_upper_terminus_variant_type( perm.segment(loop_lower(perm)).cterm_resi() );
		}
		if ( perm.pose()->residue( perm.segment(loop_upper(perm)).nterm_resi() ).is_lower_terminus() ) {
			perm.remove_lower_terminus_variant_type( perm.segment(loop_upper(perm)).nterm_resi() );
		}

		// declare bond if necessary -- this also adds cutpoint variants
		if ( polymer_connection() &&
				! perm.polymer_bond_exists( loop_lower(perm), loop_upper(perm) ) ) {
			perm.declare_polymer_bond( loop_lower(perm), loop_upper(perm) );
		}
	}
}

/// @brief finds usable/available upper termini (i.e. those for comp1)
utility::vector1< std::string >
Connection::find_available_upper_termini( components::StructureData const & perm ) const
{
	utility::vector1< std::string > local_ids = comp1_ids_;
	if ( local_ids.empty() ) {
		local_ids = perm.available_upper_termini();
	}

	utility::vector1< std::string > mod;
	for ( StringVec::const_iterator i = local_ids.begin(); i != local_ids.end(); ++i ) {
		if ( perm.has_free_upper_terminus( *i ) ) {
			mod.push_back( *i );
		} else {
			std::string const segname = add_parent_prefix( *i );
			if ( perm.has_free_upper_terminus( segname ) ) {
				mod.push_back( segname );
			}
		}
	}
	TR.Debug << "Available upper termini: " << mod << std::endl;
	return mod;
}

/// @brief finds usable/available upper termini (i.e. those for comp1)
utility::vector1< std::string >
Connection::find_available_lower_termini( components::StructureData const & perm ) const
{
	utility::vector1< std::string > local_ids = comp2_ids_;
	if ( local_ids.empty() ) {
		local_ids = perm.available_lower_termini();
	}

	utility::vector1< std::string > mod;
	for ( StringVec::const_iterator i = local_ids.begin(); i != local_ids.end(); ++i ) {
		if ( perm.has_free_lower_terminus( *i ) ) {
			mod.push_back( *i );
		} else {
			std::string const segname = add_parent_prefix( *i );
			if ( perm.has_free_lower_terminus( segname ) ) {
				mod.push_back( segname );
			}
		}
	}
	TR.Debug << "Available lower termini: " << mod << std::endl;
	return mod;
}

/// @brief sets the given pair of segments as not allowed
void
Connection::clear_disallowed_pairs()
{
	disallowed_pairs_.clear();
}

/// @brief sets the given pair of segments as not allowed
void
Connection::disallow_pair( std::string const & c1, std::string const & c2 )
{
	disallowed_pairs_.insert( std::make_pair( add_parent_prefix(c1), add_parent_prefix(c2) ) );
}

/// @brief returns true if this specific pairing of connection hasn't been explicitly disallowed by the user
bool
Connection::pair_allowed( std::string const & c1, std::string const & c2 ) const
{
	return disallowed_pairs_.find( std::make_pair( add_parent_prefix(c1), add_parent_prefix(c2) ) ) == disallowed_pairs_.end();
}

/// @brief sets up the connection mover based on the information in the permutation and a random number
protocols::moves::MoverStatus
Connection::setup_from_random( components::StructureData & perm, core::Real random ) const
{
	// if a segment exists that was built by a connection of the same name, maintain the same termini as before
	utility::vector1< std::string > local_comp1_ids;
	utility::vector1< std::string > local_comp2_ids;

	// determine alternate id -- strip off all parent names
	std::string alt_id = id();
	for ( StringList::const_iterator c = perm.segments_begin(); c != perm.segments_end(); ++c ) {
		if ( boost::ends_with( *c, id() ) ) {
			alt_id = *c;
		}
	}

	if ( perm.has_segment( alt_id ) ) {
		// if we are re-calling a connection, assume we want to rebuild it
		// copy conection points and delete the old loop
		local_comp1_ids.push_back( lower_segment_id( perm ) );
		local_comp2_ids.push_back( upper_segment_id( perm ) );
		perm.delete_segment( alt_id );
		perm.chains_from_termini();
	} else {
		if ( chain1_ ) {
			utility::vector1< std::string > const available_uppers = find_available_upper_termini( perm );
			debug_assert( chain1_ <= available_uppers.size() );
			local_comp1_ids.push_back( available_uppers[ chain1_ ] );
		} else {
			local_comp1_ids = find_available_upper_termini( perm );
		}
		if ( chain2_ ) {
			utility::vector1< std::string > const available_lowers = find_available_lower_termini( perm );
			debug_assert( chain2_ <= available_lowers.size() );
			local_comp2_ids.push_back( available_lowers[ chain2_ ] );
		} else {
			local_comp2_ids = find_available_lower_termini( perm );
		}
	}

	if ( local_comp1_ids.empty() ) {
		std::stringstream err;
		err << "Connection " << id() << ": " << " no available segment1 upper termini were found matching the user's input.";
		err << "Input ids: " << comp1_ids_ << " Perm: " << std::endl;
		err << perm << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	if ( local_comp2_ids.empty() ) {
		std::stringstream err;
		err << "Connection " << id() << ": " << " no available segment2 lower termini were found matching the user's input.";
		err << "Input ids: " << comp2_ids_ << " Perm: " << std::endl;
		err << perm << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	assert( local_comp1_ids.size() > 0 );
	assert( local_comp2_ids.size() > 0 );

	std::set< core::Size > lenset;
	if ( idealized_abego_ ) {
		for ( MotifList::const_iterator m = motifs_.begin(); m != motifs_.end(); ++m ) {
			lenset.insert( m->len );
		}
	}

	// which combinations of components and lengths are connectable
	utility::vector1< utility::vector1< core::Size > > valid_idxs;
	for ( core::Size i=1; i<=local_comp1_ids.size(); ++i ) {
		for ( core::Size j=1; j<=local_comp2_ids.size(); ++j ) {
			if ( !pair_allowed( local_comp1_ids[i], local_comp2_ids[j] ) ) {
				TR.Debug << "c1, c2, len : connectable " << local_comp1_ids[i] << ", " << local_comp2_ids[j] << " DISALLOWED by user setting." << std::endl;
				continue;
			}

			// fill in motifs
			MotifList motifs;
			if ( idealized_abego_ ) {
				motifs = calc_idealized_motifs(
					perm.abego()[ perm.segment( local_comp1_ids[i] ).stop() ],
					perm.abego()[ perm.segment( local_comp2_ids[j] ).start() ],
					lenset );
			} else {
				motifs = motifs_;
			}
			assert( motifs.size() > 0 );
			for ( core::Size k=1; k<=motifs.size(); ++k ) {
				// if we aren't doing remodel or orienting, don't use distance
				bool const use_distance = !performs_orientation() && do_remodel() && perm.pose();
				core::Real avg_dist = 0.0;
				if ( use_distance ) {
					avg_dist = calc_approx_loop_length(
						perm.abego()[ perm.segment( local_comp1_ids[i] ).stop() ] +
						motifs[k].abego +
						perm.abego()[ perm.segment( local_comp2_ids[j] ).start() ] );
					avg_dist /= static_cast< core::Real >( motifs[k].len );
				}
				if ( perm.are_connectable(
						local_comp1_ids[i],
						local_comp2_ids[j],
						motifs[k].len,
						use_distance,
						performs_orientation(),
						allow_cyclic_, connecting_bond_dist(), avg_dist ) ) {
					TR << "c1, c2, len : connectable " << local_comp1_ids[i] << ", " << local_comp2_ids[j] << ", " << motifs[k] << std::endl;
					utility::vector1< core::Size > params;
					params.push_back( i );
					params.push_back( j );
					params.push_back( k );
					valid_idxs.push_back( params );
				}
			}
		}
	}

	if ( !valid_idxs.size() ) {
		return protocols::moves::FAIL_DO_NOT_RETRY;
	}

	// now we just choose a valid combo from the list
	core::Size const combo_idx = extract_int( random, 1, valid_idxs.size() );
	assert( valid_idxs[combo_idx].size() == 3 );

	// remove anchors
	std::string const & c1 = local_comp1_ids[ valid_idxs[combo_idx][1] ];
	std::string const & c2 = local_comp2_ids[ valid_idxs[combo_idx][2] ];
	TR << "C1=" << c1 << " C2=" << c2 << std::endl;

	core::Size motifidx = valid_idxs[combo_idx][3];
	MotifList motifs;
	if ( idealized_abego_ ) {
		motifs = calc_idealized_motifs(
			perm.abego()[ perm.segment(c1).stop() ],
			perm.abego()[ perm.segment(c2).start() ],
			lenset );
	} else {
		motifs = motifs_;
	}

	assert( motifs.size() );
	assert( motifs.size() >= motifidx );
	Motif const & motif = motifs[ motifidx ];

	std::pair< std::string, std::string > comp1_segment( perm.termini(c1) );
	std::pair< std::string, std::string > comp2_segment( perm.termini(c2) );
	assert( comp1_segment.second == c1 );
	assert( comp2_segment.first == c2 );
	set_build_len( perm, motif.len );
	set_build_ss( perm, motif.ss );
	set_build_abego( perm, motif.abego );
	set_lower_segment_id( perm, c1 );
	set_upper_segment_id( perm, c2 );
	set_comp1_lower( perm, comp1_segment.first );
	set_comp2_upper( perm, comp2_segment.second );
	set_loop_lower( perm, c1 );
	set_loop_upper( perm, c2 );

	// find cutpoint (relative to built loop) if it's not set
	if ( motif.len && ( cut_resis_.size() > 0 ) && (!performs_orientation()) ) {
		if ( cut_resis_.size() == 1 ) {
			set_cut_resi( perm, cut_resis_[1] );
		} else {
			set_cut_resi( perm, cut_resis_[ extract_int( random, 1, cut_resis_.size() ) ] );
		}
	} else if ( motif.len == 0 ) {
		set_cut_resi( perm, 0 );
	} else if ( performs_orientation() ) {
		set_cut_resi( perm, 0 );
	} else {
		if ( motif.len == 1 ) {
			set_cut_resi( perm, 0 );
		} else {
			set_cut_resi( perm, extract_int( random, 1, motif.len - 1 ) );
		}
	}
	assert( ( motif.len == 0 ) || ( cut_resi( perm ) < motif.len ) );

	assert( perm.has_free_upper_terminus( loop_lower( perm ) ) );
	assert( perm.has_free_lower_terminus( loop_upper( perm ) ) );
	perm.mark_connected( loop_lower( perm ), loop_upper( perm ) );

	TR << "Going to connect " << loop_lower( perm ) << "-> " << loop_upper( perm )
		<< " len=" << build_len( perm ) << " cut=" << cut_resi( perm ) << std::endl;
	return protocols::moves::MS_SUCCESS;
}

/// @brief adds loop residues in extended conformation
void
Connection::add_loop_residues(
	components::StructureData & perm,
	std::string const & comp,
	std::string const & loop_name,
	core::Size const num_residues,
	std::string const & ss,
	utility::vector1< std::string > const & abego,
	bool const prepend = false ) const
{
	core::conformation::ResidueCOP newres;
	if ( perm.pose() ) {
		assert( perm.pose()->total_residue() > 2 );
		newres = core::conformation::ResidueFactory::create_residue(
			perm.pose()->residue(2).residue_type_set().name_map("VAL") );
	}
	if ( prepend ) {
		perm.prepend_extended_loop( comp, loop_name, num_residues, ss, abego, newres );
	} else {
		perm.append_extended_loop( comp, loop_name, num_residues, ss, abego, newres );
	}
}

/// @brief adds loop residues based on build_len and cut_resi
void
Connection::create_loop(
	components::StructureData & perm,
	std::string const & comp1_n,
	std::string const & comp1_c,
	std::string const & comp2_n,
	std::string const & comp2_c )
{
	// list of loop residues to propagate
	utility::vector1< core::Size > loop_residues;

	// detect possible cyclic connection
	std::set< std::string > chain1 = perm.connected_segments( comp1_c );
	bool cyclic_connection = ( chain1.find( comp2_n ) != chain1.end() );
	if ( ( !cyclic_connection ) && ( perm.segment(comp1_c).cterm_resi()+1 != perm.segment(comp2_n).nterm_resi() ) ) {
		TR << "Moving so that " << comp1_n << "__" << comp1_c << " precedes " << comp2_n << "__" << comp2_c << std::endl;
		perm.move_segment( comp1_n, comp1_c, comp2_n, comp2_c );
	} else {
		TR << "Not moving : " << comp1_n << "__" << comp1_c << " already precedes " << comp2_n << "__" << comp2_c << " cyclic: " << cyclic_connection << std::endl;
		perm.move_jumps_to_safety();
	}
	if ( cyclic_connection ) {
		TR << "Cyclic connection discovered: " << chain1 << " : " << perm << std::endl;
	}

	assert( cyclic_connection || ( perm.segment(comp1_c).cterm_resi()+1 == perm.segment(comp2_n).nterm_resi() ) );

	TR << "Cut resi is " << cut_resi(perm) << " out of " << build_len(perm) << std::endl;
	core::Size prepend_residues = 0, append_residues = 0;
	std::string prepend_ss = "", append_ss = "";
	utility::vector1< std::string > prepend_abego, append_abego;

	if ( cut_resi(perm) ) {
		prepend_residues = build_len(perm)-cut_resi(perm);
		append_residues = cut_resi(perm);
	} else {
		prepend_residues = 0;
		append_residues = build_len(perm);
	}

	std::string const ss = build_ss(perm);
	std::string const ab = build_abego(perm);
	for ( core::Size i=1; i<=append_residues; ++i ) {
		append_ss += ss[i-1];
		std::string achar = "";
		achar += ab[i-1];
		append_abego.push_back( achar );
	}
	for ( core::Size i=1; i<=prepend_residues; ++i ) {
		prepend_ss += ss[append_residues+i-1];
		std::string achar = "";
		achar += ab[append_residues+i-1];
		prepend_abego.push_back( achar );
	}
	assert( append_residues == append_ss.size() );
	assert( append_residues == append_abego.size() );
	assert( prepend_residues == prepend_ss.size() );
	assert( prepend_residues == prepend_abego.size() );
	assert( append_residues + prepend_residues == build_len(perm) );
	if ( prepend_residues ) {
		set_loop_upper( perm, id() + "_1" );
		add_loop_residues( perm, upper_segment_id(perm), loop_upper(perm), prepend_residues, prepend_ss, prepend_abego, true );
	}
	if ( append_residues ) {
		set_loop_lower( perm, id() );
		add_loop_residues( perm, lower_segment_id(perm), loop_lower(perm), append_residues, append_ss, append_abego, false );
	}
}

/// @brief checks the inserted region vs. the desired ss/abego.  True if it matches, false otherwise
bool
Connection::check( components::StructureData const & perm ) const
{
	if ( check_abego_ ) {
		assert( perm.pose() );
		if ( build_len(perm) > 0 ) {
			utility::vector1< std::string > abego = protocols::denovo_design::abego_vector( build_abego(perm) );
			TR << "about to check for " << abego << std::endl;
			return check_insert(
				*(perm.pose()),
				build_ss(perm),
				abego,
				build_left(perm) );
		}
	}
	return true;
}

/// @brief checks the inserted region vs. the desired ss/abego.  True if it matches, false otherwise
bool
Connection::check_insert(
	core::pose::Pose const & pose,
	std::string const & build_ss,
	utility::vector1< std::string > const & build_abego,
	core::Size const left ) const
{
	std::string const & pose_ss( pose.secstruct() );
	utility::vector1< std::string > pose_abego = core::sequence::get_abego( pose, 2 );
	return check_insert_ss_and_abego( pose_ss, pose_abego, build_ss, build_abego, left );
}

/// @brief get "left" residue of loop region taking overlap into account
core::Size
Connection::build_left( components::StructureData const & perm ) const
{
	components::Segment const & c1 = perm.segment(lower_segment_id(perm));
	core::Size const loopstart = c1.cterm_resi()+1;

	// setup left and right overlap regions
	core::Size left = loopstart;
	core::Size inc_count = 0;
	assert( perm.pose() );
	while ( ( inc_count < overlap_ ) &&
			( left > 1 ) &&
			( left <= perm.pose()->total_residue() ) &&
			( perm.pose()->chain(left) == perm.pose()->chain(c1.safe()) ) ) {
		--left;
		++inc_count;
	}
	return left;
}

/// @brief get "right" residue of loop region taking overlap into account
core::Size
Connection::build_right( components::StructureData const & perm ) const
{
	components::Segment const & c2 = perm.segment(upper_segment_id(perm));
	core::Size const loopend = c2.nterm_resi()-1;
	core::Size right = loopend;
	core::Size inc_count = 0;
	assert( perm.pose() );
	while ( ( inc_count < overlap_ ) &&
			( right > 1 ) &&
			( right <= perm.pose()->total_residue() ) &&
			( perm.pose()->chain(right) == perm.pose()->chain(c2.safe()) ) ) {
		++right;
		++inc_count;
	}
	return right;
}

/// @brief "Dummy" connection used when we wnat to track connection information but not build anything
GenericConnection::GenericConnection():
	Connection()
{
	set_performs_orientation( false );
}

GenericConnection::~GenericConnection()
{
}

std::string
GenericConnection::get_name() const
{
	return "GenericConnection";
}

protocols::moves::MoverOP
GenericConnection::fresh_instance() const
{
	return protocols::moves::MoverOP( new GenericConnection() );
}

protocols::moves::MoverOP
GenericConnection::clone() const
{
	return protocols::moves::MoverOP( new GenericConnection( *this ) );
}

/// @brief applies the loop building, in this case does nothing
void
GenericConnection::apply_connection( components::StructureData & sd )
{
	// strip off terminal residues
	sd.delete_trailing_residues( loop_lower( sd ) );
	sd.delete_leading_residues( loop_upper( sd ) );

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// StapleTomponents

// creator functions
protocols::moves::MoverOP
StapleTomponentsCreator::create_mover() const
{
	return protocols::moves::MoverOP( new StapleTomponents() );
}

std::string
StapleTomponentsCreator::keyname() const
{
	return StapleTomponentsCreator::mover_name();
}

std::string
StapleTomponentsCreator::mover_name()
{
	return "StapleTomponents";
}

// default constructor
StapleTomponents::StapleTomponents() :
	Connection()
{
	set_performs_orientation( true );
}

// destructor
StapleTomponents::~StapleTomponents()
{}

/// @brief name of the mover (i.e. type )
std::string
StapleTomponents::get_name() const
{
	return StapleTomponentsCreator::mover_name();
}

/// @brief return a fresh instance of this class in an owning pointer
protocols::moves::MoverOP
StapleTomponents::fresh_instance() const
{
	return protocols::moves::MoverOP( new StapleTomponents() );
}

/// @brief return a cloned copy of this class in an owning pointer
protocols::moves::MoverOP
StapleTomponents::clone() const
{
	return protocols::moves::MoverOP( new StapleTomponents( *this ) );
}

/// @brief setup the parameters via an xml tag
void
StapleTomponents::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	Connection::parse_my_tag( tag, data, filters, movers, pose );
}

/// @brief sets up the connection mover based on the information in the permutation
protocols::moves::MoverStatus
StapleTomponents::setup_permutation( components::StructureData & perm ) const
{
	protocols::moves::MoverStatus retval = setup_from_random( perm, numeric::random::rg().uniform() );
	// no cutpoint ever for staple connections
	set_cut_resi( perm, 0 );
	if ( perm.segment( lower_segment_id( perm ) ).movable_group == perm.segment( upper_segment_id( perm ) ).movable_group ) {
		TR.Error << "Neither of the chains specified to the connection " << id() << " are movable with respect to one another. If you are using TomponentAssembly, this may be because both components specified are subcomponents of components of the assembly. Group for " << lower_segment_id(perm) << " = " << perm.segment(lower_segment_id(perm)).movable_group << " and " << upper_segment_id(perm) << " = " << perm.segment(upper_segment_id(perm)).movable_group << std::endl;
		runtime_assert( false );
	}
	return retval;
}

/// @brief applies the loop building
void
StapleTomponents::apply_connection( components::StructureData & perm )
{
	if ( !perm.pose() ) {
		perm.delete_trailing_residues( loop_lower(perm) );
		perm.delete_leading_residues( loop_upper(perm) );
		return;
	}

	core::pose::Pose const & pose = *(perm.pose());
	int jump = perm.find_jump( upper_segment_id(perm) );
	assert( ( jump != 0 ) || allow_cyclic() );

	TR << "found jump " << std::endl;
	core::Real c1e_psi = pose.psi( perm.segment(loop_lower(perm)).stop() );
	core::Real c1e_omega = pose.omega( perm.segment(loop_lower(perm)).stop() );
	core::Real c2s_phi = pose.phi( perm.segment(loop_upper(perm)).start() );

	// move segment 2 and align with end of comp1
	perm.align_segments( loop_lower(perm), upper_segment_id(perm) );

	TR << "aligned segments " << std::endl;
	// delete residue(s) at c-terminus of N component and renumber
	perm.delete_trailing_residues( loop_lower(perm) );
	perm.delete_leading_residues( loop_upper(perm) );

	// physically join the termini
	perm.delete_jump_and_intervening_cutpoint( loop_lower(perm), loop_upper(perm) );

	perm.set_psi( perm.segment(loop_lower(perm)).stop(), c1e_psi );
	perm.set_omega( perm.segment(loop_lower(perm)).stop(), c1e_omega );
	perm.set_phi( perm.segment(loop_upper(perm)).start(), c2s_phi );

	for ( core::Size i=perm.segment(lower_segment_id(perm)).stop(); i<=perm.segment(upper_segment_id(perm)).start(); ++i ) {
		TR.Debug << "res " << i << " phi=" << perm.pose()->phi(i) << " psi=" << perm.pose()->psi(i) << " omega=" << perm.pose()->omega(i) << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// BridgeTomponents

// creator functions
protocols::moves::MoverOP
BridgeTomponentsCreator::create_mover() const
{
	return protocols::moves::MoverOP( new BridgeTomponents() );
}

std::string
BridgeTomponentsCreator::keyname() const
{
	return BridgeTomponentsCreator::mover_name();
}

std::string
BridgeTomponentsCreator::mover_name()
{
	return "BridgeByKIC";
}

// default constructor
BridgeTomponents::BridgeTomponents() :
	Connection(),
	kic_template_()
{
}

// destructor
BridgeTomponents::~BridgeTomponents()
{}

/// @brief name of the mover (i.e. type )
std::string
BridgeTomponents::get_name() const
{
	return BridgeTomponentsCreator::mover_name();
}

/// @brief return a fresh instance of this class in an owning pointer
protocols::moves::MoverOP
BridgeTomponents::fresh_instance() const
{
	return protocols::moves::MoverOP( new BridgeTomponents( *this ) );
}

/// @brief return a cloned copy of this class in an owning pointer
protocols::moves::MoverOP
BridgeTomponents::clone() const
{
	return protocols::moves::MoverOP( new BridgeTomponents( *this ) );
}

/// @brief setup the parameters via an xml tag
void
BridgeTomponents::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	protocols::generalized_kinematic_closure::GeneralizedKICOP kic =
		protocols::generalized_kinematic_closure::GeneralizedKICOP( new protocols::generalized_kinematic_closure::GeneralizedKIC() );
	kic->parse_my_tag( tag, data, filters, movers, pose );
	set_kic_mover( kic );
	Connection::parse_my_tag( tag, data, filters, movers, pose );
	if ( !tag->hasOption( "allow_cyclic" ) ) {
		set_allow_cyclic( true );
	}
	// loop for loop length of 0 -- this is illegal for KIC
	for ( core::Size i=1, s=motifs().size(); i<=s; ++i ) {
		if ( !motifs()[i].len ) {
			throw utility::excn::EXCN_RosettaScriptsOption("No loop lengths are specified to ConnectByKIC named " + id());
		}
	}
}

/// @brief sets kic selector scorefunction
void
BridgeTomponents::set_selector_scorefxn( core::scoring::ScoreFunctionCOP scorefxn )
{
	assert( kic_template_ );
	protocols::generalized_kinematic_closure::GeneralizedKICOP newkic(
		new protocols::generalized_kinematic_closure::GeneralizedKIC( *kic_template_ ) );
	newkic->set_selector_scorefunction( scorefxn->clone() );
	set_kic_mover( newkic );
}

/// @brief generates a list of loop residues for the given permutation
/// includes N- and C- terminal residues as anchors
/// second number is the overlap before the loop
std::pair< utility::vector1< core::Size >, core::Size >
BridgeTomponents::compute_loop_residues( components::StructureData & perm ) const
{
	perm.delete_trailing_residues(loop_lower(perm));
	perm.delete_leading_residues(loop_upper(perm));

	utility::vector1< core::Size > new_loop_residues;
	core::Size s = perm.segment(loop_lower(perm)).nterm_resi();
	core::Size pre_overlap = 0;
	// add overlap so that start/end anchors are not in the loop
	for ( core::Size i=1; i<=overlap(); ++i ) {
		if ( is_lower_terminus( *(perm.pose()), s ) ) {
			break;
		}
		--s;
		++pre_overlap;
	}
	core::Size e = perm.segment(loop_lower(perm)).cterm_resi();
	for ( core::Size i=s; i<=e; ++i ) {
		new_loop_residues.push_back( i );
	}

	s = perm.segment(loop_upper(perm)).nterm_resi();
	e = perm.segment(loop_upper(perm)).cterm_resi();
	// add overlap to end of loop
	for ( core::Size i=1; i<=overlap(); ++i ) {
		if ( is_upper_terminus( *(perm.pose()), e ) ) {
			break;
		}
		++e;
	}
	for ( core::Size i=s; i<=e; ++i ) {
		new_loop_residues.push_back( i );
	}
	return std::make_pair( new_loop_residues, pre_overlap );
}

/// @brief connects components via a fragment-based loop
void
BridgeTomponents::apply_connection( components::StructureData & perm )
{
	assert( build_len(perm) );

	std::pair< utility::vector1< core::Size >, core::Size > new_loop_residues =
		compute_loop_residues( perm );
	TR.Debug << "Loop residues are " << new_loop_residues.first << std::endl;

	// create and run kic protocol to close the loop
	protocols::generalized_kinematic_closure::GeneralizedKICOP kic = create_kic_mover( perm, new_loop_residues.first, new_loop_residues.second );
	assert( kic );
	perm.apply_mover( kic );
	if ( !kic->last_run_successful() || ( kic->get_last_move_status() != protocols::moves::MS_SUCCESS ) ) {
		TR.Debug << "KIC failed." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}

	// clean up fold tree if necessary
	TR << "Terminal residues are " << perm.segment(loop_lower(perm)).cterm_resi() << " and " << std::endl;
	perm.delete_jump_and_intervening_cutpoint( loop_lower(perm), loop_upper(perm) );
}

/// @brief setup kic for simple closure
void
BridgeTomponents::setup_kic_closure(
	components::StructureData const & perm,
	protocols::generalized_kinematic_closure::GeneralizedKICOP kic,
	utility::vector1< core::Size > const & lres,
	core::Size const pre_overlap ) const
{
	// add pivots to KIC
	core::Size const loop_midpoint = lres[(lres.size()/2)+1];
	kic->set_pivot_atoms( lres[1], "CA", loop_midpoint, "CA", lres[lres.size()], "CA" );
	assert( cut_resi(perm) );
	core::Size const cut_idx = cut_resi(perm)+pre_overlap;
	assert( cut_idx );
	assert( cut_idx < lres.size() );
	kic->close_bond( lres[cut_idx], "C", lres[cut_idx+1], "N",
		lres[cut_idx], "CA", //prioratom
		lres[cut_idx+1], "CA", //followingatom
		1.32, //bondlength
		123, 114, //angle1, angle2
		180, // torsion
		false,
		false );
}

/// @brief creates and configures kic mover -- assumes loop_residues vector has been set up
protocols::generalized_kinematic_closure::GeneralizedKICOP
BridgeTomponents::create_kic_mover(
	components::StructureData const & perm,
	utility::vector1< core::Size > const & lres,
	core::Size const pre_overlap ) const
{
	assert( kic_template_ );
	assert( lres.size() );
	protocols::generalized_kinematic_closure::GeneralizedKICOP kic =
		protocols::generalized_kinematic_closure::GeneralizedKICOP(
		new protocols::generalized_kinematic_closure::GeneralizedKIC( *kic_template_ ) );
	kic->clear_loop_residues();
	kic->clear_perturber_residue_lists();
	for ( core::Size i=1; i<=lres.size(); ++i ) {
		kic->add_loop_residue( lres[i] );
		kic->add_residue_to_perturber_residue_list( lres[i] );
	}
	setup_kic_closure( perm, kic, lres, pre_overlap );
	return kic;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ConnectTerminiWithDisulfide

// creator functions
protocols::moves::MoverOP
ConnectTerminiWithDisulfideCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ConnectTerminiWithDisulfide() );
}

std::string
ConnectTerminiWithDisulfideCreator::keyname() const
{
	return ConnectTerminiWithDisulfideCreator::mover_name();
}

std::string
ConnectTerminiWithDisulfideCreator::mover_name()
{
	return "ConnectTerminiWithDisulfide";
}

// default constructor
ConnectTerminiWithDisulfide::ConnectTerminiWithDisulfide() :
	BridgeTomponents()
{
	set_connecting_bond_dist( 6 );
}

// destructor
ConnectTerminiWithDisulfide::~ConnectTerminiWithDisulfide()
{}

/// @brief name of the mover (i.e. type )
std::string
ConnectTerminiWithDisulfide::get_name() const
{
	return ConnectTerminiWithDisulfideCreator::mover_name();
}

/// @brief return a fresh instance of this class in an owning pointer
protocols::moves::MoverOP
ConnectTerminiWithDisulfide::fresh_instance() const
{
	return protocols::moves::MoverOP( new ConnectTerminiWithDisulfide() );
}

/// @brief return a cloned copy of this class in an owning pointer
protocols::moves::MoverOP
ConnectTerminiWithDisulfide::clone() const
{
	return protocols::moves::MoverOP( new ConnectTerminiWithDisulfide( *this ) );
}

/// @brief given a pose and list of loop residues, creates a CYD pair
std::pair< core::Size, core::Size >
ConnectTerminiWithDisulfide::create_cyd_pair( components::StructureData & perm, utility::vector1< core::Size > const & loop_residues ) const
{
	assert( cut_resi(perm) );
	assert( loop_residues.size() >= cut_resi(perm)+1 );
	TR << "LOop residues are: " << loop_residues << std::endl;
	core::Size const pos1 = loop_residues[cut_resi(perm)];
	core::Size const pos2 = loop_residues[cut_resi(perm)+1];
	TR << "Making CYD at positions " << pos1 << " and " << pos2 << std::endl;
	using namespace protocols::simple_moves;
	protocols::simple_moves::MutateResidueOP make_cyd1 = MutateResidueOP( new MutateResidue( pos1, "CYD" ) );
	perm.apply_mover( make_cyd1 );
	protocols::simple_moves::MutateResidueOP make_cyd2 = MutateResidueOP( new MutateResidue( pos2, "CYD" ) );
	perm.apply_mover( make_cyd2 );
	utility::pointer::shared_ptr< protocols::cyclic_peptide::DeclareBond > disulf =
		utility::pointer::shared_ptr< protocols::cyclic_peptide::DeclareBond >(
		new protocols::cyclic_peptide::DeclareBond() );
	disulf->set(
		pos1, "SG", // res1, atom1
		pos2, "SG", // res2, atom2
		true, // add termini
		false, // run kic
		0, 0, // kic_res1, kic_res2
		false ); //rebuild_fold_tree
	perm.apply_mover( disulf );
	return std::pair< core::Size, core::Size >( pos1, pos2 );
}

/// @brief setup kic for simple closure
void
ConnectTerminiWithDisulfide::setup_disulf_kic_closure(
	protocols::generalized_kinematic_closure::GeneralizedKICOP kic,
	utility::vector1< core::Size > const & loop_residues,
	std::pair< core::Size, core::Size > const & disulf_pos ) const
{
	TR << "Disulfide residues are " << disulf_pos.first << " and " << disulf_pos.second << " loop=" << loop_residues << std::endl;
	// add pivots to KIC
	kic->set_pivot_atoms( loop_residues[1], "CA", disulf_pos.first, "SG", loop_residues[loop_residues.size()], "CA" );
	kic->close_bond( disulf_pos.first, "SG", disulf_pos.second, "SG",
		disulf_pos.first, "CB", //prioratom
		disulf_pos.second, "CB", //followingatom
		2.05, //bondlength
		103, 103, //angle1, angle2
		0.0, // torsion
		false,
		false );
	kic->add_perturber( "randomize_dihedral" );

	utility::vector1 < core::id::NamedAtomID > atomset;
	atomset.push_back( core::id::NamedAtomID( "CA", disulf_pos.first ) );
	atomset.push_back( core::id::NamedAtomID( "CB", disulf_pos.first ) );
	kic->add_atomset_to_perturber_atomset_list( atomset );

	atomset.clear();
	atomset.push_back( core::id::NamedAtomID( "CA", disulf_pos.second ) );
	atomset.push_back( core::id::NamedAtomID( "CB", disulf_pos.second ) );
	kic->add_atomset_to_perturber_atomset_list( atomset );

	atomset.clear();
	atomset.push_back( core::id::NamedAtomID( "SG", disulf_pos.first ) );
	atomset.push_back( core::id::NamedAtomID( "CB", disulf_pos.first ) );
	kic->add_atomset_to_perturber_atomset_list( atomset );

	atomset.clear();
	atomset.push_back( core::id::NamedAtomID( "SG", disulf_pos.second ) );
	atomset.push_back( core::id::NamedAtomID( "CB", disulf_pos.second ) );
	kic->add_atomset_to_perturber_atomset_list( atomset );

	atomset.clear();
	atomset.push_back( core::id::NamedAtomID( "SG", disulf_pos.first ) );
	atomset.push_back( core::id::NamedAtomID( "SG", disulf_pos.second ) );
	kic->add_atomset_to_perturber_atomset_list( atomset );

}

/// @brief connects components via a fragment-based loop
void
ConnectTerminiWithDisulfide::apply_connection( components::StructureData & perm )
{
	// reset build length so that last pose residues are the end pivots
	assert( build_len(perm) );
	set_build_len( perm, build_len(perm)+2 );
	set_cut_resi( perm, cut_resi(perm)+1 );

	std::pair< utility::vector1< core::Size >, core::Size > new_loop_residues =
		compute_loop_residues( perm );

	// create CYD pair
	std::pair< core::Size, core::Size > cyds = create_cyd_pair( perm, new_loop_residues.first );

	// create and run kic protocol
	protocols::generalized_kinematic_closure::GeneralizedKICOP kic = create_kic_mover( perm, new_loop_residues.first, new_loop_residues.second );
	assert( kic );

	// setup kic mover for disulfides
	setup_disulf_kic_closure( kic, new_loop_residues.first, cyds );

	// run kic
	perm.apply_mover( kic );

	set_build_len( perm, build_len(perm)-2 );
	set_cut_resi( perm, cut_resi(perm)-1 );
}

/// @brief setup kic for simple closure
void
ConnectTerminiWithDisulfide::setup_kic_closure(
	components::StructureData const &,
	protocols::generalized_kinematic_closure::GeneralizedKICOP,
	utility::vector1< core::Size > const &,
	core::Size const ) const
{
}


core::Real
calc_approx_loop_length( std::string const & abego )
{
	// these values are based on SIN(ANGLECHANGE/2)*3.8
	static std::map< std::string, core::Real > const distmap =
		boost::assign::map_list_of
		("AA",2.18)
		("AB",1.90)
		("AE",3.79)
		("AG",3.67)
		("BA",3.44)
		("BB",3.57)
		("BE",0.98)
		("BG",0.33)
		("EA",0.33)
		("EB",0.66)
		("EE",2.91)
		("EG",3.57)
		("GA",3.74)
		("GB",3.67)
		("GE",2.91)
		("GG",1.90);
	core::Real max_dist = 0.0;
	for ( core::Size i=1, endi=abego.size(); i<endi; ++i ) {
		std::string const dyad = abego.substr(i-1,2);
		std::map< std::string, core::Real >::const_iterator d = distmap.find(dyad);
		if ( d == distmap.end() ) {
			max_dist += 3.8;
		} else {
			max_dist += d->second;
		}
	}
	return max_dist;
}

/// @brief compares desired insert ss and abego to pose values, returns number of mismatches
core::Size compare_insert_ss_and_abego(
	std::string const & pose_ss,
	utility::vector1< std::string > const & pose_abego,
	std::string const & build_ss,
	utility::vector1< std::string > const & build_abego,
	core::Size const left )
{
	core::Size mismatches = 0;
	TR.Debug << "Going to compare " << pose_ss << " left=" << left << " with " << build_ss << std::endl;
	for ( core::Size i=1; i<=build_ss.size(); ++i ) {
		TR.Debug << "Comparing position " << i-1 << " to " << i+left-2 << std::endl;
		bool broken = pose_ss[i+left-2] != build_ss[i-1];
		// H is acceptable in loops almost always, so I allow it here
		/*if ( ( pose_ss[i+left-2] == 'H' ) && ( build_ss[i-1] == 'L' ) ) {
		broken = false;
		}*/
		// E is also acceptable in loops almost always, so it's also allowed.
		if ( ( pose_ss[i+left-2] == 'E' ) && ( build_ss[i-1] == 'L' ) ) {
			broken = false;
		}
		if ( broken ) {
			TR << "Connection in pose doesn't match the desired secondary structure. Wanted= " << build_ss << "; Actual= " << pose_ss.substr(left-1, build_ss.size()) << " pose=" << pose_ss << std::endl;
			++mismatches;
		} else if ( ( build_abego[i] != "X" ) && ( pose_abego[i+left-1] != build_abego[i] ) ) {
			TR << "ABEGO code in pose doesn't match the desired abego at position " << i+left-1 << ". Wanted=" << build_abego << "; Actual=" << pose_abego << std::endl;
			++mismatches;
		}
	}
	return mismatches;
}

bool check_insert_ss_and_abego(
	std::string const & pose_ss,
	utility::vector1< std::string > const & pose_abego,
	std::string const & build_ss,
	utility::vector1< std::string > const & build_abego,
	core::Size const left )
{
	if ( compare_insert_ss_and_abego( pose_ss, pose_abego, build_ss, build_abego, left ) ) {
		return false;
	} else {
		return true;
	}
}

} // namespace connection
} // namespace denovo_design
} // namespace protocols
