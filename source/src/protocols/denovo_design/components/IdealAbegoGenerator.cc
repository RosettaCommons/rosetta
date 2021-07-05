// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/IdealAbegoGenerator.cc
/// @brief Logic for selection of abego values
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/components/IdealAbegoGenerator.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
//#include <protocols/denovo_design/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility>
#include <utility/excn/Exceptions.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>

// Boost headers

#include <boost/assign/list_of.hpp> // AUTO IWYU For generic_list, assign_decay<>::type

static basic::Tracer TR( "protocols.denovo_design.components.IdealAbegoGeneratorV2" );

namespace protocols {
namespace denovo_design {
namespace components {

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::ab_no_extension = { "GB", "GBA", "BAA", "GBB", "GBAB", "BAAB" };

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::ba_no_extension = { "AB", "GBB", "BAB" };

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::bb_no_extension = { "GG", "EA", "AA", "AAG" };

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::aa_no_extension = {
"B", "G",
"BB", "GB",
"BAB", "GBB", "BBG",
"BAAB",
};

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::aa_no_extension_all = {
"B", "G",
"BB", "BG", "GB", "GG",
"BBB", "BBG", "BGB", "BGG", "BAB", "BAG",
"GBB", "GBG", "GGB", "GGG", "GAB", "GAG",
"BBBB", "BBBG", "BBGB", "BBGG", "BBAB", "BBAG",
"BGBB", "BGBG", "BGGB", "BGGG", "BGAB", "BGAG",
"BABB", "BABG", "BAGB", "BAGG", "BAAB", "BAAG",
"GBBB", "GBBG", "GBGB", "GBGG", "GBAB", "GBAG",
"GGBB", "GGBG", "GGGB", "GGGG", "GGAB", "GGAG",
"GABB", "GABG", "GAGB", "GAGG", "GAAB", "GAAG"
};

std::map< char, core::Size > const
IdealAbegoGenerator::extension_lengths = boost::assign::map_list_of
( '-', 0 ) ( 'B', 1 ) ( 'A', 3 );

IdealAbegoGenerator::Abegos
IdealAbegoGenerator::generate_extended_abegos( Abegos const & no_extension_loops, char const abego_n, char const abego_c )
{
	auto n_it = extension_lengths.find( abego_n );
	if ( n_it == extension_lengths.end() ) {
		std::stringstream msg;
		msg << "Invalid ABEGO_n type: " << abego_n << std::endl;
		CREATE_EXCEPTION( utility::excn::Exception, msg.str() );
	}
	auto c_it = extension_lengths.find( abego_c );
	if ( c_it == extension_lengths.end() ) {
		std::stringstream msg;
		msg << "Invalid ABEGO_c type: " << abego_c << std::endl;
		CREATE_EXCEPTION( utility::excn::Exception, msg.str() );
	}
	core::Size const max_res_n = n_it->second;
	core::Size const max_res_c = c_it->second;
	//TR.Debug << "Expanding " << no_extension_loops << " for " << abego_n << "-->" << abego_c << "..." << std::endl;
	Abegos extended;
	std::string n_abego = "";
	for ( core::Size n_ext=0; n_ext<=max_res_n; ++n_ext ) {
		std::string c_abego = "";
		for ( core::Size c_ext=0; c_ext<=max_res_c; ++c_ext ) {;
			//TR.Debug << " -> N: " << n_ext << " C: " << c_ext << " abegos: " << no_extension_loops << std::endl;
			for ( auto const & loop_abego : no_extension_loops ) {
				std::stringstream extended_loop;
				extended_loop << n_abego << loop_abego << c_abego;
				extended.push_back( extended_loop.str() );
			}
			c_abego += abego_c;
		}
		n_abego += abego_n;
	}
	//TR.Debug << " -> Expanded to " << extended << std::endl;
	//TR << "Allowed loops for " << abego_n << " -> " << abego_c << ": " << extended << std::endl;
	return extended;
}

IdealAbegoGenerator::IdealAbegoGenerator( std::string const & segment_name_val ):
	utility::VirtualBase(),
	segment_name_( segment_name_val ),
	extend_ss_( "1" ),
	match_ss_to_abego_( false ),
	use_hh_rules_2021_( true )
{
	extended_abegos_ = boost::assign::map_list_of
		// everything
		( "X-X-", IdealAbegoGenerator::generate_extended_abegos( aa_no_extension_all, '-', '-' ) )
		( "X-X+", IdealAbegoGenerator::generate_extended_abegos( aa_no_extension_all, '-', 'A' ) )
		( "X+X-", IdealAbegoGenerator::generate_extended_abegos( aa_no_extension_all, 'A', '-' ) )
		( "X+X+", IdealAbegoGenerator::generate_extended_abegos( aa_no_extension_all, 'A', 'A' ) )

		( "A-A-", IdealAbegoGenerator::generate_extended_abegos( aa_no_extension, '-', '-' ) )
		( "A-A+", IdealAbegoGenerator::generate_extended_abegos( aa_no_extension, '-', 'A' ) )
		( "A+A-", IdealAbegoGenerator::generate_extended_abegos( aa_no_extension, 'A', '-' ) )
		( "A+A+", IdealAbegoGenerator::generate_extended_abegos( aa_no_extension, 'A', 'A' ) )

		( "A-B-", IdealAbegoGenerator::generate_extended_abegos( ab_no_extension, '-', '-' ) )
		( "A-B+", IdealAbegoGenerator::generate_extended_abegos( ab_no_extension, '-', 'B' ) )
		( "A+B-", IdealAbegoGenerator::generate_extended_abegos( ab_no_extension, 'A', '-' ) )
		( "A+B+", IdealAbegoGenerator::generate_extended_abegos( ab_no_extension, 'A', 'B' ) )

		( "B-A-", IdealAbegoGenerator::generate_extended_abegos( ba_no_extension, '-', '-' ) )
		( "B-A+", IdealAbegoGenerator::generate_extended_abegos( ba_no_extension, '-', 'A' ) )
		( "B+A-", IdealAbegoGenerator::generate_extended_abegos( ba_no_extension, 'B', '-' ) )
		( "B+A+", IdealAbegoGenerator::generate_extended_abegos( ba_no_extension, 'B', 'A' ) )

		( "B-B-", IdealAbegoGenerator::generate_extended_abegos( bb_no_extension, '-', '-' ) )
		( "B-B+", IdealAbegoGenerator::generate_extended_abegos( bb_no_extension, '-', 'B' ) )
		( "B+B-", IdealAbegoGenerator::generate_extended_abegos( bb_no_extension, 'B', '-' ) )
		( "B+B+", IdealAbegoGenerator::generate_extended_abegos( bb_no_extension, 'B', 'B' ) );

}

IdealAbegoGenerator::~IdealAbegoGenerator() = default;

IdealAbegoGeneratorOP
IdealAbegoGenerator::clone() const
{
	return utility::pointer::make_shared< IdealAbegoGenerator >( *this );
}


bool
IdealAbegoGenerator::extend_ss( char const secstruct ) const
{
	if ( utility::is_true_string( extend_ss_ ) ) {
		return true;
	}
	if ( utility::is_false_string( extend_ss_ ) ) {
		return false;
	}
	return extend_ss_.find( secstruct ) != std::string::npos;
}

IdealAbegoGenerator::Abegos
IdealAbegoGenerator::retrieve_loop_abegos( char const abego1, char const abego2 ) const
{
	char extend_n = '-';
	if ( extend_ss( abego1 ) ) {
		extend_n = '+';
	}
	char extend_c = '-';
	if ( extend_ss( abego2 ) ) {
		extend_c = '+';
	}

	std::stringstream abegos;
	abegos << abego1 << extend_n << abego2 << extend_c;
	auto it = extended_abegos_.find( abegos.str() );
	if ( it == extended_abegos_.end() ) {
		// one or both abegos not found; since we don't know what to do,
		// we can't generate any ideal abegos
		TR.Warning << "Using idealized_abego, but the ABEGO combination we are connecting (" << abego1 << " --> " << abego2 << " is not AA, AB, BA, or BB. Loop will be all LX" << std::endl;
		Abegos abego_strings;
		return abego_strings;
	} else {
		return it->second;
	}
}

static std::map< char, char > const abego2secstruct = boost::assign::map_list_of
( 'A', 'H' ) ( 'B', 'E' ) ( 'X', 'L' );

char
get_secstruct_from_abego( char const abego )
{
	auto it = abego2secstruct.find( abego );
	if ( it == abego2secstruct.end() ) {
		std::stringstream msg;
		msg << "Error: ABEGO type " << abego << " could not be converted into a secondary structure type." << std::endl;
		CREATE_EXCEPTION( utility::excn::Exception, msg.str() );
	}
	return it->second;
}

/// @brief Given lengths and abegos to connect, generate a list of motifs
/// @param[in] abego1 Abego of residue immediately before the loop
/// @param[in] abego2 Abego of residue immediately after the loop
/// @param[in] lenset Set of allowed loop lengths
/// @returns MotifOPs of all ideal loops with given lengths/cutpoints connecting two residues of the given
///          abegos
//
/// @brief Given desired lengths, compute a set of idealized loop motifs via Nobu/Rie/YuRu rules
/// @param[in] abego1 Abego of residue immediately before the loop
/// @param[in] abego2 Abego of residue immediately after the loop
/// @param[in] lenset Set of allowed loop lengths
/// @param[in] cutpoint_set Set of allowed cutpoint residues.  A value of N indicates the Nth residue
///                         of the loop is a LOWER_CUTPOINT and the (N+1)th residue is an UPPER_CUTPOINT
/// @returns MotifOPs of all ideal loops with given lengths/cutpoints connecting two residues of the given
///          abegos
/// @details If extend_ss is true, "extension" of SS elements by adding residues of the same abego type
///          is allowed. For example, connecting A-->A, you might get a 4-residue loop
///          of type A-ABGA-A. If false, only different abegos are allowed, and you might
///          get a 4-residue loop of type A-BGAB-A
IdealAbegoGenerator::MotifOPs
IdealAbegoGenerator::generate(
	char abego1,
	char abego2,
	LengthSet const & lenset,
	LengthSet const & cutpoint_set ) const
{
	TR << "Generating motifs for abego1=" << abego1 << " abego2=" << abego2
		<< " lengths=" << lenset << " cutpoints=" << cutpoint_set << std::endl;

	if ( !use_hh_rules_2021_ ) {
		TR.Debug << "Not using H-H rules; all abegos are allowed" << std::endl;
		abego1 = 'X';
		abego2 = 'X';
	} else {
		TR.Debug << "Using H-H rules" << std::endl;
	}

	Abegos loop_abegos = retrieve_loop_abegos( abego1, abego2 );
	if ( loop_abegos.empty() ) {
		for ( unsigned long l : lenset ) {
			std::string abego( l, 'X' );
			loop_abegos.push_back( abego );
		}
		//retval.push_back( utility::pointer::make_shared< Motif >( segment_name_, secstruct, abego, false, false ) );
	}
	MotifOPs ideal_motifs = extract_ideal_motifs( loop_abegos, lenset );
	if ( match_ss_to_abego_ ) {
		for ( auto & motif : ideal_motifs ) {
			char const ss1 = get_secstruct_from_abego( abego1 );
			if ( ss1 == 'H' ) {
				note_fwd_element_extensions( *motif, ss1, abego1 );
			}
			char const ss2 = get_secstruct_from_abego( abego2 );
			if ( ss2 == 'H' ) {
				note_rev_element_extensions( *motif, ss2, abego2 );
			}
		}
	}

	if ( !cutpoint_set.empty() ) {
		ideal_motifs = add_cutpoints( ideal_motifs, cutpoint_set );
	}

	TR << "Generated ideal motifs: " <<  std::endl;
	for ( MotifOPs::const_iterator m=ideal_motifs.begin(); m!=ideal_motifs.end(); ++m ) {
		TR << **m << std::endl;
	}
	return ideal_motifs;
}

/// @brief takes a list of motifs without cutpoints and generates all permutations with cutpoints from cutpoint_set
/// @param[in] orig         List of motifs without cutpoints
/// @param[in] cutpoint_set Set of loop-relative indices of cutpoint residues
/// @returns MotifOPs of all motifs in orig with all possible cutpoints in cutpoint_set
IdealAbegoGenerator::MotifOPs
IdealAbegoGenerator::add_cutpoints( MotifOPs const & orig, LengthSet const & cutpoint_set ) const
{
	MotifOPs retval;
	for ( auto const & m : orig ) {
		for ( unsigned long c : cutpoint_set ) {
			if ( c > m->elem_length() ) continue;
			MotifOP newmotif( new Motif( *m ) );
			newmotif->set_cutpoint( c );
			retval.push_back( newmotif );
		}
	}
	return retval;
}

IdealAbegoGenerator::MotifOPs
IdealAbegoGenerator::extract_ideal_motifs(
	Abegos const & abegolist,
	LengthSet const & lenset ) const
{
	MotifOPs retval;
	for ( auto const & abego : abegolist ) {
		if ( lenset.find( abego.size() + 2 ) == lenset.end() ) continue;
		std::string const secstruct( 2 + abego.size(), 'L' );
		std::stringstream abego_str;
		abego_str << 'X' << abego << 'X';
		MotifOP m( new Motif( segment_name_, secstruct, abego_str.str(), false, false ) );
		retval.push_back( m );
	}
	return retval;
}

/// @brief if true, secondary structure may be extended to close the loop.  Overall size of the insert
///        will not change.  Therefore, a 4 residue loop might extend a helix by 2 residues and have
///        a 2-residue ideal loop if extend_ss is true.  If extend_ss is false, a 4 residue loop must
///        have four loop residues.
void
IdealAbegoGenerator::set_extend_ss( std::string const & extend_ss )
{
	extend_ss_ = extend_ss;
}

void
IdealAbegoGenerator::note_fwd_element_extensions(
	Motif & motif,
	char const sschar,
	char const element_abego ) const
{
	// starting at beginning, locate extending elements
	auto const length = static_cast< components::SegmentResid >( motif.elem_length() );
	for ( components::SegmentResid resid=1; resid<=length; ++resid ) {
		if ( motif.abego( resid ) != element_abego ) break;
		motif.set_ss( resid, sschar );
	}
}

void
IdealAbegoGenerator::note_rev_element_extensions(
	Motif & motif,
	char const sschar,
	char const element_abego ) const
{
	// starting at end, locate extending elements
	auto const length = static_cast< components::SegmentResid >( motif.elem_length() );
	for ( components::SegmentResid resid=length; resid>=1; --resid ) {
		if ( motif.abego( resid ) != element_abego ) break;
		motif.set_ss( resid, sschar );
	}
}

} //protocols
} //denovo_design
} //components

