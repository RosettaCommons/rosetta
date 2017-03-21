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
#include <utility/stream_util.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.IdealAbegoGenerator" );

namespace protocols {
namespace denovo_design {
namespace components {

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::ab_no_extension =
boost::assign::list_of
("GB") ("GBA") ("BAA") ("GBB") ("GBAB") ("BAAB");

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::ba_no_extension =
boost::assign::list_of
("AB") ("GBB") ("BAB");

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::bb_no_extension =
boost::assign::list_of
("GG") ("EA") ("AA") ("AAG");

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::aa_no_extension =
boost::assign::list_of
("B") ("G")
("BB") ("BG") ("GB") ("GG")
("BBB") ("BBG") ("BGB") ("BGG") ("BAB") ("BAG")
("GBB") ("GBG") ("GGB") ("GGG") ("GAB") ("GAG")
("BBBB") ("BBBG") ("BBGB") ("BBGG") ("BBAB") ("BBAG")
("BGBB") ("BGBG") ("BGGB") ("BGGG") ("BGAB") ("BGAG")
("BABB") ("BABG") ("BAGB") ("BAGG") ("BAAB") ("BAAG")
("GBBB") ("GBBG") ("GBGB") ("GBGG") ("GBAB") ("GBAG")
("GGBB") ("GGBG") ("GGGB") ("GGGG") ("GGAB") ("GGAG")
("GABB") ("GABG") ("GAGB") ("GAGG") ("GAAB") ("GAAG");

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::ab =
boost::assign::list_of
(   "GB" )(   "GBA" )(   "BAA" ) // Basic loops
(   "GBB")(   "GBAB")(   "BAAB") // Extend strand by one residue
(  "AGB" )(  "AGBA" )(  "ABAA" ) // Extend helix by one residue
(  "AGBB")(  "AGBAB")(  "ABAAB")
( "AAGB" )( "AAGBA" )( "AABAA" ) // Extend helix by two residues
( "AAGBB")( "AAGBAB")( "AABAAB")
("AAAGB" )("AAAGBA" )("AAABAA" ) // Extend helix by three residues
("AAAGBB")("AAAGBAB")("AAABAAB");

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::ba =
boost::assign::list_of
("AB"   )("GBB"   )("BAB"   )("BGBB"   )  // Basic loops plus extend strand by one residue
("ABA"  )("GBBA"  )("BABA"  )("BGBBA"  )  // Extend helix by one residue
("ABAA" )("GBBAA" )("BABAA" )("BGBBAA" )  // Extend helix by two residues
("ABAAA")("BGGAAA")("BABAAA")("BGBBAAA");

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::bb =
boost::assign::list_of
( "GG" )( "EA" )( "AA" )( "AAG" )   // Basic loops
("BGG" )("BEA" )("BAA" )("BAAG" )   // Extend first strand by one
( "GGB")( "EAB")( "AAB")( "AAGB")   // Extend second strand by one
("BGGB")("BEAB")("BAAB")("BAAGB"); // Extend both strands by one

IdealAbegoGenerator::Abegos const
IdealAbegoGenerator::alpha_alpha =
boost::assign::list_of
(   "B")("AB")("AAB")("AAAB")("BA")("ABA")("AABA")("AAABA")("BAA")("ABAA")   // one-residue B and extensions
( "AABAA")("AAABAA")("BAAA")("ABAAA")("AABAAA")("AAABAAA")
(   "G")("AG")("AAG")("AAAG")("GA")("AGA")("AAGA")("AAAGA")("GAA")("AGAA")   // one-residue G and extensions
( "AAGAA")("AAAGAA")("GAAA")("AGAAA")("AAGAAA")("AAAGAAA")
(   "BB"   )(   "BG"   )(   "GB"   )(   "BAB"   )(   "BAG"   )(   "GAB"   )  //loops not containing E or start/ending with A
(   "BBA"  )(   "BGA"  )(   "GBA"  )(   "BABA"  )(   "BAGA"  )(   "GABA"  )  // Extend h2 by one residue
(   "BBAA" )(   "BGAA" )(   "GBAA" )(   "BABAA" )(   "BAGAA" )(   "GABAA" )  // Extend h2 by two residues
(   "BBAAA")(   "BGAAA")(   "GBAAA")(   "BABAAA")(   "BAGAAA")(   "GABAAA")  // Extend h2 by three residues
(  "ABB"   )(  "ABG"   )(  "AGB"   )(  "ABAB"   )(  "ABAG"   )(  "AGAB"   )  // Extend h1 by one residue
(  "ABBA"  )(  "ABGA"  )(  "AGBA"  )(  "ABABA"  )(  "ABAGA"  )(  "AGABA"  )  //
(  "ABBAA" )(  "ABGAA" )(  "AGBAA" )(  "ABABAA" )(  "ABAGAA" )(  "AGABAA" )  //
(  "ABBAAA")(  "ABGAAA")(  "AGBAAA")(  "ABABAAA")(  "ABAGAAA")(  "AGABAAA")  //
( "AABB"   )( "AABG"   )( "AAGB"   )( "AABAB"   )( "AABAG"   )( "AAGAB"   )  // Extend h1 by two residue
( "AABBA"  )( "AABGA"  )( "AAGBA"  )( "AABABA"  )( "AABAGA"  )( "AAGABA"  )  //
( "AABBAA" )( "AABGAA" )( "AAGBAA" )( "AABABAA" )( "AABAGAA" )( "AAGABAA" )  //
( "AABBAAA")( "AABGAAA")( "AAGBAAA")( "AABABAAA")( "AABAGAAA")( "AAGABAAA")  //
("AAABB"   )("AAABG"   )("AAAGB"   )("AAABAB"   )("AAABAG"   )("AAAGAB"   )  // Extend h1 by three residue
("AAABBA"  )("AAABGA"  )("AAAGBA"  )("AAABABA"  )("AAABAGA"  )("AAAGABA"  )  //
("AAABBAA" )("AAABGAA" )("AAAGBAA" )("AAABABAA" )("AAABAGAA" )("AAAGABAA" )  //
("AAABBAAA")("AAABGAAA")("AAAGBAAA")("AAABABAAA")("AAABAGAAA")("AAAGABAAA");

IdealAbegoGenerator::IdealAbegoGenerator( std::string const & segment_name_val ):
	utility::pointer::ReferenceCount(),
	segment_name_( segment_name_val ),
	extend_ss_( true ),
	match_ss_to_abego_( false )
{
}

IdealAbegoGenerator::~IdealAbegoGenerator()
{}

IdealAbegoGeneratorOP
IdealAbegoGenerator::clone() const
{
	return IdealAbegoGeneratorOP( new IdealAbegoGenerator( *this ) );
}

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
	char const abego1,
	char const abego2,
	LengthSet const & lenset,
	LengthSet const & cutpoint_set ) const
{
	TR << "Generating motifs for abego1=" << abego1 << " abego2=" << abego2
		<< " lengths=" << lenset << " cutpoints=" << cutpoint_set << std::endl;

	MotifOPs retval;
	// determine previous residue by abego type
	if ( ( abego1 == 'A' ) && ( abego2 == 'A' ) ) {
		if ( extend_ss_ ) {
			retval = extract_ideal_motifs( alpha_alpha, lenset );
			if ( match_ss_to_abego_ ) {
				for ( MotifOPs::iterator m=retval.begin(); m!=retval.end(); ++m ) {
					note_fwd_element_extensions( **m, 'H', 'A' );
					note_rev_element_extensions( **m, 'H', 'A' );
				}
			}
		} else {
			retval = extract_ideal_motifs( aa_no_extension, lenset );
		}
	} else if ( ( abego1 == 'A' ) && ( abego2 == 'B' ) ) {
		if ( extend_ss_ ) {
			retval = extract_ideal_motifs( ab, lenset );
			if ( match_ss_to_abego_ ) {
				for ( MotifOPs::iterator m=retval.begin(); m!=retval.end(); ++m ) {
					note_fwd_element_extensions( **m, 'H', 'A' );
				}
			}
		} else {
			retval = extract_ideal_motifs( ab_no_extension, lenset );
		}
	} else if ( ( abego1 == 'B' ) && ( abego2 == 'A' ) ) {
		if ( extend_ss_ ) {
			retval = extract_ideal_motifs( ba, lenset );
			if ( match_ss_to_abego_ ) {
				for ( MotifOPs::iterator m=retval.begin(); m!=retval.end(); ++m ) {
					note_rev_element_extensions( **m, 'H', 'A' );
				}
			}
		} else {
			retval = extract_ideal_motifs( ba_no_extension, lenset );
		}
	} else if ( ( abego1 == 'B' ) && ( abego2 == 'B' ) ) {
		if ( extend_ss_ ) {
			retval = extract_ideal_motifs( bb, lenset );
		} else {
			retval = extract_ideal_motifs( bb_no_extension, lenset );
		}
	} else {
		// by default, throw 'X' abegos
		TR.Warning << "using idealized_abego, but secondary structure combination we are connecting (" << abego1 << " --> " << abego2 << " is not AA, AB, BA, or BB. Loop will be all LX" << std::endl;
		for ( LengthSet::const_iterator l=lenset.begin(); l!=lenset.end(); ++l ) {
			std::string const secstruct( *l, 'L' );
			std::string const abego( *l, 'X' );
			retval.push_back( MotifOP( new Motif( segment_name_, secstruct, abego, false, false ) ) );
		}
	}

	if ( !cutpoint_set.empty() ) {
		retval = add_cutpoints( retval, cutpoint_set );
	}

	TR << "Generated ideal motifs: " <<  std::endl;
	for ( MotifOPs::const_iterator m=retval.begin(); m!=retval.end(); ++m ) {
		TR << **m << std::endl;
	}
	return retval;
}

/// @brief takes a list of motifs without cutpoints and generates all permutations with cutpoints from cutpoint_set
/// @param[in] orig         List of motifs without cutpoints
/// @param[in] cutpoint_set Set of loop-relative indices of cutpoint residues
/// @returns MotifOPs of all motifs in orig with all possible cutpoints in cutpoint_set
IdealAbegoGenerator::MotifOPs
IdealAbegoGenerator::add_cutpoints( MotifOPs const & orig, LengthSet const & cutpoint_set ) const
{
	MotifOPs retval;
	for ( MotifOPs::const_iterator m=orig.begin(); m!=orig.end(); ++m ) {
		for ( LengthSet::const_iterator c=cutpoint_set.begin(); c!=cutpoint_set.end(); ++c ) {
			if ( *c > (*m)->elem_length() ) continue;
			MotifOP newmotif( new Motif( **m ) );
			newmotif->set_cutpoint( *c );
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
	for ( Abegos::const_iterator abego=abegolist.begin(); abego!=abegolist.end(); ++abego ) {
		if ( lenset.find( abego->size() + 2 ) == lenset.end() ) continue;
		std::string const secstruct( 2 + abego->size(), 'L' );
		std::stringstream abego_str;
		abego_str << 'X' << *abego << 'X';
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
IdealAbegoGenerator::set_extend_ss( bool const extend_ss )
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
	components::SegmentResid const length = static_cast< components::SegmentResid >( motif.elem_length() );
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
	components::SegmentResid const length = static_cast< components::SegmentResid >( motif.elem_length() );
	for ( components::SegmentResid resid=length; resid>=1; --resid ) {
		if ( motif.abego( resid ) != element_abego ) break;
		motif.set_ss( resid, sschar );
	}
}

} //protocols
} //denovo_design
} //components

