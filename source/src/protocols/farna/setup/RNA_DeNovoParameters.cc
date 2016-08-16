// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/setup/RNA_DeNovoParameters.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/setup/RNA_DeNovoParameters.hh>
#include <core/pose/rna/leontis_westhof_util.hh>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.setup.RNA_DeNovoParameters" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parameter files for FARFAR (following taken from http://rosie.rosettacommons.org/documentation/rna_denovo)
//
// You can specify the bounding Watson/Crick base pairs, strand boundaries, and more in a "params file", with lines like
//
// CUTPOINT_OPEN  <N1> [ <N2> ... ]                                 means that strands end after nucleotides N1, N2, etc.
//                                                                   Required if you have more than one strand!
//
// STEM    PAIR <N1> <M1> W W A  [  PAIR <N2> <M2>  W W A ... ]     means that residues N1 and M1 should form a
//                                                                   base pair with their Watson-Crick edges ('W') in an
//                                                                   antiparallel ('A') orientation; and N2 and M2, etc.. One
//                                                                   'STEM' line per contiguous helix. This will produce
//                                                                   constraints drawing the base-paired residues together,
//                                                                   and will also provide potential connection points between
//                                                                   strands for multi-strand motifs.
//
// OBLIGATE  PAIR <N> <M>  <E>  <F>  <O>                            means that a connection point between strands is forced between residues N and M
//                                                                  using their edges E and F [permitted values: W ('Watson-Crick'),
//                                                                  H ('Hoogsteen'), S ('sugar')] and orientation O [permitted
//                                                                  values: A (antiparallel) or P (parallel), based on normal
//                                                                  vectors on the two bases]. If modeling a single-strand motif,
//                                                                  forcing this 'obligate pair' will result in a 'temporary'
//                                                                  chainbreak, randomly placed in a non-stem residue.
//                                                                  Typically will not use OBLIGATE except for complex
//                                                                  topologies like pseudoknots.
//
// CUTPOINT_CLOSED <N>                                              Location of a 'temporary' chainbreak in strands. Typically
//                                                                  will not use except for complex topologies.
//
// CHAIN_CONNECTION SEGMENT1 <N1> <N2> SEGMENT2 <M1> <M2>            Used instead of obligate pair -- if we know two strands
//                                                                   are connected by a pair, but we don't know the residues, can ask
//                                                                   for some pairing to occur between strands N1-N2 and M1-M2,
//                                                                   sampled randomly. Typically will not use except for complex topologies.
//
// Recent update (2014):
//
// Allowing more flexible specification of residue sets which are connected by a base pair:
//
//  CHAIN_CONNECTION  SET1 <ints> SET2 <ints>
//
// where <ints> can be in the format like "5-7 90 92", etc.
//
//
// This file can be set up via rna_denovo_setup.py, available in tools/rna_tools/.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace ObjexxFCL::format;
using namespace ObjexxFCL;
using namespace core;
using namespace core::pose::rna;

namespace protocols {
namespace farna {
namespace setup {

//Constructor
RNA_DeNovoParameters::RNA_DeNovoParameters( std::string const & filename ):
	filename_( filename ),
	secstruct_defined_( false )
{
	if ( filename_.size() > 0 ) read_parameters_from_file( filename_ );
}

//Constructor
RNA_DeNovoParameters::RNA_DeNovoParameters():
	secstruct_defined_( false )
{
}

//Destructor
RNA_DeNovoParameters::~RNA_DeNovoParameters()
{}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoParameters::read_parameters_from_file( std::string const & filename ) {

	TR << "Reading RNA parameters file: " << filename << std::endl;

	Size a;

	utility::io::izstream pairing_stream( filename );
	if ( !pairing_stream ) {
		pairing_stream.close();
		utility_exit_with_message(  "Pairing file? " + filename );
	}

	std::string line, tag;

	while ( getline( pairing_stream, line ) ) {

		std::istringstream line_stream( line );
		line_stream >> tag;
		if ( line_stream.fail() ) continue; //Probably a blank line.

		if ( tag == "OBLIGATE" ) {

			line_stream >> tag;
			if ( line_stream.fail() || tag != "PAIR" )  utility_exit_with_message(  "Problem with OBLIGATE PAIR readin: " + tag );
			get_pairings_from_line( line_stream, true /*obligate jump*/ );

		} else if ( tag == "STEM"  || tag == "POSSIBLE" ) {

			line_stream >> tag;
			if ( line_stream.fail() || tag != "PAIR" )  utility_exit_with_message(  "Problem with STEM PAIR readin: " + tag );
			get_pairings_from_line( line_stream, false /*obligate jump*/ );

		} else if ( tag == "ALLOW_INSERT" ) {
			// deprecated!!! switch to ALLOW_INSERT_RES!!
			Size pos1, pos2;
			while ( !line_stream.fail() ) {
				line_stream >> pos1 >> pos2;
				runtime_assert( pos2 >= pos1 );
				for ( Size i = pos1; i <= pos2; i++ ) allow_insert_res_.push_back( i );
			}
			//utility_exit_with_message( "No longer reading in ALLOW_INSERT from command line. Try using -s <pdb> instead." );

		} else if ( tag == "ALLOW_INSERT_RES" ) {

			Size pos;
			while ( !line_stream.fail() ) {
				line_stream >> pos;
				allow_insert_res_.push_back( pos );
			}
			//utility_exit_with_message( "No longer reading in ALLOW_INSERT from command line. Try using -s <pdb> instead." );

		} else if ( tag == "VIRTUAL_ANCHOR" ) {

			while ( !line_stream.fail() ) {
				line_stream >> a;
				if ( !line_stream.fail()  ) {
					virtual_anchor_attachment_points_.push_back( a );
				}
			}
		} else if ( tag == "CUTPOINT_CLOSED" ) {
			while ( !line_stream.fail() ) {
				line_stream >> a;
				if ( !line_stream.fail()  ) {
					cutpoints_closed_.push_back( a );
				}
			}
		} else if ( tag == "CUTPOINT_OPEN" ) {
			while ( !line_stream.fail() ) {
				line_stream >> a;
				if ( !line_stream.fail()  ) {
					cutpoints_open_.push_back( a );
				}
			}
		} else if ( tag == "CHAIN_CONNECTION" ) {

			read_chain_connection( line_stream );

		} else if ( tag == "SECSTRUCT" ) {
			line_stream >> rna_secstruct_legacy_;
			secstruct_defined_ = true;

		} else if  ( tag[0] == '#' ) {
			continue;

		} else {
			utility_exit_with_message(   "Unrecognized tag in pairing file: " + tag );
		}
	}
	pairing_stream.close();

}


/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoParameters::get_pairings_from_line(
	std::istringstream & line_stream,
	bool const obligate_pair )
{

	Size a,b;
	BaseEdge e1( ANY_BASE_EDGE ),e2( ANY_BASE_EDGE );
	BaseDoubletOrientation o( ANY_BASE_DOUBLET_ORIENTATION );
	std::string tag;

	utility::vector1< Size > line_pairings;

	while ( !line_stream.fail() ) {

		line_stream >> a >> b ;

		if ( line_stream.fail() ) {
			std::cout << "Parse error!!!" << a << ' ' << b << std::endl;
		}

		line_stream >> tag;
		if ( line_stream.fail() || tag == "PAIR" ) {
			e1 = WATSON_CRICK;
			e2 = WATSON_CRICK;
			o  = ANTIPARALLEL;
		} else {
			e1 = get_edge_from_char( tag[0] );
			char e2_char, o_char;
			line_stream >> e2_char >> o_char;
			e2 = get_edge_from_char(e2_char );
			if ( o_char == 'C' || o_char == 'T' ) { // convert from leontis-westhof cis/trans to antiparallel/antiparallel
				o  = get_base_doublet_orientation_from_LW( e1, e2, get_LW_orientation_from_char( o_char ) );
			} else {
				o = get_orientation_from_char( o_char );
			}
			if ( line_stream.fail() )  utility_exit_with_message(  "Problem with PAIR readin: "+tag );
		}

		core::pose::rna::BasePair p;

		if ( a < b ) {
			p.set_res1( a );
			p.set_res2( b );
		} else {
			p.set_res1( b );
			p.set_res2( a );
		}

		if ( a == b ) {
			utility_exit_with_message(   "Can't base pair a residue with itself: "+I(3,a)+" "+I(3,b)  );
		}

		p.set_edge1( e1 );
		p.set_edge2( e2 );
		p.set_orientation( o );

		Size idx = check_in_pairing_sets( obligate_pairing_sets_, p );
		if ( idx == 0 ) {
			idx = check_in_pairing_sets( stem_pairing_sets_, p );
		}
		if ( idx == 0 ) {
			rna_pairing_list_.push_back( p );
			idx = rna_pairing_list_.size();
		}
		line_pairings.push_back( idx );

		line_stream >> tag;
		if ( !line_stream.fail() && tag != "PAIR" )  utility_exit_with_message(  "Problem with PAIR readin: " + tag );
	}

	if ( obligate_pair ) {
		obligate_pairing_sets_.push_back( line_pairings );
	} else {
		stem_pairing_sets_.push_back( line_pairings );
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoParameters::save_res_lists_to_chain_connections_and_clear( utility::vector1< Size > & res_list1,
	utility::vector1< Size > & res_list2 ) {
	if ( res_list1.size() > 0 || res_list2.size() > 0 ) {
		runtime_assert( res_list1.size() > 0 && res_list2.size() > 0 );
		chain_connections_.push_back( std::make_pair(res_list1, res_list2) );
		res_list1.clear();
		res_list2.clear();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoParameters::read_chain_connection( std::istringstream & line_stream ) {

	std::string tag;
	Size pos1( 0 ), pos2( 0 ), which_segment( 1 );

	bool legacy_segment_input( false );
	utility::vector1< Size > res_list1, res_list2;
	while ( !line_stream.fail() ) {

		line_stream >> tag;
		if ( line_stream.fail() ) break;

		if ( tag == "SEGMENT1" ) {
			which_segment = 1;
			legacy_segment_input = true;
			save_res_lists_to_chain_connections_and_clear( res_list1, res_list2 );
			continue;
		} else if ( tag == "SEGMENT2" ) {
			which_segment = 2;
			legacy_segment_input = true;
			continue;
		} else if ( tag == "SET1" ) {
			which_segment = 1;
			legacy_segment_input = false;
			save_res_lists_to_chain_connections_and_clear( res_list1, res_list2 );
			continue;
		} else if ( tag == "SET2" ) {
			which_segment = 2;
			legacy_segment_input = false;
			continue;
		}

		if ( legacy_segment_input ) {
			runtime_assert( is_int( tag ) );
			pos1 = int_of( tag );
			line_stream >> tag;
			runtime_assert( is_int( tag ) );
			pos2 = int_of( tag );
			runtime_assert( pos2 >= pos1 );
			for ( Size i = pos1; i <= pos2; i++ ) {
				if ( which_segment == 1 ) res_list1.push_back( i );
				if ( which_segment == 2 ) res_list2.push_back( i );
			}
		} else {
			bool string_is_ok( false );
			std::vector< int > ints = ints_of( tag, string_is_ok );
			runtime_assert( string_is_ok );
			for ( Size m = 0; m < ints.size(); m++ ) {
				if ( which_segment == 1 ) res_list1.push_back( ints[m] );
				if ( which_segment == 2 ) res_list2.push_back( ints[m] );
			}
		}

	}
	save_res_lists_to_chain_connections_and_clear( res_list1, res_list2 );

	if ( chain_connections_.size() == 0 ) {
		utility_exit_with_message(  "Did not specify SEGMENT1 or SEGMENT2 or SET1 or SET2 in CHAIN_CONNECTION line?" );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_DeNovoParameters::check_in_pairing_sets( utility::vector1 < utility::vector1 <core::Size > > pairing_sets,
	BasePair const & rna_pairing_check ) const {
	for ( Size n = 1; n <= pairing_sets.size(); n++ ) {
		for ( Size m = 1; m <= pairing_sets[n].size(); m++ ) {
			BasePair rna_pairing = rna_pairing_list_[ pairing_sets[n][m] ];
			if ( rna_pairing.res1() == rna_pairing_check.res1() && rna_pairing.res2() == rna_pairing_check.res2() ) return pairing_sets[n][m];
			if ( rna_pairing.res2() == rna_pairing_check.res1() && rna_pairing.res1() == rna_pairing_check.res2() ) return pairing_sets[n][m];
		}
	}
	return false;
}


} //setup
} //farna
} //protocols
