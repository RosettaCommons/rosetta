// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Rosetta Headers
#include <protocols/frags/VallData.hh>


#include <protocols/frags/TorsionFragment.hh>
#include <protocols/frags/heap.hh>


// utility headers
// AUTO-REMOVED #include <ObjexxFCL/ObjexxFCL.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <utility/io/izstream.hh>

#include <numeric/random/random.fwd.hh>
#include <utility/exit.hh>

// C++ Headers
//#include <cmath>
//#include <cstdlib>
#include <iostream>
//#include <fstream>
//#include <sstream>
#include <basic/Tracer.hh> // tracer output

#include <utility/vector1.hh>


//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


static basic::Tracer TR( "protocols.frags.VallData" );

namespace protocols {
namespace frags {

///////////////////////////////////////////////////////////////////////////////
///\brief read from vall database file "filename"
///
///support both old format and new format, i.e., with "trimmed" in filename
void
VallData::read_file(
	std::string const & filename
)
{

	utility::io::izstream data ( filename );
	if ( !data ) {
		utility_exit_with_message( "cant open vall file: " + filename );
	}

	//////////////////////////////////////////////
	// this file parsing should be very fast since
	// the vall is enormous
	bool const new_format( filename.find( "trimmed" ) != std::string::npos );

	if ( new_format ) {
		char line[250];
		float phi,psi,omega; // should match types in sscanf conversion
		char seq,ss;
		while ( data ) {
			data.getline( line, 250 );
			if ( data.eof() ) break;

			std::sscanf( line   , "%1c", &seq);
			std::sscanf( line+ 1, "%1c", &ss);


			std::sscanf( line+ 2, "%9f", &phi);
			std::sscanf( line+11, "%9f", &psi);
			std::sscanf( line+20, "%9f", &omega);

			add_line( seq, ss, phi, psi, omega );
		}
	} else {
		char line[250];
		float phi,psi,omega;
		char seq,ss;
		while ( data ) {
			data.getline( line, 250 );
			if ( data.eof() ) break;

			std::sscanf( line+6 , "%1c", &seq);
			std::sscanf( line+8 , "%1c", &ss);

			std::sscanf( line+52, "%9f", &phi);
			std::sscanf( line+61, "%9f", &psi);
			std::sscanf( line+70, "%9f", &omega);

			add_line( seq, ss, phi, psi, omega );
		}
	}
	data.close();

	// remove excess capacity
	shrink();

}

///////////////////////////////////////////////////////////////////////////////
///\brief pick fragments in a single window
///
///\li nfrags -- number of fragment pieces to be picked
///\li target_seq --  query sequence and this also defines fragment length,
///    such as 3mer or 9mer
///\li target_ss -- secondary structure of query sequence.
///\li seq_weight and ss_weight -- weight for sequence match and secstruct match
///\li exclude_* -- exclude such squence pattern when picking fragments
///\li library -- output as a collection of fragments for that residue window
///
///\detail scan through vall database, for each position within each window,
///if sequence does not match target sequence, penalize by seq_weight; if secstruct
///does not match target secstruct, penalize by ss_weight. In the end, picking top "nfrags"
///fragments with lowest penalty score (this internally uses a "heap" implementation).
///
void
VallData::get_frags(
	Size const nfrags,
	std::string const & target_seq,
	std::string const & target_ss,
	Real const seq_weight,
	Real const ss_weight,
	bool const exclude_gly,
	bool const exclude_pro,
	bool const exclude_cys_peptides,
	SingleResidueTorsionFragmentLibrary & library
) const
{
	Size const frag_size( target_seq.size() );
	runtime_assert( frag_size == target_ss.size() && ss_weight >= 0.0 && seq_weight >= 0.0 );

	// for randomizing the order of frags with identical score
	Real const min_nonzero_weight
		( ss_weight == 0.0 ? seq_weight : ( seq_weight == 0.0 ? ss_weight : std::min( ss_weight, seq_weight ) ) );

	// reset heaps
	FArray1D_int heap( nfrags + 2 );
	FArray1D_float coheap( nfrags + 2 );
	heap_init( heap, coheap, nfrags );

	Size const my_size( size() );

	for ( Size vall_pos=1; vall_pos <= my_size - frag_size + 1; ++vall_pos ) {
		// score this position

		bool bad_frag( false );

		Real score(0.0); // bigger is worse

		for ( Size k=0; k< frag_size; ++k ) {
			Real const phi  ( phi_      [ vall_pos+k ] );
			Real const psi  ( psi_      [ vall_pos+k ] );
			Real const omega( omega_    [ vall_pos+k ] );
			char const seq  ( sequence_ [ vall_pos+k ] );
			char const ss   ( secstruct_[ vall_pos+k ] );
			if ( ( std::abs( phi ) < 0.01 ) ||
					 ( std::abs( psi ) + std::abs( omega ) < 0.01 ) ||
					 ( seq == 'G' && exclude_gly && seq != target_seq[k] ) ||
					 ( seq == 'P' && exclude_pro && seq != target_seq[k] ) ||
					 ( std::abs( omega ) < 90.0 && exclude_cys_peptides ) ) {
				bad_frag = true;
				break;
			}
			if ( ss != target_ss[k] ) {
				score += ss_weight;
			}
			if ( seq != target_seq[k] ) {
				score += seq_weight;
			}
		}

		if ( bad_frag ) continue;

		// randomly reorder frags with the same score
		score += min_nonzero_weight * 0.1 * numeric::random::uniform();

		// insert into heap, with negative score! since bigger==better for heaps
		bool err;
		heap_insert( heap, coheap, vall_pos, -score, err );
	}

	// now extract top nfrags matches, copy Sizeo frag array
	// fragments come out of the heap from worst to best
	Size exact_matches(0);
	Real worst_score(999), best_score(999);

	for ( Size nn= nfrags; nn >=1; --nn ) {
		bool err;
		int vall_pos;
		float score;// heaps use float!!!
		heap_extract( heap, coheap, vall_pos, score, err);
		runtime_assert( !err );

		if ( score >= -0.1 * min_nonzero_weight ) ++exact_matches;

		if ( nn == nfrags ) worst_score = -score;
		else if ( nn == 1 ) best_score = -score;

		//
		TorsionFragmentOP fragment( new TorsionFragment( frag_size, 3 /* #bb torsions */ ) );
		for ( Size k=0; k< frag_size; ++k ) {
			fragment->set_torsion  ( k+1, 1, phi_   [ vall_pos + k ] );
			fragment->set_torsion  ( k+1, 2, psi_   [ vall_pos + k ] );
			fragment->set_torsion  ( k+1, 3, omega_ [ vall_pos + k ] );
			fragment->set_secstruct( k+1, secstruct_[ vall_pos + k ] );
		}
		// goes at the beginning
		library.insert_fragment( fragment );
	} // nn


	TR.Trace << "ss-frags: " <<
		" target_ss: " << target_ss <<
		" target_seq: " << target_seq <<
		" exact_matches: " << exact_matches <<
		" best_score: " << best_score <<
		" worst_score: " << worst_score << std::endl;

}





} // ns frags
} // ns protocols
