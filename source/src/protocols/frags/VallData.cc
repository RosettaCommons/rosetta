// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Rosetta Headers
#include <protocols/frags/VallData.hh>


#include <protocols/frags/TorsionFragment.hh>
#include <protocols/frags/heap.hh>

#include <basic/basic.hh> // tracer output
#include <basic/options/util.hh>
#include <basic/options/option.hh>

// utility headers
#include <ObjexxFCL/ObjexxFCL.hh>
#include <ObjexxFCL/string.functions.hh>
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
#include <basic/options/keys/loops.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.frags.VallData" );

namespace protocols {
namespace frags {

using core::Size;
using core::Real;

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


/// SEE ALSO THE DUPLICATE COPY BELOW WHICH INCLUDES FILTERING BY TARGET TORSION ANGLES
/// IE, IF YOU FIX A BUG HERE FIX IT DOWN THERE AS WELL (THANKS)

void
VallData::get_frags(
	Size const nfrags,
	std::string const & target_seq,
	std::string const & target_ss,
	Real const seq_weight,
	Real const ss_weight,
	bool const exclude_gly,
	bool const exclude_pro,
	bool const exclude_cis_peptides, // at non-pre-Pro positions
	utility::vector1< Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library,
	Real const bb_weight, // = 0.0
	std::string const & target_bb // = ""
) const
{
	Size const frag_size( target_seq.size() );
	runtime_assert( frag_size == target_ss.size() && ss_weight >= 0.0 && seq_weight >= 0.0 );

	bool const score_bb_matches( !target_bb.empty() );
	if ( score_bb_matches ) { runtime_assert( target_bb.size() == frag_size ); }

	// for randomizing the order of frags with identical score
	Real const min_nonzero_weight
		( ss_weight == 0.0 ? seq_weight : ( seq_weight == 0.0 ? ss_weight : std::min( ss_weight, seq_weight ) ) );

	// reset heaps
	ObjexxFCL::FArray1D_int heap( nfrags + 2 );
	ObjexxFCL::FArray1D_float coheap( nfrags + 2 );
	protocols::frags::heap_init( heap, coheap, nfrags );

	Size const my_size( size() );
	utility::vector1< bool > exclude_chain( chain_.back(), false );
	for ( core::Size it : homs_to_exclude ) {
		TR.Trace << "Excluding chain " << it << " from fragment picking" << std::endl;
		exclude_chain[ it ] = true;
	}

	for ( Size vall_pos=2; vall_pos <= my_size - frag_size; ++vall_pos ) { // skip 1st and last windows
		// score this position: fragment runs from vall_pos ... vall_pos + frag_size - 1

		if ( chain_[ vall_pos-1 ] != chain_[ vall_pos + frag_size ] ) continue; // frag contains a terminus
		if ( exclude_chain[ chain_[ vall_pos ] ] ) continue;

		bool bad_frag( false );

		Real score(0.0); // bigger is worse

		for ( Size k=0; k< frag_size; ++k ) {
			Real const phi  ( phi_      [ vall_pos+k ] );
			Real const psi  ( psi_      [ vall_pos+k ] );
			Real const omega( omega_    [ vall_pos+k ] );
			char const seq  ( sequence_ [ vall_pos+k ] );
			char const ss   ( secstruct_[ vall_pos+k ] );
			if ( ( std::fabs( phi ) < 0.001 ) ||
					( std::fabs( psi ) + std::fabs( omega ) < 0.001 ) ||
					( seq == 'G' && exclude_gly && seq != target_seq[k] ) ||
					( seq == 'P' && exclude_pro && seq != target_seq[k] ) ||
					( std::abs( omega ) < 90.0 && exclude_cis_peptides && ( k == frag_size-1 || target_seq[k+1] != 'P' ) ) ) {
				bad_frag = true;
				break;
			}
			if ( ss  != target_ss [k] ) score +=  ss_weight;
			if ( seq != target_seq[k] ) score += seq_weight;
			if ( score_bb_matches ) {
				if ( bigbin_[ vall_pos+k ] != target_bb[k] ) score += bb_weight;
			}
		}

		if ( bad_frag ) continue;

		// randomly reorder frags with the same score
		score += min_nonzero_weight * 0.1 * numeric::random::uniform();

		// insert into heap, with negative score! since bigger==better for heaps
		bool err;
		protocols::frags::heap_insert( heap, coheap, vall_pos, -score, err );
	}

	// now extract top nfrags matches, copy Sizeo frag array
	// fragments come out of the heap from worst to best
	Size exact_matches(0);
	Real worst_score(999), best_score(999);

	for ( Size nn= nfrags; nn >=1; --nn ) {
		bool err;
		int vall_pos;
		float score;// heaps use float!!!
		protocols::frags::heap_extract( heap, coheap, vall_pos, score, err);
		runtime_assert( !err );

		if ( score >= -0.1 * min_nonzero_weight ) ++exact_matches;

		if ( nn == nfrags ) worst_score = -score;
		else if ( nn == 1 ) best_score = -score;


		TorsionFragmentOP fragment( new TorsionFragment( frag_size, 3 /* #bb torsions */ ) );
		for ( Size k=0; k< frag_size; ++k ) {
			fragment->set_torsion  ( k+1, 1, phi_   [ vall_pos + k ] );
			fragment->set_torsion  ( k+1, 2, psi_   [ vall_pos + k ] );
			fragment->set_torsion  ( k+1, 3, omega_ [ vall_pos + k ] );
			fragment->set_secstruct( k+1, secstruct_[ vall_pos + k ] );
			fragment->set_sequence ( k+1, sequence_ [ vall_pos + k ] );
		}
		// goes at the beginning
		library.insert_fragment( fragment );
	} // nn


	TR.Trace << "ss-frags: " <<
		" target_ss: " << target_ss <<
		" target_seq: " << target_seq <<
		" exact_matches: " << exact_matches <<
		" best_score: " << best_score <<
		" worst_score: " << worst_score;
	if ( score_bb_matches ) TR.Trace << " target_bb: " << target_bb;
	TR.Trace << std::endl;

}

void
VallData::get_frags(
	Size const nfrags,
	std::string const & target_seq,
	utility::vector1< std::map< char, core::Real > > const & target_ss, // HEL
	Real const seq_weight,
	Real const ss_weight,
	bool const exclude_gly,
	bool const exclude_pro,
	bool const exclude_cis_peptides, // at non-pre-Pro positions
	utility::vector1< Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library,
	Real const bb_weight, // = 0.0
	std::string const & target_bb // = ""
) const
{
	Size const frag_size( target_seq.size() );
	runtime_assert( frag_size == target_ss.size() && ss_weight >= 0.0 && seq_weight >= 0.0 );

	std::string target_ss_string;
	for ( Size i=1; i<= frag_size; ++i ) {
		runtime_assert( target_ss[i].size() == 3 );
		runtime_assert( target_ss[i].find('H') != target_ss[i].end() );
		runtime_assert( target_ss[i].find('E') != target_ss[i].end() );
		runtime_assert( target_ss[i].find('L') != target_ss[i].end() );
		Real const H( target_ss[i].find('H')->second );
		Real const E( target_ss[i].find('E')->second );
		Real const L( target_ss[i].find('L')->second );
		runtime_assert( std::max( H, std::max( E, L ) ) <  1.001 ); // confirm between 0 and 1
		runtime_assert( std::min( H, std::min( E, L ) ) > -0.001 );
		if ( H > E && H > L ) {
			if ( H>0.666 ) target_ss_string += 'H';
			else target_ss_string += 'h';
		} else if ( E > H && E > L ) {
			if ( E>0.666 ) target_ss_string += 'E';
			else target_ss_string += 'e';
		} else {
			if ( L>0.666 ) target_ss_string += 'L';
			else target_ss_string += 'l';
		}
	}

	bool const score_bb_matches( !target_bb.empty() );
	if ( score_bb_matches ) { runtime_assert( target_bb.size() == frag_size ); }

	// for randomizing the order of frags with identical score
	// also for counting exact matches, which dont make sense here anyway
	Real min_nonzero_weight(1e6);
	if ( seq_weight > 1e-6 ) min_nonzero_weight = std::min( min_nonzero_weight, seq_weight );
	if (  bb_weight > 1e-6 ) min_nonzero_weight = std::min( min_nonzero_weight,  bb_weight );
	if (  ss_weight > 1e-6 ) min_nonzero_weight = std::min( min_nonzero_weight,  ss_weight * 0.01 );

	// Real const min_nonzero_weight
	//  ( ss_weight == 0.0 ? seq_weight : ( seq_weight == 0.0 ? ss_weight : std::min( ss_weight, seq_weight ) ) );

	// reset heaps
	ObjexxFCL::FArray1D_int heap( nfrags + 2 );
	ObjexxFCL::FArray1D_float coheap( nfrags + 2 );
	protocols::frags::heap_init( heap, coheap, nfrags );

	Size const my_size( size() );
	utility::vector1< bool > exclude_chain( chain_.back(), false );
	for ( core::Size it : homs_to_exclude ) {
		TR.Trace << "Excluding chain " << it << " from fragment picking" << std::endl;
		exclude_chain[ it ] = true;
	}

	for ( Size vall_pos=2; vall_pos <= my_size - frag_size; ++vall_pos ) { // skip 1st and last windows
		// score this position: fragment runs from vall_pos ... vall_pos + frag_size - 1

		if ( chain_[ vall_pos-1 ] != chain_[ vall_pos + frag_size ] ) continue; // frag contains a terminus
		if ( exclude_chain[ chain_[ vall_pos ] ] ) continue;

		bool bad_frag( false );

		Real score(0.0); // bigger is worse

		for ( Size k=0; k< frag_size; ++k ) {
			//Real const phi  ( phi_      [ vall_pos+k ] );
			//Real const psi  ( psi_      [ vall_pos+k ] );
			Real const omega( omega_    [ vall_pos+k ] );
			char const seq  ( sequence_ [ vall_pos+k ] );
			char const ss   ( secstruct_[ vall_pos+k ] );
			if ( ( seq == 'G' && exclude_gly && seq != target_seq[k] ) ||
					( seq == 'P' && exclude_pro && seq != target_seq[k] ) ||
					( std::abs( omega ) < 90.0 && exclude_cis_peptides && ( k == frag_size-1 || target_seq[k+1] != 'P' ) ) ) {
				bad_frag = true;
				break;
			}
			score +=  ss_weight * ( 1.0 - target_ss[k+1].find(ss)->second ); // bigger score is worse
			if ( seq != target_seq[k] ) score += seq_weight;
			if ( score_bb_matches ) {
				if ( bigbin_[ vall_pos+k ] != target_bb[k] ) score += bb_weight;
			}
		}

		if ( bad_frag ) continue;

		// randomly reorder frags with the same score
		score += min_nonzero_weight * 0.1 * numeric::random::uniform();

		// insert into heap, with negative score! since bigger==better for heaps
		bool err;
		protocols::frags::heap_insert( heap, coheap, vall_pos, -score, err );
	}

	// now extract top nfrags matches, copy Sizeo frag array
	// fragments come out of the heap from worst to best
	Size exact_matches(0);
	Real worst_score(999), best_score(999);

	for ( Size nn= nfrags; nn >=1; --nn ) {
		bool err;
		int vall_pos;
		float score;// heaps use float!!!
		protocols::frags::heap_extract( heap, coheap, vall_pos, score, err);
		runtime_assert( !err );

		if ( score >= -0.1 * min_nonzero_weight ) ++exact_matches;

		if ( nn == nfrags ) worst_score = -score;
		else if ( nn == 1 ) best_score = -score;


		TorsionFragmentOP fragment( new TorsionFragment( frag_size, 3 /* #bb torsions */ ) );
		for ( Size k=0; k< frag_size; ++k ) {
			fragment->set_torsion  ( k+1, 1, phi_   [ vall_pos + k ] );
			fragment->set_torsion  ( k+1, 2, psi_   [ vall_pos + k ] );
			fragment->set_torsion  ( k+1, 3, omega_ [ vall_pos + k ] );
			fragment->set_secstruct( k+1, secstruct_[ vall_pos + k ] );
			fragment->set_sequence ( k+1, sequence_ [ vall_pos + k ] );
		}
		// goes at the beginning
		library.insert_fragment( fragment );
	} // nn


	TR.Trace << "ss-frags: " <<
		" target_ss: " << target_ss_string <<
		" target_seq: " << target_seq <<
		" exact_matches: " << exact_matches <<
		" best_score: " << best_score <<
		" worst_score: " << worst_score;
	if ( score_bb_matches ) TR.Trace << " target_bb: " << target_bb;
	TR.Trace << std::endl;

}


Real
torsion_dev_score(
	Real const frag_angle,
	Real const target_angle,
	Real const min_dev,
	Real const max_dev
)
{
	if ( std::abs( target_angle )< 1e-3 ) return 0.0; // don't score terminal angles
	Real const dev( basic::subtract_degree_angles( frag_angle, target_angle ) );
	if ( dev < min_dev ) return 0.0;
	else if ( dev > max_dev ) return 1.0;
	else return ( dev - min_dev ) / ( max_dev - min_dev );
}


void
VallData::get_cheating_frags(
	Size const nfrags,
	std::string const & target_seq,
	std::string const & target_ss,
	utility::vector1< Real > const & target_phi,
	utility::vector1< Real > const & target_psi,
	utility::vector1< Real > const & target_omega,
	Real const seq_weight,
	Real const ss_weight,
	Real const torsion_weight,
	Real const min_torsion_dev,
	Real const max_torsion_dev,
	utility::vector1< Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library
) const
{
	Size const frag_size( target_seq.size() );
	bool const skip_exact_matches( frag_size > 5 ); // simple way to skip self frags (and a few others!)
	runtime_assert( frag_size == target_ss.size() && ss_weight >= 0.0 && seq_weight >= 0.0 );

	// for randomizing the order of frags with identical score
	Real min_nonzero_weight(1e6);
	if (      ss_weight > 1e-6 ) min_nonzero_weight = std::min( min_nonzero_weight,      ss_weight );
	if (     seq_weight > 1e-6 ) min_nonzero_weight = std::min( min_nonzero_weight,     seq_weight );
	if ( torsion_weight > 1e-6 ) min_nonzero_weight = std::min( min_nonzero_weight, torsion_weight );

	// reset heaps
	ObjexxFCL::FArray1D_int heap( nfrags + 2 );
	ObjexxFCL::FArray1D_float coheap( nfrags + 2 );
	protocols::frags::heap_init( heap, coheap, nfrags );

	Size const my_size( size() );
	utility::vector1< bool > exclude_chain( chain_.back(), false );
	for ( core::Size it : homs_to_exclude ) {
		TR.Trace << "Excluding chain " << it << " from fragment picking" << std::endl;
		exclude_chain[ it ] = true;
	}

	for ( Size vall_pos=2; vall_pos <= my_size - frag_size; ++vall_pos ) { // skip 1st and last windows
		// score this position: fragment runs from vall_pos ... vall_pos + frag_size - 1

		if ( chain_[ vall_pos-1 ] != chain_[ vall_pos + frag_size ] ) continue; // frag contains a terminus
		if ( exclude_chain[ chain_[ vall_pos ] ] ) continue;

		bool bad_frag( false );

		Real score(0.0); // bigger is worse
		Size n_seq_mismatches(0);
		for ( Size k=0; k< frag_size; ++k ) {
			Real const phi  ( phi_      [ vall_pos+k ] );
			Real const psi  ( psi_      [ vall_pos+k ] );
			Real const omega( omega_    [ vall_pos+k ] );
			char const seq  ( sequence_ [ vall_pos+k ] );
			char const ss   ( secstruct_[ vall_pos+k ] );
			if ( ss != target_ss[k] ) {
				score += ss_weight;
			}
			if ( seq != target_seq[k] ) {
				score += seq_weight;
				++n_seq_mismatches;
			}
			score += torsion_weight * ( torsion_dev_score( phi  , target_phi  [k+1], min_torsion_dev, max_torsion_dev ) +
				torsion_dev_score( psi  , target_psi  [k+1], min_torsion_dev, max_torsion_dev ) +
				torsion_dev_score( omega, target_omega[k+1], min_torsion_dev, max_torsion_dev ) );
		}

		if ( bad_frag ) continue;
		if ( skip_exact_matches && n_seq_mismatches==0 ) {
			TR.Trace << "skipping exact sequence matched fragment: " << vall_pos << ' ' << target_seq << std::endl;
			continue;
		}

		// randomly reorder frags with the same score
		score += min_nonzero_weight * 0.02 * numeric::random::uniform(); // must be smaller since torsion_dev_score not int

		// insert into heap, with negative score! since bigger==better for heaps
		bool err;
		protocols::frags::heap_insert( heap, coheap, vall_pos, -score, err );
	}

	// now extract top nfrags matches, copy Sizeo frag array
	// fragments come out of the heap from worst to best
	Size exact_matches(0);
	Real worst_score(999), best_score(999);

	for ( Size nn= nfrags; nn >=1; --nn ) {
		bool err;
		int vall_pos;
		float score;// heaps use float!!!
		protocols::frags::heap_extract( heap, coheap, vall_pos, score, err);
		runtime_assert( !err );

		if ( score >= -0.1 * min_nonzero_weight ) ++exact_matches;

		if ( nn == nfrags ) worst_score = -score;
		else if ( nn == 1 ) best_score = -score;


		TorsionFragmentOP fragment( new TorsionFragment( frag_size, 3 /* #bb torsions */ ) );
		for ( Size k=0; k< frag_size; ++k ) {
			fragment->set_torsion  ( k+1, 1, phi_   [ vall_pos + k ] );
			fragment->set_torsion  ( k+1, 2, psi_   [ vall_pos + k ] );
			fragment->set_torsion  ( k+1, 3, omega_ [ vall_pos + k ] );
			fragment->set_secstruct( k+1, secstruct_[ vall_pos + k ] );
			fragment->set_sequence ( k+1, sequence_ [ vall_pos + k ] );
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


/// @details  Handles loading default vall, then calls above routine
void
get_frags(
	Size const nfrags,
	std::string const & target_seq,
	std::string const & target_ss,
	Real const seq_weight,
	Real const ss_weight,
	bool const exclude_gly,
	bool const exclude_pro,
	bool const exclude_cis_peptides,
	utility::vector1< Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library,
	Real const bb_weight, // = 0.0
	std::string const & target_bb // = ""
)
{

	static VallData vall( basic::options::option[ basic::options::OptionKeys::loops::vall_file ] );

	vall.get_frags( nfrags,
		target_seq, target_ss,
		seq_weight, ss_weight,
		exclude_gly, exclude_pro, exclude_cis_peptides,
		homs_to_exclude,
		library,
		bb_weight, target_bb );
}

/// @details  Handles loading default vall, then calls above routine
void
get_frags(
	Size const nfrags,
	std::string const & target_seq,
	utility::vector1< std::map< char, core::Real > > const & target_ss, // HEL
	Real const seq_weight,
	Real const ss_weight,
	bool const exclude_gly,
	bool const exclude_pro,
	bool const exclude_cis_peptides,
	utility::vector1< Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library,
	Real const bb_weight, // = 0.0
	std::string const & target_bb // = ""
)
{

	static VallData vall( basic::options::option[ basic::options::OptionKeys::loops::vall_file ] );

	vall.get_frags( nfrags,
		target_seq, target_ss,
		seq_weight, ss_weight,
		exclude_gly, exclude_pro, exclude_cis_peptides,
		homs_to_exclude,
		library,
		bb_weight, target_bb );
}


void
get_cheating_frags(
	Size const nfrags,
	std::string const & target_seq,
	std::string const & target_ss,
	utility::vector1< Real > const & target_phi,
	utility::vector1< Real > const & target_psi,
	utility::vector1< Real > const & target_omega,
	Real const seq_weight,
	Real const ss_weight,
	Real const torsion_weight,
	Real const min_torsion_dev,
	Real const max_torsion_dev,
	utility::vector1< Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library
)
{

	static VallData vall( basic::options::option[ basic::options::OptionKeys::loops::vall_file ] );

	vall.get_cheating_frags( nfrags,
		target_seq, target_ss, target_phi, target_psi, target_omega,
		seq_weight, ss_weight, torsion_weight,
		min_torsion_dev, max_torsion_dev,
		homs_to_exclude,
		library );
}


void
dump_vall_fasta( std::string const & fasta_filename )
{
	VallData const vall( basic::options::option[ basic::options::OptionKeys::loops::vall_file ] );
	utility::vector1< char > const & vall_sequence( vall.sequence() );

	//  utility::vector1< Real > const & vall_phi  ( vall.phi  () );
	//  utility::vector1< Real > const & vall_psi  ( vall.psi  () );
	//  utility::vector1< Real > const & vall_omega( vall.omega() );

	utility::vector1< Size > const & vall_chain( vall.chain() );


	std::ofstream out( fasta_filename.c_str() );
	Size chain_begin( 1 );
	Size const vall_size( vall_sequence.size() );
	//if ( !is_lower_terminus( chain_begin ) || !is_upper_terminus( vall_size ) ) utility_exit_with_message("bad vall");
	//Size chain_counter(0);
	while ( chain_begin <= vall_size ) {
		//runtime_assert( is_lower_terminus( chain_begin ) );
		Size const chain( vall_chain[ chain_begin ] );
		runtime_assert( chain_begin == 1 || chain == vall_chain[ chain_begin-1 ] + 1 );

		// get the sequence
		std::string seq;
		Size pos( chain_begin );
		while ( true ) {
			seq.push_back( vall_sequence[ pos ] );
			if ( pos == vall_size || vall_chain[pos+1] != chain ) break;
			++pos;
		}
		runtime_assert( seq.size() == pos+1-chain_begin );

		// write the fasta
		std::string const id( "chain"+ObjexxFCL::lead_zero_string_of( chain, 6 ) );
		out << '>' << id << '\n' << seq << '\n';

		// move to the next chain
		chain_begin = pos+1;
	}
	out.close();
}

} // ns frags
} // ns protocols
