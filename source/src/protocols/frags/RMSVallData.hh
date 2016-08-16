// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_frags_RMSVallData_hh
#define INCLUDED_protocols_frags_RMSVallData_hh


// Rosetta Headers

#include <core/types.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>

#include <protocols/frags/heap.hh>

// ObjexxFCL Headers
#include <utility/vector1.hh>

// utility headers
#include <ObjexxFCL/FArray1A.hh>
#include <utility/io/izstream.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.fwd.hh>
#include <utility/exit.hh>

// C++ Headers
#include <string>

//Auto using namespaces
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


namespace protocols {
namespace frags {

using core::Real;
using core::Size;

class RMSVallData {
public:
	/// default constructor
	RMSVallData ()
	{
		Size const big_size( 100000 ); // rough guess
		sequence_.reserve( big_size );
		secstruct_.reserve( big_size );
		phi_.reserve( big_size );
		psi_.reserve( big_size );
		omega_.reserve( big_size );
		exclude_gly = false;
		exclude_pro = false;
		exclude_cys_peptides = false;
	}

	/// constructor from input vall database file
	RMSVallData ( std::string const & filename )
	{
		Size const big_size( 100000 ); // rough guess
		// prevent lots of redimensioning as we read file? does this even matter?
		sequence_.reserve( big_size );
		secstruct_.reserve( big_size );
		phi_.reserve( big_size );
		psi_.reserve( big_size );
		omega_.reserve( big_size );
		exclude_gly = false;
		exclude_pro = false;
		exclude_cys_peptides = false;

		read_file( filename );
	}

	/// removes excess storage capacity to minimize memory usage
	void
	shrink()
	{
		sequence_.shrink();
		secstruct_.shrink();
		phi_.shrink();
		psi_.shrink();
		omega_.shrink();
	}

	// read from vall database file "filename"
	void
	read_file( std::string const & filename ) {
		utility::io::izstream data ( filename );
		if ( !data ) {
			utility_exit_with_message( "cant open vall file: " + filename );
		}

		//////////////////////////////////////////////
		// this file parsing should be very fast since
		// the vall is enormous
		bool const new_format( filename.find( "trimmed" ) != std::string::npos );

		if ( new_format ) {
			std::cerr << "Unsupported format for RMSVallData!!!\n";
			exit(1);
		} else {
			char line[250];
			float phi,psi,omega, x,y,z;
			char seq,ss;
			while ( data ) {
				data.getline( line, 250 );
				if ( data.eof() ) break;

				std::sscanf( line+6 , "%1c", &seq);
				std::sscanf( line+8 , "%1c", &ss);

				std::sscanf( line+25 , "%9f", &x);
				std::sscanf( line+34 , "%9f", &y);
				std::sscanf( line+43 , "%9f", &z);

				std::sscanf( line+52, "%9f", &phi);
				std::sscanf( line+61, "%9f", &psi);
				std::sscanf( line+70, "%9f", &omega);

				add_line( seq, ss, x, y, z, phi, psi, omega );
			}
		}
		data.close();

		// remove excess capacity
		shrink();

	}

	/// read in one more line from Vall input file
	void
	add_line( const char sq, const char ss,
		const Real x,  const Real y,  const Real z,
		const Real ph, const Real ps, const Real om ) {
		sequence_.push_back( sq );
		secstruct_.push_back( ss );

		X_.push_back( numeric::xyzVector<Real>( x, y, z ) );

		phi_.push_back( ph );
		psi_.push_back( ps );
		omega_.push_back( om );
	}

	utility::vector1< char > const & sequence () const { return sequence_;  }
	utility::vector1< char > const & secstruct() const { return secstruct_; }

	utility::vector1< numeric::xyzVector< Real > > const & X() const {return X_;}

	utility::vector1< Real > const & phi  () const {return phi_;}
	utility::vector1< Real > const & psi  () const {return psi_;}
	utility::vector1< Real > const & omega() const {return omega_;}

	/// number of lines in Vall database
	int size() const { return sequence_.size(); }

	// pick fragments for a single residue position from vall database
	void
	get_frags(
		Size const nfrags,
		utility::vector1< numeric::xyzVector< core::Real > > const & templ,  // CA coords we wish to align to
		std::string const &pref_seq,
		char const force_ss,
		core::fragment::FrameOP & frame,
		core::Real randomness = 0.0,
		core::Real oversample = 5.0
	) const {
		Size const frag_size( templ.size() );
		assert( frag_size == pref_seq.length() );

		// reset heaps
		Size const my_size( size() );
		Size bucket1_size( (Size)std::ceil(oversample*nfrags) );
		if ( oversample<=0 ) { bucket1_size =  my_size - frag_size + 1; }
		FArray1D_int heap( bucket1_size + 2 );
		FArray1D_float coheap( bucket1_size + 2 );
		protocols::frags::heap_init( heap, coheap, bucket1_size );

		// center template
		ObjexxFCL::FArray2D< core::Real > tmpl_pos( 3,frag_size ), tgt_pos( 3,frag_size );
		numeric::xyzVector< core::Real > tmpl_com(0,0,0);
		for ( int i = 1; i <= (int)frag_size; ++i ) {
			numeric::xyzVector< core::Real > const &x_i = templ[i];
			tmpl_com += x_i;
			for ( int k = 0; k < 3; ++k ) tmpl_pos(k+1,i) = x_i[k];
		}
		tmpl_com /= frag_size;
		for ( int i=1; i<=(int)frag_size; ++i ) {
			for ( int k = 0; k < 3; ++k ) tmpl_pos(k+1,i) -= tmpl_com[k];
		}

		// set up other tmp store data
		ObjexxFCL::FArray1D< numeric::Real > ww( frag_size, 1.0 );
		//ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
		//numeric::Real ctx;

		for ( Size vall_pos=1; vall_pos <= my_size - frag_size + 1; ++vall_pos ) {
			// score this position

			bool bad_frag( false );

			Real score(0.0); // bigger is worse
			numeric::xyzVector< core::Real > tgt_com(0,0,0);
			int seq_score = 0;

			for ( Size k=0; k< frag_size; ++k ) {
				Real const phi  ( phi_      [ vall_pos+k ] );
				Real const psi  ( psi_      [ vall_pos+k ] );
				Real const omega( omega_    [ vall_pos+k ] );
				char const seq  ( sequence_ [ vall_pos+k ] );
				char const ss   ( secstruct_[ vall_pos+k ] );
				if ( ( std::abs( phi ) < 0.01 ) ||
						( std::abs( psi ) + std::abs( omega ) < 0.01 ) ||
						( seq == 'G' && exclude_gly ) ||
						( seq == 'P' && exclude_pro ) ||
						( std::abs( omega ) < 90.0 && exclude_cys_peptides ) ||
						( ss != force_ss && (force_ss == 'H' || force_ss == 'E' || force_ss == 'L') ) ) {
					bad_frag = true;
					break;
				}

				if ( seq != pref_seq[k] ) seq_score++;

				numeric::xyzVector< core::Real > const &xx( X_[ vall_pos+k ] );
				tgt_com += xx;
				for ( int j = 0; j < 3; ++j ) tgt_pos(j+1,k+1) = xx[j];
			}

			if ( bad_frag ) continue;

			// center tgt
			tgt_com /= frag_size;
			for ( int i=1; i<=(int)frag_size; ++i ) {
				for ( int j = 0; j < 3; ++j ) tgt_pos(j+1,i) -= tgt_com[j];
			}

			// score == rms (?)
			//float rms_out=999.0;
			//rms_out = numeric::model_quality::rms_wrapper( frag_size, tgt_pos, tmpl_pos);  set but never used ~Labonte

			// sort on seq score (?)
			score = seq_score;

			// insert into heap, with negative score! since bigger==better for heaps
			bool err;
			protocols::frags::heap_insert( heap, coheap, vall_pos, -score, err );  // ?? out-of-date
		}

		// from the top 5*N in terms of seq score, chose best N from these using RMS
		FArray1D_int rmsheap( nfrags + 2 );
		FArray1D_float rmscoheap( nfrags + 2 );
		protocols::frags::heap_init( rmsheap, rmscoheap, nfrags );  // ?? out-of-date

		Size exact_matches(0);
		Real worst_score(999), best_score(999);

		for ( Size nn = bucket1_size; nn >=1; --nn ) {
			bool err;
			int vall_pos;
			float score;// heaps use float!!!

			protocols::frags::heap_extract( heap, coheap, vall_pos, score, err);  // ?? out-of-date
			assert( !err );

			if ( score == 0 ) ++exact_matches;

			// grab XYZs
			numeric::xyzVector< core::Real > tgt_com(0,0,0);
			for ( Size k=0; k< frag_size; ++k ) {
				numeric::xyzVector< core::Real > const &xx( X_[ vall_pos+k ] );
				tgt_com += xx;
				for ( int j = 0; j < 3; ++j ) tgt_pos(j+1,k+1) = xx[j];
			}

			// center tgt
			tgt_com /= frag_size;
			for ( int i=1; i<=(int)frag_size; ++i ) {
				for ( int j = 0; j < 3; ++j ) tgt_pos(j+1,i) -= tgt_com[j];
			}

			// score == rms (?)
			float rms_out = numeric::model_quality::rms_wrapper( frag_size, tgt_pos, tmpl_pos);

			// make things (a bit) randomized
			rms_out += randomness * numeric::random::uniform();

			// place in heap B
			protocols::frags::heap_insert( rmsheap, rmscoheap, vall_pos, -rms_out, err );  // ?? out-of-date
		}


		for ( Size nn = nfrags; nn >=1; --nn ) {
			bool err;
			int vall_pos;
			float score;// heaps use float!!!

			protocols::frags::heap_extract( rmsheap, rmscoheap, vall_pos, score, err);  // ?? out-of-date
			assert( !err );

			if ( nn == nfrags ) worst_score = -score;
			else if ( nn == 1 ) best_score = -score;

			core::fragment::FragDataOP current_fragment( new core::fragment::FragData );
			for ( Size k=0; k< frag_size; ++k ) {
				core::fragment::BBTorsionSRFDOP res_torsions( new core::fragment::BBTorsionSRFD( 3 ,secstruct_[ vall_pos + k ], sequence_ [ vall_pos+k ] ) ); // 3 protein torsions
				res_torsions->set_torsion   ( 1, phi_   [ vall_pos + k ]  ); // ugly numbers 1-3, but pose.set_phi also uses explicit numbers
				res_torsions->set_torsion   ( 2, psi_   [ vall_pos + k ]  );
				res_torsions->set_torsion   ( 3, omega_ [ vall_pos + k ]  );
				res_torsions->set_secstruct ( secstruct_[ vall_pos + k ]  );
				current_fragment->add_residue( res_torsions );
			}
			// we must manually mark the fragment as valid
			// why???? this is dumb
			if ( current_fragment->size() == frag_size ) current_fragment->set_valid();

			// goes at the beginning
			if ( !frame->add_fragment( current_fragment ) ) {
				std::cerr << "Incompatible fragment!" << std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}
		} // nn


		std::cerr << "rms-frags: " <<
			nfrags << " of " << bucket1_size <<
			" exact_matches: " << exact_matches <<
			" best_score: "    << best_score <<
			" worst_score: "   << worst_score << std::endl;

	}

private:
	utility::vector1< char > sequence_;
	utility::vector1< char > secstruct_;

	utility::vector1< numeric::xyzVector< Real > > X_;
	utility::vector1< Real > phi_;
	utility::vector1< Real > psi_;
	utility::vector1< Real > omega_;

	bool exclude_gly;
	bool exclude_pro;
	bool exclude_cys_peptides;
};

} // ns frags
} // ns protocols

#endif
