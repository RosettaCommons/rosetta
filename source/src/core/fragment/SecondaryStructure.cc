// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief secondary structure will hold statistics about secondary structure predictions
/// sources can be from
///      - fragments
///      - psipred files ? other stuff
///
/// @details
///  from converting conformation_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange


// Unit Headers
#include <core/fragment/SecondaryStructure.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIteratorWorker_.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>

#include <core/fragment/FragData.hh>
#include <core/fragment/FrameIterator.hh>
#include <utility/vector1.hh>


//// C++ headers
//#include <cstdlib>
//#include <string>
//#include <vector>
static thread_local basic::Tracer tr( "core.fragment" );

namespace core {
namespace fragment {

/// @details Auto-generated virtual destructor
SecondaryStructure::~SecondaryStructure() {}

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

void SecondaryStructure::compute_fractions( core::fragment::FragSet const& frags, bool bJustCenterResidue /*default false */ ) {
  using namespace core::fragment;

  //make sure, that the fragments that are used start with position 1 (local position)
  if ( frags.global_offset() != 0 ){
    // JRP This is super lame, and should be pretty easy to fix adjusting everything by frags.global_offset()
    tr.Error << "[ERROR] SecondaryStructure computations must be carried out with local coordinates (global offset of fragments must be 0)." << std::endl;
    runtime_assert( false );
  }

  Size frag_nres = frags.max_pos();
  if ( total_residue_ < frag_nres ) total_residue_ = frag_nres;
  tr.Info << "compute strand/loop fractions for " << total_residue_ << " residues... " << std::endl;
  if ( total_residue_ == 0 ) utility_exit_with_message( "no fragment to compute secondary structure ");

  strand_fraction_.dimension( total_residue_, 0.0 );
  loop_fraction_.dimension( total_residue_, 0.0 );

  //Note: Confidence is rarely used, and defaults to zero unless set otherwise -rv
  confidence_.dimension( total_residue_, 0.0 );

  FArray1D_int count( total_residue_, 0 ); //keep track how many entries for each residue

  strand_fraction_(1) = 0.0;
  strand_fraction_(total_residue_) = 0.0;

  for ( FragID_Iterator it=frags.begin(), eit=frags.end(); it!=eit; ++it ) { //carefully checked that I don't change FrameData
    Size central_residue( static_cast< Size > (it->frame().length()/2.0+0.5) );
    Size loop_start = bJustCenterResidue ? central_residue : 1;
    Size loop_end = bJustCenterResidue ? central_residue : it->frame().length();
    for ( Size fpos = loop_start; fpos <= loop_end; ++fpos ) {
      char const ss = it->fragment().secstruct( fpos );
      Size pos = it->frame().seqpos( fpos );
      if ( ss == 'E' || ss == 'L' || ss == 'H' ) {
        ++count( pos );
        if (  ss == 'E' ) strand_fraction_( pos ) += 1.0;
        else if ( ss == 'L' ) loop_fraction_( pos ) += 1.0;
      } else {
        tr.Warning << "found invalid secondary structure assignment in fragment data: " << ss << std::endl;
      }
    }
  }

  for ( Size pos = 1; pos <= total_residue_; pos++ ) {
    if ( count( pos ) ) {
      strand_fraction_(pos) /= count ( pos );
      loop_fraction_(pos) /= count ( pos );
    } else {
      loop_fraction_( pos ) = 1.0;
      strand_fraction_( pos ) = 0.0;
    }
  }
}

/// @brief returns regions (in loop-class format) that belong to contiguous pieces of ss-structure
//loops::Loops SecondaryStructure::compute_ss_regions( core::Real max_loop_frac, core::Size min_length ) const {
//	Size start( 0 );
//	Size last( 0 );
//	Size max_gap( 2 );
//	loops::Loops ss_regions;
//	for ( Size pos = 1; pos <= total_residue(); ++pos ) {
//		if ( loop_fraction( pos ) <= max_loop_frac ) {
//			if ( !start ) {
//				start = pos;
//				last = pos - 1;
//			}
//			if ( last + max_gap < pos ) {
//				if ( last - start >= min_length ) {
//					ss_regions.add_loop( start, last );
//				}
//				start=0;
//			}
//			last = pos;
//		}
//	}
//	return ss_regions;
//}


char SecondaryStructure::secstruct( Size pos ) const {
  core::Real max = loop_fraction( pos );
  char ss = 'L';
  if ( max < strand_fraction( pos ) ) {
    max = strand_fraction_( pos );
    ss = 'E';
  }
  if ( max < helix_fraction( pos ) ) {
    ss = 'H';
  }
  return ss;
}

/// @detail read from file
void SecondaryStructure::read_from_file( std::string fn ) {
  utility::io::izstream data( fn );
  if ( !data ) {
    tr.Fatal << "can't secondary structure file!!!" << fn << std::endl;
    data.close();
    utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
  }

  std::string line;
  getline( data, line); //ignore header
  std::istringstream line_stream( line );
  std::string dummy;

  //read number of residues
  line_stream >> dummy >> dummy >> dummy >> dummy >> total_residue_;

  //dimension arrays
  loop_fraction_.dimension( total_residue_, 0.0 );
  strand_fraction_.dimension( total_residue_, 0.0 );

  while ( getline( data, line ) ) {
    std::istringstream line_stream( line );
    // a=i, b=j, c=orientation(1 or 2), d=pleating(1 or 2)
    int pos;
    core::Real ef,hf,lf;
    line_stream >> pos >> ef >> hf >> lf;

    if ( line_stream.fail() ) {
      tr.Warning << "parse error: " << line << std::endl;
      continue;
    }

    loop_fraction_( pos ) = lf;
    strand_fraction_( pos ) = ef;
    if ( std::abs( helix_fraction( pos ) - hf ) > 0.01 ) {
      tr.Warning << "inconsistency in secondary structure file at position "
                 << pos << " H ( read ) = " << hf << " 1.0-L-E ( expected ) "
                 << helix_fraction( pos ) << std::endl;
    }
  }

  tr.flush();
}

SecondaryStructure::SecondaryStructure( core::pose::Pose const& pose ) {
  //  pose::set_ss_from_phi_psi( pose );
  total_residue_ = pose.total_residue();
  strand_fraction_.dimension( total_residue_, 0.0 );
  loop_fraction_.dimension( total_residue_, 0.0 );
  for ( Size pos = 1; pos<=pose.total_residue(); pos++ ) {
    char ss = pose.secstruct( pos );
    if ( ss == 'L' ) loop_fraction_( pos ) += 1.0;
    if ( ss == 'E' ) strand_fraction_( pos ) += 1.0;
  }
}

SecondaryStructure::SecondaryStructure( utility::vector1<SecondaryStructureOP> & components,utility::vector1<Real> & weights ) {

    runtime_assert(components.size() == weights.size());
    Size size = components[1]->total_residue_;
    for(Size i=2;i<=components.size();i++)
  runtime_assert(components[i]->total_residue_ == size);
    total_residue_ = size;
    strand_fraction_.dimension( size, 0.0 );
    loop_fraction_.dimension( size, 0.0 );
    for ( Size pos = 1; pos <= size; pos++ ) {
  Real pE = 0.0;
  Real pH = 0.0;
  Real pL = 0.0;
  for(Size i=1;i<=components.size();i++) {
      pE += components[i]->strand_fraction(pos) * weights[i];
      pH += components[i]->helix_fraction(pos) * weights[i];
      pL += components[i]->loop_fraction(pos) * weights[i];
  }
  set_fractions(pos,pH,pE,pL);
    }
}


/// @detail write to stream ( opposite from read_from_file )
void SecondaryStructure::show( std::ostream& os ) const {
  using namespace format;
  int const width( 10 );
  os << A( width, "pos") << A( width, "E" ) << A( width, "H" ) << A( width, "L" ) << I( width, 4, total_residue() ) << std::endl;
  for ( Size i = 1; i<= total_residue(); i++ ) {
    os << I( width, 4, i)
       << F( width, 4, strand_fraction( i ) )
       << F( width, 4, helix_fraction( i ) )
       << F( width, 4, loop_fraction( i ) )
       << std::endl;
  }
}


void SecondaryStructure::read_psipred_ss2( std::string filename ) {
  utility::io::izstream data( filename );
  if ( !data ) {
    tr.Fatal << "can't secondary structure file!!!" << filename << std::endl;
    data.close();
    utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
  }
  read_psipred_ss2( data );
}

void SecondaryStructure::read_psipred_ss2( std::istream& data ) {

  std::string line;
  getline( data, line); //ignore header

  //dimension arrays
  Size total_reserved( 500 );
  loop_fraction_.dimension( total_reserved, 0.0 );
  strand_fraction_.dimension( total_reserved, 0.0 );
  Size last_pos(0);
  while ( getline( data, line ) ) {
    if ( line.size() == 0 ) continue;
    if ( line[ 0 ] == '#' ) continue;
    std::istringstream line_stream( line );
    // a=i, b=j, c=orientation(1 or 2), d=pleating(1 or 2)
    Size pos;
    std::string aa, secstruct_letter;
    core::Real ef,hf,lf;
    line_stream >> pos >> aa >> secstruct_letter >> lf >> hf >> ef;

    //fill up missing residues until pos with 0.33 0.33 0.33 probabilities:
    if ( line_stream.fail() ) {
      tr.Warning << "parse error: " << line << std::endl;
      continue;
    }
    // renormalize so that probabilities sum to 1.0
    Real const total( lf + hf + ef );
    lf /= total;
    hf /= total;
    ef /= total;

    if ( pos > total_reserved ) {
      total_reserved+=400;
      loop_fraction_.redimension( total_reserved );
      strand_fraction_.redimension( total_reserved );
    }
    for ( last_pos = last_pos+1; last_pos<pos; last_pos++ ) {
      loop_fraction_( last_pos ) = 1.0/3;
      strand_fraction_( last_pos ) = 1.0/3;
    }

    if ( total_residue_ < pos ) total_residue_ = pos;
    loop_fraction_( pos ) = lf;
    strand_fraction_( pos ) = ef;

    if ( std::abs( helix_fraction( pos ) - hf ) > 0.01 ) {
      tr.Warning << "inconsistency in secondary structure file at position "
                 << pos << " H ( read ) = " << hf << " 1.0-L-E ( expected ) "
                 << helix_fraction( pos ) << std::endl;
    }
    last_pos = pos;
  }

  tr.flush();
}

void SecondaryStructure::write_psipred_ss2(
  std::ostream & os,
  std::string const & sequence
) const {
  int const width( 10 );
  os << "# PSIPRED VFORMAT (PSIPRED V2.5 by David Jones)\n\n";
  for ( Size i = 1; i<= total_residue(); i++ ) {
    char ss = secstruct( i );
    if ( ss == 'L' ) ss = 'C'; //for psipred
    os << I( width, 4, i);
    os << A( 2, ( sequence.size() >= i ) ? sequence[ i - 1 ] : 'X' )
       << A( 2, ss )
       << F( width, 3, loop_fraction( i ) )
       << F( width, 3, helix_fraction( i ) )
       << F( width, 3, strand_fraction( i ) )
       << std::endl;
  }
}

void SecondaryStructure::read_talos_ss( std::string filename ) {
  utility::io::izstream data( filename );
  if ( !data ) {
    tr.Fatal << "can't secondary structure file!!!" << filename << std::endl;
    data.close();
    utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
  }
  read_talos_ss( data );
}

void SecondaryStructure::read_talos_ss( std::istream& data ) {

  std::string line;
  bool header_done(false);

  while (!header_done && !data.eof()) {
    getline( data, line); //ignore header
    std::istringstream line_stream( line );
    std::string keyword;
    line_stream >> keyword;

    if (keyword == "FORMAT") {
      header_done = true;
    }
  }

  //dimension arrays
  Size total_reserved( 500 );
  loop_fraction_.dimension( total_reserved, 0.0 );
  strand_fraction_.dimension( total_reserved, 0.0 );
  confidence_.dimension( total_reserved, 0.0 );
  Size last_pos( 0 );
  while ( getline( data, line ) ) {
    if ( line.size() == 0 ) continue;
    if ( line[ 0 ] == '#' ) continue;
    std::istringstream line_stream( line );
    // a=i, b=j, c=orientation(1 or 2), d=pleating(1 or 2)
    Size pos;
    std::string aa, junk, secstruct_letter;
    core::Real ef,hf,lf, confidence;

    //confidence is currently being thrown out, may want to bear that in mind
    line_stream >> pos >> aa >> junk >> junk >> hf >> ef >> lf >> confidence >> secstruct_letter;

    if ( line_stream.fail() ) {
      tr.Warning << "parse error: " << line << std::endl;
      continue;
    }
    // renormalize so that probabilities sum to 1.0
    Real const total( lf + hf + ef );
    lf /= total;
    hf /= total;
    ef /= total;

    if ( pos > total_reserved ) {
      total_reserved+=400;
      loop_fraction_.redimension( total_reserved );
      strand_fraction_.redimension( total_reserved );
      confidence_.redimension( total_reserved );
    }
    if ( total_residue_ < pos ) total_residue_ = pos;
    loop_fraction_( pos ) = lf;
    strand_fraction_( pos ) = ef;
    confidence_( pos ) = confidence;

    if ( std::abs( helix_fraction( pos ) - hf ) > 0.01 ) {
      tr.Warning << "inconsistency in secondary structure file at position "
                 << pos << " H ( read ) = " << hf << " 1.0-L-E ( expected ) "
                 << helix_fraction( pos ) << std::endl;
    }
    for ( last_pos = last_pos+1; last_pos<pos; last_pos++ ) {
      loop_fraction_( last_pos ) = 1.0/3;
      strand_fraction_( last_pos ) = 1.0/3;
    }
  }

  tr.flush();
}

void
SecondaryStructure::extend( core::Size nres ) {
  if ( nres > total_residue_ ) {
    loop_fraction_.dimension( nres );
    strand_fraction_.dimension( nres );
    confidence_.dimension( nres );
    for ( Size pos=total_residue_+1; pos<=nres; ++pos ) {
      loop_fraction_( pos ) = 1.0;
      strand_fraction_( pos ) = 0.0;
      confidence_(pos) = 0.0;
    }
    total_residue_=nres;
  }
}


} // core
} // fragment
