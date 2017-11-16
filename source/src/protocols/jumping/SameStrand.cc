// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief secondary structure will hold statistics about secondary structure predictions
/// sources can be from
///      - fragments
///      - psipred files ? other stuff
///
/// @details
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/jumping/SameStrand.hh>

// Package Headers
#include <core/fragment/SecondaryStructure.hh>

// Project Headers
#include <core/types.hh>


// Utility headers
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/StaticIndexRange.hh>


// numeric headers
#include <numeric/random/random.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


//// C++ headers
//#include <cstdlib>
//#include <string>
//#include <vector>
static basic::Tracer tr( "protocols.jumping" );

namespace protocols {
namespace jumping {

using namespace ObjexxFCL;

SameStrand::SameStrand( core::fragment::SecondaryStructureOP ss ) :
	secondary_structure_( ss )
{
	runtime_assert( ss != nullptr );

	// set total_residue
	total_residue_ = ss->total_residue();

	compute( *ss );
}

//SameStrand::SameStrand() : total_residue_( 0 )
//{}

// copy c'stor --
SameStrand::SameStrand( SameStrand const& other )
: ReferenceCount( other )
{
	total_residue_ = other.total_residue_;
	secondary_structure_ = other.secondary_structure_;
	same_strand_ = other.same_strand_;
	strand_sum_ = other.strand_sum_;
}

// nothing to be done in d'stor
SameStrand::~SameStrand() = default;


bool SameStrand::eval( Size i, Size j ) const {
	runtime_assert ( i >= 1 && i <= total_residue_ );
	runtime_assert ( j >= 1 && j <= total_residue_ );
	return same_strand_( i, j );
}

void
SameStrand::redo() const {
	runtime_assert( secondary_structure_ != nullptr );

	do_same_strand( );
	// compute( *secondary_structure_ );
}

void
SameStrand::compute( core::fragment::SecondaryStructure const& ss ) const {
	// do strand_sum
	do_strand_sum( ss );

	// and use that to do
	do_same_strand( );
}

void
SameStrand::do_same_strand( ) const {
	// do same_strand
	same_strand_.dimension( total_residue_, total_residue_ );

	float strand1,strand2,strand_ceiling,k_strand;
	//FArray1D_float strand_sum_( StaticIndexRange( 0, total_residue_ ) );
	float loop, r;

	for ( Size pos1 = 1; pos1 <= total_residue_; ++pos1 ) {
		for ( Size pos2 = 1; pos2 <= total_residue_; ++pos2 ) {
			strand1 = strand_sum_(pos1)-strand_sum_(pos1-1);
			strand2 = strand_sum_(pos2)-strand_sum_(pos2-1);
			strand_ceiling = std::max(0.2f, std::min(strand1,strand2));
			//   std::cout << "same_strand: "<< SS(pos1) << SS(pos2) << SS(strand_ceiling) << std::endl;
			int const i = std::min(pos1,pos2);
			int const j = std::max(pos1,pos2);

			if ( j-i > 5 ) {
				same_strand_(pos1,pos2) = false;
			} else if ( j-i < 2 ) {
				same_strand_(pos1,pos2) = true;
			} else {
				same_strand_(pos1,pos2) = true;
				for ( int k = i+1, ke = j-1; k <= ke; ++k ) {
					k_strand = ( strand_sum_(k) - strand_sum_(k-1) ) / strand_ceiling;
					loop = 1.0 - k_strand;
					if ( loop < 0.3 ) loop = 0.0;
					r = numeric::random::rg().uniform();
					if ( r < loop ) {
						same_strand_(pos1,pos2) = false;
					}
					//      std::cout << "loop cut:" << SS( i ) << SS( k ) << SS( j ) <<
					//       SS( loop ) << std::endl;

				}  // for k
			}    // if sep>5
		}      // for pos2
	}        // for pos1
}

void
SameStrand::do_strand_sum( core::fragment::SecondaryStructure const& ss ) const {
	runtime_assert( total_residue_ == ss.total_residue() );
	strand_sum_.dimension( StaticIndexRange( 0, total_residue_ ) );
	strand_sum_(0) = 0.0;
	for ( Size i = 1; i <= total_residue_; ++i ) {
		strand_sum_(i) = strand_sum_(i-1) + ss.strand_fraction(i);
		//  std::cout << "strandsum: " << SS(i) << SS(strand_sum_(i) ) << std::endl;
	}
}

#if 0
/// @detail read from file
void
SameStrand::read_from_file( std::string fn ) {

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
			std::cout << "parse error: " << line << std::endl;
			continue;
		}

		loop_fraction_( pos ) = lf;
		strand_fraction_( pos ) = ef;
		if ( std::abs( helix_fraction( pos ) - hf ) > 0.01 ) {
			tr.Warning << "inconsistency in secondary structure file at position "
				<< pos << " H ( read ) = " << hf << " 1.0-L-E ( expected ) " << helix_fraction( pos ) << std::endl;
		}
	}

}
#endif
/// @detail write to stream ( opposite from read_from_file )
void SameStrand::show( std::ostream& out ) const {
	using namespace format;
	for ( Size i = 1; i<=total_residue_; i++ ) {
		for ( Size j = 1; j<=total_residue_; j++ ) {
			out << ( eval(i,j) ? "E" : "." );
		}
		out << std::endl;
	}
}


} //protocols
} //jumping
